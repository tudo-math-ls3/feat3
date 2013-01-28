#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_ELL_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_ELL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/algorithm.hpp>
#include <kernel/util/graph.hpp>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>



namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief ELL based sparse matrix.
     *
     * \tparam Arch_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the ELL format.
     *
     * \author Dirk Ribbrock
     */
    template <typename Arch_, typename DT_>
    class SparseMatrixELL : public Container<Arch_, DT_>
    {
      private:
        /// The stride is the row count, rounded up to a multiple of the warp size
        Index _stride;
        /// Column count per row
        Index _num_cols_per_row;

        /// Column indices stored in a (cols_per_row x stride) matrix
        Index * _Aj;
        /// Nonzero values stored in a (cols_per_row x stride) matrix
        DT_ * _Ax;
        /// Length of every single row
        Index * _Arl;


        /// Row count.
        Index _rows;
        /// Column count.
        Index _columns;
        /// Our non zero element.
        DT_ _zero_element;
        /// Our non zero element count.
        Index _used_elements;

        void _read_from_ell(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          _read_from_ell(file);
          file.close();
        }

        void _read_from_ell(std::istream& file)
        {
          uint64_t size;
          uint64_t rows;
          uint64_t columns;
          uint64_t stride;
          uint64_t num_cols_per_row;
          file.read((char *)&size, sizeof(uint64_t));
          file.read((char *)&rows, sizeof(uint64_t));
          file.read((char *)&columns, sizeof(uint64_t));
          file.read((char *)&stride, sizeof(uint64_t));
          file.read((char *)&num_cols_per_row, sizeof(uint64_t));

          this->_size = Index(rows * columns);
          _rows = Index(rows);
          _columns = Index(columns);
          _stride = Index(stride);
          _num_cols_per_row = Index(num_cols_per_row);
          _zero_element = DT_(0);

          uint64_t * cAj = new uint64_t[size];
          file.read((char *)cAj, size * sizeof(uint64_t));
          _Aj = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((_num_cols_per_row * _stride) * sizeof(Index));
          for (Index i(0) ; i < size ; ++i)
            _Aj[i] = Index(cAj[i]);
          delete[] cAj;

          double * cAx = new double[size];
          file.read((char *)cAx, size * sizeof(double));

          _Ax = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((_num_cols_per_row * _stride) * sizeof(DT_));
          for (Index i(0) ; i < size ; ++i)
            _Ax[i] = DT_(cAx[i]);
          delete[] cAx;

          _Arl = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((_rows) * sizeof(Index));
          //compute row length vector
          _used_elements = 0;
          for (Index row(0) ; row < _rows ; ++row)
          {
            Index count(0);
            for (Index i(row) ; i < Index(size) ; i += Index(stride))
            {
                if (_Ax[i] == DT_(0))
                {
                  i = Index(size);
                  break;
                }
                ++count;
                ++_used_elements;
            }
            _Arl[row] = count;
          }

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(_num_cols_per_row * _stride * sizeof(DT_)));
          this->_elements_size.push_back(_num_cols_per_row * _stride);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_num_cols_per_row * _stride * sizeof(Index)));
          this->_indices_size.push_back(_num_cols_per_row * _stride);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_rows * sizeof(Index)));
          this->_indices_size.push_back(_rows);

          MemoryPool<Arch_>::upload(this->get_elements().at(0), _Ax, _num_cols_per_row * _stride * sizeof(DT_));
          MemoryPool<Arch_>::upload(this->get_indices().at(0), _Aj, _num_cols_per_row * _stride * sizeof(Index));
          MemoryPool<Arch_>::upload(this->get_indices().at(1), _Arl, _rows * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(_Ax);
          MemoryPool<Mem::Main>::instance()->release_memory(_Aj);
          MemoryPool<Mem::Main>::instance()->release_memory(_Arl);

          _Ax = this->_elements.at(0);
          _Aj = this->_indices.at(0);
          _Arl = this->_indices.at(1);
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Arch_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit SparseMatrixELL() :
          Container<Arch_, DT_> (0),
          _stride(0),
          _num_cols_per_row(0),
          _rows(0),
          _columns(0),
          _zero_element(DT_(0)),
          _used_elements(0)
        {
          CONTEXT("When creating SparseMatrixELL");
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix in CSR format.
         *
         * Creates a ELL matrix based on the CSR source matrix.
         */
        explicit SparseMatrixELL(const SparseMatrixCSR<Arch_, DT_> & other_orig) :
          Container<Arch_, DT_>(other_orig.size()),
          _rows(other_orig.rows()),
          _columns(other_orig.columns()),
          _zero_element(other_orig.zero_element()),
          _used_elements(other_orig.used_elements())
        {
          CONTEXT("When creating SparseMatrixELL");

          SparseMatrixCSR<Mem::Main, DT_> other(other_orig);

          _Arl = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((_rows) * sizeof(Index));

          _num_cols_per_row = 0;
          for (Index i(0) ; i < _rows ; ++i)
          {
            _Arl[i] = other.row_ptr_end()[i] - other.row_ptr()[i];
            if (_Arl[i] > _num_cols_per_row)
              _num_cols_per_row = _Arl[i];
          }

          Index alignment(32);
          _stride = alignment * ((_rows + alignment - 1)/ alignment);


          _Ax = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((_num_cols_per_row * _stride) * sizeof(DT_));
          MemoryPool<Mem::Main>::instance()->set_memory(_Ax, DT_(0), _num_cols_per_row * _stride);
          _Aj = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((_num_cols_per_row * _stride) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->set_memory(_Aj, Index(0), _num_cols_per_row * _stride);

          for (Index row(0); row < _rows ; ++row)
          {
            Index target(0);
            for (Index i(0) ; i < _Arl[row] ; ++i)
            {
              const Index row_start(other.row_ptr()[row]);
              //if(other.val()[row_start + i] != DT_(0))
              {
                _Aj[row + target * _stride] = (other.col_ind())[row_start + i];
                _Ax[row + target * _stride] = (other.val())[row_start + i];
                target++;
              }
            }
          }

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(_num_cols_per_row * _stride * sizeof(DT_)));
          this->_elements_size.push_back(_num_cols_per_row * _stride);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_num_cols_per_row * _stride * sizeof(Index)));
          this->_indices_size.push_back(_num_cols_per_row * _stride);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_rows * sizeof(Index)));
          this->_indices_size.push_back(_rows);

          MemoryPool<Arch_>::upload(this->get_elements().at(0), _Ax, _num_cols_per_row * _stride * sizeof(DT_));
          MemoryPool<Arch_>::upload(this->get_indices().at(0), _Aj, _num_cols_per_row * _stride * sizeof(Index));
          MemoryPool<Arch_>::upload(this->get_indices().at(1), _Arl, _rows * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(_Ax);
          MemoryPool<Mem::Main>::instance()->release_memory(_Aj);
          MemoryPool<Mem::Main>::instance()->release_memory(_Arl);

          _Ax = this->_elements.at(0);
          _Aj = this->_indices.at(0);
          _Arl = this->_indices.at(1);
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix in COO format.
         *
         * Creates a ELL matrix based on the COO source matrix.
         */
        explicit SparseMatrixELL(const SparseMatrixCOO<Mem::Main, DT_> & other_orig) :
          Container<Arch_, DT_>(other_orig.size()),
          _rows(other_orig.rows()),
          _columns(other_orig.columns()),
          _zero_element(other_orig.zero_element()),
          _used_elements(other_orig.used_elements())
        {
          CONTEXT("When creating SparseMatrixELL");

          SparseMatrixCSR<Mem::Main, DT_> other(other_orig);

          _Arl = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((_rows) * sizeof(Index));

          _num_cols_per_row = 0;
          for (Index i(0) ; i < _rows ; ++i)
          {
            _Arl[i] = other.row_ptr_end()[i] - other.row_ptr()[i];
            if (_Arl[i] > _num_cols_per_row)
              _num_cols_per_row = _Arl[i];
          }

          Index alignment(32);
          _stride = alignment * ((_rows + alignment - 1)/ alignment);


          _Ax = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((_num_cols_per_row * _stride) * sizeof(DT_));
          MemoryPool<Mem::Main>::instance()->set_memory(_Ax, DT_(0), _num_cols_per_row * _stride);
          _Aj = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((_num_cols_per_row * _stride) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->set_memory(_Aj, Index(0), _num_cols_per_row * _stride);

          for (Index row(0); row < _rows ; ++row)
          {
            Index target(0);
            for (Index i(0) ; i < _Arl[row] ; ++i)
            {
              const Index row_start(other.row_ptr()[row]);
              //if(other.val()[row_start + i] != DT_(0))
              {
                _Aj[row + target * _stride] = (other.col_ind())[row_start + i];
                _Ax[row + target * _stride] = (other.val())[row_start + i];
                target++;
              }
            }
          }

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(_num_cols_per_row * _stride * sizeof(DT_)));
          this->_elements_size.push_back(_num_cols_per_row * _stride);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_num_cols_per_row * _stride * sizeof(Index)));
          this->_indices_size.push_back(_num_cols_per_row * _stride);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_rows * sizeof(Index)));
          this->_indices_size.push_back(_rows);


          MemoryPool<Arch_>::upload(this->get_elements().at(0), _Ax, _num_cols_per_row * _stride * sizeof(DT_));
          MemoryPool<Arch_>::upload(this->get_indices().at(0), _Aj, _num_cols_per_row * _stride * sizeof(Index));
          MemoryPool<Arch_>::upload(this->get_indices().at(1), _Arl, _rows * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(_Ax);
          MemoryPool<Mem::Main>::instance()->release_memory(_Aj);
          MemoryPool<Mem::Main>::instance()->release_memory(_Arl);

          _Ax = this->_elements.at(0);
          _Aj = this->_indices.at(0);
          _Arl = this->_indices.at(1);

        }

        /**
         * \brief Constructor
         *
         * \param[in] filename The source file in HONEI ELL format.
         *
         * Creates a ELL matrix based on the source file.
         */
        explicit SparseMatrixELL(String filename) :
          Container<Arch_, DT_>(0)
        {
          CONTEXT("When creating SparseMatrixELL");

          _read_from_ell(filename);
        }

        /**
         * \brief Constructor
         *
         * \param[in] file The stream that shall be read from.
         *
         * Creates a ELL matrix based on the source file.
         */
        explicit SparseMatrixELL(std::istream& file) :
          Container<Arch_, DT_>(0)
        {
          CONTEXT("When creating SparseMatrixELL");

          _read_from_ell(file);
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a shallow copy of a given matrix.
         */
        SparseMatrixELL(const SparseMatrixELL<Arch_, DT_> & other) :
          Container<Arch_, DT_>(other),
          _stride(other._stride),
          _num_cols_per_row(other._num_cols_per_row),
          _rows(other._rows),
          _columns(other._columns),
          _zero_element(other._zero_element),
          _used_elements(other._used_elements)
        {
          CONTEXT("When copying SparseMatrixELL");

          this->_Ax = this->_elements.at(0);
          this->_Aj = this->_indices.at(0);
          this->_Arl = this->_indices.at(1);
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a copy of a given matrix from another memory architecture.
         */
        template <typename Arch2_, typename DT2_>
        SparseMatrixELL(const SparseMatrixELL<Arch2_, DT2_> & other) :
          Container<Arch_, DT_>(other),
          _stride(other.stride()),
          _num_cols_per_row(other.num_cols_per_row()),
          _rows(other.rows()),
          _columns(other.columns()),
          _zero_element(other.zero_element()),
          _used_elements(other.used_elements())
        {
          CONTEXT("When copying SparseMatrixELL");

          this->_Ax = this->_elements.at(0);
          this->_Aj = this->_indices.at(0);
          this->_Arl = this->_indices.at(1);
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        SparseMatrixELL<Arch_, DT_> clone()
        {
          CONTEXT("When cloning SparseMatrixELL");


          SparseMatrixCOO<Arch_, DT_> tcoo(*this);
          SparseMatrixELL<Arch_, DT_> t(tcoo);

          return t;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns another matrix to the target matrix.
         */
        SparseMatrixELL<Arch_, DT_> & operator= (const SparseMatrixELL<Arch_, DT_> & other)
        {
          CONTEXT("When assigning SparseMatrixELL");

          if (this == &other)
            return *this;

          this->_size = other.size();
          this->_rows = other.rows();
          this->_columns = other.columns();
          this->_used_elements = other.used_elements();
          this->_zero_element = other._zero_element;
          this->_stride = other.stride();
          this->_num_cols_per_row = other.num_cols_per_row();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();

          std::vector<DT_ *> new_elements = other.get_elements();
          std::vector<Index *> new_indices = other.get_indices();

          this->_elements.assign(new_elements.begin(), new_elements.end());
          this->_indices.assign(new_indices.begin(), new_indices.end());
          this->_elements_size.assign(other.get_elements_size().begin(), other.get_elements_size().end());
          this->_indices_size.assign(other.get_indices_size().begin(), other.get_indices_size().end());

          this->_Ax = this->_elements.at(0);
          this->_Aj = this->_indices.at(0);
          this->_Arl = this->_indices.at(1);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

          return *this;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns a matrix from another memory architecture to the target matrix.
         */
        template <typename Arch2_, typename DT2_>
        SparseMatrixELL<Arch_, DT_> & operator= (const SparseMatrixELL<Arch2_, DT2_> & other)
        {
          CONTEXT("When assigning SparseMatrixELL");

          this->_size = other.size();
          this->_rows = other.rows();
          this->_columns = other.columns();
          this->_used_elements = other.used_elements();
          this->_zero_element = other.zero_element();
          this->_stride = other.stride();
          this->_num_cols_per_row = other.num_cols_per_row();
          this->_stride = other.stride();
          this->_num_cols_per_row = other.num_cols_per_row();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(_num_cols_per_row * _stride * sizeof(DT_)));
          this->_elements_size.push_back(_num_cols_per_row * _stride);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_num_cols_per_row * _stride * sizeof(Index)));
          this->_indices_size.push_back(_num_cols_per_row * _stride);
          this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(_rows * sizeof(Index)));
          this->_indices_size.push_back(_rows);

          this->_Ax = this->_elements.at(0);
          this->_Aj = this->_indices.at(0);
          this->_Arl = this->_indices.at(1);

          Index src_size(other.get_elements_size().at(0) * sizeof(DT2_));
          Index dest_size(other.get_elements_size().at(0) * sizeof(DT_));
          void * temp(::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_elements().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_elements().at(0), temp, dest_size);
          ::free(temp);

          src_size = (other.get_indices_size().at(0) * sizeof(Index));
          dest_size = (other.get_indices_size().at(0) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_indices().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(0), temp, dest_size);
          ::free(temp);

          src_size = (other.get_indices_size().at(1) * sizeof(Index));
          dest_size = (other.get_indices_size().at(1) * sizeof(Index));
          temp = (::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_indices().at(1), src_size);
          MemoryPool<Arch_>::upload(this->get_indices().at(1), temp, dest_size);
          ::free(temp);

          return *this;
        }

        /**
         * \brief Write out matrix to file.
         *
         * \param[in] mode The used file format.
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out(FileMode mode, String filename) const
        {
          CONTEXT("When writing out SparseMatrixELL");

          switch(mode)
          {
            case fm_ell:
              write_out_ell(filename);
              break;
            case fm_m:
              write_out_m(filename);
              break;
            default:
                throw InternalError("Filemode not supported!");
          }
        }

        /**
         * \brief Write out matrix to file.
         *
         * \param[in] mode The used file format.
         * \param[in] file The stream that shall be written to.
         */
        void write_out(FileMode mode, std::ostream& file) const
        {
          CONTEXT("When writing out SparseMatrixELL");

          switch(mode)
          {
            case fm_ell:
              write_out_ell(file);
              break;
            case fm_m:
              write_out_m(file);
              break;
            default:
                throw InternalError("Filemode not supported!");
          }
        }

        /**
         * \brief Write out matrix to ell binary file.
         *
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out_ell(String filename) const
        {
          std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          write_out_ell(file);
          file.close();
        }

        /**
         * \brief Write out matrix to ell binary file.
         *
         * \param[in] file The stream that shall be written to.
         */
        void write_out_ell(std::ostream& file) const
        {
          if (typeid(DT_) != typeid(double))
            std::cout<<"Warning: You are writing out an ell matrix with less than double precission!"<<std::endl;

          Index * Aj = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((_num_cols_per_row * _stride) * sizeof(Index));
          MemoryPool<Arch_>::download(Aj, _Aj, _num_cols_per_row * _stride * sizeof(Index));
          uint64_t * cAj = new uint64_t[_num_cols_per_row * _stride];
          for (Index i(0) ; i < _num_cols_per_row * _stride ; ++i)
            cAj[i] = Aj[i];
          MemoryPool<Mem::Main>::instance()->release_memory(Aj);

          DT_ * Ax = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((_num_cols_per_row * _stride) * sizeof(DT_));
          MemoryPool<Arch_>::download(Ax, _Ax, _num_cols_per_row * _stride * sizeof(DT_));
          double * cAx = new double[_num_cols_per_row * _stride];
          for (Index i(0) ; i < _num_cols_per_row * _stride ; ++i)
            cAx[i] = Ax[i];
          MemoryPool<Mem::Main>::instance()->release_memory(Ax);

          uint64_t size(_num_cols_per_row * _stride);
          uint64_t rows(_rows);
          uint64_t columns(_columns);
          uint64_t stride(_stride);
          uint64_t num_cols_per_row(_num_cols_per_row);
          file.write((const char *)&size, sizeof(uint64_t));
          file.write((const char *)&rows, sizeof(uint64_t));
          file.write((const char *)&columns, sizeof(uint64_t));
          file.write((const char *)&stride, sizeof(uint64_t));
          file.write((const char *)&num_cols_per_row, sizeof(uint64_t));
          file.write((const char *)cAj, size * sizeof(uint64_t));
          file.write((const char *)cAx, size * sizeof(double));

          delete[] cAj;
          delete[] cAx;
        }

        /**
         * \brief Write out matrix to matlab m file.
         *
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out_m(String filename) const
        {
          std::ofstream file(filename.c_str(), std::ofstream::out);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          write_out_m(file);
          file.close();
        }

        /**
         * \brief Write out matrix to matlab m file.
         *
         * \param[in] file The stream that shall be written to.
         */
        void write_out_m(std::ostream& file) const
        {
          SparseMatrixELL<Mem::Main, DT_> temp(*this);

          file << "data = [" << std::endl;
          for (Index i(0) ; i < _num_cols_per_row * _stride ; ++i)
          {
            if (temp.Ax()[i] != DT_(0))
            {
              file << (i%_stride) + 1 << " " << temp.Aj()[i] + 1 << " " << std::scientific << (double)temp.Ax()[i] << ";" << std::endl;
            }
          }
          file << "];" << std::endl;
          file << "mat=sparse(data(:,1),data(:,2),data(:,3));";
        }

        /**
         * \brief Retrieve specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         *
         * \returns Specific matrix element.
         */
        DT_ operator()(Index row, Index col) const
        {
          CONTEXT("When retrieving SparseMatrixELL element");

          ASSERT(row < this->_rows, "Error: " + stringify(row) + " exceeds sparse matrix ell row size " + stringify(this->_rows) + " !");
          ASSERT(col < this->_columns, "Error: " + stringify(col) + " exceeds sparse matrix ell column size " + stringify(this->_columns) + " !");

          if (typeid(Arch_) == typeid(Mem::Main))
          {
            Index max(this->_Arl[row]);
            for (Index i(row), j(0) ; j < max && this->_Aj[i] <= col ; i += this->_stride, ++j)
            {
                if (this->_Aj[i] == col)
                  return this->_Ax[i];
            }
            return this->_zero_element;
          }
          else
          {
            SparseMatrixELL<Mem::Main, DT_> temp(*this);
            return temp(row, col);
          }
        }

        /**
         * \brief Retrieve matrix row count.
         *
         * \returns Matrix row count.
         */
        const Index & rows() const
        {
          return this->_rows;
        }

        /**
         * \brief Retrieve matrix column count.
         *
         * \returns Matrix column count.
         */
        const Index & columns() const
        {
          return this->_columns;
        }

        /**
         * \brief Retrieve non zero element count.
         *
         * \returns Non zero element count.
         */
        const Index & used_elements() const
        {
          return this->_used_elements;
        }

        /**
         * \brief Retrieve column indices array.
         *
         * \returns Column indices array.
         */
        Index * Aj() const
        {
          return _Aj;
        }

        /**
         * \brief Retrieve non zero element array.
         *
         * \returns Non zero element array.
         */
        DT_ * Ax() const
        {
          return _Ax;
        }

        /**
         * \brief Retrieve row length array.
         *
         * \returns Row lenght array.
         */
        Index * Arl() const
        {
          return _Arl;
        }

        /**
         * \brief Retrieve non zero element.
         *
         * \returns Non zero element.
         */
        DT_ zero_element() const
        {
          return _zero_element;
        }

        /**
         * \brief Retrieve stride, i.e. lowest common multiple of row count and warp size.
         *
         * \returns Stride.
         */
        Index stride() const
        {
          return _stride;
        }

        /**
         * \brief Retrieve the maximum amount of non zero columns in a single row.
         *
         * \returns Columns per row count.
         */
        Index num_cols_per_row() const
        {
          return _num_cols_per_row;
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String type_name()
        {
          return "SparseMatrixELL";
        }
    };

    /**
     * \brief SparseMatrixELL comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Arch_, typename Arch2_, typename DT_> bool operator== (const SparseMatrixELL<Arch_, DT_> & a, const SparseMatrixELL<Arch2_, DT_> & b)
    {
      CONTEXT("When comparing SparseMatrixELLs");

      if (a.rows() != b.rows())
        return false;
      if (a.columns() != b.columns())
        return false;
      if (a.used_elements() != b.used_elements())
        return false;
      if (a.zero_element() != b.zero_element())
        return false;

      for (Index i(0) ; i < a.rows() ; ++i)
      {
        for (Index j(0) ; j < a.columns() ; ++j)
        {
          if (a(i, j) != b(i, j))
            return false;
        }
      }

      return true;
    }

    /**
     * \brief SparseMatrixELL streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Arch_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixELL<Arch_, DT_> & b)
    {
      lhs << "[" << std::endl;
      for (Index i(0) ; i < b.rows() ; ++i)
      {
        lhs << "[";
        for (Index j(0) ; j < b.columns() ; ++j)
        {
          lhs << "  " << b(i, j);
        }
        lhs << "]" << std::endl;
      }
      lhs << "]" << std::endl;

      return lhs;
    }

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_MATRIX_ELL_HPP
