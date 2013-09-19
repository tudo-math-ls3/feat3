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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>



namespace FEAST
{
  namespace LAFEM
  {
    //forward declarations
    template <typename Mem_, typename DT_>
    class SparseMatrixCSR;

    template <typename Mem_, typename DT_>
    class SparseMatrixCOO;

    /**
     * \brief ELL based sparse matrix.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the ELL-R format.\n\n
     * Data survey: \n
     * _elements[0]: Ax - raw non zero number values, stored in a (cols per row x stride) matrix \n
     * _indices[0]: Aj - column index per non zero element, stored in a (cols per row x stride) matrix \n
     * _indices[1]: Arl - length of every single row\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: stride, aka the row count, rounded up to a multiple of the warp size \n
     * _scalar_index[4]: column count per row \n
     * _scalar_index[5]: non zero element count (used elements) \n
     * _scalar_dt[0]: zero element
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_>
    class SparseMatrixELL : public Container<Mem_, DT_>
    {
      private:
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

          this->_scalar_index.at(0) = Index(rows * columns);
          this->_scalar_index.push_back((Index)rows);
          this->_scalar_index.push_back((Index)columns);
          this->_scalar_index.push_back((Index)stride);
          this->_scalar_index.push_back((Index)num_cols_per_row);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));

          uint64_t * cAj = new uint64_t[size];
          file.read((char *)cAj, (long)(size * sizeof(uint64_t)));
          Index* tAj = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(4) * this->_scalar_index.at(3)) * sizeof(Index));
          for (Index i(0) ; i < size ; ++i)
            tAj[i] = Index(cAj[i]);
          delete[] cAj;

          double * cAx = new double[size];
          file.read((char *)cAx, (long)(size * sizeof(double)));

          DT_* tAx = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(4) * this->_scalar_index.at(3)) * sizeof(DT_));
          for (Index i(0) ; i < size ; ++i)
            tAx[i] = DT_(cAx[i]);
          delete[] cAx;

          Index* tArl = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(1)) * sizeof(Index));
          //compute row length vector
          this->_scalar_index.at(5) = 0;
          for (Index row(0) ; row < this->_scalar_index.at(1) ; ++row)
          {
            Index count(0);
            for (Index i(row) ; i < Index(size) ; i += Index(stride))
            {
                if (tAx[i] == DT_(0))
                {
                  i = Index(size);
                  break;
                }
                ++count;
                ++this->_scalar_index.at(5);
            }
            tArl[row] = count;
          }

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(1) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::upload(this->get_elements().at(0), tAx, this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(DT_));
          MemoryPool<Mem_>::upload(this->get_indices().at(0), tAj, this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(Index));
          MemoryPool<Mem_>::upload(this->get_indices().at(1), tArl, this->_scalar_index.at(1) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(tAx);
          MemoryPool<Mem::Main>::instance()->release_memory(tAj);
          MemoryPool<Mem::Main>::instance()->release_memory(tArl);
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Mem_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit SparseMatrixELL() :
          Container<Mem_, DT_> (0)
        {
          CONTEXT("When creating SparseMatrixELL");
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix in CSR format.
         *
         * Creates a ELL matrix based on the CSR source matrix.
         */
        template <typename Arch2_>
        explicit SparseMatrixELL(const SparseMatrixCSR<Arch2_, DT_> & other_orig) :
          Container<Mem_, DT_>(other_orig.size())
        {
          CONTEXT("When creating SparseMatrixELL");
          this->_scalar_index.push_back(other_orig.rows());
          this->_scalar_index.push_back(other_orig.columns());
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(other_orig.used_elements());
          this->_scalar_dt.push_back(other_orig.zero_element());

          SparseMatrixCSR<Mem::Main, DT_> other(other_orig);

          Index* tArl = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(1)) * sizeof(Index));

          this->_scalar_index.at(4) = 0;
          for (Index i(0) ; i < this->_scalar_index.at(1) ; ++i)
          {
            tArl[i] = other.row_ptr_end()[i] - other.row_ptr()[i];
            if (tArl[i] > this->_scalar_index.at(4))
              this->_scalar_index.at(4) = tArl[i];
          }

          Index alignment(32);
          this->_scalar_index.at(3) = alignment * ((this->_scalar_index.at(1) + alignment - 1)/ alignment);


          DT_* tAx = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(4) * this->_scalar_index.at(3)) * sizeof(DT_));
          MemoryPool<Mem::Main>::instance()->set_memory(tAx, DT_(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));
          Index* tAj = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(4) * this->_scalar_index.at(3)) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->set_memory(tAj, Index(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));

          for (Index row(0); row < this->_scalar_index.at(1) ; ++row)
          {
            Index target(0);
            for (Index i(0) ; i < tArl[row] ; ++i)
            {
              const Index row_start(other.row_ptr()[row]);
              //if(other.val()[row_start + i] != DT_(0))
              {
                tAj[row + target * this->_scalar_index.at(3)] = (other.col_ind())[row_start + i];
                tAx[row + target * this->_scalar_index.at(3)] = (other.val())[row_start + i];
                target++;
              }
            }
          }

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(1) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::upload(this->get_elements().at(0), tAx, this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(DT_));
          MemoryPool<Mem_>::upload(this->get_indices().at(0), tAj, this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(Index));
          MemoryPool<Mem_>::upload(this->get_indices().at(1), tArl, this->_scalar_index.at(1) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(tAx);
          MemoryPool<Mem::Main>::instance()->release_memory(tAj);
          MemoryPool<Mem::Main>::instance()->release_memory(tArl);
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix in COO format.
         *
         * Creates a ELL matrix based on the COO source matrix.
         */
        template <typename Mem2_>
        explicit SparseMatrixELL(const SparseMatrixCOO<Mem2_, DT_> & other) :
          Container<Mem_, DT_>(other.size())
        {
          CONTEXT("When creating SparseMatrixELL");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_dt.push_back(other.zero_element());

          SparseMatrixCOO<Mem::Main, DT_> cother(other);

          Index* tArl = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(1)) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->set_memory(tArl, Index(0), this->_scalar_index.at(1));

          this->_scalar_index.at(4) = 0;
          for (Index i(0) ; i < this->_scalar_index.at(5) ; ++i)
          {
            Index cur_row(cother.row()[i]);
            ++tArl[cur_row];
            if (tArl[cur_row] > this->_scalar_index.at(4))
              this->_scalar_index.at(4) = tArl[cur_row];
          }

          Index alignment(32);
          this->_scalar_index.at(3) = alignment * ((this->_scalar_index.at(1) + alignment - 1)/ alignment);


          DT_* tAx = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(4) * this->_scalar_index.at(3)) * sizeof(DT_));
          MemoryPool<Mem::Main>::instance()->set_memory(tAx, DT_(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));
          Index* tAj = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(4) * this->_scalar_index.at(3)) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->set_memory(tAj, Index(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));

          Index last_row(cother.row()[0]);
          Index target(0);
          for (Index i(0) ; i < this->_scalar_index.at(5) ; ++i)
          {
            Index row(cother.row()[i]);
            if (row != last_row)
            {
              target = 0;
              last_row = row;
            }
            tAj[row + target * this->_scalar_index.at(3)] = (cother.column())[i];
            tAx[row + target * this->_scalar_index.at(3)] = (cother.val())[i];
            target++;
          }

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(1) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(1));


          MemoryPool<Mem_>::upload(this->get_elements().at(0), tAx, this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(DT_));
          MemoryPool<Mem_>::upload(this->get_indices().at(0), tAj, this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(Index));
          MemoryPool<Mem_>::upload(this->get_indices().at(1), tArl, this->_scalar_index.at(1) * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(tAx);
          MemoryPool<Mem::Main>::instance()->release_memory(tAj);
          MemoryPool<Mem::Main>::instance()->release_memory(tArl);
        }

        /**
         * \brief Constructor
         *
         * \param[in] filename The source file in HONEI ELL format.
         *
         * Creates a ELL matrix based on the source file.
         */
        explicit SparseMatrixELL(String filename) :
          Container<Mem_, DT_>(0)
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
          Container<Mem_, DT_>(0)
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
        SparseMatrixELL(const SparseMatrixELL<Mem_, DT_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseMatrixELL");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.stride());
          this->_scalar_index.push_back(other.num_cols_per_row());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_dt.push_back(other.zero_element());
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
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseMatrixELL");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.stride());
          this->_scalar_index.push_back(other.num_cols_per_row());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_dt.push_back(other.zero_element());
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        SparseMatrixELL<Mem_, DT_> clone()
        {
          CONTEXT("When cloning SparseMatrixELL");

          SparseMatrixELL<Mem_, DT_> t;
          ((Container<Mem_, DT_>&)t).clone(*this);
          return t;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns another matrix to the target matrix.
         */
        SparseMatrixELL<Mem_, DT_> & operator= (const SparseMatrixELL<Mem_, DT_> & other)
        {
          CONTEXT("When assigning SparseMatrixELL");

         this->assign(other);

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
        SparseMatrixELL<Mem_, DT_> & operator= (const SparseMatrixELL<Arch2_, DT2_> & other)
        {
          CONTEXT("When assigning SparseMatrixELL");

         this->assign(other);

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

          Index * Aj = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(4) * this->_scalar_index.at(3)) * sizeof(Index));
          MemoryPool<Mem_>::download(Aj, this->_indices.at(0), this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(Index));
          uint64_t * cAj = new uint64_t[this->_scalar_index.at(4) * this->_scalar_index.at(3)];
          for (Index i(0) ; i < this->_scalar_index.at(4) * this->_scalar_index.at(3) ; ++i)
            cAj[i] = Aj[i];
          MemoryPool<Mem::Main>::instance()->release_memory(Aj);

          DT_ * Ax = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(4) * this->_scalar_index.at(3)) * sizeof(DT_));
          MemoryPool<Mem_>::download(Ax, this->_elements.at(0), this->_scalar_index.at(4) * this->_scalar_index.at(3) * sizeof(DT_));
          double * cAx = new double[this->_scalar_index.at(4) * this->_scalar_index.at(3)];
          for (Index i(0) ; i < this->_scalar_index.at(4) * this->_scalar_index.at(3) ; ++i)
            cAx[i] = Ax[i];
          MemoryPool<Mem::Main>::instance()->release_memory(Ax);

          uint64_t size(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          uint64_t rows(this->_scalar_index.at(1));
          uint64_t columns(this->_scalar_index.at(2));
          uint64_t stride(this->_scalar_index.at(3));
          uint64_t num_cols_per_row(this->_scalar_index.at(4));
          file.write((const char *)&size, sizeof(uint64_t));
          file.write((const char *)&rows, sizeof(uint64_t));
          file.write((const char *)&columns, sizeof(uint64_t));
          file.write((const char *)&stride, sizeof(uint64_t));
          file.write((const char *)&num_cols_per_row, sizeof(uint64_t));
          file.write((const char *)cAj, (long)(size * sizeof(uint64_t)));
          file.write((const char *)cAx, (long)(size * sizeof(double)));

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
          for (Index i(0) ; i < this->_scalar_index.at(4) * this->_scalar_index.at(3) ; ++i)
          {
            if (temp.Ax()[i] != DT_(0))
            {
              file << (i%this->_scalar_index.at(3)) + 1 << " " << temp.Aj()[i] + 1 << " " << std::scientific << (double)temp.Ax()[i] << ";" << std::endl;
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

          ASSERT(row < this->_scalar_index.at(1), "Error: " + stringify(row) + " exceeds sparse matrix ell row size " + stringify(this->_scalar_index.at(1)) + " !");
          ASSERT(col < this->_scalar_index.at(2), "Error: " + stringify(col) + " exceeds sparse matrix ell column size " + stringify(this->_scalar_index.at(2)) + " !");

          if (typeid(Mem_) == typeid(Mem::Main))
          {
            Index max(this->_indices.at(1)[row]);
            for (Index i(row), j(0) ; j < max && this->_indices.at(0)[i] <= col ; i += this->_scalar_index.at(3), ++j)
            {
                if (this->_indices.at(0)[i] == col)
                  return this->_elements.at(0)[i];
            }
            return this->_scalar_dt.at(0);
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
          return this->_scalar_index.at(1);
        }

        /**
         * \brief Retrieve matrix column count.
         *
         * \returns Matrix column count.
         */
        const Index & columns() const
        {
          return this->_scalar_index.at(2);
        }

        /**
         * \brief Retrieve non zero element count.
         *
         * \returns Non zero element count.
         */
        virtual Index used_elements() const
        {
          return this->_scalar_index.at(5);
        }

        /**
         * \brief Retrieve column indices array.
         *
         * \returns Column indices array.
         */
        Index * Aj() const
        {
          return this->_indices.at(0);
        }

        /**
         * \brief Retrieve non zero element array.
         *
         * \returns Non zero element array.
         */
        DT_ * Ax() const
        {
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve row length array.
         *
         * \returns Row lenght array.
         */
        Index * Arl() const
        {
          return this->_indices.at(1);
        }

        /**
         * \brief Retrieve non zero element.
         *
         * \returns Non zero element.
         */
        DT_ zero_element() const
        {
          return this->_scalar_dt.at(0);
        }

        /**
         * \brief Retrieve stride, i.e. lowest common multiple of row count and warp size.
         *
         * \returns Stride.
         */
        Index stride() const
        {
          return this->_scalar_index.at(3);
        }

        /**
         * \brief Retrieve the maximum amount of non zero columns in a single row.
         *
         * \returns Columns per row count.
         */
        Index num_cols_per_row() const
        {
          return this->_scalar_index.at(4);
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
    template <typename Mem_, typename Arch2_, typename DT_> bool operator== (const SparseMatrixELL<Mem_, DT_> & a, const SparseMatrixELL<Arch2_, DT_> & b)
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
      if (a.stride() != b.stride())
        return false;
      if (a.num_cols_per_row() != b.num_cols_per_row())
        return false;

      for (Index i(0) ; i < a.stride() * a.num_cols_per_row() ; ++i)
      {
        if (MemoryPool<Mem_>::get_element(a.Ax(), i) != MemoryPool<Arch2_>::get_element(b.Ax(), i))
          return false;
        if (MemoryPool<Mem_>::get_element(a.Aj(), i) != MemoryPool<Arch2_>::get_element(b.Aj(), i))
          return false;
      }

      return true;
    }

    /**
     * \brief SparseMatrixELL streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Mem_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixELL<Mem_, DT_> & b)
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
