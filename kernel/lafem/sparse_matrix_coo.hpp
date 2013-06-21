#pragma once
#ifndef KERNEL_LAFEM_DENSE_MATRIX_HPP
#define KERNEL_LAFEM_DENSE_MATRIX_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/matrix_base.hpp>
#include <kernel/util/graph.hpp>

#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>


namespace FEAST
{
  namespace LAFEM
  {
    //forward declarations
    template <typename Mem_, typename DT_>
    class SparseMatrixELL;

    template <typename Mem_, typename DT_>
    class SparseMatrixCSR;

    /**
     * \brief Coordinate based sparse matrix.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements alongside with its coordinates explicitly. \n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: row index \n
     * _indices[1]: column index \n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_dt[0]: zero element
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_>
    class SparseMatrixCOO : public Container<Mem_, DT_>, public MatrixBase
    {
      private:
        void _read_from_m(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          _read_from_m(file);
          file.close();
        }

        void _read_from_m(std::istream& file)
        {
          std::vector<Index> rowsv;
          std::vector<Index> colsv;
          std::vector<DT_> valsv;
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));

          String line;
          std::getline(file, line);
          while(!file.eof())
          {
            std::getline(file, line);

            if(line.find("]", 0) < line.npos)
              break;

            if(line[line.size()-1] == ';')
              line.resize(line.size()-1);

            String::size_type begin(line.find_first_not_of(" "));
            line.erase(0, begin);
            String::size_type end(line.find_first_of(" "));
            String srow(line, 0, end);
            Index row(atol(srow.c_str()));
            --row;
            this->_scalar_index.at(1) = std::max(row+1, this->rows());
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scol(line, 0, end);
            Index col(atol(scol.c_str()));
            --col;
            this->_scalar_index.at(2) = std::max(col+1, this->columns());
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String sval(line, 0, end);
            DT_ val((DT_)atof(sval.c_str()));

            rowsv.push_back(row);
            colsv.push_back(col);
            valsv.push_back(val);

          }
          this->_scalar_index.at(0) = this->rows() * this->columns();

          for (Index i(0) ; i < rowsv.size() ; ++i)
          {
            (*this)(rowsv.at(i), colsv.at(i), valsv.at(i));
          }
        }

        void _read_from_mtx(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          _read_from_mtx(file);
          file.close();
        }

        void _read_from_mtx(std::istream& file)
        {
          std::vector<Index> rowsv;
          std::vector<Index> colsv;
          std::vector<DT_> valsv;
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));

          String line;
          std::getline(file, line);
          {
            std::getline(file, line);

            String::size_type begin(line.find_first_not_of(" "));
            line.erase(0, begin);
            String::size_type end(line.find_first_of(" "));
            String srow(line, 0, end);
            Index row(atol(srow.c_str()));
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scol(line, 0, end);
            Index col(atol(scol.c_str()));
            line.erase(0, end);
            this->_scalar_index.at(1) = row;
            this->_scalar_index.at(2) = col;
            this->_scalar_index.at(0) = this->rows() * this->columns();
          }
          while(!file.eof())
          {
            std::getline(file, line);
            if (file.eof())
              break;

            String::size_type begin(line.find_first_not_of(" "));
            line.erase(0, begin);
            String::size_type end(line.find_first_of(" "));
            String srow(line, 0, end);
            Index row(atol(srow.c_str()));
            --row;
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scol(line, 0, end);
            Index col(atol(scol.c_str()));
            --col;
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String sval(line, 0, end);
            DT_ val((DT_)atof(sval.c_str()));

            rowsv.push_back(row);
            colsv.push_back(col);
            valsv.push_back(val);
          }

          for (Index i(0) ; i < rowsv.size() ; ++i)
          {
            (*this)(rowsv.at(i), colsv.at(i), valsv.at(i));
          }
        }

        void _read_from_coo(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          _read_from_coo(file);
          file.close();
        }

        void _read_from_coo(std::istream& file)
        {
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));

          uint64_t rows;
          uint64_t columns;
          uint64_t elements64;
          file.read((char *)&rows, sizeof(uint64_t));
          file.read((char *)&columns, sizeof(uint64_t));
          file.read((char *)&elements64, sizeof(uint64_t));
          Index elements = (Index)elements64;

          this->_scalar_index.at(0) = Index(rows * columns);
          this->_scalar_index.at(1) = Index(rows);
          this->_scalar_index.at(2) = Index(columns);
          this->_scalar_index.at(3) = (Index)elements;

          uint64_t * crow_ptr = new uint64_t[elements];
          file.read((char *)crow_ptr, (elements) * sizeof(uint64_t));
          Index * trow_ptr = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((elements) * sizeof(Index));
          for (Index i(0) ; i < elements ; ++i)
            trow_ptr[i] = Index(crow_ptr[i]);
          delete[] crow_ptr;

          uint64_t * ccol_ptr = new uint64_t[elements];
          file.read((char *)ccol_ptr, (elements) * sizeof(uint64_t));
          Index * tcol_ptr = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory((elements) * sizeof(Index));
          for (Index i(0) ; i < elements ; ++i)
            tcol_ptr[i] = Index(ccol_ptr[i]);
          delete[] ccol_ptr;

          double * cval = new double[elements];
          file.read((char *)cval, elements * sizeof(double));
           DT_ * tval_ptr = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((elements) * sizeof(DT_));
          for (Index i(0) ; i < elements ; ++i)
            tval_ptr[i] = DT_(cval[i]);
          delete[] cval;

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(elements * sizeof(DT_)));
          this->_elements_size.push_back(elements);
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory((elements)* sizeof(Index)));
          this->_indices_size.push_back(elements);
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory((elements)* sizeof(Index)));
          this->_indices_size.push_back(elements);

          MemoryPool<Mem_>::upload(this->get_elements().at(0), tval_ptr, elements * sizeof(DT_));
          MemoryPool<Mem_>::upload(this->get_indices().at(0), trow_ptr, elements * sizeof(Index));
          MemoryPool<Mem_>::upload(this->get_indices().at(1), tcol_ptr, elements * sizeof(Index));
          MemoryPool<Mem::Main>::instance()->release_memory(tval_ptr);
          MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr);
          MemoryPool<Mem::Main>::instance()->release_memory(tcol_ptr);
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
        explicit SparseMatrixCOO() :
          Container<Mem_, DT_> (0)
        {
          CONTEXT("When creating SparseMatrixCOO");
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] dimensions The row/column count of the created matrix.
         *
         * Creates a matrix with given dimensions.
         */
        explicit SparseMatrixCOO(Index dimensions) :
          Container<Mem_, DT_>(dimensions * dimensions)
        {
          CONTEXT("When creating SparseMatrixCOO");
          this->_scalar_index.push_back(dimensions);
          this->_scalar_index.push_back(dimensions);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] rows The row count of the created matrix.
         * \param[in] columns The column count of the created matrix.
         *
         * Creates a matrix with given dimensions.
         */
        explicit SparseMatrixCOO(Index rows, Index columns) :
          Container<Mem_, DT_>(rows * columns)
        {
          CONTEXT("When creating SparseMatrixCOO");
          this->_scalar_index.push_back(rows);
          this->_scalar_index.push_back(columns);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source ell matrix.
         *
         * Creates a matrix from the given source matrix.
         */
        template <typename Mem2_>
        explicit SparseMatrixCOO(const SparseMatrixELL<Mem2_, DT_> & other) :
          Container<Mem_, DT_>(other.rows() * other.columns())
        {
          CONTEXT("When creating SparseMatrixCOO from SparseMatrixELL");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));

          SparseMatrixELL<Mem::Main, DT_> cother(other);

          SparseMatrixCOO<Mem::Main, DT_> coo(other.rows(), other.columns());

          for (Index i(0) ; i < other.num_cols_per_row() * other.stride() ; ++i)
          {
            if (cother.Ax()[i] != DT_(0))
            {
              coo(i%other.stride(), cother.Aj()[i], cother.Ax()[i]);
              ++this->_scalar_index.at(3);
            }
          }

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(3) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(3));

          MemoryPool<Mem_>::upload(this->_elements.at(0), coo.val(), this->_scalar_index.at(3) * sizeof(DT_));
          MemoryPool<Mem_>::upload(this->_indices.at(0), coo.row(), this->_scalar_index.at(3) * sizeof(Index));
          MemoryPool<Mem_>::upload(this->_indices.at(1), coo.column(), this->_scalar_index.at(3) * sizeof(Index));
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source csr matrix.
         *
         * Creates a matrix from the given source matrix.
         */
        template <typename Mem2_>
        explicit SparseMatrixCOO(const SparseMatrixCSR<Mem2_, DT_> & other) :
          Container<Mem_, DT_>(other.rows() * other.columns())
        {
          CONTEXT("When creating SparseMatrixCOO from SparseMatrixCSR");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));

          SparseMatrixCSR<Mem::Main, DT_> cother(other);

          SparseMatrixCOO<Mem::Main, DT_> coo(other.rows(), other.columns());

          for (Index row(0) ; row < other.rows() ; ++row)
          {
            const Index end(cother.row_ptr_end()[row]);
            for (Index i(cother.row_ptr()[row]) ; i < end ; ++i)
            {
              if (cother.val()[i] != DT_(0))
              {
                coo(row, cother.col_ind()[i], cother.val()[i]);
                ++this->_scalar_index.at(3);
              }
            }
          }

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(3) * sizeof(DT_)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(3) * sizeof(Index)));
          this->_indices_size.push_back(this->_scalar_index.at(3));

          MemoryPool<Mem_>::upload(this->_elements.at(0), coo.val(), this->_scalar_index.at(3) * sizeof(DT_));
          MemoryPool<Mem_>::upload(this->_indices.at(0), coo.row(), this->_scalar_index.at(3) * sizeof(Index));
          MemoryPool<Mem_>::upload(this->_indices.at(1), coo.column(), this->_scalar_index.at(3) * sizeof(Index));
        }

        /**
         * \brief Constructor
         *
         * \param[in] graph The Graph, the matrix will be created from.
         *
         * Creates a matrix from a given graph.
         */
        explicit SparseMatrixCOO(const Graph & graph) :
          Container<Mem_, DT_>(graph.get_num_nodes_domain() * graph.get_num_nodes_image())
        {
          CONTEXT("When creating SparseMatrixCOO");

          this->_scalar_index.push_back(graph.get_num_nodes_domain());
          this->_scalar_index.push_back(graph.get_num_nodes_image());
          this->_scalar_index.push_back(graph.get_num_indices());
          this->_scalar_dt.push_back(DT_(0));

          for (Index i(0) ; i < this->rows() ; ++i)
          {
            typename Graph::ImageIterator it(graph.image_begin(i));
            typename Graph::ImageIterator jt(graph.image_end(i));
            for( ; it != jt ; ++it)
            {
              (*this)(i, *it, DT_(0));
            }
          }
        }

        /**
         * \brief Constructor
         *
         * \param[in] mode The used file format.
         * \param[in] filename The source file to be read in.
         *
         * Creates a matrix based on the source file.
         */
        explicit SparseMatrixCOO(FileMode mode, String filename) :
          Container<Mem_, DT_>(0)
        {
          CONTEXT("When creating SparseMatrixCOO");

          switch(mode)
          {
            case fm_m:
              _read_from_m(filename);
              break;
            case fm_mtx:
              _read_from_mtx(filename);
              break;
            case fm_coo:
              _read_from_coo(filename);
              break;
            default:
              throw InternalError("Filemode not supported!");
          }
        }

        /**
         * \brief Constructor
         *
         * \param[in] mode The used file format.
         * \param[in] file The stream that is to be read from.
         *
         * Creates a matrix based on the source filestream.
         */
        explicit SparseMatrixCOO(FileMode mode, std::istream& file) :
          Container<Mem_, DT_>(0)
        {
          CONTEXT("When creating SparseMatrixCOO");

          switch(mode)
          {
            case fm_m:
              _read_from_m(file);
              break;
            case fm_mtx:
              _read_from_mtx(file);
              break;
            case fm_coo:
              _read_from_coo(file);
              break;
            default:
              throw InternalError("Filemode not supported!");
          }
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a shallow copy of a given matrix.
         */
        SparseMatrixCOO(const SparseMatrixCOO<Mem_, DT_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseMatrixCOO");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
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
        SparseMatrixCOO(const SparseMatrixCOO<Arch2_, DT2_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseMatrixCOO");
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_dt.push_back(other.zero_element());
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        SparseMatrixCOO<Mem_, DT_> clone()
        {
          CONTEXT("When cloning SparseMatrixCOO");

          SparseMatrixCOO<Mem_, DT_> t;
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
        SparseMatrixCOO<Mem_, DT_> & operator= (const SparseMatrixCOO<Mem_, DT_> & other)
        {
          CONTEXT("When assigning SparseMatrixCOO");

          assign(other);

          return *this;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source matrix.
         *
         * Assigns a matrix from another memory architecture to the target matrix.
         */
        template <typename Mem2_, typename DT2_>
        SparseMatrixCOO<Mem_, DT_> & operator= (const SparseMatrixCOO<Mem2_, DT2_> & other)
        {
          CONTEXT("When assigning SparseMatrixCOO");

          assign(other);

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
          CONTEXT("When writing out SparseMatrixCOO");

          switch(mode)
          {
            case fm_coo:
              write_out_coo(filename);
              break;
            case fm_m:
              write_out_m(filename);
              break;
            case fm_mtx:
              write_out_mtx(filename);
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
          CONTEXT("When writing out SparseMatrixCOO");

          switch(mode)
          {
            case fm_coo:
              write_out_coo(file);
              break;
            case fm_m:
              write_out_m(file);
              break;
            case fm_mtx:
              write_out_mtx(file);
              break;
            default:
                throw InternalError("Filemode not supported!");
          }
        }

        /**
         * \brief Write out matrix to coo binary file.
         *
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out_coo(String filename) const
        {
          std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          write_out_coo(file);
          file.close();
        }

        /**
         * \brief Write out matrix to coo binary file.
         *
         * \param[in] file The stream that shall be written to.
         */
        void write_out_coo(std::ostream& file) const
        {
          if (typeid(DT_) != typeid(double))
            std::cout<<"Warning: You are writing out an coo matrix with less than double precission!"<<std::endl;

          Index * row_ptr = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory(this->_indices_size.at(0) * sizeof(Index));
          MemoryPool<Mem_>::download(row_ptr, this->_indices.at(0), this->_indices_size.at(0) * sizeof(Index));
          uint64_t * crow_ptr = new uint64_t[this->_indices_size.at(0)];
          for (Index i(0) ; i < this->_indices_size.at(0) ; ++i)
            crow_ptr[i] = row_ptr[i];
          MemoryPool<Mem::Main>::instance()->release_memory(row_ptr);

          Index * col_ptr = (Index*)MemoryPool<Mem::Main>::instance()->allocate_memory(this->_indices_size.at(1) * sizeof(Index));
          MemoryPool<Mem_>::download(col_ptr, this->_indices.at(1), this->_indices_size.at(1) * sizeof(Index));
          uint64_t * ccol_ptr = new uint64_t[this->_indices_size.at(1)];
          for (Index i(0) ; i < this->_indices_size.at(1) ; ++i)
            ccol_ptr[i] = col_ptr[i];
          MemoryPool<Mem::Main>::instance()->release_memory(col_ptr);

          DT_ * val = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory(this->_elements_size.at(0) * sizeof(DT_));
          MemoryPool<Mem_>::download(val, this->_elements.at(0), this->_elements_size.at(0) * sizeof(DT_));
          double * cval = new double[this->_elements_size.at(0)];
          for (Index i(0) ; i < this->_elements_size.at(0) ; ++i)
            cval[i] = val[i];
          MemoryPool<Mem::Main>::instance()->release_memory(val);

          uint64_t rows(this->_scalar_index.at(1));
          uint64_t columns(this->_scalar_index.at(2));
          uint64_t elements(this->_indices_size.at(0));
          file.write((const char *)&rows, sizeof(uint64_t));
          file.write((const char *)&columns, sizeof(uint64_t));
          file.write((const char *)&elements, sizeof(uint64_t));
          file.write((const char *)crow_ptr, elements * sizeof(uint64_t));
          file.write((const char *)ccol_ptr, elements * sizeof(uint64_t));
          file.write((const char *)cval, elements * sizeof(double));

          delete[] crow_ptr;
          delete[] ccol_ptr;
          delete[] cval;
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
          SparseMatrixCOO<Mem::Main, DT_> temp(*this);

          file << "data = [" << std::endl;
          for (Index i(0) ; i < this->_scalar_index.at(3) ; ++i)
          {
            file << temp.row()[i] + 1 << " " << temp.column()[i] + 1 << " " << std::scientific << (double)temp.val()[i] << ";" << std::endl;
          }
          file << "];" << std::endl;
          file << "mat=sparse(data(:,1),data(:,2),data(:,3));";
        }

        /**
         * \brief Write out matrix to matrix market mtx file.
         *
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out_mtx(String filename) const
        {
          std::ofstream file(filename.c_str(), std::ofstream::out);
          if (! file.is_open())
            throw InternalError("Unable to open Matrix file " + filename);
          write_out_mtx(file);
          file.close();
        }

        /**
         * \brief Write out matrix to matrix market mtx file.
         *
         * \param[in] file The stream that shall be written to.
         */
        void write_out_mtx(std::ostream& file) const
        {
          SparseMatrixCOO<Mem::Main, DT_> temp(*this);

          file << "%%MatrixMarket matrix coordinate real general" << std::endl;
          file << temp.rows() << " " << temp.columns() << " " << temp.used_elements() << std::endl;

          for (Index i(0) ; i < this->_scalar_index.at(3) ; ++i)
          {
            file << temp.row()[i] + 1 << " " << temp.column()[i] + 1 << " " << std::scientific << (double)temp.val()[i] << std::endl;
          }
        }

        /**
         * \brief Set specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         * \param[in] value The value to be set.
         */
        void operator()(Index row, Index col, DT_ val)
        {
          CONTEXT("When setting SparseMatrixCOO element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds sparse matrix coo row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds sparse matrix coo column size " + stringify(this->columns()) + " !");

          if (this->_elements.size() == 0)
          {
            this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(sizeof(DT_)));
            this->_elements_size.push_back(1);
            this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(sizeof(Index)));
            this->_indices_size.push_back(1);
            this->_indices.push_back((Index*)MemoryPool<Mem_>::instance()->allocate_memory(sizeof(Index)));
            this->_indices_size.push_back(1);

            MemoryPool<Mem_>::set_memory(this->_elements.at(0), val);
            MemoryPool<Mem_>::set_memory(this->_indices.at(0), row);
            MemoryPool<Mem_>::set_memory(this->_indices.at(1), col);

            this->_scalar_index.at(3) = 1;
            return;
          }

          Index i(0);
          while (i < this->_scalar_index.at(3))
          {
            if (MemoryPool<Mem_>::get_element(this->_indices.at(0), i) >= row)
              break;
            ++i;
          }

          while (i < this->_scalar_index.at(3))
          {
            if (MemoryPool<Mem_>::get_element(this->_indices.at(1),i) >= col || MemoryPool<Mem_>::get_element(this->_indices.at(0), i) > row)
              break;
            ++i;
          }

          if(i < this->_scalar_index.at(3) && MemoryPool<Mem_>::get_element(this->_indices.at(0), i) == row && MemoryPool<Mem_>::get_element(this->_indices.at(1), i) == col)
          {
            MemoryPool<Mem_>::set_memory(this->_elements.at(0) + i, val);
            return;
          }
          else
          {
            DT_ * t_val((DT_*)MemoryPool<Mem_>::instance()->allocate_memory((this->_scalar_index.at(3) + 1) * sizeof(DT_)));
            Index * t_row((Index*)MemoryPool<Mem_>::instance()->allocate_memory((this->_scalar_index.at(3) + 1) * sizeof(Index)));
            Index * t_col((Index*)MemoryPool<Mem_>::instance()->allocate_memory((this->_scalar_index.at(3) + 1) * sizeof(Index)));

            if (i > 0)
            {
              MemoryPool<Mem_>::copy(t_val, this->_elements.at(0), (i) * sizeof(DT_));
              MemoryPool<Mem_>::copy(t_row, this->_indices.at(0), (i) * sizeof(Index));
              MemoryPool<Mem_>::copy(t_col, this->_indices.at(1), (i) * sizeof(Index));
            }

            MemoryPool<Mem_>::set_memory(t_val + i, val);
            MemoryPool<Mem_>::set_memory(t_row + i, row);
            MemoryPool<Mem_>::set_memory(t_col + i, col);

            MemoryPool<Mem_>::copy(t_val + i + 1, this->_elements.at(0) + i, (this->_scalar_index.at(3) - i) * sizeof(DT_));
            MemoryPool<Mem_>::copy(t_row + i + 1, this->_indices.at(0) + i, (this->_scalar_index.at(3) - i) * sizeof(Index));
            MemoryPool<Mem_>::copy(t_col+ i + 1, this->_indices.at(1) + i, (this->_scalar_index.at(3) - i) * sizeof(Index));

            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(0));
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(0));
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(1));

            this->_elements.at(0) = t_val;
            this->_indices.at(0) = t_row;
            this->_indices.at(1) = t_col;
            ++this->_scalar_index.at(3);

            this->_elements_size.at(0) = this->_scalar_index.at(3);
            this->_indices_size.at(0) = this->_scalar_index.at(3);
            this->_indices_size.at(1) = this->_scalar_index.at(3);
          }
        }

        /**
         * \brief Retrieve specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         *
         * \returns Specific matrix element.
         */
        const DT_ operator()(Index row, Index col) const
        {
          CONTEXT("When retrieving SparseMatrixCOO element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds sparse matrix coo row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds sparse matrix coo column size " + stringify(this->columns()) + " !");

          if (this->_elements.size() == 0)
            return this->_scalar_dt.at(0);

          Index i(0);
          while (i < this->_scalar_index.at(3))
          {
            if (MemoryPool<Mem_>::get_element(this->_indices.at(0), i) >= row)
              break;
            ++i;
          }

          while (i < this->_scalar_index.at(3))
          {
            if (MemoryPool<Mem_>::get_element(this->_indices.at(1),i) >= col || MemoryPool<Mem_>::get_element(this->_indices.at(0), i) > row)
              break;
            ++i;
          }
          if (i == this->_scalar_index.at(3))
            return this->_scalar_dt.at(0);

          if(MemoryPool<Mem_>::get_element(this->_indices.at(0), i) == row && MemoryPool<Mem_>::get_element(this->_indices.at(1), i) == col)
          {
            return MemoryPool<Mem_>::get_element(this->_elements.at(0), i);
          }
          else
            return this->_scalar_dt.at(0);
        }

        /**
         * \brief Reset all elements to zero.
         */
        void clear()
        {
          CONTEXT("When clearing SparseMatrixCOO");
          MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(0));
          MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(0));
          MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(1));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();
          this->_scalar_index.at(3) = 0;
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
        const Index & used_elements() const
        {
          return this->_scalar_index.at(3);
        }

        /**
         * \brief Retrieve non zero element array.
         *
         * \returns Non zero element array.
         */
        DT_ * val() const
        {
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve row coordinate array.
         *
         * \returns Row coordinate array.
         */
        Index * row() const
        {
          return this->_indices.at(0);
        }

        /**
         * \brief Retrieve columns coordinate array.
         *
         * \returns Column coordinate array.
         */
        Index * column() const
        {
          return this->_indices.at(1);
        }

        /**
         * \brief Retrieve non zero element.
         *
         * \returns Non zero element.
         */
        const DT_ zero_element() const
        {
          return this->_scalar_dt.at(0);
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String type_name()
        {
          return "SparseMatrixCOO";
        }
    };

    /**
     * \brief SparseMatrixCOO comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Mem_, typename Mem2_, typename DT_> bool operator== (const SparseMatrixCOO<Mem_, DT_> & a, const SparseMatrixCOO<Mem2_, DT_> & b)
    {
      CONTEXT("When comparing SparseMatrixCOOs");

      if (a.rows() != b.rows())
        return false;
      if (a.columns() != b.columns())
        return false;
      if (a.used_elements() != b.used_elements())
        return false;
      if (a.zero_element() != b.zero_element())
        return false;

      for (Index i(0) ; i < a.used_elements() ; ++i)
      {
        if (MemoryPool<Mem_>::get_element(a.val(), i) != MemoryPool<Mem2_>::get_element(b.val(), i))
            return false;
        if (MemoryPool<Mem_>::get_element(a.row(), i) != MemoryPool<Mem2_>::get_element(b.row(), i))
            return false;
        if (MemoryPool<Mem_>::get_element(a.column(), i) != MemoryPool<Mem2_>::get_element(b.column(), i))
            return false;
      }

      return true;
    }

    /**
     * \brief SparseMatrixCOO streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The matrix to be streamed.
     */
    template <typename Mem_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixCOO<Mem_, DT_> & b)
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

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
