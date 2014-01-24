#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_COO_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_COO_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/matrix_base.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/arch/scale.hpp>

#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdint.h>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Coordinate based sparse matrix.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements alongside with its coordinates explicitly. \n\n
     * Note, that the elements are sorted in a row major order.
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: row index \n
     * _indices[1]: column index \n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_index[4]: allocated elements \n
     * _scalar_index[5]: allocation size increment \n
     * _scalar_index[6]: boolean flag, if container is sorted \n
     * _scalar_dt[0]: zero element
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_>
    class SparseMatrixCOO : public Container<Mem_, DT_>, public MatrixBase
    {
      private:
        template <typename T1_, typename T2_, typename T3_>
        static void _insertion_sort(T1_ * key, T2_ * val1, T3_ * val2, Index size)
        {
          T1_ swap_key;
          T2_ swap1;
          T3_ swap2;
          for (Index i(1), j ; i < size ; ++i)
          {
            swap_key = MemoryPool<Mem_>::get_element(key, i);
            swap1 = MemoryPool<Mem_>::get_element(val1, i);
            swap2 = MemoryPool<Mem_>::get_element(val2, i);
            j = i;
            while (j > 0 && MemoryPool<Mem_>::get_element(key, j - 1) > swap_key)
            {
              MemoryPool<Mem_>::copy(key + j, key + j - 1, 1);
              MemoryPool<Mem_>::copy(val1 + j, val1 + j - 1, 1);
              MemoryPool<Mem_>::copy(val2 + j, val2 + j - 1, 1);
              --j;
            }
            MemoryPool<Mem_>::set_memory(key + j, swap_key);
            MemoryPool<Mem_>::set_memory(val1 + j, swap1);
            MemoryPool<Mem_>::set_memory(val2 + j, swap2);
          }
        }

        /*template <typename T1_, typename T2_, typename T3_>
        static void _comb_sort(T1_ * key, T2_ * val1, T3_ * val2, Index size)
        {
          const float shrink = 1.3f;
          T1_ swap_key;
          T2_ swap1;
          T3_ swap2;
          Index gap = size;
          bool swapped = false;
          while ((gap > 1) || swapped)
          {
            if (gap > 1)
            {
              gap = (Index)((float)gap / shrink);
            }

            swapped = false;

            for (Index i = 0 ; gap + i < size ; ++i)
            {
              if (MemoryPool<Mem_>::get_element(key, i) > MemoryPool<Mem_>::get_element(key, i + gap))
              {
                swap_key = MemoryPool<Mem_>::get_element(key, i);
                MemoryPool<Mem_>::copy(key + i, key + i + gap, 1);
                MemoryPool<Mem_>::set_memory(key + i + gap, swap_key);

                swap1 = MemoryPool<Mem_>::get_element(val1, i);
                MemoryPool<Mem_>::copy(val1 + i, val1 + i + gap, 1);
                MemoryPool<Mem_>::set_memory(val1 + i + gap, swap1);

                swap2 = MemoryPool<Mem_>::get_element(val2, i);
                MemoryPool<Mem_>::copy(val2 + i, val2 + i + gap, 1);
                MemoryPool<Mem_>::set_memory(val2 + i + gap, swap2);

                swapped = true;
              }
            }
          }
        }*/

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
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
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
            Index row((Index)atol(srow.c_str()));
            --row;
            this->_scalar_index.at(1) = std::max(row+1, this->rows());
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scol(line, 0, end);
            Index col((Index)atol(scol.c_str()));
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

          this->_scalar_index.at(5) = rowsv.size();
          for (Index i(0) ; i < rowsv.size() ; ++i)
          {
            (*this)(rowsv.at(i), colsv.at(i), valsv.at(i));
          }
          this->_scalar_index.at(5) = 1000;
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
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));

          String line;
          std::getline(file, line);
          {
            std::getline(file, line);

            String::size_type begin(line.find_first_not_of(" "));
            line.erase(0, begin);
            String::size_type end(line.find_first_of(" "));
            String srow(line, 0, end);
            Index row((Index)atol(srow.c_str()));
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scol(line, 0, end);
            Index col((Index)atol(scol.c_str()));
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
            Index row((Index)atol(srow.c_str()));
            --row;
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scol(line, 0, end);
            Index col((Index)atol(scol.c_str()));
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

          this->_scalar_index.at(5) = rowsv.size();
          for (Index i(0) ; i < rowsv.size() ; ++i)
          {
            (*this)(rowsv.at(i), colsv.at(i), valsv.at(i));
          }
          this->_scalar_index.at(5) = 1000;
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
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
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
          this->_scalar_index.at(4) = (Index)elements;

          uint64_t * crow_ptr = new uint64_t[elements];
          file.read((char *)crow_ptr, (long)((elements) * sizeof(uint64_t)));
          Index * trow_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(elements);
          for (Index i(0) ; i < elements ; ++i)
            trow_ptr[i] = Index(crow_ptr[i]);
          delete[] crow_ptr;

          uint64_t * ccol_ptr = new uint64_t[elements];
          file.read((char *)ccol_ptr, (long)((elements) * sizeof(uint64_t)));
          Index * tcol_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(elements);
          for (Index i(0) ; i < elements ; ++i)
            tcol_ptr[i] = Index(ccol_ptr[i]);
          delete[] ccol_ptr;

          double * cval = new double[elements];
          file.read((char *)cval, (long)(elements * sizeof(double)));
           DT_ * tval_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(elements);
          for (Index i(0) ; i < elements ; ++i)
            tval_ptr[i] = DT_(cval[i]);
          delete[] cval;

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(elements));
          this->_elements_size.push_back(elements);
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(elements));
          this->_indices_size.push_back(elements);
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(elements));
          this->_indices_size.push_back(elements);

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval_ptr, elements);
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(0), trow_ptr, elements);
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(1), tcol_ptr, elements);
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
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
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
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
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
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
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
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));

          SparseMatrixELL<Mem::Main, DT_> cother(other);
          DT_ * tval = new DT_[other.used_elements()];
          Index * trow = new Index[other.used_elements()];
          Index * tcolumn = new Index[other.used_elements()];

          const Index stride(cother.stride());
          for (Index row(0), ue(0); row < cother.rows() ; ++row)
          {
            const Index * tAj(cother.Aj());
            const DT_ * tAx(cother.Ax());
            tAj += row;
            tAx += row;

            const Index max(cother.Arl()[row]);
            for(Index n(0); n < max ; n++)
            {
              tval[ue] = *tAx;
              trow[ue] = row;
              tcolumn[ue] = *tAj;

              tAj += stride;
              tAx += stride;
              ++ue;
            }
          }

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));

          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(0), trow, this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(1), tcolumn, this->_scalar_index.at(3));

          delete[] tval;
          delete[] trow;
          delete[] tcolumn;
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
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));

          SparseMatrixCSR<Mem::Main, DT_> cother(other);
          DT_ * tval = new DT_[other.used_elements()];
          Index * trow = new Index[other.used_elements()];
          Index * tcolumn = new Index[other.used_elements()];

          for (Index row(0), ue(0) ; row < other.rows() ; ++row)
          {
            const Index end(cother.row_ptr()[row + 1]);
            for (Index i(cother.row_ptr()[row]) ; i < end ; ++i)
            {
              tval[ue] = cother.val()[i];
              trow[ue] = row;
              tcolumn[ue] = cother.col_ind()[i];
              ++ue;
            }
          }

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));

          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(0), trow, this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->_indices.at(1), tcolumn, this->_scalar_index.at(3));

          delete[] tval;
          delete[] trow;
          delete[] tcolumn;
        }

        /**
         * \brief Constructor
         *
         * \param[in] rows The row count of the created matrix.
         * \param[in] columns The column count of the created matrix.
         * \param[in] row_ind Vector with row indices.
         * \param[in] col_ind Vector with column indices.
         * \param[in] val Vector with non zero elements.
         *
         * Creates a matrix with given dimensions and content.
         */
        explicit SparseMatrixCOO(const Index rows, const Index columns, DenseVector<Mem_, Index> & row_ind,
            DenseVector<Mem_, Index> & col_ind, DenseVector<Mem_, DT_> & val) :
          Container<Mem_, DT_>(rows * columns)
        {
          CONTEXT("When creating SparseMatrixCOO");
          this->_scalar_index.push_back(rows);
          this->_scalar_index.push_back(columns);
          this->_scalar_index.push_back(val.size());
          this->_scalar_index.push_back(val.size());
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(0);
          this->_scalar_dt.push_back(DT_(0));

          this->_elements.push_back(val.elements());
          this->_elements_size.push_back(val.size());
          this->_indices.push_back(row_ind.elements());
          this->_indices_size.push_back(row_ind.size());
          this->_indices.push_back(col_ind.elements());
          this->_indices_size.push_back(col_ind.size());

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
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
          this->_scalar_index.push_back(other.allocated_elements());
          this->_scalar_index.push_back(other.alloc_increment());
          this->_scalar_index.push_back(other.sorted());
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
          this->_scalar_index.push_back(other.allocated_elements());
          this->_scalar_index.push_back(other.alloc_increment());
          this->_scalar_index.push_back(other.sorted());
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
        template <typename Mem2_, typename DT2_>
        SparseMatrixCOO<Mem_, DT_> & operator= (const SparseMatrixCOO<Mem2_, DT2_> & other)
        {
          CONTEXT("When assigning SparseMatrixCOO");

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
          CONTEXT("When writing out SparseMatrixCOO");

          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();

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

          const Index ue(used_elements());
          Index * row_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(ue);
          MemoryPool<Mem_>::template download<Index>(row_ptr, this->_indices.at(0), ue);
          uint64_t * crow_ptr = new uint64_t[ue];
          for (Index i(0) ; i < ue ; ++i)
            crow_ptr[i] = row_ptr[i];
          MemoryPool<Mem::Main>::instance()->release_memory(row_ptr);

          Index * col_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(ue);
          MemoryPool<Mem_>::template download<Index>(col_ptr, this->_indices.at(1), ue);
          uint64_t * ccol_ptr = new uint64_t[ue];
          for (Index i(0) ; i < ue ; ++i)
            ccol_ptr[i] = col_ptr[i];
          MemoryPool<Mem::Main>::instance()->release_memory(col_ptr);

          DT_ * val = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(ue);
          MemoryPool<Mem_>::template download<DT_>(val, this->_elements.at(0), ue);
          double * cval = new double[ue];
          for (Index i(0) ; i < ue ; ++i)
            cval[i] = Type::Traits<DT_>::to_double(val[i]);
          MemoryPool<Mem::Main>::instance()->release_memory(val);

          uint64_t rows(this->_scalar_index.at(1));
          uint64_t columns(this->_scalar_index.at(2));
          uint64_t elements(ue);
          file.write((const char *)&rows, sizeof(uint64_t));
          file.write((const char *)&columns, sizeof(uint64_t));
          file.write((const char *)&elements, sizeof(uint64_t));
          file.write((const char *)crow_ptr, (long)(elements * sizeof(uint64_t)));
          file.write((const char *)ccol_ptr, (long)(elements * sizeof(uint64_t)));
          file.write((const char *)cval, (long)(elements * sizeof(double)));

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
          for (Index i(0) ; i < used_elements() ; ++i)
          {
            file << temp.row()[i] + 1 << " " << temp.column()[i] + 1 << " " << std::scientific << Type::Traits<DT_>::to_double(temp.val()[i]) << ";" << std::endl;
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

          for (Index i(0) ; i < used_elements() ; ++i)
          {
            file << temp.row()[i] + 1 << " " << temp.column()[i] + 1 << " " << std::scientific << Type::Traits<DT_>::to_double(temp.val()[i]) << std::endl;
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

          this->_scalar_index.at(6) = 0;

          if (this->_elements.size() == 0)
          {
            this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(alloc_increment()));
            this->_elements_size.push_back(alloc_increment());
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(alloc_increment()));
            this->_indices_size.push_back(alloc_increment());
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(alloc_increment()));
            this->_indices_size.push_back(alloc_increment());

            MemoryPool<Mem_>::set_memory(this->_elements.at(0), val);
            MemoryPool<Mem_>::set_memory(this->_indices.at(0), row);
            MemoryPool<Mem_>::set_memory(this->_indices.at(1), col);

            this->_scalar_index.at(3) = 1;
            this->_scalar_index.at(4) = alloc_increment();
            return;
          }

          //append element
          if (this->_scalar_index.at(3) < this->_scalar_index.at(4))
          {
            MemoryPool<Mem_>::set_memory(this->_elements.at(0) + this->_scalar_index.at(3), val);
            MemoryPool<Mem_>::set_memory(this->_indices.at(0) + this->_scalar_index.at(3), row);
            MemoryPool<Mem_>::set_memory(this->_indices.at(1) + this->_scalar_index.at(3), col);

            ++this->_scalar_index.at(3);
          }
          //reallocate
          else
          {
            this->_scalar_index.at(4) += alloc_increment();

            DT_ * elements_new(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(allocated_elements()));
            Index * rows_new(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(allocated_elements()));
            Index * cols_new(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(allocated_elements()));

            MemoryPool<Mem_>::copy(elements_new, this->_elements.at(0), used_elements());
            MemoryPool<Mem_>::copy(rows_new, this->_indices.at(0), used_elements());
            MemoryPool<Mem_>::copy(cols_new, this->_indices.at(1), used_elements());

            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(0));
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(0));
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(1));

            this->_elements.at(0) = elements_new;
            this->_indices.at(0) = rows_new;
            this->_indices.at(1) = cols_new;

            MemoryPool<Mem_>::set_memory(this->_elements.at(0) + used_elements(), val);
            MemoryPool<Mem_>::set_memory(this->_indices.at(0) + used_elements(), row);
            MemoryPool<Mem_>::set_memory(this->_indices.at(1) + used_elements(), col);

            ++this->_scalar_index.at(3);
            this->_elements_size.at(0) = allocated_elements();
            this->_indices_size.at(0) = allocated_elements();
            this->_indices_size.at(1) = allocated_elements();
          }
        }

        void sort()
        {
          if (this->_scalar_index.at(6) == 0)
          {
            //first of all, mark matrix as sorted, because otherwise we would call ourselves inifite times
            this->_scalar_index.at(6) = 1;

            // sort elements by row index
            _insertion_sort(row(), val(), column(), used_elements());

            Index row_start(0);
            for (Index rowi(0) ; rowi < rows() ; ++rowi)
            {
              Index offset(0);
              // explore range of elements in given row
              for ( ; row_start + offset < used_elements() && MemoryPool<Mem_>::get_element(row(), row_start + offset) == rowi ; ++offset) ;
              // sort range of elements in given row by column index
              _insertion_sort(column() + row_start, val() + row_start, row() + row_start, offset);
              // find and mark duplicate entries
              for (Index i(row_start + 1) ; i < row_start + offset ; ++i)
              {
                if (MemoryPool<Mem_>::get_element(column(), i - 1) == MemoryPool<Mem_>::get_element(column(), i))
                {
                  MemoryPool<Mem_>::set_memory(row() + i - 1, std::numeric_limits<Index>::max());
                }
              }
              row_start += offset;
            }

            // sort out marked duplicated elements
            _insertion_sort(row(), val(), column(), used_elements());
            Index junk(0);
            while (MemoryPool<Mem_>::get_element(row(), used_elements() - 1 - junk) == std::numeric_limits<Index>::max()
                && junk < used_elements())
              ++junk;
            this->_scalar_index.at(3) -= junk;
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
        DT_ operator()(Index row, Index col) const
        {
          CONTEXT("When retrieving SparseMatrixCOO element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds sparse matrix coo row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds sparse matrix coo column size " + stringify(this->columns()) + " !");

          if (this->_elements.size() == 0)
            return this->_scalar_dt.at(0);

          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();

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
          this->_scalar_index.at(4) = 0;
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
        Index used_elements() const
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();
          return this->_scalar_index.at(3);
        }

        /**
         * \brief Retrieve non zero element array.
         *
         * \returns Non zero element array.
         */
        DT_ * val()
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();
          return this->_elements.at(0);
        }

        DT_ const * val() const
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve row coordinate array.
         *
         * \returns Row coordinate array.
         */
        Index * row()
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();
          return this->_indices.at(0);
        }

        Index const * row() const
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();
          return this->_indices.at(0);
        }

        /**
         * \brief Retrieve columns coordinate array.
         *
         * \returns Column coordinate array.
         */
        Index * column()
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();
          return this->_indices.at(1);
        }

        Index const * column() const
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO<Mem_, DT_> *>(this)->sort();
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
         * \brief Retrieve amount of allocated elements.
         *
         * \return Allocated element count.
         */
        Index allocated_elements() const
        {
          return this->_scalar_index.at(4);
        }

        /**
         * \brief Retrieve allocation incrementation value.
         *
         * \return Allocation increment.
         */
        Index alloc_increment() const
        {
          return this->_scalar_index.at(5);
        }

        /**
         * \brief Retrieve status of element row-wise sorting.
         *
         * \return Sorting status.
         */
        Index sorted() const
        {
          return this->_scalar_index.at(6);
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

        /**
         * \brief Calculate \f$this \leftarrow x \cdot s\f$
         *
         * \param[in] x The matrix to be scaled.
         * \param[in] s A scalar to scale x with.
         */
        template <typename Algo_>
        void scale(const SparseMatrixCOO<Mem_, DT_> & x, const DT_ s)
        {
          if (x.rows() != this->rows())
            throw InternalError("Row count does not match!");
          if (x.columns() != this->columns())
            throw InternalError("Column count does not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError("Nonzero count does not match!");

          Arch::Scale<Mem_, Algo_>::value(this->val(), x.val(), s, this->used_elements());
        }

        template <typename Algo_>
        void product(const SparseMatrixCOO<Mem_, DT_> & x, const DT_ s)
        {
          this->template scale<Algo_>(x, s);
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

#endif // KERNEL_LAFEM_SPARSE_MATRIX_COO_HPP
