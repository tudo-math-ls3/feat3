#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_COO_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_COO_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/matrix_base.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>

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
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseMatrixCOO : public Container<Mem_, DT_, IT_>, public MatrixBase
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
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
          _read_from_m(file);
          file.close();
        }

        void _read_from_m(std::istream& file)
        {
          std::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));

          Index ue(0);
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

            entries[row].insert(std::pair<IT_, DT_>(col, val));
            ++ue;
          }
          this->_scalar_index.at(0) = this->rows() * this->columns();
          this->_scalar_index.at(3) = ue;
          this->_scalar_index.at(4) = ue;

          DT_ * tval = new DT_[ue];
          IT_ * trow = new IT_[ue];
          IT_ * tcolumn = new IT_[ue];

          Index idx(0);
          for (auto row : entries)
          {
            for (auto col : row.second )
            {
              trow[idx] = row.first;
              tcolumn[idx] = col.first;
              tval[idx] = col.second;
              ++idx;
            }
            row.second.clear();
          }
          entries.clear();

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));

          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), trow, this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), tcolumn, this->_scalar_index.at(3));

          delete[] tval;
          delete[] trow;
          delete[] tcolumn;
        }

        void _read_from_mtx(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
          _read_from_mtx(file);
          file.close();
        }

        void _read_from_mtx(std::istream& file)
        {
          std::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));

          Index ue(0);
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

            entries[row].insert(std::pair<IT_, DT_>(col, val));
            ++ue;
          }

          this->_scalar_index.at(0) = this->rows() * this->columns();
          this->_scalar_index.at(3) = ue;
          this->_scalar_index.at(4) = ue;

          DT_ * tval = new DT_[ue];
          IT_ * trow = new IT_[ue];
          IT_ * tcolumn = new IT_[ue];

          Index idx(0);
          for (auto row : entries)
          {
            for (auto col : row.second )
            {
              trow[idx] = row.first;
              tcolumn[idx] = col.first;
              tval[idx] = col.second;
              ++idx;
            }
            row.second.clear();
          }
          entries.clear();

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));

          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), trow, this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), tcolumn, this->_scalar_index.at(3));

          delete[] tval;
          delete[] trow;
          delete[] tcolumn;
        }

        void _read_from_coo(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
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
          IT_ * trow_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(elements);
          for (Index i(0) ; i < elements ; ++i)
            trow_ptr[i] = IT_(crow_ptr[i]);
          delete[] crow_ptr;

          uint64_t * ccol_ptr = new uint64_t[elements];
          file.read((char *)ccol_ptr, (long)((elements) * sizeof(uint64_t)));
          IT_ * tcol_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(elements);
          for (Index i(0) ; i < elements ; ++i)
            tcol_ptr[i] = IT_(ccol_ptr[i]);
          delete[] ccol_ptr;

          double * cval = new double[elements];
          file.read((char *)cval, (long)(elements * sizeof(double)));
           DT_ * tval_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(elements);
          for (Index i(0) ; i < elements ; ++i)
            tval_ptr[i] = DT_(cval[i]);
          delete[] cval;

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(elements));
          this->_elements_size.push_back(elements);
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(elements));
          this->_indices_size.push_back(elements);
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(elements));
          this->_indices_size.push_back(elements);

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tval_ptr, elements);
          MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(0), trow_ptr, elements);
          MemoryPool<Mem_>::template upload<IT_>(this->get_indices().at(1), tcol_ptr, elements);
          MemoryPool<Mem::Main>::instance()->release_memory(tval_ptr);
          MemoryPool<Mem::Main>::instance()->release_memory(trow_ptr);
          MemoryPool<Mem::Main>::instance()->release_memory(tcol_ptr);
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our indextype
        typedef IT_ IndexType;
        /// Our memory architecture type
        typedef Mem_ MemType;
        /// Compatible L-vector type
        typedef DenseVector<MemType, DT_, IT_> VectorTypeL;
        /// Compatible R-vector type
        typedef DenseVector<MemType, DataType, IT_> VectorTypeR;
        /// Our used layout type
        const static SparseLayoutType LayoutType = SparseLayoutType::lt_coo;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional matrix.
         */
        explicit SparseMatrixCOO() :
          Container<Mem_, DT_, IT_> (0)
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
          Container<Mem_, DT_, IT_>(dimensions * dimensions)
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
          Container<Mem_, DT_, IT_>(rows * columns)
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
         * \param[in] layout The layout to be used.
         *
         * Creates an empty matrix with given layout.
         */
        explicit SparseMatrixCOO(const SparseLayout<Mem_, LayoutType> & layout) :
          Container<Mem_, DT_> (layout._scalar_index.at(0))
        {
          CONTEXT("When creating SparseMatrixCOO");
          this->_indices.assign(layout._indices.begin(), layout._indices.end());
          this->_indices_size.assign(layout._indices_size.begin(), layout._indices_size.end());
          this->_scalar_index.assign(layout._scalar_index.begin(), layout._scalar_index.end());
          this->_scalar_dt.push_back(DT_(0));

          for (auto i : this->_indices)
            MemoryPool<Mem_>::instance()->increase_memory(i);

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4)));
          this->_elements_size.push_back(this->_scalar_index.at(4));
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source ell matrix.
         *
         * Creates a matrix from the given source matrix.
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        explicit SparseMatrixCOO(const SparseMatrixELL<Mem2_, DT2_, IT2_> & other) :
          Container<Mem_, DT_, IT_>(other.rows() * other.columns())
        {
          CONTEXT("When creating SparseMatrixCOO from SparseMatrixELL");

          convert(other);
        }

        /**
         * \brief Constructor
         *
         * \param[in] other The source csr matrix.
         *
         * Creates a matrix from the given source matrix.
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        explicit SparseMatrixCOO(const SparseMatrixCSR<Mem2_, DT2_, IT2_> & other) :
          Container<Mem_, DT_, IT_>(other.rows() * other.columns())
        {
          CONTEXT("When creating SparseMatrixCOO from SparseMatrixCSR");

          convert(other);
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
        explicit SparseMatrixCOO(const Index rows, const Index columns, DenseVector<Mem_, IT_> & row_ind,
            DenseVector<Mem_, IT_> & col_ind, DenseVector<Mem_, DT_, IT_> & val) :
          Container<Mem_, DT_, IT_>(rows * columns)
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
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseMatrixCOO");

          switch(mode)
          {
            case FileMode::fm_m:
              _read_from_m(filename);
              break;
            case FileMode::fm_mtx:
              _read_from_mtx(filename);
              break;
            case FileMode::fm_coo:
              _read_from_coo(filename);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
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
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseMatrixCOO");

          switch(mode)
          {
            case FileMode::fm_m:
              _read_from_m(file);
              break;
            case FileMode::fm_mtx:
              _read_from_mtx(file);
              break;
            case FileMode::fm_coo:
              _read_from_coo(file);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source matrix.
         *
         * Moves a given matrix to this matrix.
         */
        SparseMatrixCOO(SparseMatrixCOO && other) :
          Container<Mem_, DT_, IT_>(std::move(other))
        {
          CONTEXT("When moving SparseMatrixCOO");
        }

        /**
         * \brief Move operator
         *
         * \param[in] other The source matrix.
         *
         * Moves another matrix to the target matrix.
         */
        SparseMatrixCOO & operator= (SparseMatrixCOO && other)
        {
          CONTEXT("When moving SparseMatrixCOO");

          this->move(std::move(other));

          return *this;
        }

        /**
         * \brief Convertion method
         *
         * Use source matrix content as content of current matrix
         *
         * \param[in] other The source Matrix.
         *
         * \note This creates a deep copy in any case!
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const SparseMatrixCOO<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting SparseMatrixCOO");
          this->clone(other);
        }

        /**
         * \brief Convertion method
         *
         * \param[in] other The source Matrix.
         *
         * Use source matrix content as content of current matrix
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const SparseMatrixCSR<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting SparseMatrixCOO");

          this->clear();

          this->_scalar_index.push_back(other.size());
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));

          SparseMatrixCSR<Mem::Main, DT_, IT_> cother;
          cother.convert(other);

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));

          DT_ * tval(nullptr);
          IT_ * trow(nullptr);
          IT_ * tcolumn(nullptr);
          if (std::is_same<Mem_, Mem::Main>::value)
          {
            tval = this->_elements.at(0);
            trow = this->_indices.at(0);
            tcolumn = this->_indices.at(1);
          }
          else
          {
            tval = new DT_[other.used_elements()];
            trow = new IT_[other.used_elements()];
            tcolumn = new IT_[other.used_elements()];
          }

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

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, this->_scalar_index.at(3));
            MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), trow, this->_scalar_index.at(3));
            MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), tcolumn, this->_scalar_index.at(3));
            delete[] tval;
            delete[] trow;
            delete[] tcolumn;
          }
        }

        /**
         * \brief Convertion method
         *
         * \param[in] other The source Matrix.
         *
         * Use source matrix content as content of current matrix
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const SparseMatrixELL<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting SparseMatrixCOO");

          this->clear();

          this->_scalar_index.push_back(other.size());
          this->_scalar_index.push_back(other.rows());
          this->_scalar_index.push_back(other.columns());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(other.used_elements());
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));

          SparseMatrixELL<Mem::Main, DT_, IT_> cother;
          cother.convert(other);

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(3));

          DT_ * tval(nullptr);
          IT_ * trow(nullptr);
          IT_ * tcolumn(nullptr);
          if (std::is_same<Mem_, Mem::Main>::value)
          {
            tval = this->_elements.at(0);
            trow = this->_indices.at(0);
            tcolumn = this->_indices.at(1);
          }
          else
          {
            tval = new DT_[other.used_elements()];
            trow = new IT_[other.used_elements()];
            tcolumn = new IT_[other.used_elements()];
          }

          const Index stride(cother.stride());
          for (Index row(0), ue(0); row < cother.rows() ; ++row)
          {
            const IT_ * tAj(cother.Aj());
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

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, this->_scalar_index.at(3));
            MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), trow, this->_scalar_index.at(3));
            MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), tcolumn, this->_scalar_index.at(3));
            delete[] tval;
            delete[] trow;
            delete[] tcolumn;
          }
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
            const_cast<SparseMatrixCOO *>(this)->sort();

          switch(mode)
          {
            case FileMode::fm_coo:
              write_out_coo(filename);
              break;
            case FileMode::fm_m:
              write_out_m(filename);
              break;
            case FileMode::fm_mtx:
              write_out_mtx(filename);
              break;
            default:
                throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
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
            case FileMode::fm_coo:
              write_out_coo(file);
              break;
            case FileMode::fm_m:
              write_out_m(file);
              break;
            case FileMode::fm_mtx:
              write_out_mtx(file);
              break;
            default:
                throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
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
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
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
          if (! std::is_same<DT_, double>::value)
            std::cout<<"Warning: You are writing out an coo matrix with less than double precission!"<<std::endl;

          const Index ue(used_elements());
          IT_ * row_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(ue);
          MemoryPool<Mem_>::template download<IT_>(row_ptr, this->_indices.at(0), ue);
          uint64_t * crow_ptr = new uint64_t[ue];
          for (Index i(0) ; i < ue ; ++i)
            crow_ptr[i] = row_ptr[i];
          MemoryPool<Mem::Main>::instance()->release_memory(row_ptr);

          IT_ * col_ptr = MemoryPool<Mem::Main>::instance()->template allocate_memory<IT_>(ue);
          MemoryPool<Mem_>::template download<IT_>(col_ptr, this->_indices.at(1), ue);
          uint64_t * ccol_ptr = new uint64_t[ue];
          for (Index i(0) ; i < ue ; ++i)
            ccol_ptr[i] = col_ptr[i];
          MemoryPool<Mem::Main>::instance()->release_memory(col_ptr);

          DT_ * val = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(ue);
          MemoryPool<Mem_>::template download<DT_>(val, this->_elements.at(0), ue);
          double * cval = new double[ue];
          for (Index i(0) ; i < ue ; ++i)
            cval[i] = (double)val[i];
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
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
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
          SparseMatrixCOO<Mem::Main, DT_, IT_> temp;
          temp.convert(*this);

          file << "data = [" << std::endl;
          for (Index i(0) ; i < used_elements() ; ++i)
          {
            file << temp.row()[i] + 1 << " " << temp.column()[i] + 1 << " " << std::scientific << temp.val()[i] << ";" << std::endl;
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
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
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
          SparseMatrixCOO<Mem::Main, DT_, IT_> temp;
          temp.convert(*this);

          file << "%%MatrixMarket matrix coordinate real general" << std::endl;
          file << temp.rows() << " " << temp.columns() << " " << temp.used_elements() << std::endl;

          for (Index i(0) ; i < used_elements() ; ++i)
          {
            file << temp.row()[i] + 1 << " " << temp.column()[i] + 1 << " " << std::scientific << temp.val()[i] << std::endl;
          }
        }

        /**
         * \brief Set specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         * \param[in] value The value to be set.
         */
        void operator()(IT_ row, IT_ col, DT_ val)
        {
          CONTEXT("When setting SparseMatrixCOO element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds sparse matrix coo row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds sparse matrix coo column size " + stringify(this->columns()) + " !");

          this->_scalar_index.at(6) = 0;

          if (this->_elements.size() == 0)
          {
            this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(alloc_increment()));
            this->_elements_size.push_back(alloc_increment());
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(alloc_increment()));
            this->_indices_size.push_back(alloc_increment());
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(alloc_increment()));
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
            IT_ * rows_new(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(allocated_elements()));
            IT_ * cols_new(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(allocated_elements()));

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

            // check if there is anything to be sorted
            if(used_elements() <= Index(0))
              return;

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
                  MemoryPool<Mem_>::set_memory(row() + i - 1, std::numeric_limits<IT_>::max());
                }
              }
              row_start += offset;
            }

            // sort out marked duplicated elements
            _insertion_sort(row(), val(), column(), used_elements());
            Index junk(0);
            while (MemoryPool<Mem_>::get_element(row(), used_elements() - 1 - junk) == std::numeric_limits<IT_>::max()
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
            const_cast<SparseMatrixCOO *>(this)->sort();

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
         * \brief Reset all elements of the container to a given value or zero if missing.
         *
         * \param[in] value The value to be set (defaults to 0)
         *
         */
        void format(DT_ value = DT_(0))
        {
          CONTEXT("When clearing SparseMatrixCOO");

          if (value == DT_(0))
          {
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
          else
          {
            ((Container<Mem_, DT_, IT_> &)*this).format(value);
          }
        }

        /**
         * \brief Retrieve convenient sparse matrix layout object.
         *
         * \return An object containing the sparse matrix layout.
         *
         * \note This methods creates a deep copy of its own layout and returns it.
         * This is necessary because coo layouts may change after creation and thus cannot be used by two different SparseMatrix Objects.
         * Nevertheless it is usefull to extract a matrix' layout, to create another matrix with the same matrix (same as 'same' at the moment of creation).
         */
        SparseLayout<Mem_, LayoutType> layout() const
        {
          SparseMatrixCOO t;
          t = this->clone();
          SparseLayout<Mem_, LayoutType> layout(t._indices, t._indices_size, t._scalar_index);
          return layout;
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
        const Index & used_elements() const override
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
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
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_elements.at(0);
        }

        DT_ const * val() const
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve row coordinate array.
         *
         * \returns Row coordinate array.
         */
        IT_ * row()
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_indices.at(0);
        }

        IT_ const * row() const
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_indices.at(0);
        }

        /**
         * \brief Retrieve columns coordinate array.
         *
         * \returns Column coordinate array.
         */
        IT_ * column()
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_indices.at(1);
        }

        IT_ const * column() const
        {
          if (this->_scalar_index.at(6) == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
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
        const Index & allocated_elements() const
        {
          return this->_scalar_index.at(4);
        }

        /**
         * \brief Retrieve allocation incrementation value.
         *
         * \return Allocation increment.
         */
        const Index & alloc_increment() const
        {
          return this->_scalar_index.at(5);
        }

        /**
         * \brief Retrieve status of element row-wise sorting.
         *
         * \return Sorting status.
         */
        const Index & sorted() const
        {
          return this->_scalar_index.at(6);
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String name()
        {
          return "SparseMatrixCOO";
        }

        /**
         * \brief Performs \f$this \leftarrow x\f$.
         *
         * \param[in] x The Matrix to be copied.
         */
        void copy(const SparseMatrixCOO & x)
        {
          this->sort();
          const_cast<SparseMatrixCOO *>(&x)->sort();
          this->_copy_content(x);
        }

        /**
         * \brief Performs \f$this \leftarrow x\f$.
         *
         * \param[in] x The Matrix to be copied.
         */
        template <typename Mem2_>
        void copy(const SparseMatrixCOO<Mem2_, DT_, IT_> & x)
        {
          this->sort();
          const_cast<SparseMatrixCOO<Mem2_, DT_, IT_> *>(&x)->sort();
          this->_copy_content(x);
        }

        /**
         * \brief Calculate \f$this \leftarrow y + \alpha x\f$
         *
         * \param[in] x The first summand matrix to be scaled.
         * \param[in] y The second summand matrix
         * \param[in] alpha A scalar to multiply x with.
         */
        template <typename Algo_>
        void axpy(
          const SparseMatrixCOO & x,
          const SparseMatrixCOO & y,
          const DT_ alpha = DT_(1))
        {
          if (x.rows() != y.rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
          if (x.rows() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix rows do not match!");
          if (x.columns() != y.columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
          if (x.columns() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix columns do not match!");
          if (x.used_elements() != y.used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix used_elements do not match!");

          // check for special cases
          // r <- x + y
          if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
            Arch::Sum<Mem_, Algo_>::value(this->val(), x.val(), y.val(), this->used_elements());
          // r <- y - x
          else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
            Arch::Difference<Mem_, Algo_>::value(this->val(), y.val(), x.val(), this->used_elements());
          // r <- y
          else if (Math::abs(alpha) < Math::eps<DT_>())
            this->copy(y);
          // r <- y + alpha*x
          else
            Arch::Axpy<Mem_, Algo_>::dv(this->val(), alpha, x.val(), y.val(), this->used_elements());
        }

        /**
         * \brief Calculate \f$this \leftarrow \alpha x \f$
         *
         * \param[in] x The matrix to be scaled.
         * \param[in] alpha A scalar to scale x with.
         */
        template <typename Algo_>
        void scale(const SparseMatrixCOO & x, const DT_ alpha)
        {
          if (x.rows() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
          if (x.columns() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");

          Arch::Scale<Mem_, Algo_>::value(this->val(), x.val(), alpha, this->used_elements());
        }

        /**
         * \brief Calculate \f$r \leftarrow this\cdot x \f$
         *
         * \param[out] r The vector that recieves the result.
         * \param[in] x The vector to be multiplied by this matrix.
         */
        template<typename Algo_>
        void apply(DenseVector<Mem_,DT_, IT_>& r, const DenseVector<Mem_, DT_, IT_>& x) const
        {
          if (r.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
          if (x.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

          Arch::ProductMatVec<Mem_, Algo_>::coo(r.elements(), this->val(), this->row(),
            this->column(), x.elements(), this->rows(), this->used_elements());
        }

        /**
         * \brief Calculate \f$r \leftarrow y + \alpha this\cdot x \f$
         *
         * \param[out] r The vector that recieves the result.
         * \param[in] x The vector to be multiplied by this matrix.
         * \param[in] y The summand vector.
         * \param[in] alpha A scalar to scale the product with.
         */
        template<typename Algo_>
        void apply(
          DenseVector<Mem_,DT_, IT_>& r,
          const DenseVector<Mem_, DT_, IT_>& x,
          const DenseVector<Mem_, DT_, IT_>& y,
          const DT_ alpha = DT_(1)) const
        {
          if (r.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
          if (x.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
          if (y.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

          // check for special cases
          // r <- y - A*x
          if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          {
            Arch::Defect<Mem_, Algo_>::coo(r.elements(), y.elements(), this->val(),
                this->row(), this->column(), x.elements(), this->rows(), this->used_elements());
          }
          // r <- y
          else if(Math::abs(alpha) < Math::eps<DT_>())
            r.copy(y);
          // r <- y + alpha*x
          else
          {
            Arch::Axpy<Mem_, Algo_>::coo(r.elements(), alpha, x.elements(), y.elements(),
              this->val(), this->row(), this->column(), this->rows(), this->used_elements());
          }
        }

        /// Returns a new compatible L-Vector.
        VectorTypeL create_vector_l() const
        {
          return VectorTypeL(this->rows());
        }

        /// Returns a new compatible R-Vector.
        VectorTypeR create_vector_r() const
        {
          return VectorTypeR(this->columns());
        }
    };

    /**
     * \brief SparseMatrixCOO comparison operator
     *
     * \param[in] a A matrix to compare with.
     * \param[in] b A matrix to compare with.
     */
    template <typename Mem_, typename Mem2_, typename DT_, typename IT_> bool operator== (const SparseMatrixCOO<Mem_, DT_, IT_> & a, const SparseMatrixCOO<Mem2_, DT_, IT_> & b)
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
    template <typename Mem_, typename DT_, typename IT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixCOO<Mem_, DT_, IT_> & b)
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
