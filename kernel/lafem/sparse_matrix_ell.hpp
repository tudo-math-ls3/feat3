#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_ELL_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_ELL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>

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
          std::map<Index, std::map<Index, DT_> > entries; // map<row, map<column, value> >
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
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

            entries[row].insert(std::pair<Index, DT_>(col, val));
            ++ue;
          }
          this->_scalar_index.at(0) = this->rows() * this->columns();
          this->_scalar_index.at(5) = ue;
          this->_scalar_index.at(4) = 0;

          Index* tArl = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1));

          Index idx(0);
          for (auto row : entries)
          {
            tArl[idx] = Index(row.second.size());
            this->_scalar_index.at(4) = std::max(this->_scalar_index.at(4), tArl[idx]);
            ++idx;
          }

          Index alignment(32);
          this->_scalar_index.at(3) = alignment * ((this->_scalar_index.at(1) + alignment - 1)/ alignment);

          DT_* tAx = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem::Main>::instance()->set_memory(tAx, DT_(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));
          Index* tAj = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem::Main>::instance()->set_memory(tAj, Index(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));

          Index row_idx(0);
          for (auto row : entries)
          {
            Index target(0);
            for (auto col : row.second )
            {
              tAj[row_idx + target * stride()] = col.first;
              tAx[row_idx + target * stride()] = col.second;
              ++target;
            }
            row.second.clear();
            ++row_idx;
          }
          entries.clear();

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tAx, this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(0), tAj, this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(1), tArl, this->_scalar_index.at(1));

          MemoryPool<Mem::Main>::instance()->release_memory(tAx);
          MemoryPool<Mem::Main>::instance()->release_memory(tAj);
          MemoryPool<Mem::Main>::instance()->release_memory(tArl);
        }

        void _read_from_mtx(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
          _read_from_m(file);
          file.close();
        }

        void _read_from_mtx(std::istream& file)
        {
          std::map<Index, std::map<Index, DT_> > entries; // map<row, map<column, value> >
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
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

            entries[row].insert(std::pair<Index, DT_>(col, val));
            ++ue;
          }
          this->_scalar_index.at(0) = this->rows() * this->columns();
          this->_scalar_index.at(5) = ue;
          this->_scalar_index.at(4) = 0;

          Index* tArl = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1));

          Index idx(0);
          for (auto row : entries)
          {
            tArl[idx] = Index(row.second.size());
            this->_scalar_index.at(4) = std::max(this->_scalar_index.at(4), tArl[idx]);
            ++idx;
          }

          Index alignment(32);
          this->_scalar_index.at(3) = alignment * ((this->_scalar_index.at(1) + alignment - 1)/ alignment);

          DT_* tAx = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem::Main>::instance()->set_memory(tAx, DT_(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));
          Index* tAj = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem::Main>::instance()->set_memory(tAj, Index(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));

          Index row_idx(0);
          for (auto row : entries)
          {
            Index target(0);
            for (auto col : row.second )
            {
              tAj[row_idx + target * stride()] = col.first;
              tAx[row_idx + target * stride()] = col.second;
              ++target;
            }
            row.second.clear();
            ++row_idx;
          }
          entries.clear();

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tAx, this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(0), tAj, this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(1), tArl, this->_scalar_index.at(1));

          MemoryPool<Mem::Main>::instance()->release_memory(tAx);
          MemoryPool<Mem::Main>::instance()->release_memory(tAj);
          MemoryPool<Mem::Main>::instance()->release_memory(tArl);
        }

        void _read_from_ell(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
          _read_from_ell(file);
          file.close();
        }

        void _read_from_ell(std::istream& file)
        {
          uint64_t dim;
          uint64_t rows;
          uint64_t columns;
          uint64_t stride;
          uint64_t num_cols_per_row;
          file.read((char *)&dim, sizeof(uint64_t));
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

          uint64_t * cAj = new uint64_t[std::size_t(dim)];
          file.read((char *)cAj, (long)(dim * sizeof(uint64_t)));
          Index* tAj = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(Index(dim));
          for (Index i(0) ; i < dim ; ++i)
            tAj[i] = Index(cAj[i]);
          delete[] cAj;

          double * cAx = new double[std::size_t(dim)];
          file.read((char *)cAx, (long)(dim * sizeof(double)));

          DT_* tAx = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(Index(dim));
          for (Index i(0) ; i < dim ; ++i)
            tAx[i] = DT_(cAx[i]);
          delete[] cAx;

          Index* tArl = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1));
          //compute row length vector
          this->_scalar_index.at(5) = 0;
          for (Index row(0) ; row < this->_scalar_index.at(1) ; ++row)
          {
            Index count(0);
            for (Index i(row) ; i < Index(dim) ; i += Index(stride))
            {
                if (tAx[i] == DT_(0))
                {
                  i = Index(dim);
                  break;
                }
                ++count;
                ++this->_scalar_index.at(5);
            }
            tArl[row] = count;
          }

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(Index(dim)));
          this->_elements_size.push_back(Index(dim));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_indices_size.push_back(Index(dim));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tAx, Index(dim));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(0), tAj, Index(dim));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(1), tArl, this->_scalar_index.at(1));
          MemoryPool<Mem::Main>::instance()->release_memory(tAx);
          MemoryPool<Mem::Main>::instance()->release_memory(tAj);
          MemoryPool<Mem::Main>::instance()->release_memory(tArl);
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Mem_ MemType;
        /// Compatible L-vector type
        typedef DenseVector<MemType, DataType> VectorTypeL;
        /// Compatible R-vector type
        typedef DenseVector<MemType, DataType> VectorTypeR;
        /// Our used layout type
        const static SparseLayoutType LayoutType = SparseLayoutType::lt_ell;

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
         * \param[in] layout The layout to be used.
         *
         * Creates an empty matrix with given layout.
         */
        explicit SparseMatrixELL(const SparseLayout<Mem_, LayoutType> & layout) :
          Container<Mem_, DT_> (layout._scalar_index.at(0))
        {
          CONTEXT("When creating SparseMatrixELL");
          this->_indices.assign(layout._indices.begin(), layout._indices.end());
          this->_indices_size.assign(layout._indices_size.begin(), layout._indices_size.end());
          this->_scalar_index.assign(layout._scalar_index.begin(), layout._scalar_index.end());
          this->_scalar_dt.push_back(DT_(0));

          for (auto i : this->_indices)
            MemoryPool<Mem_>::instance()->increase_memory(i);

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
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

          Index* tArl = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1));

          this->_scalar_index.at(4) = 0;
          for (Index i(0) ; i < this->_scalar_index.at(1) ; ++i)
          {
            tArl[i] = other.row_ptr()[i + 1] - other.row_ptr()[i];
            if (tArl[i] > this->_scalar_index.at(4))
              this->_scalar_index.at(4) = tArl[i];
          }

          Index alignment(32);
          this->_scalar_index.at(3) = alignment * ((this->_scalar_index.at(1) + alignment - 1)/ alignment);


          DT_* tAx = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem::Main>::instance()->set_memory(tAx, DT_(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));
          Index* tAj = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3));
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

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1)));
          this->_indices_size.push_back(this->_scalar_index.at(1));

          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tAx, this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(0), tAj, this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(1), tArl, this->_scalar_index.at(1));
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

          Index* tArl = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1));
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


          DT_* tAx = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem::Main>::instance()->set_memory(tAx, DT_(0), this->_scalar_index.at(4) * this->_scalar_index.at(3));
          Index* tAj = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3));
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

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_indices_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(this->_scalar_index.at(1)));
          this->_indices_size.push_back(this->_scalar_index.at(1));


          MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(0), tAx, this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(0), tAj, this->_scalar_index.at(4) * this->_scalar_index.at(3));
          MemoryPool<Mem_>::template upload<Index>(this->get_indices().at(1), tArl, this->_scalar_index.at(1));
          MemoryPool<Mem::Main>::instance()->release_memory(tAx);
          MemoryPool<Mem::Main>::instance()->release_memory(tAj);
          MemoryPool<Mem::Main>::instance()->release_memory(tArl);
        }

        /**
         * \brief Constructor
         *
         * \param[in] rows The row count of the created matrix.
         * \param[in] columns The column count of the created matrix.
         * \param[in] Ax Vector with non zero elements.
         * \param[in] Aj Vector with column indices.
         * \param[in] Arl Vector with row length.
         *
         * Creates a matrix with given dimensions and content.
         */
        explicit SparseMatrixELL(const Index rows, const Index columns,
            const Index stride, const Index num_cols_per_row, const Index used_elements,
            DenseVector<Mem_, DT_> & Ax,
            DenseVector<Mem_, Index> & Aj,
            DenseVector<Mem_, Index> & Arl) :
          Container<Mem_, DT_>(rows * columns)
        {
          CONTEXT("When creating SparseMatrixELL");
          this->_scalar_index.push_back(rows);
          this->_scalar_index.push_back(columns);
          this->_scalar_index.push_back(stride);
          this->_scalar_index.push_back(num_cols_per_row);
          this->_scalar_index.push_back(used_elements);
          this->_scalar_dt.push_back(DT_(0));

          this->_elements.push_back(Ax.elements());
          this->_elements_size.push_back(Ax.size());
          this->_indices.push_back(Aj.elements());
          this->_indices_size.push_back(Aj.size());
          this->_indices.push_back(Arl.elements());
          this->_indices_size.push_back(Arl.size());

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
        }

        /**
         * \brief Constructor
         *
         * \param[in] mode The used file format.
         * \param[in] filename The source file.
         *
         * Creates a ELL matrix based on the source file.
         */
        explicit SparseMatrixELL(FileMode mode, String filename) :
          Container<Mem_, DT_>(0)
        {
          CONTEXT("When creating SparseMatrixELL");

          switch(mode)
          {
            case FileMode::fm_m:
              _read_from_m(filename);
              break;
            case FileMode::fm_mtx:
              _read_from_mtx(filename);
              break;
            case FileMode::fm_ell:
              _read_from_ell(filename);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
        }

        /**
         * \brief Constructor
         *
         * \param[in] mode The used file format.
         * \param[in] file The source filestream.
         *
         * Creates a ELL matrix based on the source filestream.
         */
        explicit SparseMatrixELL(FileMode mode, std::istream& file) :
          Container<Mem_, DT_>(0)
        {
          CONTEXT("When creating SparseMatrixELL");

          switch(mode)
          {
            case FileMode::fm_m:
              _read_from_m(file);
              break;
            case FileMode::fm_mtx:
              _read_from_mtx(file);
              break;
            case FileMode::fm_ell:
              _read_from_ell(file);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
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
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source matrix.
         *
         * Moves a given matrix to this matrix.
         */
        SparseMatrixELL(SparseMatrixELL<Mem_, DT_> && other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When moving SparseMatrixELL");
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a copy of a given matrix from another memory architecture.
         */
        template <typename Arch2_, typename DT2_>
        explicit SparseMatrixELL(const SparseMatrixELL<Arch2_, DT2_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseMatrixELL");
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this matrix.
         */
        SparseMatrixELL<Mem_, DT_> clone() const
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
         * \brief Move operator
         *
         * \param[in] other The source matrix.
         *
         * Moves another matrix to the target matrix.
         */
        SparseMatrixELL<Mem_, DT_> & operator= (SparseMatrixELL<Mem_, DT_> && other)
        {
          CONTEXT("When moving SparseMatrixELL");

          this->move(std::move(other));

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
         * \brief Assignment operator
         *
         * \param[in] layout A sparse matrix layout.
         *
         * Assigns a new matrix layout, discarding all old data
         */
        SparseMatrixELL<Mem_, DT_> & operator= (const SparseLayout<Mem_, LayoutType> & layout)
        {
          CONTEXT("When assigning SparseMatrixELL");

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();
          this->_scalar_index.clear();
          this->_scalar_dt.clear();

          this->_indices.assign(layout._indices.begin(), layout._indices.end());
          this->_indices_size.assign(layout._indices_size.begin(), layout._indices_size.end());
          this->_scalar_index.assign(layout._scalar_index.begin(), layout._scalar_index.end());
          this->_scalar_dt.push_back(DT_(0));

          for (auto i : this->_indices)
            MemoryPool<Mem_>::instance()->increase_memory(i);

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(4) * this->_scalar_index.at(3)));
          this->_elements_size.push_back(this->_scalar_index.at(4) * this->_scalar_index.at(3));

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
            case FileMode::fm_ell:
              write_out_ell(filename);
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
          CONTEXT("When writing out SparseMatrixELL");

          switch(mode)
          {
            case FileMode::fm_ell:
              write_out_ell(file);
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
         * \brief Write out matrix to ell binary file.
         *
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out_ell(String filename) const
        {
          std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
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
          if (! std::is_same<DT_, double>::value)
            std::cout<<"Warning: You are writing out an ell matrix with less than double precission!"<<std::endl;

          const Index dim(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          Index * Aj = MemoryPool<Mem::Main>::instance()->template allocate_memory<Index>(dim);
          MemoryPool<Mem_>::template download<Index>(Aj, this->_indices.at(0), dim);
          uint64_t * cAj = new uint64_t[dim];
          for (Index i(0) ; i < dim ; ++i)
            cAj[i] = Aj[i];
          MemoryPool<Mem::Main>::instance()->release_memory(Aj);

          DT_ * Ax = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(dim);
          MemoryPool<Mem_>::template download<DT_>(Ax, this->_elements.at(0), dim);
          double * cAx = new double[dim];
          for (Index i(0) ; i < dim ; ++i)
            cAx[i] = Type::Traits<DT_>::to_double(Ax[i]);
          MemoryPool<Mem::Main>::instance()->release_memory(Ax);

          uint64_t tdim(dim);
          uint64_t rows(this->_scalar_index.at(1));
          uint64_t columns(this->_scalar_index.at(2));
          uint64_t stride(this->_scalar_index.at(3));
          uint64_t num_cols_per_row(this->_scalar_index.at(4));
          file.write((const char *)&tdim, sizeof(uint64_t));
          file.write((const char *)&rows, sizeof(uint64_t));
          file.write((const char *)&columns, sizeof(uint64_t));
          file.write((const char *)&stride, sizeof(uint64_t));
          file.write((const char *)&num_cols_per_row, sizeof(uint64_t));
          file.write((const char *)cAj, (long)(tdim * sizeof(uint64_t)));
          file.write((const char *)cAx, (long)(tdim * sizeof(double)));

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
          SparseMatrixELL<Mem::Main, DT_> temp(*this);

          file << "data = [" << std::endl;
          const Index dim(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          for (Index i(0) ; i < dim ; ++i)
          {
            if (temp.Ax()[i] != DT_(0))
            {
              file << (i%this->_scalar_index.at(3)) + 1 << " " << temp.Aj()[i] + 1 << " " << std::scientific << Type::Traits<DT_>::to_double(temp.Ax()[i]) << ";" << std::endl;
            }
          }
          file << "];" << std::endl;
          file << "mat=sparse(data(:,1),data(:,2),data(:,3));";
        }

        /**
         * \brief Write out matrix to MatrixMarktet mtx file.
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
         * \brief Write out matrix to MatrixMarktet mtx file.
         *
         * \param[in] file The stream that shall be written to.
         */
        void write_out_mtx(std::ostream& file) const
        {
          SparseMatrixELL<Mem::Main, DT_> temp(*this);

          file << "%%MatrixMarket matrix coordinate real general" << std::endl;
          file << temp.rows() << " " << temp.columns() << " " << temp.used_elements() << std::endl;

          const Index dim(this->_scalar_index.at(4) * this->_scalar_index.at(3));
          for (Index i(0) ; i < dim ; ++i)
          {
            if (temp.Ax()[i] != DT_(0))
            {
              file << (i%this->_scalar_index.at(3)) + 1 << " " << temp.Aj()[i] + 1 << " " << std::scientific << Type::Traits<DT_>::to_double(temp.Ax()[i]) << ";" << std::endl;
            }
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
          CONTEXT("When retrieving SparseMatrixELL element");

          ASSERT(row < this->_scalar_index.at(1), "Error: " + stringify(row) + " exceeds sparse matrix ell row size " + stringify(this->_scalar_index.at(1)) + " !");
          ASSERT(col < this->_scalar_index.at(2), "Error: " + stringify(col) + " exceeds sparse matrix ell column size " + stringify(this->_scalar_index.at(2)) + " !");

          if (std::is_same<Mem_, Mem::Main>::value)
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
         * \brief Retrieve convenient sparse matrix layout object.
         *
         * \return An object containing the sparse matrix layout.
         */
        SparseLayout<Mem_, LayoutType> layout() const
        {
          SparseLayout<Mem_, LayoutType> layout(this->_indices, this->_indices_size, this->_scalar_index);
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
        virtual Index used_elements() const
        {
          return this->_scalar_index.at(5);
        }

        /**
         * \brief Retrieve column indices array.
         *
         * \returns Column indices array.
         */
        Index const * Aj() const
        {
          return this->_indices.at(0);
        }

        /**
         * \brief Retrieve non zero element array.
         *
         * \returns Non zero element array.
         */
        DT_ * Ax()
        {
          return this->_elements.at(0);
        }

        DT_ const * Ax() const
        {
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve row length array.
         *
         * \returns Row lenght array.
         */
        Index const * Arl() const
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
        static String name()
        {
          return "SparseMatrixELL";
        }

        /**
         * \brief Performs \f$this \leftarrow x\f$.
         *
         * \param[in] x The Matrix to be copied.
         */
        void copy(const SparseMatrixELL<Mem_, DT_> & x)
        {
          this->_copy_content(x);
        }

        /**
         * \brief Performs \f$this \leftarrow x\f$.
         *
         * \param[in] x The Matrix to be copied.
         */
        template <typename Mem2_>
        void copy(const SparseMatrixELL<Mem2_, DT_> & x)
        {
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
          const SparseMatrixELL<Mem_, DT_> & x,
          const SparseMatrixELL<Mem_, DT_> & y,
          DT_ alpha = DT_(1))
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
          if (x.stride() != y.stride())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix stride do not match!");
          if (x.stride() != this->stride())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix stride do not match!");
          if (x.num_cols_per_row() != y.num_cols_per_row())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_cols_per_row do not match!");
          if (x.num_cols_per_row() != this->num_cols_per_row())
            throw InternalError(__func__, __FILE__, __LINE__, "Matrix num_cols_per_row do not match!");

          // check for special cases
          // r <- x + y
          if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
            Arch::Sum<Mem_, Algo_>::value(this->Ax(), x.Ax(), y.Ax(), this->stride() * this->num_cols_per_row());
          // r <- y - x
          else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
            Arch::Difference<Mem_, Algo_>::value(this->Ax(), y.Ax(), x.Ax(), this->stride() * this->num_cols_per_row());
          // r <- y
          else if (Math::abs(alpha) < Math::eps<DT_>())
            this->copy(y);
          // r <- y + alpha*x
          else
            Arch::Axpy<Mem_, Algo_>::dv(this->Ax(), alpha, x.Ax(), y.Ax(), this->stride() * this->num_cols_per_row());
        }

        /**
         * \brief Calculate \f$this \leftarrow \alpha x\f$
         *
         * \param[in] x The matrix to be scaled.
         * \param[in] alpha A scalar to scale x with.
         */
        template <typename Algo_>
        void scale(const SparseMatrixELL<Mem_, DT_> & x, const DT_ alpha)
        {
          if (x.rows() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
          if (x.columns() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");

          Arch::Scale<Mem_, Algo_>::value(this->Ax(), x.Ax(), alpha, this->stride() * this->num_cols_per_row());
        }

        /**
         * \brief Calculate \f$ r \leftarrow this\cdot x \f$
         *
         * \param[out] r The vector that recieves the result.
         * \param[in] x The vector to be multiplied by this matrix.
         */
        template<typename Algo_>
        void apply(DenseVector<Mem_,DT_>& r, const DenseVector<Mem_, DT_>& x) const
        {
          if (r.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
          if (x.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

          Arch::ProductMatVec<Mem_, Algo_>::ell(r.elements(), this->Ax(), this->Aj(), this->Arl(),
            x.elements(), this->stride(), this->rows());
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
          DenseVector<Mem_,DT_>& r,
          const DenseVector<Mem_, DT_>& x,
          const DenseVector<Mem_, DT_>& y,
          const DT_ alpha = DT_(1)) const
        {
          if (r.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
          if (x.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
          if (y.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

          // check for special cases
          // r <- y - A*x
          if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          {
            Arch::Defect<Mem_, Algo_>::ell(r.elements(), y.elements(), this->Ax(), this->Aj(),
              this->Arl(), x.elements(), this->stride(), this->rows());
          }
          // r <- y
          else if (Math::abs(alpha) < Math::eps<DT_>())
            r.copy(y);
          // r <- y + alpha*x
          else
          {
            Arch::Axpy<Mem_, Algo_>::ell(r.elements(), alpha, x.elements(), y.elements(),
              this->Ax(), this->Aj(), this->Arl(), this->stride(), this->rows());
          }
        }

        /// Returns a new compatible L-Vector.
        VectorTypeL create_vector_l() const
        {
          return VectorTypeL(this->columns());
        }

        /// Returns a new compatible R-Vector.
        VectorTypeR create_vector_r() const
        {
          return VectorTypeR(this->rows());
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
