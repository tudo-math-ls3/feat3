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
#include <kernel/lafem/arch/scale_row_col.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/adjacency/graph.hpp>

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
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
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
     * Refer to \ref lafem_design for general usage informations.
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
            swap_key = Util::MemoryPool<Mem_>::get_element(key, i);
            swap1 = Util::MemoryPool<Mem_>::get_element(val1, i);
            swap2 = Util::MemoryPool<Mem_>::get_element(val2, i);
            j = i;
            while (j > 0 && Util::MemoryPool<Mem_>::get_element(key, j - 1) > swap_key)
            {
              Util::MemoryPool<Mem_>::copy(key + j, key + j - 1, 1);
              Util::MemoryPool<Mem_>::copy(val1 + j, val1 + j - 1, 1);
              Util::MemoryPool<Mem_>::copy(val2 + j, val2 + j - 1, 1);
              --j;
            }
            Util::MemoryPool<Mem_>::set_memory(key + j, swap_key);
            Util::MemoryPool<Mem_>::set_memory(val1 + j, swap1);
            Util::MemoryPool<Mem_>::set_memory(val2 + j, swap2);
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
          if (Util::MemoryPool<Mem_>::get_element(key, i) > Util::MemoryPool<Mem_>::get_element(key, i + gap))
          {
          swap_key = Util::MemoryPool<Mem_>::get_element(key, i);
          Util::MemoryPool<Mem_>::copy(key + i, key + i + gap, 1);
          Util::MemoryPool<Mem_>::set_memory(key + i + gap, swap_key);

          swap1 = Util::MemoryPool<Mem_>::get_element(val1, i);
          Util::MemoryPool<Mem_>::copy(val1 + i, val1 + i + gap, 1);
          Util::MemoryPool<Mem_>::set_memory(val1 + i + gap, swap1);

          swap2 = Util::MemoryPool<Mem_>::get_element(val2, i);
          Util::MemoryPool<Mem_>::copy(val2 + i, val2 + i + gap, 1);
          Util::MemoryPool<Mem_>::set_memory(val2 + i + gap, swap2);

          swapped = true;
          }
          }
          }
          }*/

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
          const bool general((line.find("%%MatrixMarket matrix coordinate real general") != String::npos) ? true : false);
          const bool symmetric((line.find("%%MatrixMarket matrix coordinate real symmetric") != String::npos) ? true : false);

          if (symmetric == false && general == false)
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a compatible mtx-file");
          }

          while(!file.eof())
          {
            std::getline(file,line);
            if (file.eof())
              throw InternalError(__func__, __FILE__, __LINE__, "Input-file is empty");

            String::size_type begin(line.find_first_not_of(" "));
            if (line.at(begin) != '%')
              break;
          }
          {
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
            _rows() = row;
            _columns() = col;
            _size() = this->rows() * this->columns();
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
            IT_ row((IT_)atol(srow.c_str()));
            --row;
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scol(line, 0, end);
            IT_ col((IT_)atol(scol.c_str()));
            --col;
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String sval(line, 0, end);
            DT_ tval((DT_)atof(sval.c_str()));

            entries[(unsigned int)row].insert(std::pair<IT_, DT_>(col, tval));
            ++ue;
            if (symmetric == true && row != col)
            {
              entries[(unsigned int)col].insert(std::pair<IT_, DT_>(row, tval));
              ++ue;
            }
          }

          _size() = this->rows() * this->columns();
          _used_elements() = ue;
          _allocated_elements() = ue;

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

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());

          Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), trow, _used_elements());
          Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), tcolumn, _used_elements());

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
          this->template _deserialise<double, uint64_t>(FileMode::fm_coo, file);
        }

        Index & _size()
        {
          return this->_scalar_index.at(0);
        }

        Index & _rows()
        {
          return this->_scalar_index.at(1);
        }

        Index & _columns()
        {
          return this->_scalar_index.at(2);
        }

        Index & _used_elements()
        {
          return this->_scalar_index.at(3);
        }

        Index & _allocated_elements()
        {
          return this->_scalar_index.at(4);
        }

        Index & _alloc_increment()
        {
          return this->_scalar_index.at(5);
        }

        Index & _sorted()
        {
          return this->_scalar_index.at(6);
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our indextype
        typedef IT_ IndexType;
        /// Our memory architecture type
        typedef Mem_ MemType;
        /// Compatible L-vector type
        typedef DenseVector<Mem_, DT_, IT_> VectorTypeL;
        /// Compatible R-vector type
        typedef DenseVector<Mem_, DT_, IT_> VectorTypeR;
        /// Our used layout type
        static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_coo;
        /// Our 'base' class type
        template <typename Mem2_, typename DT2_, typename IT2_ = IT_>
        using ContainerType = class SparseMatrixCOO<Mem2_, DT2_, IT2_>;

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
        explicit SparseMatrixCOO(Index rows_in, Index columns_in) :
          Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        CONTEXT("When creating SparseMatrixCOO");
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
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
        explicit SparseMatrixCOO(const SparseLayout<Mem_, IT_, layout_id> & layout_in) :
          Container<Mem_, DT_, IT_> (layout_in._scalar_index.at(0))
      {
        CONTEXT("When creating SparseMatrixCOO");
        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          Util::MemoryPool<Mem_>::instance()->increase_memory(i);

        this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_allocated_elements()));
        this->_elements_size.push_back(_allocated_elements());
      }

        /**
         * \brief Constructor
         *
         * \param[in] other The source matrix.
         *
         * Creates a COO matrix based on the source matrix.
         */
        template <typename MT_>
        explicit SparseMatrixCOO(const MT_ & other) :
          Container<Mem_, DT_, IT_>(other.size())
      {
        CONTEXT("When creating SparseMatrixCOO");

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
        explicit SparseMatrixCOO(const Index rows_in, const Index columns_in, DenseVector<Mem_, IT_, IT_> & row_ind,
        DenseVector<Mem_, IT_, IT_> & col_ind, DenseVector<Mem_, DT_, IT_> & val_in) :
          Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        CONTEXT("When creating SparseMatrixCOO");
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(val_in.size());
        this->_scalar_index.push_back(val_in.size());
        this->_scalar_index.push_back(1000);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(row_ind.elements());
        this->_indices_size.push_back(row_ind.size());
        this->_indices.push_back(col_ind.elements());
        this->_indices_size.push_back(col_ind.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
      }

        /**
         * \brief Constructor
         *
         * \param[in] graph The.graph to create the matrix from
         *
         * Creates a COO matrix based on a given adjacency graph, representing the sparsity pattern.
         */
        explicit SparseMatrixCOO(const Adjacency::Graph & graph) :
          Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating SparseMatrixCOO");

        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nnze = graph.get_num_indices();

        // create temporary vectors
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vrow_idx(num_nnze);
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vcol_idx(num_nnze);
        LAFEM::DenseVector<Mem::Main, DT_, IT_> vdata(num_nnze, DT_(0));

        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());
        IT_ * prow_idx(vrow_idx.elements());
        IT_ * pcol_idx(vcol_idx.elements());

        // build col-idx and row-idx
        for (Index i(0); i < num_rows; ++i)
        {
          for (Index j(dom_ptr[i]); j < dom_ptr[i+1]; ++j)
          {
            prow_idx[j] = IT_(i);
            pcol_idx[j] = IT_(img_idx[j]);
          }
        }

        // build the matrix
        this->assign(SparseMatrixCOO<Mem::Main, DT_, IT_>(num_rows, num_cols, vrow_idx, vcol_idx, vdata));
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
         * \brief Constructor
         *
         * \param[in] std::pair<Index, char *> A std::pair, containing byte array size and byte array pointer.
         *
         * Creates a matrix from the given byte array.
         */
        template <typename DT2_ = DT_, typename IT2_ = IT_>
        explicit SparseMatrixCOO(std::pair<Index, char *> input) :
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseMatrixCOO");
          deserialise<DT2_, IT2_>(input);
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source matrix.
         *
         * Moves a given matrix to this matrix.
         */
        SparseMatrixCOO(SparseMatrixCOO && other) :
          Container<Mem_, DT_, IT_>(std::forward<SparseMatrixCOO>(other))
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

          this->move(std::forward<SparseMatrixCOO>(other));

          return *this;
        }

        /** \brief Clone operation
         *
         * Create a deep copy of itself.
         *
         * \param[in] clone_indices Should we create a deep copy of the index arrays, too ?
         *
         * \return A deep copy of itself.
         *
         */
        SparseMatrixCOO clone(bool clone_indices = true) const
        {
          CONTEXT("When cloning SparseMatrixCOO");
          SparseMatrixCOO t;
          t.clone(*this, clone_indices);
          return t;
        }

        /** \brief Clone operation
         *
         * Become a deep copy of a given matrix.
         *
         * \param[in] other The source matrix.
         * \param[in] clone_indices Should we create a deep copy of the index arrays, too ?
         *
         */
        void clone(const SparseMatrixCOO & other, bool clone_indices = true)
        {
          CONTEXT("When cloning SparseMatrixCOO");
          Container<Mem_, DT_, IT_>::clone(other, clone_indices);
        }

        /** \brief Clone operation
         *
         * Become a deep copy of a given matrix.
         *
         * \param[in] other The source matrix.
         * \param[in] clone_indices Should we create a deep copy of the index arrays, too ?
         *
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void clone(const SparseMatrixCOO<Mem2_, DT2_, IT2_> & other, bool clone_indices = true)
        {
          CONTEXT("When cloning SparseMatrixCOO");
          SparseMatrixCOO t;
          t.assign(other);
          Container<Mem_, DT_, IT_>::clone(t, clone_indices);
        }

        /**
         * \brief Conversion method
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
         * \brief Conversion method
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

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());

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
            const IT_ end(cother.row_ptr()[row + 1]);
            for (IT_ i(cother.row_ptr()[row]) ; i < end ; ++i)
            {
              tval[ue] = cother.val()[i];
              trow[ue] = IT_(row);
              tcolumn[ue] = cother.col_ind()[i];
              ++ue;
            }
          }

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), trow, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), tcolumn, _used_elements());
            delete[] tval;
            delete[] trow;
            delete[] tcolumn;
          }
        }

        /**
         * \brief Conversion method
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

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());

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
              trow[ue] = IT_(row);
              tcolumn[ue] = *tAj;

              tAj += stride;
              tAx += stride;
              ++ue;
            }
          }

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), trow, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), tcolumn, _used_elements());
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
        void convert(const SparseMatrixBanded<Mem2_, DT2_, IT2_> & other)
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

          SparseMatrixBanded<Mem::Main, DT_, IT_> cother;
          cother.convert(other);

          this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_used_elements()));
          this->_elements_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());
          this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_used_elements()));
          this->_indices_size.push_back(_used_elements());

          DT_ * tval(nullptr);
          IT_ * tcol_ind(nullptr);
          IT_ * trow_ind(nullptr);
          if (std::is_same<Mem_, Mem::Main>::value)
          {
            tval = this->_elements.at(0);
            trow_ind = this->_indices.at(0);
            tcol_ind = this->_indices.at(1);
          }
          else
          {
            tval = new DT_[other.used_elements()];
            trow_ind = new IT_[other.used_elements()];
            tcol_ind = new IT_[other.used_elements()];
          }

          const DT_ * cval(cother.val());
          const IT2_ * coffsets(cother.offsets());
          const Index cnum_of_offsets(cother.num_of_offsets());
          const Index crows(cother.rows());

          // Search first offset of the upper triangular matrix
          Index k(0);
          while (k < cnum_of_offsets && coffsets[k] + 1 < crows)
          {
            ++k;
          }

          // iteration over all offsets of the lower triangular matrix
          Index ue(0);
          for (Index i(k + 1); i > 0;)
          {
            --i;

            // iteration over all offsets of the upper triangular matrix
            for (Index j(cnum_of_offsets + 1); j > 0;)
            {
              --j;

              // iteration over all rows which contain the offsets between offset i and offset j
              const Index start(Math::max(cother.start_offset(i),
                                          cother.end_offset(j) + 1));
              const Index end  (Math::min(cother.start_offset(i-1),
                                          cother.end_offset(j-1) + 1));
              for (Index l(start); l < end; ++l)
              {
                for (Index a(i); a < j; ++a)
                {
                  tval[ue] = cval[a * crows + l];
                  tcol_ind[ue] = IT_(l + coffsets[a] + 1 - crows);
                  trow_ind[ue] = IT_(l);
                  ++ue;
                }
              }
            }
          }

          if (! std::is_same<Mem_, Mem::Main>::value)
          {
            Util::MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), trow_ind, _used_elements());
            Util::MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), tcol_ind, _used_elements());
            delete[] tval;
            delete[] tcol_ind;
            delete[] trow_ind;
          }
        }

        /**
         * \brief Conversion method
         *
         * \param[in] a The input matrix.
         *
         * Converts any matrix to SparseMatrixCOO-format
         */
        template <typename MT_>
        void convert(const MT_ & a)
        {
          CONTEXT("When converting SparseMatrixCOO");

          typename MT_::template ContainerType<Mem::Main, DT_, IT_> ta;
          ta.convert(a);

          const Index arows(ta.rows());
          const Index acolumns(ta.columns());
          const Index aused_elements(ta.used_elements());

          DenseVector<Mem::Main, DT_, IT_> tval(aused_elements);
          DenseVector<Mem::Main, IT_, IT_> tcol_ind(aused_elements);
          DenseVector<Mem::Main, IT_, IT_> trow_ind(aused_elements);

          DT_ * pval(tval.elements());
          IT_ * pcol_ind(tcol_ind.elements());
          IT_ * prow_ind(trow_ind.elements());

          DenseVector<Mem::Main, IT_, IT_> trow_ptr(arows + 1);
          IT_ * prow_ptr(trow_ptr.elements());

          for (Index i(0); i < arows; ++i)
          {
            prow_ptr[i + 1] = IT_(ta.get_length_of_line(i));
          }

          prow_ptr[0] = IT_(0);

          for (Index i(1); i < arows + 1; ++i)
          {
            prow_ptr[i] += prow_ptr[i - 1];
          }

          for (Index i(0); i < arows; ++i)
          {
            ta.set_line(i, pval + prow_ptr[i], pcol_ind + prow_ptr[i], 0);
          }

          for (Index i(0); i < arows; ++i)
          {
            for (Index k(prow_ptr[i]); k < prow_ptr[i+1]; ++k)
            {
              prow_ind[k] = IT_(i);
            }
          }

          SparseMatrixCOO<Mem::Main, DT_, IT_> ta_coo(arows, acolumns, trow_ind, tcol_ind, tval);
          SparseMatrixCOO<Mem_, DT_, IT_> a_coo;
          a_coo.convert(ta_coo);

          this->assign(a_coo);
        }

        /**
         * \brief Deserialisation of complete container entity.
         *
         * \param[in] std::pair<Index, char *> A std::pair, containing byte array size and byte array pointer.
         *
         * Recreate a complete container entity by a single binary array.
         */
        template <typename DT2_ = DT_, typename IT2_ = IT_>
        void deserialise(std::pair<Index, char *> input)
        {
          this->template _deserialise<DT2_, IT2_>(FileMode::fm_coo, input);
        }

        /**
         * \brief Serialisation of complete container entity.
         *
         * \param[in] mode FileMode enum, describing the actual container specialisation.
         * \param[out] std::pair<Index, char *> A std::pair, containing byte array size and byte array pointer.
         *
         * Serialise a complete container entity into a single binary array.
         *
         * \warning The allocated array must be freed by the user!
         *
         * See \ref FEAST::LAFEM::Container::_serialise for details.
         */
        template <typename DT2_ = DT_, typename IT2_ = IT_>
        std::pair<Index, char *> serialise()
        {
          return this->template _serialise<DT2_, IT2_>(FileMode::fm_coo);
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

          if (_sorted() == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();

          switch(mode)
          {
            case FileMode::fm_mtx:
              write_out_mtx(filename);
              break;
            case FileMode::fm_coo:
              write_out_coo(filename);
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
            case FileMode::fm_mtx:
              write_out_mtx(file);
              break;
            case FileMode::fm_coo:
              write_out_coo(file);
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

          this->template _serialise<double, uint64_t>(FileMode::fm_coo, file);
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
            file << temp.row_indices()[i] + 1 << " " << temp.column_indices()[i] + 1 << " " << std::scientific << temp.val()[i] << std::endl;
          }
        }

        /**
         * \brief Set specific matrix element.
         *
         * \param[in] row The row of the matrix element.
         * \param[in] col The column of the matrix element.
         * \param[in] value The value to be set.
         */
        void operator()(Index row, Index col, DT_ value)
        {
          CONTEXT("When setting SparseMatrixCOO element");

          ASSERT(row < this->rows(), "Error: " + stringify(row) + " exceeds sparse matrix coo row size " + stringify(this->rows()) + " !");
          ASSERT(col < this->columns(), "Error: " + stringify(col) + " exceeds sparse matrix coo column size " + stringify(this->columns()) + " !");

          _sorted() = 0;

          if (this->_elements.size() == 0)
          {
            this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(alloc_increment()));
            this->_elements_size.push_back(alloc_increment());
            Util::MemoryPool<Mem_>::instance()->template set_memory<DT_>(this->_elements.back(), DT_(4711), alloc_increment());
            this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(alloc_increment()));
            this->_indices_size.push_back(alloc_increment());
            Util::MemoryPool<Mem_>::instance()->template set_memory<IT_>(this->_indices.back(), IT_(4711), alloc_increment());
            this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(alloc_increment()));
            this->_indices_size.push_back(alloc_increment());
            Util::MemoryPool<Mem_>::instance()->template set_memory<IT_>(this->_indices.back(), IT_(4711), alloc_increment());

            Util::MemoryPool<Mem_>::set_memory(this->_elements.at(0), value);
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(0), IT_(row));
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(1), IT_(col));

            _used_elements() = 1;
            _allocated_elements() = alloc_increment();
            return;
          }

          //append element
          if (_used_elements() < _allocated_elements())
          {
            Util::MemoryPool<Mem_>::set_memory(this->_elements.at(0) + _used_elements(), value);
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(0) + _used_elements(), IT_(row));
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(1) + _used_elements(), IT_(col));

            ++_used_elements();
          }
          //reallocate
          else
          {
            _allocated_elements() += alloc_increment();

            DT_ * elements_new(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(_allocated_elements()));
            Util::MemoryPool<Mem_>::instance()->template set_memory<DT_>(elements_new, DT_(4711), _allocated_elements());
            IT_ * rows_new(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_allocated_elements()));
            Util::MemoryPool<Mem_>::instance()->template set_memory<IT_>(rows_new, IT_(4711), _allocated_elements());
            IT_ * cols_new(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(_allocated_elements()));
            Util::MemoryPool<Mem_>::instance()->template set_memory<IT_>(cols_new, IT_(4711), _allocated_elements());

            Util::MemoryPool<Mem_>::copy(elements_new, this->_elements.at(0), _used_elements());
            Util::MemoryPool<Mem_>::copy(rows_new, this->_indices.at(0), _used_elements());
            Util::MemoryPool<Mem_>::copy(cols_new, this->_indices.at(1), _used_elements());

            Util::MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(0));
            Util::MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(0));
            Util::MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(1));

            this->_elements.at(0) = elements_new;
            this->_indices.at(0) = rows_new;
            this->_indices.at(1) = cols_new;

            Util::MemoryPool<Mem_>::set_memory(this->_elements.at(0) + _used_elements(), value);
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(0) + _used_elements(), IT_(row));
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(1) + _used_elements(), IT_(col));

            ++_used_elements();
            this->_elements_size.at(0) = _allocated_elements();
            this->_indices_size.at(0) = _allocated_elements();
            this->_indices_size.at(1) = _allocated_elements();
          }
        }

        void sort()
        {
          if (_sorted() == 0)
          {
            //first of all, mark matrix as sorted, because otherwise we would call ourselves inifite times
            _sorted() = 1;

            // check if there is anything to be sorted
            if(_used_elements() <= Index(0))
              return;

            // sort elements by row index
            _insertion_sort(this->_indices.at(0), val(), this->_indices.at(1), _used_elements());

            Index row_start(0);
            for (IT_ rowi(0) ; rowi < IT_(rows()) ; ++rowi)
            {
              Index offset(0);
              // explore range of elements in given row
              for ( ; row_start + offset < _used_elements() && Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), row_start + offset) == rowi ; ++offset) ;
              // sort range of elements in given row by column index
              _insertion_sort(this->_indices.at(1) + row_start, val() + row_start, this->_indices.at(0) + row_start, offset);
              // find and mark duplicate entries
              for (Index i(row_start + 1) ; i < row_start + offset ; ++i)
              {
                if (Util::MemoryPool<Mem_>::get_element(this->_indices.at(1), i - 1) == Util::MemoryPool<Mem_>::get_element(this->_indices.at(1), i))
                {
                  Util::MemoryPool<Mem_>::set_memory(this->_indices.at(0) + i - 1, std::numeric_limits<IT_>::max());
                }
              }
              row_start += offset;
            }

            // sort out marked duplicated elements
            _insertion_sort(this->_indices.at(0), val(), this->_indices.at(1), _used_elements());
            Index junk(0);
            while (Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), _used_elements() - 1 - junk) == std::numeric_limits<IT_>::max()
                && junk < _used_elements())
              ++junk;
            _used_elements() -= junk;
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
            return zero_element();

          if (sorted() == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();

          Index i(0);
          while (i < used_elements())
          {
            if (Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i) >= IT_(row))
              break;
            ++i;
          }

          while (i < used_elements())
          {
            if (Util::MemoryPool<Mem_>::get_element(this->_indices.at(1),i) >= IT_(col) || Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i) > IT_(row))
              break;
            ++i;
          }

          if(i < used_elements() && Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i) == IT_(row) && Util::MemoryPool<Mem_>::get_element(this->_indices.at(1), i) == IT_(col))
          {
            return Util::MemoryPool<Mem_>::get_element(this->_elements.at(0), i);
          }
          else
            return zero_element();
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
        SparseLayout<Mem_, IT_, layout_id> layout() const
        {
          SparseMatrixCOO t;
          t.clone(*this);
          return SparseLayout<Mem_, IT_, layout_id>(t._indices, t._indices_size, t._scalar_index);
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
          if (sorted() == 0)
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
          if (sorted() == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_elements.at(0);
        }

        DT_ const * val() const
        {
          if (sorted() == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_elements.at(0);
        }

        /**
         * \brief Retrieve row coordinate array.
         *
         * \returns Row coordinate array.
         */
        IT_ * row_indices()
        {
          if (sorted() == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_indices.at(0);
        }

        IT_ const * row_indices() const
        {
          if (sorted() == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_indices.at(0);
        }

        /**
         * \brief Retrieve columns coordinate array.
         *
         * \returns Column coordinate array.
         */
        IT_ * column_indices()
        {
          if (sorted() == 0)
            const_cast<SparseMatrixCOO *>(this)->sort();
          return this->_indices.at(1);
        }

        IT_ const * column_indices() const
        {
          if (sorted() == 0)
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

        ///@name Linear algebra operations
        ///@{
        /**
         * \brief Calculate \f$this \leftarrow y + \alpha x\f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
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
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
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
         * \brief Calculates the Frobenius norm of this matrix.
         *
         * \returns The Frobenius norm of this matrix.
         */
        template <typename Algo_>
        DT_ norm_frobenius() const
        {
          return Arch::Norm2<Mem_, Algo_>::value(this->val(), this->used_elements());
        }

        /**
         * \brief Calculate \f$this^\top \f$
         *
         * \return The transposed matrix
         */
        SparseMatrixCOO transpose()
        {
          SparseMatrixCOO x_t;
          x_t.transpose(*this);
          return x_t;
        }

        /**
         * \brief Calculate \f$this \leftarrow x^\top \f$
         *
         * \param[in] x The matrix to be transposed.
         */
        void transpose(const SparseMatrixCOO & x)
        {
          SparseMatrixCOO<Mem::Main, DT_, IT_> tx;
          tx.convert(x);

          const Index txrows(tx.rows());
          const Index txcolumns(tx.columns());
          const Index txused_elements(tx.used_elements());

          const IT_ * ptxrow_ind(tx.row_indices());
          const IT_ * ptxcol_ind(tx.column_indices());
          const DT_ * ptxval(tx.val());

          DenseVector<Mem::Main, IT_, IT_> tcol_ind(txused_elements);
          DenseVector<Mem::Main, DT_, IT_> tval(txused_elements);
          DenseVector<Mem::Main, IT_, IT_> trow_ind(txused_elements);

          IT_ * ptcol_ind(tcol_ind.elements());
          DT_ * ptval(tval.elements());
          IT_ * ptrow_ind(trow_ind.elements());

          DenseVector<Mem::Main, IT_, IT_> trow_ptr(txcolumns + 1, IT_(0));
          IT_ * ptrow_ptr(trow_ptr.elements());

          for (Index i(0); i < txused_elements; ++i)
          {
            ++ptrow_ptr[ptxcol_ind[i] + 1];
          }

          ptrow_ptr[0] = 0;

          for (Index i(1); i < txcolumns - 1; ++i)
          {
            ptrow_ptr[i + 1] += ptrow_ptr[i];
          }

          for (Index i(0); i < txused_elements; ++i)
          {
            const Index l(ptxcol_ind[i]);
            const Index j(ptrow_ptr[l]);
            ptval[j] = ptxval[i];
            ptcol_ind[j] = ptxrow_ind[i];
            ptrow_ind[j] = ptxcol_ind[i];
            ++ptrow_ptr[l];
          }

          SparseMatrixCOO<Mem::Main, DT_, IT_> tx_t(txcolumns, txrows, trow_ind, tcol_ind, tval);

          SparseMatrixCOO<Mem_, DT_, IT_> x_t;
          x_t.convert(tx_t);
          this->assign(x_t);
        }

        /**
         * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_i\f$
         *
         * \param[in] x The matrix whose rows are to be scaled.
         * \param[in] s The vector to the scale the rows by.
         */
        template<typename Algo_>
        void scale_rows(const SparseMatrixCOO & x, const DenseVector<Mem_,DT_,IT_> & s)
        {
          if (x.rows() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
          if (x.columns() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
          if (s.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

          Arch::ScaleRows<Mem_, Algo_>::coo(this->val(), x.val(), this->row_indices(), this->column_indices(),
          s.elements(), this->rows(), this->columns(), this->used_elements());
        }

        /**
         * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_j\f$
         *
         * \param[in] x The matrix whose columns are to be scaled.
         * \param[in] s The vector to the scale the columns by.
         */
        template<typename Algo_>
        void scale_cols(const SparseMatrixCOO & x, const DenseVector<Mem_,DT_,IT_> & s)
        {
          if (x.rows() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
          if (x.columns() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
          if (x.used_elements() != this->used_elements())
            throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
          if (s.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

          Arch::ScaleCols<Mem_, Algo_>::coo(this->val(), x.val(), this->row_indices(), this->column_indices(),
          s.elements(), this->rows(), this->columns(), this->used_elements());
        }

        /**
         * \brief Calculate \f$r \leftarrow this\cdot x \f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
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

          Arch::ProductMatVec<Mem_, Algo_>::coo(r.elements(), this->val(), this->row_indices(),
          this->column_indices(), x.elements(), this->rows(), this->used_elements());
        }

        /**
         * \brief Calculate \f$ r \leftarrow this\cdot x \f$, global version.
         *
         * \param[out] r The vector that recieves the result.
         * \param[in] x The vector to be multiplied by this matrix.
         * \param[in] gate The gate base pointer
         */
        template<typename Algo_>
        void apply(DenseVector<Mem_,DT_, IT_>& r,
                   const DenseVector<Mem_, DT_, IT_>& x,
                   Arch::ProductMat0Vec1GatewayBase<Mem_, Algo_, DenseVector<Mem_, DT_, IT_>, SparseMatrixCOO<Mem_, DT_, IT_> >* gate
                  )
        {
          if (r.size() != this->rows())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
          if (x.size() != this->columns())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

          gate->value(r, *this, x);
        }

        /**
         * \brief Calculate \f$r \leftarrow y + \alpha this\cdot x \f$
         *
         * \tparam Algo_ The \ref FEAST::Algo "algorithm" to be used.
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
            this->row_indices(), this->column_indices(), x.elements(), this->rows(), this->used_elements());
          }
          // r <- y
          else if(Math::abs(alpha) < Math::eps<DT_>())
            r.copy(y);
          // r <- y + alpha*x
          else
          {
            Arch::Axpy<Mem_, Algo_>::coo(r.elements(), alpha, x.elements(), y.elements(),
            this->val(), this->row_indices(), this->column_indices(), this->rows(), this->used_elements());
          }
        }
        ///@}

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

        /// Returns the number of NNZ-elements of the selected row
        Index get_length_of_line(const Index row) const
        {
          const IT_ * prow(this->row_indices());
          const Index tused_elements(this->used_elements());

          Index i(0);
          while (prow[i] < row)
          {
            ++i;
          }

          Index j(i);
          while (j < tused_elements && prow[j] == row)
          {
            ++j;
          }
          return (j - i);
        }

        /// \cond internal

        /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
        void set_line(const Index row, DT_ * const pval_set, IT_ * const pcol_set,
        const Index col_start, const Index stride = 1) const
        {
          const IT_ * prow(this->row_indices());
          const IT_ * pcol(this->column_indices());
          const DT_ * pval(this->val());

          const Index tused_elements(this->used_elements());

          Index start(0);
          while (prow[start] < row)
          {
            ++start;
          }

          for (Index i(0); start + i < tused_elements; ++i)
          {
            if (prow[start + i] != row)
            {
              return;
            }
            pval_set[i * stride] = pval[start + i];
            pcol_set[i * stride] = pcol[start + i] + IT_(col_start);
          }
        }
        /// \endcond

        /**
         * \brief SparseMatrixCOO comparison operator
         *
         * \param[in] a A matrix to compare with.
         * \param[in] b A matrix to compare with.
         */
        template <typename Mem2_> friend bool operator== (const SparseMatrixCOO & a, const SparseMatrixCOO<Mem2_, DT_, IT_> & b)
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
            if (Util::MemoryPool<Mem_>::get_element(a.val(), i) != Util::MemoryPool<Mem2_>::get_element(b.val(), i))
              return false;
            if (Util::MemoryPool<Mem_>::get_element(a.row_indices(), i) != Util::MemoryPool<Mem2_>::get_element(b.row_indices(), i))
              return false;
            if (Util::MemoryPool<Mem_>::get_element(a.column_indices(), i) != Util::MemoryPool<Mem2_>::get_element(b.column_indices(), i))
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
        friend std::ostream & operator<< (std::ostream & lhs, const SparseMatrixCOO & b)
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
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_MATRIX_COO_HPP
