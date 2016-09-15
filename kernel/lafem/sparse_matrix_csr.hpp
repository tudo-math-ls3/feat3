#pragma once
#ifndef KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP
#define KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_layout.hpp>
#include <kernel/lafem/arch/scale_row_col.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/lafem/arch/diagonal.hpp>
#include <kernel/lafem/arch/lumping.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>

#include <fstream>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief CSR based sparse matrix.
     *
     * \tparam Mem_ The \ref FEAT::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a sparse matrix, that stores its non zero elements in the compressed sparse row format.\n\n
     * Data survey: \n
     * _elements[0]: raw non zero number values \n
     * _indices[0]: column index per non zero element \n
     * _indices[1]: row start index (including matrix end index)\n
     *
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: row count \n
     * _scalar_index[2]: column count \n
     * _scalar_index[3]: non zero element count (used elements) \n
     * _scalar_dt[0]: zero element
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseMatrixCSR : public Container<Mem_, DT_, IT_>
    {
    public:
      /**
       * \brief Scatter-Axpy operation for SparseMatrixCSR
       *
       * \author Peter Zajac
       */
      class ScatterAxpy
      {
      public:
        typedef LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_> MatrixType;
        typedef Mem::Main MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
#ifdef DEBUG
        const IT_ _deadcode;
#endif
        Index _num_rows;
        Index _num_cols;
        IT_* _row_ptr;
        IT_* _col_idx;
        IT_* _col_ptr;
        DT_ *_data;

      public:
        explicit ScatterAxpy(MatrixType& matrix) :
#ifdef DEBUG
          _deadcode(~IT_(0)),
#endif
          _num_rows(matrix.rows()),
          _num_cols(matrix.columns()),
          _row_ptr(matrix.row_ptr()),
          _col_idx(matrix.col_ind()),
          _col_ptr(nullptr),
          _data(matrix.val())
        {
          // allocate column-pointer array
          _col_ptr = new IT_[matrix.columns()];
#ifdef DEBUG
          for(Index i(0); i < _num_cols; ++i)
          {
            _col_ptr[i] = _deadcode;
          }
#endif
        }

        virtual ~ScatterAxpy()
        {
          if(_col_ptr != nullptr)
          {
            delete [] _col_ptr;
          }
        }

        template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
        void operator()(const LocalMatrix_& loc_mat, const RowMapping_& row_map,
                        const ColMapping_& col_map, DT_ alpha = DT_(1))
        {
          // loop over all local row entries
          for(int i(0); i < row_map.get_num_local_dofs(); ++i)
          {
            // loop over all row entry contributations
            for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
            {
              // fetch row entry weight and pre-multiply by alpha
              DT_ iw = alpha * DT_(row_map.get_weight(i, ic));

              // fetch row index
              Index ix = row_map.get_index(i, ic);

              // build column pointer for this row entry contribution
              for(IT_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
              {
                _col_ptr[_col_idx[k]] = k;
              }

              // loop over all local column entries
              for(int j(0); j < col_map.get_num_local_dofs(); ++j)
              {
                // loop over all column entry contributions
                for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
                {
                  // fetch trial function dof weight
                  DT_ jw = DT_(col_map.get_weight(j, jc));

                  // fetch column index
                  Index jx = col_map.get_index(j, jc);

                  // ensure that the column pointer is valid for this index
                  ASSERTM(_col_ptr[jx] != _deadcode, "invalid column index");

                  // incorporate data into global matrix
                  _data[_col_ptr[jx]] += (iw * jw) * loc_mat[i][j];

                  // continue with next column contribution
                }
                // continue with next column entry
              }

#ifdef DEBUG
              // reformat column-pointer array
              for(IT_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
              {
                _col_ptr[_col_idx[k]] = _deadcode;
              }
#endif
              // continue with next row contribution
            }
            // continue with next row entry
          }
        }
      }; // class ScatterAxpy

      /**
       * \brief Gather-Axpy operation for SparseMatrixCSR
       *
       * \author Peter Zajac
       */
      class GatherAxpy
      {
      public:
        typedef LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_> MatrixType;
        typedef Mem::Main MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
#ifdef DEBUG
        const IT_ _deadcode;
#endif
        Index _num_rows;
        Index _num_cols;
        const IT_* _row_ptr;
        const IT_* _col_idx;
        IT_* _col_ptr;
        const DT_ *_data;

      public:
        explicit GatherAxpy(const MatrixType& matrix) :
#ifdef DEBUG
          _deadcode(~IT_(0)),
#endif
          _num_rows(matrix.rows()),
          _num_cols(matrix.columns()),
          _row_ptr(matrix.row_ptr()),
          _col_idx(matrix.col_ind()),
          _col_ptr(nullptr),
          _data(matrix.val())
        {
          // allocate column-pointer array
          _col_ptr = new IT_[matrix.columns()];
#ifdef DEBUG
          for(Index i(0); i < _num_cols; ++i)
          {
            _col_ptr[i] = _deadcode;
          }
#endif
        }

        virtual ~GatherAxpy()
        {
          if(_col_ptr != nullptr)
          {
            delete [] _col_ptr;
          }
        }

        template<typename LocalMatrix_, typename RowMapping_, typename ColMapping_>
        void operator()(LocalMatrix_& loc_mat, const RowMapping_& row_map,
                        const ColMapping_& col_map, DT_ alpha = DT_(1))
        {
          // loop over all local row entries
          for(int i(0); i < row_map.get_num_local_dofs(); ++i)
          {
            // loop over all row entry contributations
            for(int ic(0); ic < row_map.get_num_contribs(i); ++ic)
            {
              // fetch row index
              Index ix = row_map.get_index(i, ic);

              // build column pointer for this row entry contribution
              for(IT_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
              {
                _col_ptr[_col_idx[k]] = k;
              }

              // loop over all local column entries
              for(int j(0); j < col_map.get_num_local_dofs(); ++j)
              {
                // clear  accumulation entry
                DT_ dx(DT_(0));

                // loop over all column entry contributions
                for(int jc(0); jc < col_map.get_num_contribs(j); ++jc)
                {
                  // fetch column index
                  Index jx = col_map.get_index(j, jc);

                  // ensure that the column pointer is valid for this index
                  ASSERTM(_col_ptr[jx] != _deadcode, "invalid column index");

                  // update accumulator
                  dx += DT_(col_map.get_weight(j, jc)) * _data[_col_ptr[jx]];

                  // continue with next column contribution
                }

                // update local matrix data
                loc_mat[i][j] += (alpha * DT_(row_map.get_weight(i, ic))) * dx;

                // continue with next column entry
              }

#ifdef DEBUG
              // reformat column-pointer array
              for(IT_ k(_row_ptr[ix]); k < _row_ptr[ix + 1]; ++k)
              {
                _col_ptr[_col_idx[k]] = _deadcode;
              }
#endif

              // continue with next row contribution
            }
            // continue with next row entry
          }
        }
      }; // class GatherAxpy

    private:
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
      static constexpr SparseLayoutId layout_id = SparseLayoutId::lt_csr;
      /// ImageIterator typedef for Adjactor interface implementation
      typedef const IT_* ImageIterator;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = class SparseMatrixCSR<Mem2_, DT2_, IT2_>;

      /// this typedef lets you create a matrix container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;


      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional matrix.
       */
      explicit SparseMatrixCSR() :
        Container<Mem_, DT_, IT_> (0)
      {
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       *
       * Creates an empty matrix.
       * Because SparseMatrixCSR is a read-only container, it stays empty.
       *
       * \note This matrix does not allocate any memory
       */
      explicit SparseMatrixCSR(Index rows_in, Index columns_in) :
        Container<Mem_, DT_, IT_> (rows_in * columns_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] used_elements_in The amount of non zero elements of the created matrix.
       *
       * Creates an empty (but allocated) matrix.
       *
       * \note The allocated memory will not be initialised.
       */
      explicit SparseMatrixCSR(Index rows_in, Index columns_in, Index used_elements_in) :
        Container<Mem_, DT_, IT_> (rows_in * columns_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(used_elements_in);
        this->_scalar_dt.push_back(DT_(0));

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());

        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_rows() + 1));
        this->_indices_size.push_back(_rows() + 1);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
      }

      /**
       * \brief Constructor
       *
       * \param[in] layout_in The layout to be used.
       *
       * Creates an empty matrix with given layout.
       */
      explicit SparseMatrixCSR(const SparseLayout<Mem_, IT_, layout_id> & layout_in) :
        Container<Mem_, DT_, IT_> (layout_in._scalar_index.at(0))
      {
        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          MemoryPool<Mem_>::increase_memory(i);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
      }

      /**
       * \brief Constructor
       *
       * \param[in] other The source matrix.
       *
       * Creates a CSR matrix based on the source matrix.
       */
      template <typename MT_>
      explicit SparseMatrixCSR(const MT_ & other) :
        Container<Mem_, DT_, IT_>(other.size())
      {
        convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph The.graph to create the matrix from
       *
       * Creates a CSR matrix based on a given adjacency graph, representing the sparsity pattern.
       */
      explicit SparseMatrixCSR(const Adjacency::Graph & graph) :
        Container<Mem_, DT_, IT_>(0)
      {
        Index num_rows = graph.get_num_nodes_domain();
        Index num_cols = graph.get_num_nodes_image();
        Index num_nnze = graph.get_num_indices();

        // create temporary vectors
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vrow_ptr(num_rows+1);
        LAFEM::DenseVector<Mem::Main, IT_, IT_> vcol_idx(num_nnze);
        LAFEM::DenseVector<Mem::Main, DT_, IT_> vdata(num_nnze, DT_(0));

        const Index * dom_ptr(graph.get_domain_ptr());
        const Index * img_idx(graph.get_image_idx());
        IT_ * prow_ptr(vrow_ptr.elements());
        IT_ * pcol_idx(vcol_idx.elements());

        // build row-end
        prow_ptr[0] = IT_(dom_ptr[0]);
        for(Index i(0); i < num_rows; ++i)
          prow_ptr[i+1] = IT_(dom_ptr[i+1]);

        // build col-idx
        for(Index i(0); i < num_nnze; ++i)
          pcol_idx[i] = IT_(img_idx[i]);

        // build the matrix
        this->assign(SparseMatrixCSR<Mem::Main, DT_, IT_>(num_rows, num_cols, vcol_idx, vdata, vrow_ptr));
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a CSR matrix based on the source file.
       */
      explicit SparseMatrixCSR(FileMode mode, String filename) :
        Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a CSR matrix based on the source filestream.
       */
      explicit SparseMatrixCSR(FileMode mode, std::istream& file) :
        Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, file);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rows_in The row count of the created matrix.
       * \param[in] columns_in The column count of the created matrix.
       * \param[in] col_ind_in Vector with column indices.
       * \param[in] val_in Vector with non zero elements.
       * \param[in] row_ptr_in Vector with start indices of all rows into the val/col_ind arrays.
       * Note that this vector must also contain the end index of the last row and thus has a size of row_count + 1.
       *
       * Creates a matrix with given dimensions and content.
       */
      explicit SparseMatrixCSR(const Index rows_in, const Index columns_in,
                               DenseVector<Mem_, IT_, IT_> & col_ind_in, DenseVector<Mem_, DT_, IT_> & val_in, DenseVector<Mem_, IT_, IT_> & row_ptr_in) :
        Container<Mem_, DT_, IT_>(rows_in * columns_in)
      {
        this->_scalar_index.push_back(rows_in);
        this->_scalar_index.push_back(columns_in);
        this->_scalar_index.push_back(val_in.size());
        this->_scalar_dt.push_back(DT_(0));

        this->_elements.push_back(val_in.elements());
        this->_elements_size.push_back(val_in.size());
        this->_indices.push_back(col_ind_in.elements());
        this->_indices_size.push_back(col_ind_in.size());
        this->_indices.push_back(row_ptr_in.elements());
        this->_indices_size.push_back(row_ptr_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a matrix from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseMatrixCSR(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        deserialise<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source matrix.
       *
       * Moves a given matrix to this matrix.
       */
      SparseMatrixCSR(SparseMatrixCSR && other) :
        Container<Mem_, DT_, IT_>(std::forward<SparseMatrixCSR>(other))
      {
      }

      /**
       * \brief Move operator
       *
       * \param[in] other The source matrix.
       *
       * Moves another matrix to the target matrix.
       */
      SparseMatrixCSR & operator= (SparseMatrixCSR && other)
      {
        this->move(std::forward<SparseMatrixCSR>(other));

        return *this;
      }

      InsertWeakClone( SparseMatrixCSR )

      /** \brief Shallow copy operation
       *
       * Create a shallow copy of itself.
       *
       */
      SparseMatrixCSR shared() const
      {
        SparseMatrixCSR r;
        r.assign(*this);
        return r;
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
        this->assign(other);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const SparseMatrixCOO<Mem2_, DT2_, IT2_> & other)
      {
        this->clear();

        this->_scalar_index.push_back(other.size());
        this->_scalar_index.push_back(other.rows());
        this->_scalar_index.push_back(other.columns());
        this->_scalar_index.push_back(other.used_elements());

        SparseMatrixCOO<Mem::Main, DT_, IT_> cother;
        cother.convert(other);
        this->_scalar_dt.push_back(cother.zero_element());

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_rows() + 1));
        this->_indices_size.push_back(_rows() + 1);

        DT_ * tval(nullptr);
        IT_ * tcol_ind(nullptr);
        IT_ * trow_ptr(nullptr);
        if (std::is_same<Mem_, Mem::Main>::value)
        {
          tval = this->_elements.at(0);
          tcol_ind = this->_indices.at(0);
          trow_ptr = this->_indices.at(1);
        }
        else
        {
          tval = new DT_[other.used_elements()];
          tcol_ind = new IT_[other.used_elements()];
          trow_ptr = new IT_[other.rows() + 1];
        }

        IT_ ait(0);
        IT_ current_row(IT_(0));
        trow_ptr[current_row] = IT_(0);
        for (Index it(0) ; it < cother.used_elements() ; ++it)
        {
          IT_ row(cother.row_indices()[it]);
          IT_ column(cother.column_indices()[it]);

          if (current_row < row)
          {
            for (IT_ i(current_row + IT_(1)) ; i < row ; ++i)
            {
              trow_ptr[i] = ait;
            }
            current_row = row;
            trow_ptr[current_row] = ait;
          }
          tval[ait] = cother.val()[it];
          tcol_ind[ait] = column;
          ++ait;
        }
        for (IT_ i(current_row + IT_(1)) ; i < _rows() ; ++i)
        {
          trow_ptr[i] = ait;
        }
        trow_ptr[_rows()] = ait;

        if (! std::is_same<Mem_, Mem::Main>::value)
        {
          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, _rows() + 1);
          delete[] tval;
          delete[] tcol_ind;
          delete[] trow_ptr;
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
        this->clear();

        this->_scalar_index.push_back(other.size());
        this->_scalar_index.push_back(other.rows());
        this->_scalar_index.push_back(other.columns());
        this->_scalar_index.push_back(other.used_elements());

        SparseMatrixELL<Mem::Main, DT_, IT_> cother;
        cother.convert(other);
        this->_scalar_dt.push_back(cother.zero_element());

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_rows() + 1));
        this->_indices_size.push_back(_rows() + 1);

        DT_ * tval(nullptr);
        IT_ * tcol_ind(nullptr);
        IT_ * trow_ptr(nullptr);
        if (std::is_same<Mem_, Mem::Main>::value)
        {
          tval = this->_elements.at(0);
          tcol_ind = this->_indices.at(0);
          trow_ptr = this->_indices.at(1);
        }
        else
        {
          tval = new DT_[used_elements()];
          tcol_ind = new IT_[used_elements()];
          trow_ptr = new IT_[rows() + 1];
        }

        trow_ptr[0] = IT_(0);

        const Index cC(cother.C());
        const DT_ * cval(cother.val());
        const IT_ * ccol(cother.col_ind());
        const IT_ * crl(cother.rl());
        const IT_ * ccs(cother.cs());

        IT_ ue(0);
        for (Index row(0) ; row < cother.rows() ; ++row)
        {
          for (Index i(0); i < crl[row]; ++i)
          {
            tval[ue] = cval[ccs[row/cC] + row%cC + i*cC];
            tcol_ind[ue] = ccol[ccs[row/cC] + row%cC + i*cC];
            ++ue;
          }
          trow_ptr[row + 1] = ue;
        }

        if (! std::is_same<Mem_, Mem::Main>::value)
        {
          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, _rows() + 1);
          delete[] tval;
          delete[] tcol_ind;
          delete[] trow_ptr;
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
        this->clear();

        this->_scalar_index.push_back(other.size());
        this->_scalar_index.push_back(other.rows());
        this->_scalar_index.push_back(other.columns());
        this->_scalar_index.push_back(other.used_elements());

        SparseMatrixBanded<Mem::Main, DT_, IT_> cother;
        cother.convert(other);
        this->_scalar_dt.push_back(cother.zero_element());

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_rows() + 1));
        this->_indices_size.push_back(_rows() + 1);

        DT_ * tval(nullptr);
        IT_ * tcol_ind(nullptr);
        IT_ * trow_ptr(nullptr);
        if (std::is_same<Mem_, Mem::Main>::value)
        {
          tval = this->_elements.at(0);
          tcol_ind = this->_indices.at(0);
          trow_ptr = this->_indices.at(1);
        }
        else
        {
          tval = new DT_[other.used_elements()];
          tcol_ind = new IT_[other.used_elements()];
          trow_ptr = new IT_[other.rows() + 1];
        }

        trow_ptr[0] = 0;

        const DT_ * cval(cother.val());
        const IT_ * coffsets(cother.offsets());
        const Index cnum_of_offsets(cother.num_of_offsets());
        const Index crows(cother.rows());

        // Search first offset of the upper triangular matrix
        Index k(0);
        while (k < cnum_of_offsets && coffsets[k] + 1 < crows)
        {
          ++k;
        }

        IT_ ue(0);
        // iteration over all offsets of the lower triangular matrix
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
                ++ue;
              }
              trow_ptr[l + 1] = ue;
            }
          }
        }

        if (! std::is_same<Mem_, Mem::Main>::value)
        {
          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, _rows() + 1);
          delete[] tval;
          delete[] tcol_ind;
          delete[] trow_ptr;
        }
      }

      /**
       * \brief Convertion method
       *
       * \param[in] other The source Matrix.
       *
       * Use source matrix content as content of current matrix
       */
      template <typename Mem2_, typename DT2_, typename IT2_, int BlockHeight_, int BlockWidth_>
      void convert(const SparseMatrixBCSR<Mem2_, DT2_, IT2_, BlockHeight_, BlockWidth_> & other)
      {
        this->clear();

        this->_scalar_index.push_back(other.template rows<Perspective::pod>() * other.template columns<Perspective::pod>());
        this->_scalar_index.push_back(other.template rows<Perspective::pod>());
        this->_scalar_index.push_back(other.template columns<Perspective::pod>());
        this->_scalar_index.push_back(other.template used_elements<Perspective::pod>());

        SparseMatrixBCSR<Mem::Main, DT_, IT_, BlockHeight_, BlockWidth_> cother;
        cother.convert(other);
        this->_scalar_dt.push_back(cother.zero_element());

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_rows() + 1));
        this->_indices_size.push_back(_rows() + 1);

        DT_ * tval(nullptr);
        IT_ * tcol_ind(nullptr);
        IT_ * trow_ptr(nullptr);
        if (std::is_same<Mem_, Mem::Main>::value)
        {
          tval = this->_elements.at(0);
          tcol_ind = this->_indices.at(0);
          trow_ptr = this->_indices.at(1);
        }
        else
        {
          tval = new DT_[used_elements()];
          tcol_ind = new IT_[used_elements()];
          trow_ptr = new IT_[rows() + 1];
        }

        Index ait(0);
        trow_ptr[0] = IT_(0);
        Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *mval(reinterpret_cast<Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *>(cother.val()));
        for (Index orow(0) ; orow < cother.rows() ; ++orow)
        {
          for (int row(0) ; row < BlockHeight_ ; ++row)
          {
            for(Index ocol(cother.row_ptr()[orow]) ; ocol < cother.row_ptr()[orow + 1] ; ++ocol)
            {
              Tiny::Matrix<DT_, BlockHeight_, BlockWidth_> *tbm(mval + ocol);
              for (int col(0) ; col < BlockWidth_ ; ++col)
              {
                tval[ait] = (*tbm)(row,col);
                tcol_ind[ait] = cother.col_ind()[ocol] * (IT_)BlockWidth_ + (IT_)col;
                ++ait;
              }
            }
            trow_ptr[orow * Index(BlockHeight_) + Index(row) + 1] = (IT_)ait;
          }
        }

        if (! std::is_same<Mem_, Mem::Main>::value)
        {
          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
          MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, _rows() + 1);
          delete[] tval;
          delete[] tcol_ind;
          delete[] trow_ptr;
        }
      }

      /**
       * \brief Conversion method
       *
       * \param[in] a The input matrix.
       *
       * Converts any matrix to SparseMatrixCSR-format
       */
      template <typename MT_>
      void convert(const MT_ & a)
      {
        typename MT_::template ContainerType<Mem::Main, DT_, IT_> ta;
        ta.convert(a);

        const Index arows(ta.template rows<Perspective::pod>());
        const Index acolumns(ta.template columns<Perspective::pod>());
        const Index aused_elements(ta.template used_elements<Perspective::pod>());

        DenseVector<Mem::Main, DT_, IT_> tval(aused_elements);
        DenseVector<Mem::Main, IT_, IT_> tcol_ind(aused_elements);
        DenseVector<Mem::Main, IT_, IT_> trow_ptr(arows + 1);

        DT_ * pval(tval.elements());
        IT_ * pcol_ind(tcol_ind.elements());
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

        SparseMatrixCSR<Mem::Main, DT_, IT_> ta_csr(arows, acolumns, tcol_ind, tval, trow_ptr);
        SparseMatrixCSR<Mem_, DT_, IT_> a_csr;
        a_csr.convert(ta_csr);

        this->assign(a_csr);
      }

      /**
       * \brief Assignment operator
       *
       * \param[in] layout_in A sparse matrix layout.
       *
       * Assigns a new matrix layout, discarding all old data
       */
      SparseMatrixCSR & operator= (const SparseLayout<Mem_, IT_, layout_id> & layout_in)
      {
        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_indices.at(i));

        this->_elements.clear();
        this->_indices.clear();
        this->_elements_size.clear();
        this->_indices_size.clear();
        this->_scalar_index.clear();
        this->_scalar_dt.clear();

        this->_indices.assign(layout_in._indices.begin(), layout_in._indices.end());
        this->_indices_size.assign(layout_in._indices_size.begin(), layout_in._indices_size.end());
        this->_scalar_index.assign(layout_in._scalar_index.begin(), layout_in._scalar_index.end());
        this->_scalar_dt.push_back(DT_(0));

        for (auto i : this->_indices)
          MemoryPool<Mem_>::increase_memory(i);

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());

        return *this;
      }

      /**
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialise(std::vector<char> input)
      {
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_csr, input);
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialise for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialise()
      {
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_csr);
      }

      /**
       * \brief Read in matrix from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, String filename)
      {
        switch(mode)
        {
        case FileMode::fm_mtx:
          read_from_mtx(filename);
          break;
        case FileMode::fm_csr:
          read_from_csr(filename);
          break;
        case FileMode::fm_binary:
          read_from_csr(filename);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Read in matrix from stream.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be read in.
       */
      void read_from(FileMode mode, std::istream& file)
      {
        switch(mode)
        {
        case FileMode::fm_mtx:
          read_from_mtx(file);
          break;
        case FileMode::fm_csr:
          read_from_csr(file);
          break;
        case FileMode::fm_binary:
          read_from_csr(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Read in matrix from MatrixMarket mtx file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_mtx(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        read_from_mtx(file);
        file.close();
      }

      /**
       * \brief Read in matrix from MatrixMarket mtx stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_mtx(std::istream& file)
      {
        this->clear();
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_dt.push_back(DT_(0));

        std::map<IT_, std::map<IT_, DT_> > entries; // map<row, map<column, value> >

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

          entries[IT_(row)].insert(std::pair<IT_, DT_>(col, tval));
          ++ue;
          if (symmetric == true && row != col)
          {
            entries[IT_(col)].insert(std::pair<IT_, DT_>(row, tval));
            ++ue;
          }
        }
        _size() = this->rows() * this->columns();
        _used_elements() = ue;

        DT_ * tval = new DT_[ue];
        IT_ * tcol_ind = new IT_[ue];
        IT_ * trow_ptr = new IT_[rows() + 1];

        IT_ idx(0);
        Index row_idx(0);
        for (auto row : entries)
        {
          trow_ptr[row_idx] = idx;
          for (auto col : row.second )
          {
            tcol_ind[idx] = col.first;
            tval[idx] = col.second;
            ++idx;
          }
          row.second.clear();
          ++row_idx;
        }
        trow_ptr[rows()] = IT_(ue);
        entries.clear();

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(_used_elements()));
        this->_elements_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(_used_elements()));
        this->_indices_size.push_back(_used_elements());
        this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(rows() + 1));
        this->_indices_size.push_back(rows() + 1);

        MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), tval, _used_elements());
        MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(0), tcol_ind, _used_elements());
        MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(1), trow_ptr, rows() + 1);

        delete[] tval;
        delete[] tcol_ind;
        delete[] trow_ptr;
      }

      /**
       * \brief Read in matrix from binary file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_csr(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        read_from_csr(file);
        file.close();
      }

      /**
       * \brief Read in matrix from binary stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_csr(std::istream& file)
      {
        this->template _deserialise<double, uint64_t>(FileMode::fm_csr, file);
      }


      /**
       * \brief Write out matrix to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        switch(mode)
        {
        case FileMode::fm_csr:
          write_out_csr(filename);
          break;
        case FileMode::fm_binary:
          write_out_csr(filename);
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
        switch(mode)
        {
        case FileMode::fm_csr:
          write_out_csr(file);
          break;
        case FileMode::fm_binary:
          write_out_csr(file);
          break;
        case FileMode::fm_mtx:
          write_out_mtx(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Write out matrix to csr binary file.
       *
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out_csr(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        write_out_csr(file);
        file.close();
      }

      /**
       * \brief Write out matrix to csr binary file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_csr(std::ostream& file) const
      {
        if (! std::is_same<DT_, double>::value)
          std::cout<<"Warning: You are writing out a csr matrix that is not double precision!"<<std::endl;

        this->template _serialise<double, uint64_t>(FileMode::fm_csr, file);
      }

      /**
       * \brief Write out matrix to MatrixMarktet mtx file.
       *
       * \param[in] filename The file where the matrix shall be stored.
       * \param[in] symmetric Should we store only the lower half of the matrix in symmetric format?
       *
       * \warning This routine does no check on symmetric properties of the source matrix!
       */
      void write_out_mtx(String filename, bool symmetric = false) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Matrix file " + filename);
        write_out_mtx(file, symmetric);
        file.close();
      }

      /**
       * \brief Write out matrix to MatrixMarktet mtx file.
       *
       * \param[in] file The stream that shall be written to.
       * \param[in] symmetric Should we store only the LD part of the matrix in symmetric format?
       *
       * \warning This routine does no check on symmetric properties of the source matrix!
       */
      void write_out_mtx(std::ostream& file, bool symmetric = false) const
      {
        SparseMatrixCSR<Mem::Main, DT_, IT_> temp;
        temp.convert(*this);

        if (symmetric)
        {
          file << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
          std::vector<IT_> rowv;
          std::vector<IT_> colv;
          std::vector<DT_> valv;
          for (Index row(0) ; row < rows() ; ++row)
          {
            const IT_ end(temp.row_ptr()[row + 1]);
            for (IT_ i(temp.row_ptr()[row]) ; i < end ; ++i)
            {
              const IT_ col(temp.col_ind()[i]);
              if (row >= col)
              {
                rowv.push_back(IT_(row + 1));
                colv.push_back(col + 1);
                valv.push_back(temp.val()[i]);
              }
            }
          }
          file << temp.rows() << " " << temp.columns() << " " << valv.size() << std::endl;
          for (Index i(0) ; i < valv.size() ; ++i)
          {
            file << rowv.at(i) << " " << colv.at(i) << " " << std::scientific << valv.at(i) << std::endl;
          }
        }
        else
        {
          file << "%%MatrixMarket matrix coordinate real general" << std::endl;
          file << temp.rows() << " " << temp.columns() << " " << temp.used_elements() << std::endl;

          for (Index row(0) ; row < rows() ; ++row)
          {
            const IT_ end(temp.row_ptr()[row + 1]);
            for (IT_ i(temp.row_ptr()[row]) ; i < end ; ++i)
            {
              file << row + 1 << " " << temp.col_ind()[i] + 1 << " " << std::scientific << temp.val()[i] << std::endl;
            }
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
        ASSERT(row < rows());
        ASSERT(col < columns());

        for (Index i(Index(MemoryPool<Mem_>::get_element(this->_indices.at(1), row))) ; i < Index(MemoryPool<Mem_>::get_element(this->_indices.at(1), row + 1)) ; ++i)
        {
          if (Index(MemoryPool<Mem_>::get_element(this->_indices.at(0), i)) == col)
            return MemoryPool<Mem_>::get_element(this->_elements.at(0), i);
          if (Index(MemoryPool<Mem_>::get_element(this->_indices.at(0), i)) > col)
            return zero_element();
        }
        return zero_element();
      }

      /**
       * \brief Retrieve convenient sparse matrix layout object.
       *
       * \return An object containing the sparse matrix layout.
       */
      SparseLayout<Mem_, IT_, layout_id> layout() const
      {
        return SparseLayout<Mem_, IT_, layout_id>(this->_indices, this->_indices_size, this->_scalar_index);
      }

      /**
       * \brief Retrieve matrix row count.
       *
       * \returns Matrix row count.
       */
      template <Perspective = Perspective::native>
      const Index & rows() const
      {
        return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve matrix column count.
       *
       * \returns Matrix column count.
       */
      template <Perspective = Perspective::native>
      const Index & columns() const
      {
        return this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      template <Perspective = Perspective::native>
      Index used_elements() const
      {
        return this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve column indices array.
       *
       * \returns Column indices array.
       */
      IT_ * col_ind()
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(0);
      }

      IT_ const * col_ind() const
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(0);
      }

      /**
       * \brief Retrieve non zero element array.
       *
       * \returns Non zero element array.
       */
      DT_ * val()
      {
        if (this->_elements.size() == 0)
          return nullptr;

        return this->_elements.at(0);
      }

      DT_ const * val() const
      {
        if (this->_elements.size() == 0)
          return nullptr;

        return this->_elements.at(0);
      }

      /**
       * \brief Retrieve row start index array.
       *
       * \returns Row start index array.
       */
      IT_ * row_ptr()
      {
        if (this->_indices.size() == 0)
          return nullptr;

        return this->_indices.at(1);
      }

      IT_ const * row_ptr() const
      {
        if (this->_indices.size() == 0)
          return nullptr;

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
       * \brief Retrieve maximum bandwidth among all rows.
       *
       * \param[out] bandw The maximum bandwidth.
       * \param[out] bandw_i The row, where the bandwidth is maximal.
       */
      void bandwidth_row(Index & bandw, Index & bandw_i)
      {
        SparseMatrixCSR<Mem::Main, DT_, IT_> tm;
        tm.convert(*this);
        bandw = 0;
        bandw_i = 0;
        for (Index row(0) ; row < rows() ; ++row)
        {
          if (tm.row_ptr()[row+1] == tm.row_ptr()[row])
            continue;

          Index temp = tm.col_ind()[tm.row_ptr()[row+1]-1] - tm.col_ind()[tm.row_ptr()[row]] + 1;
          if(temp > bandw)
          {
            bandw = temp;
            bandw_i = row;
          }
        }
      }

      /**
       * \brief Retrieve maximum bandwidth among all columns.
       *
       * \param[out] bandw The maximum bandwidth.
       * \param[out] bandw_i The column, where the bandwidth is maximal.
       */
      void bandwidth_column(Index & bandw, Index & bandw_i)
      {
        SparseMatrixCSR<Mem::Main, DT_, IT_> tm;
        tm.convert(*this);
        tm.transpose(tm);
        tm.bandwidth_row(bandw, bandw_i);
      }

      /**
       * \brief Retrieve maximum radius among all rows.
       *
       * The radius is defined as the maximal distance to the main diagonal
       *
       * \param[out] radius The maximum radius.
       * \param[out] radius_i The row, where the radius is maximal.
       */
      void radius_row(Index & radius, Index & radius_i)
      {
        SparseMatrixCSR<Mem::Main, DT_, IT_> tm;
        tm.convert(*this);
        radius = 0;
        radius_i = 0;

        for (Index row(0) ; row < rows() ; ++row)
        {
          if (tm.row_ptr()[row+1] == tm.row_ptr()[row])
            continue;

          if (tm.row_ptr()[row+1] > 0)
          {
            Index temp1 = tm.col_ind()[tm.row_ptr()[row+1]-1];
            if(Math::max(temp1,row) - Math::min(temp1, row) > radius)
            {
              radius = Math::max(temp1,row) - Math::min(temp1, row);
              radius_i = row;
            }
          }
          Index temp2 = tm.col_ind()[tm.row_ptr()[row]];
          if(Math::max(temp2,row) - Math::min(temp2, row) > radius)
          {
            radius = Math::max(temp2,row) - Math::min(temp2, row);
            radius_i = row;
          }
        }
      }

      /**
       * \brief Retrieve maximum radius among all column.
       *
       * The radius is defined as the maximal distance to the main diagonal
       *
       * \param[out] radius The maximum radius.
       * \param[out] radius_i The column, where the radius is maximal.
       */
      void radius_column(Index & radius, Index & radius_i)
      {
        SparseMatrixCSR<Mem::Main, DT_, IT_> tm;
        tm.convert(*this);
        tm.transpose(tm);
        tm.radius_row(radius, radius_i);
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SparseMatrixCSR";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const SparseMatrixCSR & x, bool full = false)
      {
        this->_copy_content(x, full);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The Matrix to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      template <typename Mem2_>
      void copy(const SparseMatrixCSR<Mem2_, DT_, IT_> & x, bool full = false)
      {
        this->_copy_content(x, full);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow y + \alpha~ x\f$
       *
       * \param[in] x The first summand matrix to be scaled.
       * \param[in] y The second summand matrix
       * \param[in] alpha A scalar to multiply x with.
       *
       * \warning All three matrices must have the same non zero layout. This operation assumes this silently and does not check this on its own!
       */
      void axpy(
                const SparseMatrixCSR & x,
                const SparseMatrixCSR & y,
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

        TimeStamp ts_start;

        // check for special cases
        // r <- x + y
        if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements());
          Arch::Sum<Mem_>::value(this->val(), x.val(), y.val(), this->used_elements());
        }
        // r <- y - x
        else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements());
          Arch::Difference<Mem_>::value(this->val(), y.val(), x.val(), this->used_elements());
        }
        // r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          this->copy(y);
        // r <- y + alpha*x
        else
        {
          Statistics::add_flops(this->used_elements() * 2);
          Arch::Axpy<Mem_>::dv(this->val(), alpha, x.val(), y.val(), this->used_elements());
        }

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x \f$
       *
       * \param[in] x The matrix to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const SparseMatrixCSR & x, const DT_ alpha)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements());
        Arch::Scale<Mem_>::value(this->val(), x.val(), alpha, this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculates the Frobenius norm of this matrix.
       *
       * \returns The Frobenius norm of this matrix.
       */
      DT_ norm_frobenius() const
      {
        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements() * 2);
        DT_ result = Arch::Norm2<Mem_>::value(this->val(), this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculate \f$this^\top \f$
       *
       * \return The transposed matrix
       */
      SparseMatrixCSR transpose() const
      {
        SparseMatrixCSR x_t;
        x_t.transpose(*this);
        return x_t;
      }

      /**
       * \brief Calculate \f$this \leftarrow x^\top \f$
       *
       * \param[in] x The matrix to be transposed.
       */
      void transpose(const SparseMatrixCSR & x)
      {
        SparseMatrixCSR<Mem::Main, DT_, IT_> tx;
        tx.convert(x);

        const Index txrows(tx.rows());
        const Index txcolumns(tx.columns());
        const Index txused_elements(tx.used_elements());

        const IT_ * ptxcol_ind(tx.col_ind());
        const IT_ * ptxrow_ptr(tx.row_ptr());
        const DT_ * ptxval(tx.val());

        DenseVector<Mem::Main, IT_, IT_> tcol_ind(txused_elements);
        DenseVector<Mem::Main, DT_, IT_> tval(txused_elements);
        DenseVector<Mem::Main, IT_, IT_> trow_ptr(txcolumns + 1, IT_(0));

        IT_ * ptcol_ind(tcol_ind.elements());
        DT_ * ptval(tval.elements());
        IT_ * ptrow_ptr(trow_ptr.elements());

        ptrow_ptr[0] = 0;

        for (Index i(0); i < txused_elements; ++i)
        {
          ++ptrow_ptr[ptxcol_ind[i] + 1];
        }

        for (Index i(1); i < txcolumns - 1; ++i)
        {
          ptrow_ptr[i + 1] += ptrow_ptr[i];
        }

        for (Index i(0); i < txrows; ++i)
        {
          for (IT_ k(ptxrow_ptr[i]); k < ptxrow_ptr[i+1]; ++k)
          {
            const IT_ l(ptxcol_ind[k]);
            const IT_ j(ptrow_ptr[l]);
            ptval[j] = ptxval[k];
            ptcol_ind[j] = IT_(i);
            ++ptrow_ptr[l];
          }
        }

        for (Index i(txcolumns); i > 0; --i)
        {
          ptrow_ptr[i] = ptrow_ptr[i - 1];
        }
        ptrow_ptr[0] = 0;

        SparseMatrixCSR<Mem::Main, DT_, IT_> tx_t(txcolumns, txrows, tcol_ind, tval, trow_ptr);

        SparseMatrixCSR<Mem_, DT_, IT_> x_t;
        x_t.convert(tx_t);
        this->assign(x_t);
      }

      /**
       * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_i\f$
       *
       * \param[in] x The matrix whose rows are to be scaled.
       * \param[in] s The vector to the scale the rows by.
       */
      void scale_rows(const SparseMatrixCSR & x, const DenseVector<Mem_,DT_,IT_> & s)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
        if (s.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements());
        Arch::ScaleRows<Mem_>::csr(this->val(), x.val(), this->col_ind(), this->row_ptr(),
                                          s.elements(), this->rows(), this->columns(), this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ this_{ij} \leftarrow x_{ij}\cdot s_j\f$
       *
       * \param[in] x The matrix whose columns are to be scaled.
       * \param[in] s The vector to the scale the columns by.
       */
      void scale_cols(const SparseMatrixCSR & x, const DenseVector<Mem_,DT_,IT_> & s)
      {
        if (x.rows() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Row count does not match!");
        if (x.columns() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Column count does not match!");
        if (x.used_elements() != this->used_elements())
          throw InternalError(__func__, __FILE__, __LINE__, "Nonzero count does not match!");
        if (s.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->used_elements());
        Arch::ScaleCols<Mem_>::csr(this->val(), x.val(), this->col_ind(), this->row_ptr(),
                                          s.elements(), this->rows(), this->columns(), this->used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(DenseVector<Mem_,DT_, IT_> & r, const DenseVector<Mem_, DT_, IT_> & x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        TimeStamp ts_start;

        if (this->used_elements() == 0)
        {
          r.format();
          return;
        }

        if (r.template elements<Perspective::pod>() == x.template elements<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector x and r must not share the same memory!");

        Statistics::add_flops(this->used_elements() * 2);
        Arch::ProductMatVec<Mem_>::csr(r.elements(), this->val(), this->col_ind(), this->row_ptr(),
                                              x.elements(), this->rows(), columns(), used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The block vector that recieves the result.
       * \param[in] x The block vector to be multiplied by this matrix.
       *
       * \note Every element of each block in the vector x is multiplied with the corresponding single scalar entry in the matrix.
       */
      template<int BlockSize_>
      void apply(DenseVectorBlocked<Mem_,DT_, IT_, BlockSize_> & r, const DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_> & x) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");

        TimeStamp ts_start;

        if (this->used_elements() == 0)
        {
          r.format();
          return;
        }

        if (r.template elements<Perspective::pod>() == x.template elements<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector x and r must not share the same memory!");


        Statistics::add_flops(this->used_elements() * 2 * BlockSize_);
        Arch::ProductMatVec<Mem_>::template csrsb<DT_, IT_, BlockSize_>(
          r.template elements<Perspective::pod>(), this->val(), this->col_ind(), this->row_ptr(),
          x.template elements<Perspective::pod>(), this->rows(), columns(), used_elements());

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that recieves the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(
                 DenseVector<Mem_,DT_, IT_> & r,
                 const DenseVector<Mem_, DT_, IT_> & x,
                 const DenseVector<Mem_, DT_, IT_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        TimeStamp ts_start;

        if (this->used_elements() == 0)
        {
          r.copy(y);
          return;
        }

        if (r.template elements<Perspective::pod>() == x.template elements<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector x and r must not share the same memory!");

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements() * 3);
          Arch::Defect<Mem_>::csr(r.elements(), y.elements(), this->val(), this->col_ind(),
                                         this->row_ptr(), x.elements(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*A*x
        else
        {
          Statistics::add_flops(this->used_elements() * 3);
          Arch::Axpy<Mem_>::csr(r.elements(), alpha, x.elements(), y.elements(),
                                       this->val(), this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The block vector that recieves the result.
       * \param[in] x The block vector to be multiplied by this matrix.
       * \param[in] y The summand block vector.
       * \param[in] alpha A scalar to scale the product with.
       *
       * \note Every element of each block in the vector x is multiplied with the corresponding single scalar entry in the matrix.
       */
      template<int BlockSize_>
      void apply(
                 DenseVectorBlocked<Mem_,DT_, IT_, BlockSize_> & r,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_> & x,
                 const DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_> & y,
                 const DT_ alpha = DT_(1)) const
      {
        if (r.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of r does not match!");
        if (x.size() != this->columns())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of x does not match!");
        if (y.size() != this->rows())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size of y does not match!");

        TimeStamp ts_start;

        if (this->used_elements() == 0)
        {
          r.copy(y);
          return;
        }

        if (r.template elements<Perspective::pod>() == x.template elements<Perspective::pod>())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector x and r must not share the same memory!");

        // check for special cases
        // r <- y - A*x
        if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
        {
          Statistics::add_flops(this->used_elements() * 3 * BlockSize_);
          Arch::Defect<Mem_>::template csrsb<DT_, IT_, BlockSize_>(
              r.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->val(), this->col_ind(),
              this->row_ptr(), x.template elements<Perspective::pod>(), this->rows(), this->columns(), this->used_elements());
        }
        //r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          r.copy(y);
        // r <- y + alpha*A*x
        else
        {
          Statistics::add_flops(this->used_elements() * 3 * BlockSize_);
          Arch::Axpy<Mem_>::template csrsb<DT_, IT_, BlockSize_>(
              r.template elements<Perspective::pod>(), alpha, x.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(),
              this->val(), this->col_ind(), this->row_ptr(), this->rows(), this->columns(), this->used_elements());
        }

        TimeStamp ts_stop;
        Statistics::add_time_spmv(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Adds a double-matrix product onto this matrix
       *
       * This function performs the following computation:
       * \f[ X \leftarrow X + \alpha D\cdot A\cdot B\f]
       *
       * where
       * - \e X denotes this m-by-n matrix
       * - \e D denotes a m-by-k matrix
       * - \e A denotes a k-by-l matrix
       * - \e B denotes a l-by-n matrix
       *
       * \attention
       * This function assumes that the output matrix already contains the
       * required sparsity pattern. This function will throw an exception
       * if the sparsity pattern of the output matrix is incomplete.
       *
       * \note
       * This function currently only supports data in main memory.
       *
       * \param[in] d, a, b
       * The three matrices to be multiplied
       *
       * \param[in] alpha
       * The scaling factor for the product.
       */
      void add_double_mat_mult(
        const LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>& d,
        const LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>& a,
        const LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>& b,
        const DT_ alpha = DT_(1))
      {
        // validate matrix dimensions
        XASSERT(this->rows() == d.rows());
        XASSERT(d.columns() == a.rows());
        XASSERT(a.columns() == b.rows());
        XASSERT(b.columns() == this->columns());

        // fetch matrix arrays:
        DT_* data_x = this->val();
        const DT_* data_d = d.val();
        const DT_* data_a = a.val();
        const DT_* data_b = b.val();
        const IT_* row_ptr_x = this->row_ptr();
        const IT_* col_idx_x = this->col_ind();
        const IT_* row_ptr_d = d.row_ptr();
        const IT_* col_idx_d = d.col_ind();
        const IT_* row_ptr_a = a.row_ptr();
        const IT_* col_idx_a = a.col_ind();
        const IT_* row_ptr_b = b.row_ptr();
        const IT_* col_idx_b = b.col_ind();

        // loop over all rows of D and X, resp.
        for(IT_ i(0); i < IT_(this->rows()); ++i)
        {
          // loop over all non-zeros D_ik in row i of D
          for(IT_ ik(row_ptr_d[i]); ik  < row_ptr_d[i+1]; ++ik)
          {
            // get column index k
            const IT_ k = col_idx_d[ik];

            // loop over all non-zeros A_kl in row k of A
            for(IT_ kl(row_ptr_a[k]); kl < row_ptr_a[k+1]; ++kl)
            {
              // get column index l
              const IT_ l = col_idx_a[kl];

              // pre-compute factor (alpha * D_ik * A_kl)
              const DT_ omega = alpha *  data_d[ik] * data_a[kl];

              // loop over all non-zeros B_lj in row j of B and
              // loop over all non-zeros X_ij in row i of X and
              // perform a "sparse axpy" of B_l onto X_i, i.e.:
              //   X_i. += (alpha * D_ik * A_kl) * B_l.
              IT_ ij = row_ptr_x[i];
              IT_ lj = row_ptr_b[l];
              while((ij < row_ptr_x[i+1]) && (lj < row_ptr_b[l+1]))
              {
                if(col_idx_x[ij] == col_idx_b[lj])
                {
                  // okay: B_lj contributes to X_ij
                  data_x[ij] += omega * data_b[lj];
                  ++ij;
                  ++lj;
                }
                else if(col_idx_x[ij] < col_idx_b[lj])
                {
                  // entry X_ij exists, but B_lj is missing:
                  // this is a perfectly valid case, so continue with the next non-zero of X_i
                  ++ij;
                }
                else //if(col_idx_x[ij] > col_idx_b[lj])
                {
                  // If we come out here, then the sparsity pattern of X is incomplete:
                  // B_lj is meant to be added onto X_ij, but the entry X_ij is missing
                  throw InternalError(__func__, __FILE__, __LINE__, "Incomplete output matrix structure");
                  //++lj;
                }
              }
            }
          }
        }
      }

      /**
       * \brief Adds a double-matrix product onto this matrix
       *
       * This function performs the following computation:
       * \f[ X \leftarrow X + \alpha D\cdot \textnormal{diag}(A)\cdot B\f]
       *
       * where
       * - \e X denotes this m-by-n matrix
       * - \e D denotes a m-by-l matrix
       * - \e A denotes a vector representing a l-by-l diagonal matrix
       * - \e B denotes a l-by-n matrix
       *
       * \attention
       * This function assumes that the output matrix already contains the
       * required sparsity pattern. This function will throw an exception
       * if the sparsity pattern of the output matrix is incomplete.
       *
       * \note
       * This function currently only supports data in main memory.
       *
       * \param[in] a
       * The vector representing the diagonal matrix A.
       *
       * \param[in] d, b
       * The left and right multiplicant matrices
       *
       * \param[in] alpha
       * The scaling factor for the product.
       */
      void add_double_mat_mult(
        const LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>& d,
        const LAFEM::DenseVector<Mem::Main, DT_, IT_>& a,
        const LAFEM::SparseMatrixCSR<Mem::Main, DT_, IT_>& b,
        const DT_ alpha = DT_(1))
      {
        // validate matrix dimensions
        XASSERT(this->rows() == d.rows());
        XASSERT(d.columns() == a.size());
        XASSERT(a.size() == b.rows());
        XASSERT(b.columns() == this->columns());

        // fetch matrix arrays:
        DT_* data_x = this->val();
        const DT_* data_d = d.val();
        const DT_* data_a = a.elements();
        const DT_* data_b = b.val();
        const IT_* row_ptr_x = this->row_ptr();
        const IT_* col_idx_x = this->col_ind();
        const IT_* row_ptr_d = d.row_ptr();
        const IT_* col_idx_d = d.col_ind();
        const IT_* row_ptr_b = b.row_ptr();
        const IT_* col_idx_b = b.col_ind();

        // loop over all rows of D and X, resp.
        for(IT_ i(0); i < IT_(this->rows()); ++i)
        {
          // loop over all non-zeros D_ik in row i of D
          for(IT_ ik(row_ptr_d[i]); ik  < row_ptr_d[i+1]; ++ik)
          {
            // get column index k
            const IT_ k = col_idx_d[ik];

            // pre-compute factor (alpha * D_ik * A_kk)
            const DT_ omega = alpha *  data_d[ik] * data_a[k];

            // loop over all non-zeros B_kj in row j of B and
            // loop over all non-zeros X_ij in row i of X and
            // perform a "sparse axpy" of B_l onto X_i, i.e.:
            //   X_i. += (alpha * D_ik * A_kk) * B_k.
            IT_ ij = row_ptr_x[i];
            IT_ kj = row_ptr_b[k];
            while((ij < row_ptr_x[i+1]) && (kj < row_ptr_b[k+1]))
            {
              if(col_idx_x[ij] == col_idx_b[kj])
              {
                // okay: B_kj contributes to X_ij
                data_x[ij] += omega * data_b[kj];
                ++ij;
                ++kj;
              }
              else if(col_idx_x[ij] < col_idx_b[kj])
              {
                // entry X_ij exists, but B_kj is missing:
                // this is a perfectly valid case, so continue with the next non-zero of X_i
                ++ij;
              }
              else //if(col_idx_x[ij] > col_idx_b[kj])
              {
                // If we come out here, then the sparsity pattern of X is incomplete:
                // B_kj is meant to be added onto X_ij, but the entry X_ij is missing
                throw InternalError(__func__, __FILE__, __LINE__, "Incomplete output matrix structure");
                //++kj;
              }
            }
          }
        }
      }
      ///@}

      /// \copydoc lump_rows()
      void lump_rows(VectorTypeL& lump) const
      {
        XASSERTM(lump.size() == rows(), "lump vector size does not match matrix row count!");

        Arch::Lumping<Mem_>::csr(lump.elements(), val(), col_ind(), row_ptr(), rows());
      }

      /**
       * \brief Returns the lumped rows vector
       *
       * Each entry in the returned lumped rows vector contains the
       * the sum of all matrix elements in the corresponding row.
       *
       * \returns
       * The lumped vector.
       */
      VectorTypeL lump_rows() const
      {
        VectorTypeL lump = create_vector_l();
        lump_rows(lump);
        return lump;
      }

      /// \copydoc extract_diag()
      void extract_diag(VectorTypeL & diag) const
      {
        XASSERTM(diag.size() == rows(), "diag size does not match matrix row count!");
        XASSERTM(rows() == columns(), "matrix is not square!");

        Arch::Diagonal<Mem_>::csr(diag.elements(), val(), col_ind(), row_ptr(), rows());
      }

      /// extract main diagonal vector from matrix
      VectorTypeL extract_diag() const
      {
        VectorTypeL diag = create_vector_l();
        extract_diag(diag);
        return diag;
      }

      /// Permutate matrix rows and columns according to the given Permutations
      void permute(Adjacency::Permutation & perm_row, Adjacency::Permutation & perm_col)
      {
        if (perm_row.size() == 0 && perm_col.size() == 0)
          return;

        XASSERTM(perm_row.size() == this->rows(), "Container rows does not match permutation size");
        XASSERTM(perm_col.size() == this->columns(), "Container columns does not match permutation size");

        // http://de.mathworks.com/help/matlab/math/sparse-matrix-operations.html#f6-13070
        SparseMatrixCSR<Mem::Main, DT_, IT_> local;
        local.convert(*this);
        IT_ * temp_row_ptr = new IT_[rows() + 1];
        IT_ * temp_col_ind = new IT_[used_elements()];
        DT_ * temp_val = new DT_[used_elements()];

        Index * perm_pos;
        perm_pos = perm_row.get_perm_pos();

        //permute rows from local to temp_*
        Index new_start(0);
        temp_row_ptr[0] = 0;
        for (Index row(0) ; row < local.rows() ; ++row)
        {
          Index row_size(local.row_ptr()[perm_pos[row] + 1] - local.row_ptr()[perm_pos[row]]);

          //iterate over all elements in single one new and old row
          for (Index i(new_start), j(local.row_ptr()[perm_pos[row]]) ; i < new_start + row_size ; ++i, ++j)
          {
            temp_col_ind[i] = local.col_ind()[j];
            temp_val[i] = local.val()[j];
          }

          new_start += row_size;
          temp_row_ptr[row+1] = (IT_)new_start;
        }

        //use inverse col permutation as lookup table: i -> new location of i
        Adjacency::Permutation perm_col_inv = perm_col.inverse();
        perm_pos = perm_col_inv.get_perm_pos();

        //permute columns from temp_* to local
        ::memcpy(local.row_ptr(), temp_row_ptr, (rows() + 1) * sizeof(IT_));
        ::memcpy(local.val(), temp_val, used_elements() * sizeof(DT_));
        for (Index i(0) ; i < used_elements() ; ++i)
        {
          local.col_ind()[i] = (IT_)perm_pos[temp_col_ind[i]];
        }

        delete[] temp_row_ptr;
        delete[] temp_col_ind;
        delete[] temp_val;

        //sort columns in every row by column index
        IT_ swap_key;
        DT_ swap_val;
        for (Index row(0) ; row < rows() ; ++row)
        {
          Index offset(local.row_ptr()[row]);
          Index row_size(local.row_ptr()[row+1] - local.row_ptr()[row]);
          for (Index i(1), j ; i < row_size ; ++i)
          {
            swap_key = local.col_ind()[i + offset];
            swap_val = local.val()[i + offset];
            j = i;
            while (j > 0 && local.col_ind()[j - 1 + offset] > swap_key)
            {
              local.col_ind()[j + offset] = local.col_ind()[j - 1 + offset];
              local.val()[j + offset] = local.val()[j - 1 + offset];
              --j;
            }
            local.col_ind()[j + offset] = swap_key;
            local.val()[j + offset] = swap_val;
          }
        }

        this->assign(local);
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

      /// Returns the number of NNZ-elements of the selected row
      Index get_length_of_line(const Index row) const
      {
        const IT_ * prow_ptr(this->row_ptr());
        return Index(prow_ptr[row + 1] - prow_ptr[row]);
      }

      /// \cond internal

      /// Writes the non-zero-values and matching col-indices of the selected row in allocated arrays
      void set_line(const Index row, DT_ * const pval_set, IT_ * const pcol_set,
                    const Index col_start, const Index stride = 1) const
      {
        const IT_ * prow_ptr(this->row_ptr());
        const IT_ * pcol_ind(this->col_ind());
        const DT_ * pval(this->val());

        const Index start((Index(prow_ptr[row])));
        const Index end((Index(prow_ptr[row + 1] - prow_ptr[row])));
        for (Index i(0); i < end; ++i)
        {
          pval_set[i * stride] = pval[start + i];
          pcol_set[i * stride] = pcol_ind[start + i] + IT_(col_start);
        }
      }
      /// \endcond

      /* ******************************************************************* */
      /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
      /* ******************************************************************* */
    public:
      /** \copydoc Adjactor::get_num_nodes_domain() */
      const inline Index & get_num_nodes_domain() const
      {
        return rows();
      }

      /** \copydoc Adjactor::get_num_nodes_image() */
      const inline Index & get_num_nodes_image() const
      {
        return columns();
      }

      /** \copydoc Adjactor::image_begin() */
      inline ImageIterator image_begin(Index domain_node) const
      {
        XASSERTM(domain_node < rows(), "Domain node index out of range");
        return &this->_indices.at(0)[this->_indices.at(1)[domain_node]];
      }

      /** \copydoc Adjactor::image_end() */
      inline ImageIterator image_end(Index domain_node) const
      {
        XASSERTM(domain_node < rows(), "Domain node index out of range");
        return &this->_indices.at(0)[this->_indices.at(1)[domain_node + 1]];
      }


      /**
       * \brief SparseMatrixCSR comparison operator
       *
       * \param[in] a A matrix to compare with.
       * \param[in] b A matrix to compare with.
       */
      template <typename Mem2_> friend bool operator== (const SparseMatrixCSR & a, const SparseMatrixCSR<Mem2_, DT_, IT_> & b)
      {
        if (a.rows() != b.rows())
          return false;
        if (a.columns() != b.columns())
          return false;
        if (a.used_elements() != b.used_elements())
          return false;
        if (a.zero_element() != b.zero_element())
          return false;

        if(a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && a.get_indices().size() == 0 && b.get_elements().size() == 0 && b.get_indices().size() == 0)
          return true;

        IT_ * col_ind_a;
        IT_ * col_ind_b;
        DT_ * val_a;
        DT_ * val_b;
        IT_ * row_ptr_a;
        IT_ * row_ptr_b;

        bool ret(true);

        if(std::is_same<Mem::Main, Mem_>::value)
        {
          col_ind_a = (IT_*)a.col_ind();
          val_a = (DT_*)a.val();
          row_ptr_a = (IT_*)a.row_ptr();
        }
        else
        {
          col_ind_a = new IT_[a.used_elements()];
          MemoryPool<Mem_>::template download<IT_>(col_ind_a, a.col_ind(), a.used_elements());
          val_a = new DT_[a.used_elements()];
          MemoryPool<Mem_>::template download<DT_>(val_a, a.val(), a.used_elements());
          row_ptr_a = new IT_[a.rows() + 1];
          MemoryPool<Mem_>::template download<IT_>(row_ptr_a, a.row_ptr(), a.rows() + 1);
        }
        if(std::is_same<Mem::Main, Mem2_>::value)
        {
          col_ind_b = (IT_*)b.col_ind();
          val_b = (DT_*)b.val();
          row_ptr_b = (IT_*)b.row_ptr();
        }
        else
        {
          col_ind_b = new IT_[b.used_elements()];
          MemoryPool<Mem2_>::template download<IT_>(col_ind_b, b.col_ind(), b.used_elements());
          val_b = new DT_[b.used_elements()];
          MemoryPool<Mem2_>::template download<DT_>(val_b, b.val(), b.used_elements());
          row_ptr_b = new IT_[b.rows() + 1];
          MemoryPool<Mem2_>::template download<IT_>(row_ptr_b, b.row_ptr(), b.rows() + 1);
        }

        for (Index i(0) ; i < a.used_elements() ; ++i)
        {
          if (col_ind_a[i] != col_ind_b[i])
          {
            ret = false;
            break;
          }
          if (val_a[i] != val_b[i])
          {
            ret = false;
            break;
          }
        }
        if (ret)
        {
          for (Index i(0) ; i < a.rows() + 1; ++i)
          {
            if (row_ptr_a[i] != row_ptr_b[i])
            {
              ret = false;
              break;
            }
          }
        }

        if(! std::is_same<Mem::Main, Mem_>::value)
        {
          delete[] col_ind_a;
          delete[] val_a;
          delete[] row_ptr_a;
        }
        if(! std::is_same<Mem::Main, Mem2_>::value)
        {
          delete[] col_ind_b;
          delete[] val_b;
          delete[] row_ptr_b;
        }

        return ret;
      }

      /**
       * \brief SparseMatrixCSR streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The matrix to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const SparseMatrixCSR & b)
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
    }; //SparseMatrixCSR

    extern template class SparseMatrixCSR<Mem::Main, float, unsigned int>;
    extern template class SparseMatrixCSR<Mem::Main, double, unsigned int>;
    extern template class SparseMatrixCSR<Mem::Main, float, unsigned long>;
    extern template class SparseMatrixCSR<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
    extern template class SparseMatrixCSR<Mem::CUDA, float, unsigned int>;
    extern template class SparseMatrixCSR<Mem::CUDA, double, unsigned int>;
    extern template class SparseMatrixCSR<Mem::CUDA, float, unsigned long>;
    extern template class SparseMatrixCSR<Mem::CUDA, double, unsigned long>;
#endif
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_SPARSE_MATRIX_CSR_HPP
