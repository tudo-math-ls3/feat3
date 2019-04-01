// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_MATRIX_MIRROR_HPP
#define KERNEL_LAFEM_MATRIX_MIRROR_HPP 1

// includes, FEAT
#include <kernel/adjacency/graph.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/matrix_mirror_buffer.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>

// includes, system
#include <vector>
#include <map>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Matrix-Mirror class template.
     *
     * \todo reimplement operations for SparseMatrixELL and SparseMatrixBanded
     *
     * \author Peter Zajac
     */
    template<typename Mem_, typename DT_, typename IT_>
    class MatrixMirror
    {
    public:
      /// arch typedef
      typedef Mem_ MemType;
      /// data-type typedef
      typedef DT_ DataType;
      /// index-type typedef
      typedef IT_ IndexType;
      /// vector mirror type
      typedef LAFEM::VectorMirror<Mem_, DT_, IT_> VectorMirrorType;

    protected:
      /// row-mirror reference
      const VectorMirrorType& _row_mirror;
      /// col-mirror reference
      const VectorMirrorType& _col_mirror;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] row_mirror
       * A reference to the vector-mirror for the rows, assembled on the test space.
       *
       * \param[in] col_mirror
       * A reference to the vector-mirror for the columns, assembled on the trial space.
       */
      explicit MatrixMirror(const VectorMirrorType& row_mirror, const VectorMirrorType& col_mirror) :
        _row_mirror(row_mirror),
        _col_mirror(col_mirror)
      {
      }

      /// move constructor
      MatrixMirror(MatrixMirror&& other) :
        _row_mirror(other._row_mirror),
        _col_mirror(other._col_mirror)
      {
      }

      /// virtual destructor
      virtual ~MatrixMirror()
      {
      }

      /// \returns A reference to the internal row vector-mirror.
      const VectorMirrorType& get_row_mirror() const
      {
        return _row_mirror;
      }

      /// \returns A reference to the internal column vector-mirror.
      const VectorMirrorType& get_col_mirror() const
      {
        return _col_mirror;
      }

      /**
       * \brief Creates a buffer matrix based on a BCSR template matrix.
       *
       * \param[in] matrix
       * A reference to the template matrix. Its structure must be initialised,
       * but the numerical content is ignored
       *
       * \returns
       * A matrix buffer of the same data- and index-type as the template matrix.
       */
      template<int bw_, int bh_>
      LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_> create_buffer(
        const LAFEM::SparseMatrixBCSR<MemType, DT_, IT_, bw_, bh_>& matrix) const
      {
        return LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>(_create_buffer_graph(matrix), Index(bw_*bh_));
      }

      /**
       * \brief Creates a buffer matrix based on a CSR template matrix.
       *
       * \param[in] matrix
       * A reference to the template matrix. Its structure must be initialised,
       * but the numerical content is ignored
       *
       * \returns
       * A matrix buffer of the same data- and index-type as the template matrix.
       */
      LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_> create_buffer(
        const LAFEM::SparseMatrixCSR<MemType, DT_, IT_>& matrix) const
      {
        return LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>(_create_buffer_graph(matrix), Index(1));
      }

      /**
       * \brief Creates a buffer matrix based on an ELL template matrix.
       *
       * \param[in] matrix
       * A reference to the template matrix. Its structure must be initialised,
       * but the numerical content is ignored
       *
       * \returns
       * A matrix buffer of the same data- and index-type as the template matrix.
       */
      LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_> create_buffer(
        const LAFEM::SparseMatrixELL<MemType, DT_, IT_>& matrix) const
      {
        return LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>(_create_buffer_graph(matrix), Index(1));
      }

      /**
       * \brief Creates a buffer matrix based on a banded template matrix.
       *
       * \param[in] matrix
       * A reference to the template matrix. Its structure must be initialised,
       * but the numerical content is ignored
       *
       * \returns
       * A matrix buffer of the same data- and index-type as the template matrix.
       */
      LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_> create_buffer(
        const LAFEM::SparseMatrixBanded<MemType, DT_, IT_>& matrix) const
      {
        return LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>(_create_buffer_graph(matrix), Index(1));
      }

      /**
       * \brief Performs a gather-operation on an operator matrix.
       *
       * \param[in,out] buffer
       * A reference to the buffer matrix.
       *
       * \param[in] matrix
       * A reference to the operator matrix whose entries are to be gathered.
       */
      void gather(
        LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>& buffer,
        const LAFEM::SparseMatrixCSR<MemType, DT_, IT_>& matrix) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));
        XASSERT(buffer.rows() == this->_row_mirror.num_indices());
        XASSERT(buffer.columns() == this->_col_mirror.num_indices());
        XASSERT(matrix.rows() == this->_row_mirror.size());
        XASSERT(matrix.columns() == this->_col_mirror.size());

        // fetch system matrix arrays
        const IT_* row_ptr_a(matrix.row_ptr());
        const IT_* col_idx_a(matrix.col_ind());
        const DT_* val_a(matrix.val());

        // fetch buffer arrays
        const IT_* row_ptr_b(buffer.row_ptr());
        const IT_* col_idx_b(buffer.col_ind());
        DT_* val_b(buffer.val());

        // fetch row/column mirror indices
        const IT_* mir_idx_r(this->_row_mirror.indices());
        const IT_* mir_idx_c(this->_col_mirror.indices());

        // loop over all buffer matrix rows
        for(IT_ i(0); i < IT_(buffer.rows()); ++i)
        {
          // get the row-index
          const IT_ ridx = mir_idx_r[i];

          // loop over all non-zeros in current row
          for(IT_ j(row_ptr_b[i]); j < row_ptr_b[i+1]; ++j)
          {
            // format buffer entry
            val_b[j] = DT_(0);

            // get the column index
            const IT_ cidx = mir_idx_c[col_idx_b[j]];

            // try to find this entry in the input matrix
            for(IT_ k(row_ptr_a[ridx]); k < row_ptr_a[ridx+1]; ++k)
            {
              if(col_idx_a[k] == cidx)
              {
                val_b[j] = val_a[k];
                break;
              }
            }
          }
        }
      }

      /**
       * \brief Performs a scatter-axpy-operation on an operator matrix.
       *
       * \param[in,out] matrix
       * A reference to the operator matrix.
       *
       * \param[in] buffer
       * A reference to the buffer matrix whose entries are to be scattered.
       *
       * \param[in] alpha
       * The scaling factor for the operation.
       */
      void scatter_axpy(
        LAFEM::SparseMatrixCSR<MemType, DT_, IT_>& matrix,
        const LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>& buffer,
        const DT_ alpha = DT_(1)) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));
        XASSERT(buffer.rows() == this->_row_mirror.num_indices());
        XASSERT(buffer.columns() == this->_col_mirror.num_indices());
        XASSERT(matrix.rows() == this->_row_mirror.size());
        XASSERT(matrix.columns() == this->_col_mirror.size());

        // fetch system matrix arrays
        const IT_* row_ptr_a(matrix.row_ptr());
        const IT_* col_idx_a(matrix.col_ind());
        DT_* val_a(matrix.val());

        // fetch buffer arrays
        const IT_* row_ptr_b(buffer.row_ptr());
        const IT_* col_idx_b(buffer.col_ind());
        const DT_* val_b(buffer.val());

        // fetch row/column mirror indices
        const IT_* mir_idx_r(this->_row_mirror.indices());
        const IT_* mir_idx_c(this->_col_mirror.indices());

        // loop over all buffer matrix rows
        for(IT_ i(0); i < IT_(buffer.rows()); ++i)
        {
          // get the row-index
          const IT_ ridx = mir_idx_r[i];

          // loop over all non-zeros in current row
          for(IT_ j(row_ptr_b[i]); j < row_ptr_b[i+1]; ++j)
          {
            // get the column index
            const IT_ cidx = mir_idx_c[col_idx_b[j]];

            // try to find this entry in the input matrix
            for(IT_ k(row_ptr_a[ridx]); k < row_ptr_a[ridx+1]; ++k)
            {
              if(col_idx_a[k] == cidx)
              {
                val_a[k] += alpha*val_b[j];
                break;
              }
            }
          }
        }
      }

      /**
       * \brief Performs a gather-operation on an operator matrix.
       *
       * \param[in,out] buffer
       * A reference to the buffer matrix.
       *
       * \param[in] matrix
       * A reference to the operator matrix whose entries are to be gathered.
       */
      template<int bw_, int bh_>
      void gather(
        LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>& buffer,
        const LAFEM::SparseMatrixBCSR<MemType, DT_, IT_, bw_, bh_>& matrix) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(bw_*bh_));
        XASSERT(buffer.rows() == this->_row_mirror.num_indices());
        XASSERT(buffer.columns() == this->_col_mirror.num_indices());
        XASSERT(matrix.rows() == this->_row_mirror.size());
        XASSERT(matrix.columns() == this->_col_mirror.size());

        typedef Tiny::Matrix<DT_, bw_, bh_> ValueType;

        // fetch system matrix arrays
        const IT_* row_ptr_a(matrix.row_ptr());
        const IT_* col_idx_a(matrix.col_ind());
        const ValueType* val_a(matrix.val());

        // fetch buffer arrays
        const IT_* row_ptr_b(buffer.row_ptr());
        const IT_* col_idx_b(buffer.col_ind());
        ValueType* val_b = reinterpret_cast<ValueType*>(buffer.val());

        // fetch row/column mirror indices
        const IT_* mir_idx_r(this->_row_mirror.indices());
        const IT_* mir_idx_c(this->_col_mirror.indices());

        // loop over all buffer matrix rows
        for(IT_ i(0); i < IT_(buffer.rows()); ++i)
        {
          // get the row-index
          const IT_ ridx = mir_idx_r[i];

          // loop over all non-zeros in current row
          for(IT_ j(row_ptr_b[i]); j < row_ptr_b[i+1]; ++j)
          {
            // format buffer entry
            val_b[j].format();

            // get the column index
            const IT_ cidx = mir_idx_c[col_idx_b[j]];

            // try to find this entry in the input matrix
            for(IT_ k(row_ptr_a[ridx]); k < row_ptr_a[ridx+1]; ++k)
            {
              if(col_idx_a[k] == cidx)
              {
                val_b[j] = val_a[k];
                break;
              }
            }
          }
        }
      }

      /**
       * \brief Performs a scatter-axpy-operation on an operator matrix.
       *
       * \param[in,out] matrix
       * A reference to the operator matrix.
       *
       * \param[in] buffer
       * A reference to the buffer matrix whose entries are to be scattered.
       *
       * \param[in] alpha
       * The scaling factor for the operation.
       */
      template<int bw_, int bh_>
      void scatter_axpy(
        LAFEM::SparseMatrixBCSR<MemType, DT_, IT_, bw_, bh_>& matrix,
        const LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>& buffer,
        const DT_ alpha = DT_(1)) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(bw_*bh_));
        XASSERT(buffer.rows() == this->_row_mirror.num_indices());
        XASSERT(buffer.columns() == this->_col_mirror.num_indices());
        XASSERT(matrix.rows() == this->_row_mirror.size());
        XASSERT(matrix.columns() == this->_col_mirror.size());

        typedef Tiny::Matrix<DT_, bw_, bh_> ValueType;

        // fetch system matrix arrays
        const IT_* row_ptr_a(matrix.row_ptr());
        const IT_* col_idx_a(matrix.col_ind());
        ValueType* val_a(matrix.val());

        // fetch buffer arrays
        const IT_* row_ptr_b(buffer.row_ptr());
        const IT_* col_idx_b(buffer.col_ind());
        const ValueType* val_b = reinterpret_cast<const ValueType*>(buffer.val());

        // fetch row/column mirror indices
        const IT_* mir_idx_r(this->_row_mirror.indices());
        const IT_* mir_idx_c(this->_col_mirror.indices());

        // loop over all buffer matrix rows
        for(IT_ i(0); i < IT_(buffer.rows()); ++i)
        {
          // get the row-index
          const IT_ ridx = mir_idx_r[i];

          // loop over all non-zeros in current row
          for(IT_ j(row_ptr_b[i]); j < row_ptr_b[i+1]; ++j)
          {
            // get the column index
            const IT_ cidx = mir_idx_c[col_idx_b[j]];

            // try to find this entry in the input matrix
            for(IT_ k(row_ptr_a[ridx]); k < row_ptr_a[ridx+1]; ++k)
            {
              if(col_idx_a[k] == cidx)
              {
                val_a[k].axpy(alpha, val_b[j]);
                break;
              }
            }
          }
        }
      }

      /**
       * \brief Performs a gather-operation on an operator matrix.
       *
       * \param[in,out] buffer
       * A reference to the buffer matrix.
       *
       * \param[in] matrix
       * A reference to the operator matrix whose entries are to be gathered.
       */
      void gather(
        LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>& buffer,
        const LAFEM::SparseMatrixELL<MemType, DT_, IT_>& matrix) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));
        XASSERT(buffer.rows() == this->_row_mirror.num_indices());
        XASSERT(buffer.columns() == this->_col_mirror.num_indices());
        XASSERT(matrix.rows() == this->_row_mirror.size());
        XASSERT(matrix.columns() == this->_col_mirror.size());

        // fetch system matrix arrays
        const IT_* cs_a(matrix.cs());
        const IT_* rl_a(matrix.rl());
        const IT_* col_idx_a(matrix.col_ind());
        const DT_* val_a(matrix.val());
        const IT_ C_a(matrix.C());

        // fetch buffer arrays
        const IT_* row_ptr_b(buffer.row_ptr());
        const IT_* col_idx_b(buffer.col_ind());
        DT_* val_b(buffer.val());

        // fetch row/column mirror indices
        const IT_* mir_idx_r(this->_row_mirror.indices());
        const IT_* mir_idx_c(this->_col_mirror.indices());

        // loop over all buffer matrix rows
        for(IT_ i(0); i < IT_(buffer.rows()); ++i)
        {
          // get the row-index
          const IT_ ridx = mir_idx_r[i];

          // get the corresponding ell chunk
          const Index nchunk(Index(Math::floor(ridx / float(C_a))));
          // get the corresponding ell chunk start index
          const Index start(cs_a[nchunk]);
          // get the corresponding ell max row length
          const Index max(rl_a[ridx]);

          // loop over all non-zeros in current row
          for(IT_ j(row_ptr_b[i]); j < row_ptr_b[i+1]; ++j)
          {
            // format buffer entry
            val_b[j] = DT_(0);

            // get the column index
            const IT_ cidx = mir_idx_c[col_idx_b[j]];

            // try to find this entry in the input matrix
            for (Index k(start + ridx - nchunk * C_a), l(0) ; l < max && col_idx_a[k] <= cidx ; k += C_a, ++l)
            {
              if(col_idx_a[k] == cidx)
              {
                val_b[j] = val_a[k];
                break;
              }
            }
          }
        }
      }

      /**
       * \brief Performs a scatter-axpy-operation on an operator matrix.
       *
       * \param[in,out] matrix
       * A reference to the operator matrix.
       *
       * \param[in] buffer
       * A reference to the buffer matrix whose entries are to be scattered.
       *
       * \param[in] alpha
       * The scaling factor for the operation.
       */
      void scatter_axpy(
        LAFEM::SparseMatrixELL<MemType, DT_, IT_>& matrix,
        const LAFEM::MatrixMirrorBuffer<MemType, DT_, IT_>& buffer,
        const DT_ alpha = DT_(1)) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));
        XASSERT(buffer.rows() == this->_row_mirror.num_indices());
        XASSERT(buffer.columns() == this->_col_mirror.num_indices());
        XASSERT(matrix.rows() == this->_row_mirror.size());
        XASSERT(matrix.columns() == this->_col_mirror.size());

        // fetch system matrix arrays
        const IT_* cs_a(matrix.cs());
        const IT_* rl_a(matrix.rl());
        const IT_* col_idx_a(matrix.col_ind());
        DT_* val_a(matrix.val());
        const IT_ C_a(matrix.C());

        // fetch buffer arrays
        const IT_* row_ptr_b(buffer.row_ptr());
        const IT_* col_idx_b(buffer.col_ind());
        const DT_* val_b(buffer.val());

        // fetch row/column mirror indices
        const IT_* mir_idx_r(this->_row_mirror.indices());
        const IT_* mir_idx_c(this->_col_mirror.indices());

        // loop over all buffer matrix rows
        for(IT_ i(0); i < IT_(buffer.rows()); ++i)
        {
          // get the row-index
          const IT_ ridx = mir_idx_r[i];

          // get the corresponding ell chunk
          const Index nchunk(Index(Math::floor(ridx / float(C_a))));
          // get the corresponding ell chunk start index
          const Index start(cs_a[nchunk]);
          // get the corresponding ell max row length
          const Index max(rl_a[ridx]);

          // loop over all non-zeros in current row
          for(IT_ j(row_ptr_b[i]); j < row_ptr_b[i+1]; ++j)
          {
            // get the column index
            const IT_ cidx = mir_idx_c[col_idx_b[j]];

            // try to find this entry in the input matrix
            for (Index k(start + ridx - nchunk * C_a), l(0) ; l < max && col_idx_a[k] <= cidx ; k += C_a, ++l)
            {
              if(col_idx_a[k] == cidx)
              {
                val_a[k] += alpha*val_b[j];
                break;
              }
            }
          }
        }
      }

      /*
       * \brief Performs a gather-operation on an operator matrix.
       *
       * \param[in,out] buffer
       * A reference to the buffer matrix.
       *
       * \param[in] matrix
       * A reference to the operator matrix whose entries are to be gathered.
       */
      /*template<typename Tx_, typename Ix_, typename Ty_, typename Iy_>
      void gather(
        LAFEM::MatrixMirrorBuffer<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::SparseMatrixBanded<Mem::Main, Ty_, Iy_>& matrix) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));

        //const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_gather());
        //const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_gather());
        const auto& row_mir_mat = _row_gather;
        const auto& col_mir_mat = _col_gather;
        DataType* _work = _vec_work.data();

        // fetch row mirror arrays
        const IndexType* row_ptr_a(row_mir_mat.row_ptr());
        const IndexType* col_idx_a(row_mir_mat.col_ind());
        const DataType* av(row_mir_mat.val());

        // fetch col mirror arrays
        const IndexType* row_ptr_b(col_mir_mat.row_ptr());
        const IndexType* col_idx_b(col_mir_mat.col_ind());
        const DataType* bv(col_mir_mat.val());

        // fetch system matrix arrays
        const Iy_ num_rows_y(Iy_(matrix.rows()));
        const Iy_ num_of_offsets_y(Iy_(matrix.num_of_offsets()));
        const Iy_* offsets_y(matrix.offsets());
        const Ty_* yv(matrix.val());

        // fetch buffer arrays
        const Ix_* row_ptr_x(buffer.row_ptr());
        const Ix_* col_idx_x(buffer.col_ind());
        Tx_* xv(buffer.val());

        // In the following, we have to compute:
        //    X := A * Z := A * Y * B^T,
        // where:
        // X is the buffer matrix
        // Y is the system matrix
        // A is the row-mirror gather matrix
        // B is the col-mirror gather matrix

        // loop over all buffer rows (X)
        Index nrows_buf(buffer.rows());
        for(Index irow_x(0); irow_x < nrows_buf; ++irow_x)
        {
          Index irow_a(irow_x); // row of a := row of x

          // loop over all non-zeroes in the buffer row (X_i.)
          for(Ix_ ix(row_ptr_x[irow_x]); ix < row_ptr_x[irow_x + 1]; ++ix)
          {
            // init result
            Tx_ x_ij(Tx_(0));

            // fetch the column index
            Ix_ irow_b(col_idx_x[ix]); // row of b := col of x

            // loop over all non-zeroes of the col-mirror (B_j.)
            for(IndexType ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
            {
              // and densify the sparse row B_j.
              _work[col_idx_b[ib]] = bv[ib];
            }

            // loop over all non-zeroes of the row-mirror (A_i.)
            for(IndexType ia(row_ptr_a[irow_a]); ia < row_ptr_a[irow_a + 1]; ++ia)
            {
              // fetch the column index
              Iy_ irow_y(col_idx_a[ia]); // row of y := col of a

              // temporary entry: Z_kj := Y_k. * B_j.
              Tx_ z_kj(Tx_(0));

              // loop over all non-zeroes of the system matrix (Y_k.)
              for(IndexType k(0); k < num_of_offsets_y; ++k)
              {
                if(offsets_y[k] + irow_y + 1 >= num_rows_y && offsets_y[k] + irow_y + 1 < 2 * num_rows_y)
                {
                  z_kj += Tx_(yv[k * num_rows_y + irow_y]) * Tx_(_work[offsets_y[k] + irow_y + 1 - num_rows_y]);
                }
              }

              // update x_ij += a_ik * z_kj
              x_ij += Tx_(av[ia]) * z_kj;
            }

            // store X_ij
            xv[ix] = x_ij;

            // reset temporary data
            for(IndexType ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
            {
              _work[col_idx_b[ib]] = DataType(0);
            }
          }
        }
      }*/

      /*
       * \brief Performs a scatter-axpy-operation on an operator matrix.
       *
       * \param[in,out] matrix
       * A reference to the operator matrix.
       *
       * \param[in] buffer
       * A reference to the buffer matrix whose entries are to be scattered.
       *
       * \param[in] alpha
       * The scaling factor for the operation.
       */
      /*template<typename Ty_, typename Iy_, typename Tx_, typename Ix_>
      void scatter_axpy(
        LAFEM::SparseMatrixBanded<Mem::Main, Ty_, Iy_>& matrix,
        const LAFEM::MatrixMirrorBuffer<Mem::Main, Tx_, Ix_>& buffer,
        const Ty_ alpha = Ty_(1)) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));

        //const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_scatter());
        //const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_scatter());
        const auto& row_mir_mat = _row_scatter;
        const auto& col_mir_mat = _col_scatter;
        DataType* _work = _vec_work.data();

        // fetch row-mirror arrays
        const IndexType* row_ptr_a(row_mir_mat.row_ptr());
        const IndexType* col_idx_a(row_mir_mat.col_ind());
        const DataType* av(row_mir_mat.val());

        // fetch col-mirror arrays
        const IndexType* row_ptr_b(col_mir_mat.row_ptr());
        const IndexType* col_idx_b(col_mir_mat.col_ind());
        const DataType* bv(col_mir_mat.val());

        // fetch system matrix arrays
        const Iy_ num_rows_y(Iy_(matrix.rows()));
        const Iy_ num_of_offsets_y(Iy_(matrix.num_of_offsets()));
        const Iy_* offsets_y(matrix.offsets());
        Ty_* yv(matrix.val());

        // fetch buffer arrays
        const Ix_* row_ptr_x(buffer.row_ptr());
        const Ix_* col_idx_x(buffer.col_ind());
        const Tx_* xv(buffer.val());

        // In the following, we have to compute:
        //    Y := B * Z := B * X * A^T,
        // where:
        // Y is the system matrix
        // X is the buffer matrix
        // A is the row-mirror scatter matrix
        // B is the col-mirror scatter matrix

        // loop over all system matrix rows (Y)
        Index nrows_sys(matrix.rows());
        for(Index irow_y(0); irow_y < nrows_sys; ++irow_y)
        {
          Index irow_a(irow_y); // row of a := row of y

          // skip if the row of a is empty
          if(row_ptr_a[irow_a] >= row_ptr_a[irow_a + 1])
            continue;

          // loop over all non-zeroes in the system row (Y_i.)
          for(IndexType k(0); k < num_of_offsets_y; ++k)
          {
            if(offsets_y[k] + irow_y + 1 >= num_rows_y && offsets_y[k] + irow_y + 1 < 2 * num_rows_y)
            {
              Iy_ iy(k * num_rows_y + irow_y);
              // init result
              Ty_ y_ij(Ty_(0));

              // fetch the column index
              IndexType irow_b(offsets_y[k] + irow_y + 1 - num_rows_y); // row of b := col of y

              // skip if the row of b is empty
              if(row_ptr_b[irow_b] >= row_ptr_b[irow_b + 1])
                continue;

              // loop over all non-zeroes of the col-mirror (B_j.)
              for(IndexType ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
              {
                // and densify the sparse row B_j.
                _work[col_idx_b[ib]] = bv[ib];
              }

              // loop over all non-zeroes of the row-mirror (A_i.)
              for(IndexType ia(row_ptr_a[irow_a]); ia < row_ptr_a[irow_a + 1]; ++ia)
              {
                // fetch the column index
                IndexType irow_x(col_idx_a[ia]); // row of x := col of a

                // temporary entry: Z_kj := X_k. * B_j.
                Ty_ z_kj(Ty_(0));

                // loop over all non-zeroes of the buffer matrix (X_k.)
                for(Ix_ ix(row_ptr_x[irow_x]); ix < row_ptr_x[irow_x + 1]; ++ix)
                {
                  z_kj += Ty_(xv[ix]) * Ty_(_work[col_idx_x[ix]]);
                }

                // update Y_ij += A_ik * Z_kj
                y_ij += Ty_(av[ia]) * z_kj;
              }

              // store Y_ij
              yv[iy] += alpha*y_ij;

              // reset temporary data
              for(IndexType ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
              {
                _work[col_idx_b[ib]] = DataType(0);
              }
            }
          }
        }
      }*/

    protected:
      template<typename MT_>
      Adjacency::Graph _create_buffer_graph(const MT_& tmpl_mat) const
      {
        // compute the sparsity pattern of
        //  X := A * Y * B^T,
        // where:
        // X is the buffer matrix
        // Y is the system matrix
        // A is the row-mirror gather matrix
        // B is the col-mirror gather matrix

        // render the matrix structure to a graph to obtain
        // the celebrated ptr/idx array pair
        Adjacency::Graph mat_graph(Adjacency::RenderType::injectify, tmpl_mat);
        const Index* dom_ptr = mat_graph.get_domain_ptr();
        const Index* img_idx = mat_graph.get_image_idx();

        // get row/column mirror indices
        const Index nrows = _row_mirror.num_indices();
        const Index ncols = _col_mirror.num_indices();
        const auto* cidx = _col_mirror.indices();
        const auto* ridx = _row_mirror.indices();

        // build map of column mirror indices
        std::map<Index,Index> col_map;
        for(Index i(0); i < ncols; ++i)
          col_map.emplace(cidx[i], i);

        // count number of non-zeros in matrix buffer
        Index count(0);
        for(Index i(0); i < nrows; ++i)
        {
          const Index irow = ridx[i];
          for(Index j(dom_ptr[irow]); j < dom_ptr[irow+1]; ++j)
          {
            if(col_map.find(img_idx[j]) != col_map.end())
              ++count;
          }
        }

        // allocate output graph
        Adjacency::Graph graph(nrows, ncols, count);
        Index* row_ptr = graph.get_domain_ptr();
        Index* col_idx = graph.get_image_idx();

        // compute output graph
        row_ptr[0] = Index(0);
        for(Index i(0); i < nrows; ++i)
        {
          const Index irow = ridx[i];
          Index k = row_ptr[i];
          for(Index j(dom_ptr[irow]); j < dom_ptr[irow+1]; ++j)
          {
            auto it = col_map.find(img_idx[j]);
            if(it != col_map.end())
            {
              col_idx[k] = it->second;
              ++k;
            }
          }
          row_ptr[i+1] = k;
        }
        XASSERT(row_ptr[nrows] == count);

        // sort column indices
        graph.sort_indices();

        return graph;
      }
    }; // class MatrixMirror<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_MATRIX_MIRROR_HPP
