#pragma once
#ifndef KERNEL_LAFEM_MATRIX_MIRROR_HPP
#define KERNEL_LAFEM_MATRIX_MIRROR_HPP 1

// includes, FEAST
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/util/dynamic_graph.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Matrix-Mirror class template.
     *
     * \author Peter Zajac
     */
    template<
      typename Arch_,
      typename DataType_>
    class MatrixMirror
    {
    public:
      /// arch typedef
      typedef Arch_ MemType;
      /// data-type typedef
      typedef DataType_ DataType;

      /// corresponding vector-mirror type
      typedef VectorMirror<MemType, DataType> VectorMirrorType;

    protected:
      /// row-mirror reference
      const VectorMirrorType& _row_mirror;
      /// col-mirror reference
      const VectorMirrorType& _col_mirror;
      /// \cond internal
      /// mutable work array
      mutable DataType* _work;
      /// \endcond

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
        _col_mirror(col_mirror),
        _work(nullptr)
      {
        Index n = std::max(
          col_mirror.get_gather_dual().columns(),
          col_mirror.get_scatter_dual().columns());
        _work = new DataType[n];
        for(Index i(0); i < n; ++i)
        {
          _work[i] = DataType(0);
        }
      }

      /// virtual destructor
      virtual ~MatrixMirror()
      {
        if(_work != nullptr)
        {
          delete [] _work;
        }
      }

      /**
       * \brief Creates a buffer matrix based on a template matrix.
       *
       * \param[in] tmpl_mat
       * A reference to the template matrix.
       *
       * \returns
       * A new matrix of the same data-type and arch as the template matrix.
       */
      template<typename DataType2_>
      LAFEM::SparseMatrixCSR<Mem::Main, DataType2_>
        create_buffer(const LAFEM::SparseMatrixCSR<Mem::Main, DataType2_>& tmpl_mat) const
      {
        // compute the sparsity pattern of
        //  X := A * Y * B^T,
        // where:
        // X is the buffer matrix
        // Y is the system matrix
        // A is the row-mirror dual-gather matrix
        // B is the col-mirror dual-gather matrix
        DynamicGraph graph(DynamicGraph::rt_as_is, _row_mirror.get_gather_dual());
        graph.compose(tmpl_mat);
        graph.compose(_col_mirror.get_scatter_dual());
        return LAFEM::SparseMatrixCSR<Mem::Main, DataType2_>(Graph(Graph::rt_as_is, graph));
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
      template<
        typename Tx_,
        typename Ty_>
      void gather_op(
        LAFEM::SparseMatrixCSR<Mem::Main, Tx_>& buffer,
        const LAFEM::SparseMatrixCSR<Mem::Main, Ty_>& matrix) const
      {
        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_gather_dual());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_gather_dual());

        // fetch row mirror arrays
        const Index* row_ptr_a(row_mir_mat.row_ptr());
        const Index* row_end_a(row_mir_mat.row_ptr_end());
        const Index* col_idx_a(row_mir_mat.col_ind());
        const DataType* av(row_mir_mat.val());

        // fetch col mirror arrays
        const Index* row_ptr_b(col_mir_mat.row_ptr());
        const Index* row_end_b(col_mir_mat.row_ptr_end());
        const Index* col_idx_b(col_mir_mat.col_ind());
        const DataType* bv(col_mir_mat.val());

        // fetch system matrix arrays
        const Index* row_ptr_y(matrix.row_ptr());
        const Index* row_end_y(matrix.row_ptr_end());
        const Index* col_idx_y(matrix.col_ind());
        const Ty_* yv(matrix.val());

        // fetch buffer arrays
        const Index* row_ptr_x(buffer.row_ptr());
        const Index* row_end_x(buffer.row_ptr_end());
        const Index* col_idx_x(buffer.col_ind());
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
          for(Index ix(row_ptr_x[irow_x]); ix < row_end_x[irow_x]; ++ix)
          {
            // init result
            Tx_ x_ij(0);

            // fetch the column index
            Index irow_b(col_idx_x[ix]); // row of b := col of x

            // loop over all non-zeroes of the col-mirror (B_j.)
            for(Index ib(row_ptr_b[irow_b]); ib < row_end_b[irow_b]; ++ib)
            {
              // and densify the sparse row B_j.
              _work[col_idx_b[ib]] = bv[ib];
            }

            // loop over all non-zeroes of the row-mirror (A_i.)
            for(Index ia(row_ptr_a[irow_a]); ia < row_end_a[irow_a]; ++ia)
            {
              // fetch the column index
              Index irow_y(col_idx_a[ia]); // row of y := col of a

              // temporary entry: Z_kj := Y_k. * B_j.
              Tx_ z_kj(0);

              // loop over all non-zeroes of the system matrix (Y_k.)
              for(Index iy(row_ptr_y[irow_y]); iy < row_end_y[irow_y]; ++iy)
              {
                z_kj += Tx_(yv[iy]) * Tx_(_work[col_idx_y[iy]]);
              }

              // update x_ij += a_ik * z_kj
              x_ij += Tx_(av[ia]) * z_kj;
            }

            // store X_ij
            xv[ix] = x_ij;

            // reset temporary data
            for(Index ib(row_ptr_b[irow_b]); ib < row_end_b[irow_b]; ++ib)
            {
              _work[col_idx_b[ib]] = DataType(0);
            }
          }
        }
      }

      /**
       * \brief Performs a scatter-operation on an operator matrix.
       *
       * \param[in,out] matrix
       * A reference to the operator matrix.
       *
       * \param[in]
       * A reference to the buffer matrix whose entries are to be scattered.
       */
      template<
        typename Ty_,
        typename Tx_>
      void scatter_op(
        LAFEM::SparseMatrixCSR<Mem::Main, Ty_>& matrix,
        const LAFEM::SparseMatrixCSR<Mem::Main, Tx_>& buffer) const
      {
        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_scatter_dual());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_scatter_dual());

        // fetch row-mirror arrays
        const Index* row_ptr_a(row_mir_mat.row_ptr());
        const Index* row_end_a(row_mir_mat.row_ptr_end());
        const Index* col_idx_a(row_mir_mat.col_ind());
        const DataType* av(row_mir_mat.val());

        // fetch col-mirror arrays
        const Index* row_ptr_b(col_mir_mat.row_ptr());
        const Index* row_end_b(col_mir_mat.row_ptr_end());
        const Index* col_idx_b(col_mir_mat.col_ind());
        const DataType* bv(col_mir_mat.val());

        // fetch system matrix arrays
        const Index* row_ptr_y(matrix.row_ptr());
        const Index* row_end_y(matrix.row_ptr_end());
        const Index* col_idx_y(matrix.col_ind());
        Ty_* yv(matrix.val());

        // fetch buffer arrays
        const Index* row_ptr_x(buffer.row_ptr());
        const Index* row_end_x(buffer.row_ptr_end());
        const Index* col_idx_x(buffer.col_ind());
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
          if(row_ptr_a[irow_a] >= row_end_a[irow_a])
            continue;

          // loop over all non-zeroes in the system row (Y_i.)
          for(Index iy(row_ptr_y[irow_y]); iy < row_end_y[irow_y]; ++iy)
          {
            // init result
            Ty_ y_ij(0);

            // fetch the column index
            Index irow_b(col_idx_y[iy]); // row of b := col of y

            // skip if the row of b is empty
            if(row_ptr_b[irow_b] >= row_end_b[irow_b])
              continue;

            // loop over all non-zeroes of the col-mirror (B_j.)
            for(Index ib(row_ptr_b[irow_b]); ib < row_end_b[irow_b]; ++ib)
            {
              // and densify the sparse row B_j.
              _work[col_idx_b[ib]] = bv[ib];
            }

            // loop over all non-zeroes of the row-mirror (A_i.)
            for(Index ia(row_ptr_a[irow_a]); ia < row_end_a[irow_a]; ++ia)
            {
              // fetch the column index
              Index irow_x(col_idx_a[ia]); // row of x := col of a

              // temporary entry: Z_kj := X_k. * B_j.
              Ty_ z_kj(0);

              // loop over all non-zeroes of the buffer matrix (X_k.)
              for(Index ix(row_ptr_x[irow_x]); ix < row_end_x[irow_x]; ++ix)
              {
                z_kj += Ty_(xv[ix]) * Ty_(_work[col_idx_x[ix]]);
              }

              // update Y_ij += A_ik * Z_kj
              y_ij += Ty_(av[ia]) * z_kj;
            }

            // store Y_ij
            yv[iy] = y_ij;

            // reset temporary data
            for(Index ib(row_ptr_b[irow_b]); ib < row_end_b[irow_b]; ++ib)
            {
              _work[col_idx_b[ib]] = DataType(0);
            }
          }
        }
      }
    }; // class MatrixMirror<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_MATRIX_MIRROR_HPP
