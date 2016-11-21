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

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Matrix-Mirror class template.
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

      template<typename DT2_, typename IT2_>
      using BufferType = LAFEM::MatrixMirrorBuffer<Mem::Main, DT2_, IT2_>;

    protected:
      /// row-mirror reference
      const VectorMirrorType& _row_mirror;
      /// col-mirror reference
      const VectorMirrorType& _col_mirror;
      /// \cond internal
      // mutable work array
      mutable std::vector<DataType> _vec_work;
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
        _vec_work(std::max(col_mirror.get_gather().columns(), col_mirror.get_scatter().columns()), DataType(0))
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
      template<typename DT2_, typename IT2_, int bw_, int bh_>
      LAFEM::MatrixMirrorBuffer<Mem::Main, DT2_, IT2_> create_buffer(
        const LAFEM::SparseMatrixBCSR<Mem::Main, DT2_, IT2_, bw_, bh_>& matrix) const
      {
        Adjacency::Graph graph = _create_buffer_graph(matrix);
        return LAFEM::MatrixMirrorBuffer<Mem::Main, DT_, IT_>(graph, Index(bw_*bh_));
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
      template<typename DT2_, typename IT2_>
      LAFEM::MatrixMirrorBuffer<Mem::Main, DT2_, IT2_> create_buffer(
        const LAFEM::SparseMatrixCSR<Mem::Main, DT2_, IT2_>& matrix) const
      {
        Adjacency::Graph graph = _create_buffer_graph(matrix);
        return LAFEM::MatrixMirrorBuffer<Mem::Main, DT_, IT_>(graph, Index(1));
      }

      /**
       * \brief Creates a buffer matrix based on a ELL template matrix.
       *
       * \param[in] matrix
       * A reference to the template matrix. Its structure must be initialised,
       * but the numerical content is ignored
       *
       * \returns
       * A matrix buffer of the same data- and index-type as the template matrix.
       */
      template<typename DT2_, typename IT2_>
      LAFEM::MatrixMirrorBuffer<Mem::Main, DT2_, IT2_> create_buffer(
        const LAFEM::SparseMatrixELL<Mem::Main, DT2_, IT2_>& matrix) const
      {
        Adjacency::Graph graph = _create_buffer_graph(matrix);
        return LAFEM::MatrixMirrorBuffer<Mem::Main, DT_, IT_>(graph, Index(1));
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
      template<typename DT2_, typename IT2_>
      LAFEM::MatrixMirrorBuffer<Mem::Main, DT2_, IT2_> create_buffer(
        const LAFEM::SparseMatrixBanded<Mem::Main, DT2_, IT2_>& matrix) const
      {
        Adjacency::Graph graph = _create_buffer_graph(matrix);
        return LAFEM::MatrixMirrorBuffer<Mem::Main, DT_, IT_>(graph, Index(1));
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
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_>
      void gather(
        LAFEM::MatrixMirrorBuffer<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::SparseMatrixCSR<Mem::Main, Ty_, Iy_>& matrix) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));

        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_gather());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_gather());
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
        const Iy_* row_ptr_y(matrix.row_ptr());
        const Iy_* col_idx_y(matrix.col_ind());
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
              for(Iy_ iy(row_ptr_y[irow_y]); iy < row_ptr_y[irow_y + 1]; ++iy)
              {
                z_kj += Tx_(yv[iy]) * Tx_(_work[col_idx_y[iy]]);
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
      template<typename Ty_, typename Iy_, typename Tx_, typename Ix_>
      void scatter_axpy(
        LAFEM::SparseMatrixCSR<Mem::Main, Ty_, Iy_>& matrix,
        const LAFEM::MatrixMirrorBuffer<Mem::Main, Tx_, Ix_>& buffer,
        const Ty_ alpha = Ty_(1)) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));

        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_scatter());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_scatter());
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
        const Iy_* row_ptr_y(matrix.row_ptr());
        const Iy_* col_idx_y(matrix.col_ind());
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
          for(Iy_ iy(row_ptr_y[irow_y]); iy < row_ptr_y[irow_y + 1]; ++iy)
          {
            // init result
            Ty_ y_ij(Ty_(0));

            // fetch the column index
            IndexType irow_b(col_idx_y[iy]); // row of b := col of y

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

            // update Y_ij
            yv[iy] += alpha*y_ij;

            // reset temporary data
            for(IndexType ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
            {
              _work[col_idx_b[ib]] = DataType(0);
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
      template</*typename DT_, typename IT_,*/ int bw_, int bh_>
      void gather(
        LAFEM::MatrixMirrorBuffer<Mem::Main, DT_, IT_>& buffer,
        const LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, bw_, bh_>& matrix) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(bw_*bh_));

        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_gather());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_gather());
        DataType* _work = _vec_work.data();

        typedef Tiny::Matrix<DT_, bw_, bh_> ValueType;

        // fetch row mirror arrays
        const IndexType* row_ptr_a(row_mir_mat.row_ptr());
        const IndexType* col_idx_a(row_mir_mat.col_ind());
        const DataType* av(row_mir_mat.val());

        // fetch col mirror arrays
        const IndexType* row_ptr_b(col_mir_mat.row_ptr());
        const IndexType* col_idx_b(col_mir_mat.col_ind());
        const DataType* bv(col_mir_mat.val());

        // fetch system matrix arrays
        const IT_* row_ptr_y(matrix.row_ptr());
        const IT_* col_idx_y(matrix.col_ind());
        const ValueType* yv(matrix.val());

        // fetch buffer arrays
        const IT_* row_ptr_x(buffer.row_ptr());
        const IT_* col_idx_x(buffer.col_ind());
        ValueType* xv = reinterpret_cast<ValueType*>(buffer.val());

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
          for(IT_ ix(row_ptr_x[irow_x]); ix < row_ptr_x[irow_x + 1]; ++ix)
          {
            // init result
            ValueType& x_ij = xv[ix];
            x_ij.format();

            // fetch the column index
            IT_ irow_b(col_idx_x[ix]); // row of b := col of x

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
              IT_ irow_y(col_idx_a[ia]); // row of y := col of a

              // temporary entry: Z_kj := Y_k. * B_j.
              ValueType z_kj(DT_(0));

              // loop over all non-zeroes of the system matrix (Y_k.)
              for(IT_ iy(row_ptr_y[irow_y]); iy < row_ptr_y[irow_y + 1]; ++iy)
              {
                z_kj.axpy(_work[col_idx_y[iy]], yv[iy]);
              }

              // update x_ij += a_ik * z_kj
              x_ij.axpy(av[ia], z_kj);
            }

            // reset temporary data
            for(IndexType ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
            {
              _work[col_idx_b[ib]] = DataType(0);
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
      template</*typename DT_, typename IT_,*/ int bw_, int bh_>
      void scatter_axpy(
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, bw_, bh_>& matrix,
        const LAFEM::MatrixMirrorBuffer<Mem::Main, DT_, IT_>& buffer,
        const DT_ alpha = DT_(1)) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(bw_*bh_));

        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_scatter());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_scatter());
        DataType* _work = _vec_work.data();

        typedef Tiny::Matrix<DT_, bw_, bh_> ValueType;

        // fetch row-mirror arrays
        const IndexType* row_ptr_a(row_mir_mat.row_ptr());
        const IndexType* col_idx_a(row_mir_mat.col_ind());
        const DataType* av(row_mir_mat.val());

        // fetch col-mirror arrays
        const IndexType* row_ptr_b(col_mir_mat.row_ptr());
        const IndexType* col_idx_b(col_mir_mat.col_ind());
        const DataType* bv(col_mir_mat.val());

        // fetch system matrix arrays
        const IT_* row_ptr_y(matrix.row_ptr());
        const IT_* col_idx_y(matrix.col_ind());
        ValueType* yv(matrix.val());

        // fetch buffer arrays
        const IT_* row_ptr_x(buffer.row_ptr());
        const IT_* col_idx_x(buffer.col_ind());
        const ValueType* xv = reinterpret_cast<const ValueType*>(buffer.val());

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
          for(IT_ iy(row_ptr_y[irow_y]); iy < row_ptr_y[irow_y + 1]; ++iy)
          {
            // init result
            ValueType& y_ij = yv[iy];

            // fetch the column index
            IndexType irow_b(col_idx_y[iy]); // row of b := col of y

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
              ValueType z_kj(DT_(0));
              z_kj.format();

              // loop over all non-zeroes of the buffer matrix (X_k.)
              for(IT_ ix(row_ptr_x[irow_x]); ix < row_ptr_x[irow_x + 1]; ++ix)
              {
                z_kj.axpy(_work[col_idx_x[ix]], xv[ix]);
              }

              // update Y_ij += A_ik * Z_kj
              y_ij.axpy(alpha*av[ia], z_kj);
            }

            // reset temporary data
            for(IndexType ib(row_ptr_b[irow_b]); ib < row_ptr_b[irow_b + 1]; ++ib)
            {
              _work[col_idx_b[ib]] = DataType(0);
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
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_>
      void gather(
        LAFEM::MatrixMirrorBuffer<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::SparseMatrixELL<Mem::Main, Ty_, Iy_>& matrix) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));

        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_gather());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_gather());
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
        const Iy_ C_y(Iy_(matrix.C()));
        const Iy_* cs_y(matrix.cs());
        const Iy_* rl_y(matrix.rl());
        const Iy_* col_idx_y(matrix.col_ind());
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
              for(Iy_ iy(cs_y[irow_y/C_y] + irow_y%C_y); iy < cs_y[irow_y/C_y] + irow_y%C_y + rl_y[irow_y] * C_y; iy += C_y)
              {
                z_kj += Tx_(yv[iy]) * Tx_(_work[col_idx_y[iy]]);
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
      template<typename Ty_, typename Iy_, typename Tx_, typename Ix_>
      void scatter_axpy(
        LAFEM::SparseMatrixELL<Mem::Main, Ty_, Iy_>& matrix,
        const LAFEM::MatrixMirrorBuffer<Mem::Main, Tx_, Ix_>& buffer,
        const Ty_ alpha = Ty_(1)) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));

        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_scatter());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_scatter());
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
        const Iy_ C_y(Iy_(matrix.C()));
        const Iy_* cs_y(matrix.cs());
        const Iy_* rl_y(matrix.rl());
        const Iy_* col_idx_y(matrix.col_ind());
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
          for(Iy_ iy(cs_y[irow_y/C_y] + irow_y%C_y); iy < cs_y[irow_y/C_y] + irow_y%C_y + rl_y[irow_y] * C_y; iy += C_y)
          {
            // init result
            Ty_ y_ij(Ty_(0));

            // fetch the column index
            IndexType irow_b(col_idx_y[iy]); // row of b := col of y

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

      /**
       * \brief Performs a gather-operation on an operator matrix.
       *
       * \param[in,out] buffer
       * A reference to the buffer matrix.
       *
       * \param[in] matrix
       * A reference to the operator matrix whose entries are to be gathered.
       */
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_>
      void gather(
        LAFEM::MatrixMirrorBuffer<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::SparseMatrixBanded<Mem::Main, Ty_, Iy_>& matrix) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));

        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_gather());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_gather());
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
      template<typename Ty_, typename Iy_, typename Tx_, typename Ix_>
      void scatter_axpy(
        LAFEM::SparseMatrixBanded<Mem::Main, Ty_, Iy_>& matrix,
        const LAFEM::MatrixMirrorBuffer<Mem::Main, Tx_, Ix_>& buffer,
        const Ty_ alpha = Ty_(1)) const
      {
        XASSERT(buffer.entries_per_nonzero() == Index(1));

        const typename VectorMirrorType::MirrorMatrixType& row_mir_mat(_row_mirror.get_scatter());
        const typename VectorMirrorType::MirrorMatrixType& col_mir_mat(_col_mirror.get_scatter());
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
      }

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
        Adjacency::Graph tmp1(Adjacency::rt_injectify, _row_mirror.get_gather(), tmpl_mat);
        Adjacency::Graph tmp2(Adjacency::rt_injectify, tmp1, _col_mirror.get_scatter());
        tmp2.sort_indices();
        return tmp2;
      }
    }; // class MatrixMirror<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_MATRIX_MIRROR_HPP
