#pragma once
#ifndef KERNEL_LAFEM_VECTOR_MIRROR_HPP
#define KERNEL_LAFEM_VECTOR_MIRROR_HPP 1

// includes, FEAST
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_vector.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/arch/gather_prim.hpp>
#include <kernel/lafem/arch/scatter_prim.hpp>
#include <kernel/lafem/arch/scatter_axpy_prim.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Handles vector prolongation, restriction and serialisation
     *
     * \tparam Mem_
     * Memory architecture
     *
     * \tparam DataType_
     * Data type in the vector
     *
     * \tparam IndexType_
     * Type for indexing
     *
     * A Mirror handles the restriction of a given vector to a subvector of itself and the serialisation into
     * buffer vector (which always is a LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>), as well as the
     * prolongation of such buffer vectors. The buffer vectors then can be communicated via MPI etc.
     *
     * This can be used for restricting a vector associated with a Mesh to a subset of that Mesh, i.e. identified by
     * a MeshPart. The standard use case is that we have an FE coefficient vector living on a patch that we want to
     * communicate. For this, we need all coefficients (= vector entries) that contribute to other patches, and
     * these are identified by a MeshPart called Halo. The mirror is constructed from the restriction (gather) and
     * prolongation (scatter) matrices representing the adjacency structure of the underlying FE space and thus the
     * coefficients.
     *
     * \author Peter Zajac
     * \author Markus Geveler
     */
    template<
      typename Mem_,
      typename DataType_,
      typename IndexType_ = Index>
    class VectorMirror
    {
    public:
      /// arch typedef
      typedef Mem_ MemType;
      /// data-type typedef
      typedef DataType_ DataType;
      /// index-type typedef
      typedef IndexType_ IndexType;

      /// mirror matrix typedef
      typedef SparseMatrixCSR<Mem_, DataType_, IndexType_> MirrorMatrixType;

      /// corresponding vector
      typedef DenseVector<MemType, DataType, IndexType> VectorType;

      /// gather-mirror matrix
      MirrorMatrixType _mirror_gather;
      /// scatter-mirror matrix
      MirrorMatrixType _mirror_scatter;

      /// default constructor
      VectorMirror() :
        _mirror_gather(),
        _mirror_scatter()
      {
      }

      /// explicit constructor
      explicit VectorMirror(MirrorMatrixType&& mirror_gather, MirrorMatrixType&& mirror_scatter) :
        _mirror_gather(std::move(mirror_gather)),
        _mirror_scatter(std::move(mirror_scatter))
      {
      }

      /// move-ctor
      VectorMirror(VectorMirror&& other) :
        _mirror_gather(std::move(other._mirror_gather)),
        _mirror_scatter(std::move(other._mirror_scatter))
      {
      }

      /// move-assign operator
      VectorMirror& operator=(VectorMirror&& other)
      {
        if(this == &other)
          return *this;
        _mirror_gather = std::move(other._mirror_gather);
        _mirror_scatter = std::move(other._mirror_scatter);
        return *this;
      }

      /**
       * \brief Clone operation
       *
       * \returns A deep copy of this vector mirror.
       */
      VectorMirror clone() const
      {
        return VectorMirror(_mirror_gather.clone(), _mirror_scatter.clone());
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other
       * The source mirror.
       */
      template<typename Mem2_, typename DT2_, typename IT2_>
      void convert(const VectorMirror<Mem2_, DT2_, IT2_>& other)
      {
        this->_mirror_gather.convert(other._mirror_gather);
        this->_mirror_scatter.convert(other._mirror_scatter);
      }

      /// \cond internal
      const MirrorMatrixType& get_gather_prim() const
      {
        return _mirror_gather;
      }

      const MirrorMatrixType& get_gather_dual() const
      {
        return _mirror_gather;
      }

      const MirrorMatrixType& get_scatter_prim() const
      {
        return _mirror_scatter;
      }

      const MirrorMatrixType& get_scatter_dual() const
      {
        return _mirror_scatter;
      }

      MirrorMatrixType& get_gather_prim()
      {
        return _mirror_gather;
      }

      MirrorMatrixType& get_gather_dual()
      {
        return _mirror_gather;
      }

      MirrorMatrixType& get_scatter_prim()
      {
        return _mirror_scatter;
      }

      MirrorMatrixType& get_scatter_dual()
      {
        return _mirror_scatter;
      }
      /// \endcond

      /**
       * \brief Returns the number of entries in the mirror.
       */
      Index size() const
      {
        return _mirror_gather.rows();
      }

      /**
       * \brief Creates a new buffer vector.
       */
      DenseVector<MemType, DataType, IndexType> create_buffer_vector() const
      {
        return DenseVector<MemType, DataType, IndexType>(size());
      }

      /**
       * \brief Creates a new (local) vector.
       */
      VectorType create_vector() const
      {
        return _mirror_gather.create_vector_r();
      }

      /**
       * \brief Performs a gather-operation on a primal vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A primal vector whose entries are to be gathered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_prim(
                       LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                       const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& vector,
                       const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        //temp-->
        SparseMatrixCSR<Mem::Main, DataType_, IndexType_> tmp_mirror_gather;
        tmp_mirror_gather.convert(_mirror_gather);
        //<--temp

        Tx_ * x(buffer.elements());
        const Ty_ * y(vector.elements());
        const Index * col_idx(tmp_mirror_gather.col_ind());
        const DataType_* val(tmp_mirror_gather.val());
        const Index * row_ptr(tmp_mirror_gather.row_ptr());
        Index num_rows(tmp_mirror_gather.rows());

        ASSERT_(num_rows + buffer_offset <= buffer.size());

        Arch::GatherPrim<Mem::Main>::dv_csr(x, y, col_idx, val, row_ptr, num_rows, buffer_offset);

      }

      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_prim(
                       LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                       const LAFEM::DenseVector<Mem::CUDA, Ty_, Iy_>& cuda_vector,
                       const Index buffer_offset = Index(0)) const
      {
        if(_mirror_gather.empty())
          return;

        LAFEM::DenseVector<Mem::CUDA, Tx_, Ix_> cuda_buffer(buffer.size());

        Tx_ * x(cuda_buffer.elements());
        const Ty_ * y(cuda_vector.elements());
        const Index * col_idx(_mirror_gather.col_ind());
        const DataType_* val(_mirror_gather.val());
        const Index * row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

        ASSERT_(num_rows + buffer_offset <= cuda_buffer.size());

        Arch::GatherPrim<Mem::CUDA>::dv_csr(x, y, col_idx, val, row_ptr, num_rows, buffer_offset);

        //download
        buffer.convert(cuda_buffer);
      }

      /**
       * \brief Performs a gather-axpy-operation on a primal vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A primal vector whose entries are to be gathered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_axpy_prim(
                            LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                            const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& vector,
                            const Tx_ alpha = Tx_(1),
                            const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        //temp-->
        SparseMatrixCSR<Mem::Main, DataType_, IndexType_> tmp_mirror_gather;
        tmp_mirror_gather.convert(_mirror_gather);
        //<--temp

        Tx_ * x(buffer.elements());
        const Ty_ * y(vector.elements());
        const Index * col_idx(tmp_mirror_gather.col_ind());
        const DataType_* val(tmp_mirror_gather.val());
        const Index * row_ptr(tmp_mirror_gather.row_ptr());
        Index num_rows(tmp_mirror_gather.rows());

        ASSERT_(num_rows + buffer_offset <= buffer.size());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          Tx_ sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            sum += Tx_(val[i]) * Tx_(y[col_idx[i]]);
          }
          x[buffer_offset + row] += alpha*sum;
        }
      }

      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_axpy_prim(
                            LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                            const LAFEM::DenseVector<Mem::CUDA, Ty_, Iy_>& cuda_vector,
                            const Tx_ alpha = Tx_(1),
                            const Index buffer_offset = Index(0)) const
      {
        ///TODO
        LAFEM::DenseVector<Mem::Main, Ty_, Iy_> vector;
        vector.convert(cuda_vector);

        gather_axpy_prim(buffer, vector, alpha, buffer_offset);
      }

      /**
       * \brief Performs a scatter-operation on a primal vector.
       *
       * \param[in,out] vector
       * A reference to a primal vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_prim(
                        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& vector,
                        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
                        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        //temp-->
        SparseMatrixCSR<Mem::Main, DataType_, IndexType_> tmp_mirror_scatter;
        tmp_mirror_scatter.convert(_mirror_scatter);
        //<--temp

        Tx_ * x(vector.elements());
        const Ty_ * y(buffer.elements());
        const Index * col_idx(tmp_mirror_scatter.col_ind());
        const DataType_* val(tmp_mirror_scatter.val());
        const Index * row_ptr(tmp_mirror_scatter.row_ptr());
        const Index num_rows(tmp_mirror_scatter.rows());
        const Index num_cols(tmp_mirror_scatter.columns());

        ASSERT_(num_cols + buffer_offset <= buffer.size());
#ifndef DEBUG
        (void)num_cols;
#endif

        Arch::ScatterPrim<Mem::Main>::dv_csr(x, y, col_idx, val, row_ptr, num_rows, buffer_offset);

        // loop over all scatter-matrix rows
        /*for (Index row(0) ; row < num_rows ; ++row)
        {
          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          Tx_ sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            sum += Tx_(val[i]) * Tx_(y[buffer_offset + col_idx[i]]);
          }
          x[row] = sum;
        }*/
      }

      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_prim(
                        LAFEM::DenseVector<Mem::CUDA, Tx_, Ix_>& cuda_vector,
                        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
                        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        LAFEM::DenseVector<Mem::CUDA, Tx_, Ix_> cuda_buffer;
        cuda_buffer.copy(buffer);

        Tx_ * x(cuda_vector.elements());
        const Ty_ * y(cuda_buffer.elements());
        const Index * col_idx(_mirror_scatter.col_ind());
        const DataType_* val(_mirror_scatter.val());
        const Index * row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(num_cols + buffer_offset <= buffer.size());

        Arch::ScatterPrim<Mem::CUDA>::dv_csr(x, y, col_idx, val, row_ptr, num_rows, buffer_offset);
      }

      /**
       * \brief Performs a scatter-axpy-operation on a primal vector.
       *
       * \param[in,out] vector
       * A reference to a primal vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_axpy_prim(
                             LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& vector,
                             const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
                             const Tx_ alpha = Tx_(1),
                             const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        //temp-->
        SparseMatrixCSR<Mem::Main, DataType_, IndexType_> tmp_mirror_scatter;
        tmp_mirror_scatter.convert(_mirror_scatter);
        //<--temp

        Tx_ * x(vector.elements());
        const Ty_ * y(buffer.elements());
        const Index * col_idx(tmp_mirror_scatter.col_ind());
        const DataType_* val(tmp_mirror_scatter.val());
        const Index * row_ptr(tmp_mirror_scatter.row_ptr());
        const Index num_rows(tmp_mirror_scatter.rows());
        const Index num_cols(tmp_mirror_scatter.columns());

        ASSERT_(num_cols + buffer_offset <= buffer.size());
#ifndef DEBUG
        (void)num_cols;
#endif

        Arch::ScatterAxpyPrim<Mem::Main>::dv_csr(x, y, col_idx, val, row_ptr, alpha, num_rows, buffer_offset);
      }

      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_axpy_prim(
                             LAFEM::DenseVector<Mem::CUDA, Tx_, Ix_>& cuda_vector,
                             const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
                             const Tx_ alpha = Tx_(1),
                             const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        LAFEM::DenseVector<Mem::CUDA, Tx_, Ix_> cuda_buffer;
        cuda_buffer.copy(buffer);

        Tx_ * x(cuda_vector.elements());
        const Ty_ * y(cuda_buffer.elements());
        const Index * col_idx(_mirror_scatter.col_ind());
        const DataType_* val(_mirror_scatter.val());
        const Index * row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(num_cols + buffer_offset <= buffer.size());

        Arch::ScatterAxpyPrim<Mem::Main>::dv_csr(x, y, col_idx, val, row_ptr, alpha, num_rows, buffer_offset);
      }

      /**
       * \brief Performs a gather-operation on a dual vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A dual vector whose entries are to be gathered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<
        typename Mx_,
        typename Tx_,
        typename Ix_,
        typename My_,
        typename Ty_,
        typename Iy_>
      void gather_dual(
                       LAFEM::DenseVector<Mx_, Tx_, Ix_>& buffer,
                       const LAFEM::DenseVector<My_, Ty_, Iy_>& vector,
                       const Index buffer_offset = Index(0)) const
      {
        this->gather_prim<Tx_, Ix_, Ty_, Iy_>(buffer, vector, buffer_offset);
      }

      /**
       * \brief Performs a gather-axpy-operation on a dual vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A dual vector whose entries are to be gathered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<
        typename Mx_,
        typename Tx_,
        typename Ix_,
        typename My_,
        typename Ty_,
        typename Iy_>
      void gather_axpy_dual(
                            LAFEM::DenseVector<Mx_, Tx_, Ix_>& buffer,
                            const LAFEM::DenseVector<My_, Ty_, Iy_>& vector,
                            const Tx_ alpha = Tx_(1),
                            const Index buffer_offset = Index(0)) const
      {
        this->gather_axpy_prim<Tx_, Ix_, Ty_, Iy_>(buffer, vector, alpha, buffer_offset);
      }

      /**
       * \brief Performs a scatter-operation on a dual vector.
       *
       * \param[in,out] vector
       * A reference to a dual vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<
        typename Mx_,
        typename Tx_,
        typename Ix_,
        typename My_,
        typename Ty_,
        typename Iy_>
      void scatter_dual(
                        LAFEM::DenseVector<Mx_, Tx_, Ix_>& vector,
                        const LAFEM::DenseVector<My_, Ty_, Iy_>& buffer,
                        const Index buffer_offset = Index(0)) const
      {
        this->scatter_prim<Tx_, Ix_, Ty_, Iy_>(vector, buffer, buffer_offset);
      }

      /**
       * \brief Performs a scatter-axpy-operation on a dual vector.
       *
       * \param[in,out] vector
       * A reference to a dual vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<
        typename Mx_,
        typename Tx_,
        typename Ix_,
        typename My_,
        typename Ty_,
        typename Iy_>
      void scatter_axpy_dual(
                             LAFEM::DenseVector<Mx_, Tx_, Ix_>& vector,
                             const LAFEM::DenseVector<My_, Ty_, Iy_>& buffer,
                             const Tx_ alpha = Tx_(1),
                             const Index buffer_offset = Index(0)) const
      {
        this->scatter_axpy_prim<Tx_, Ix_, Ty_, Iy_>(vector, buffer, alpha, buffer_offset);
      }

      /// \copydoc gather_prim()
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_prim(
                       LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                       const LAFEM::SparseVector<Mem::Main, Ty_, Iy_>& vector,
                       const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        Tx_ * x(buffer.elements());
        const Ty_ * y(vector.elements());
        const IndexType_ * col_idx(_mirror_gather.col_ind());
        const DataType_* val(_mirror_gather.val());
        const IndexType_ * row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

        ASSERT_(num_rows + buffer_offset <= buffer.size());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          // Because the index array of the SparseVector AND the row_ptr are sorted, we remember the last position
          // in the SparseVector which contained a required entry OR told us that the required entry is not there.
          Index istart(0);

          Tx_ sum(Tx_(0));
          // Loop over all columns. For every column, we need to find the corresponding entry in the SparseVector
          // (if any).
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            for (Index isparse(istart); isparse < vector.used_elements(); ++isparse)
            {
              // This is only safe because vector is in Mem::Main
              Index iv(vector.indices()[isparse]);

              if(iv < col_idx[i])
                continue;

              if(iv == col_idx[i])
                sum += Tx_(val[i]) * Tx_(y[isparse]);

              // If we come to here, we know we found something, so set istart and break
              istart = isparse;
              break;
            }
          }
          x[buffer_offset + row] = sum;
        }
      }

      /// \copydoc gather_axpy_prim()
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_axpy_prim(
                            LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                            const LAFEM::SparseVector<Mem::Main, Ty_, Iy_>& vector,
                            const Tx_ alpha = Tx_(1),
                            const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        Tx_ * x(buffer.elements());
        const Ty_ * y(vector.elements());
        const IndexType_ * col_idx(_mirror_gather.col_ind());
        const DataType_* val(_mirror_gather.val());
        const IndexType_* row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

        ASSERT_(num_rows + buffer_offset <= buffer.size());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          // Because the index array of the SparseVector AND the row_ptr are sorted, we remember the last position
          // in the SparseVector which contained a required entry OR told us that the required entry is not there.
          Index istart(0);

          Tx_ sum(Tx_(0));
          // Loop over all columns. For every column, we need to find the corresponding entry in the SparseVector
          // (if any).
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            for (Index isparse(istart); isparse < vector.used_elements(); ++isparse)
            {
              // This is only safe because vector is in Mem::Main
              Index iv(vector.indices()[isparse]);

              if(iv < col_idx[i])
                continue;

              if(iv == col_idx[i])
                sum += Tx_(val[i]) * Tx_(y[isparse]);

              // If we come to here, we know we found something, so set istart and break
              istart = isparse;
              break;
            }

          }
          x[buffer_offset + row] += alpha*sum;
        }
      }

      /// \copydoc scatter_prim()
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_prim(
                        LAFEM::SparseVector<Mem::Main, Tx_, Ix_>& vector,
                        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
                        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        Tx_ * x(vector.elements());
        const Ty_ * y(buffer.elements());
        const IndexType_ * col_idx(_mirror_scatter.col_ind());
        const DataType_* val(_mirror_scatter.val());
        const IndexType_ * row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(num_cols + buffer_offset <= buffer.size());
        ASSERT_(num_rows >= vector.indices()[vector.used_elements()-1]);
#ifndef DEBUG
        (void)num_cols;
        (void)num_rows;
#endif

        for (Index isparse(0); isparse < vector.used_elements(); ++isparse)
        {
          Index row(vector.indices()[isparse]);
          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          Tx_ sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            sum += Tx_(val[i]) * Tx_(y[buffer_offset + col_idx[i]]);
          }
          x[isparse] = sum;
        }
      }

      /// \copydoc scatter_axpy_prim()
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_axpy_prim(
                             LAFEM::SparseVector<Mem::Main, Tx_, Ix_>& vector,
                             const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
                             const Tx_ alpha = Tx_(1),
                             const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        Tx_ * x(vector.elements());
        const Ty_ * y(buffer.elements());
        const IndexType_ * col_idx(_mirror_scatter.col_ind());
        const DataType_* val(_mirror_scatter.val());
        const IndexType_ * row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(num_cols + buffer_offset <= buffer.size());
        ASSERT_(num_rows >= vector.indices()[vector.used_elements()-1]);
#ifndef DEBUG
        (void)num_cols;
        (void)num_rows;
#endif

        for (Index isparse(0); isparse < vector.used_elements(); ++isparse)
        {
          Index row(vector.indices()[isparse]);

          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          Tx_ sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            sum += Tx_(val[i]) * Tx_(y[buffer_offset + col_idx[i]]);
          }
          x[isparse] += alpha*sum;
        }
      }

      /// \copydoc gather_dual()
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_dual(
                       LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                       const LAFEM::SparseVector<Mem::Main, Ty_, Iy_>& vector,
                       const Index buffer_offset = Index(0)) const
      {
        this->gather_prim<Tx_, Ix_, Ty_, Iy_>(buffer, vector, buffer_offset);
      }

      /// \copydoc gather_axpy_dual()
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_axpy_dual(
                            LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                            const LAFEM::SparseVector<Mem::Main, Ty_, Iy_>& vector,
                            const Tx_ alpha = Tx_(1),
                            const Index buffer_offset = Index(0)) const
      {
        this->gather_axpy_prim<Tx_, Ix_, Ty_, Iy_>(buffer, vector, alpha, buffer_offset);
      }

      /// \copydoc scatter_dual()
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_dual(
                        LAFEM::SparseVector<Mem::Main, Tx_, Ix_>& vector,
                        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
                        const Index buffer_offset = Index(0)) const
      {
        this->scatter_prim<Tx_, Ix_, Ty_, Iy_>(vector, buffer, buffer_offset);
      }

      /// \copydoc scatter_axpy_dual()
      template<
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_axpy_dual(
                             LAFEM::SparseVector<Mem::Main, Tx_, Ix_>& vector,
                             const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
                             const Tx_ alpha = Tx_(1),
                             const Index buffer_offset = Index(0)) const
      {
        this->scatter_axpy_prim<Tx_, Ix_, Ty_, Iy_>(vector, buffer, alpha, buffer_offset);
      }

    }; // class VectorMirror<...>

    /**
     * \brief Handles DenseVectorBlocked prolongation, restriction and serialisation
     *
     * \tparam Mem_
     * Memory architecture
     *
     * \tparam DataType_
     * Data type in the vector
     *
     * \tparam IndexType_
     * Type for indexing
     *
     * \tparam BlockSize_
     * Size of the DenseVectorBlocked blocks.
     *
     * \see VectorMirror for more details on Mirrors in general.
     *
     * \author Jordi Paul
     */
    template<typename Mem_, typename DT_, typename IT_, int BlockSize_>
    class VectorMirrorBlocked
    {
    public:
      /// arch typedef
      typedef Mem_ MemType;
      /// data-type typedef
      typedef DT_ DataType;
      /// index-type typedef
      typedef IT_ IndexType;
      /// Our block size
      static constexpr int BlockSize = BlockSize_;

      /// mirror matrix typedef
      typedef SparseMatrixCSR<Mem_, DT_, IT_> MirrorMatrixType;

      /// corresponding vector
      typedef DenseVector<Mem_, DT_, IT_> BufferVectorType;

    protected:
      /// gather-mirror matrix
      MirrorMatrixType _mirror_gather;
      /// scatter-mirror matrix
      MirrorMatrixType _mirror_scatter;

    public:
      /// default constructor
      VectorMirrorBlocked() :
        _mirror_gather(),
        _mirror_scatter()
      {
      }

      /// explicit constructor
      explicit VectorMirrorBlocked(MirrorMatrixType&& mirror_gather, MirrorMatrixType&& mirror_scatter) :
        _mirror_gather(std::move(mirror_gather)),
        _mirror_scatter(std::move(mirror_scatter))
      {
      }

      /// move-ctor
      VectorMirrorBlocked(VectorMirrorBlocked&& other) :
        _mirror_gather(std::move(other._mirror_gather)),
        _mirror_scatter(std::move(other._mirror_scatter))
      {
      }

      /// move-assign operator
      VectorMirrorBlocked& operator=(VectorMirrorBlocked&& other)
      {
        if(this == &other)
          return *this;
        _mirror_gather = std::move(other._mirror_gather);
        _mirror_scatter = std::move(other._mirror_scatter);
        return *this;
      }

      /**
       * \brief Clone operation
       *
       * \returns A deep copy of this vector mirror.
       */
      VectorMirrorBlocked clone() const
      {
        return VectorMirrorBlocked(_mirror_gather.clone(), _mirror_scatter.clone());
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other
       * The source mirror.
       */
      template<typename Mem2_, typename DT2_, typename IT2_>
      void convert(const VectorMirrorBlocked<Mem2_, DT2_, IT2_, BlockSize>& other)
      {
        this->_mirror_gather.convert(other._mirror_gather);
        this->_mirror_scatter.convert(other._mirror_scatter);
      }

      /// \cond internal
      const MirrorMatrixType& get_gather_prim() const
      {
        return _mirror_gather;
      }

      const MirrorMatrixType& get_gather_dual() const
      {
        return _mirror_gather;
      }

      const MirrorMatrixType& get_scatter_prim() const
      {
        return _mirror_scatter;
      }

      const MirrorMatrixType& get_scatter_dual() const
      {
        return _mirror_scatter;
      }

      MirrorMatrixType& get_gather_prim()
      {
        return _mirror_gather;
      }

      MirrorMatrixType& get_gather_dual()
      {
        return _mirror_gather;
      }

      MirrorMatrixType& get_scatter_prim()
      {
        return _mirror_scatter;
      }

      MirrorMatrixType& get_scatter_dual()
      {
        return _mirror_scatter;
      }
      /// \endcond

      /**
       * \brief Returns the number of entries in the mirror.
       */
      Index size() const
      {
        return _mirror_gather.rows();
      }

      /**
       * \brief Creates a new buffer vector.
       */
       BufferVectorType create_buffer_vector() const
      {
        return BufferVectorType(Index(BlockSize)*_mirror_gather.rows());
      }

      /**
       * \brief Creates a new (local) vector.
       */
      DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_>create_vector() const
      {
        return DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_>(_mirror_gather.columns());
      }

      /**
       * \brief Performs a gather-operation on a primal vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A primal vector whose entries are to be gathered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template< typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void gather_prim(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::DenseVectorBlocked<Mem::Main, Ty_, Iy_, BS_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        Tx_* x(buffer.elements());
        const Ty_ * y(vector.template elements<Perspective::pod>());
        const IndexType* col_idx(_mirror_gather.col_ind());
        const DataType* val(_mirror_gather.val());
        const IndexType* row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

        ASSERT_(Index(BS_)*num_rows + buffer_offset <= buffer.size());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          for(Index j(0); j < Index(BS_); ++j)
          {
            Tx_ sum(Tx_(0));
            for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
            {
              sum += Tx_(val[i]) * Tx_(y[Index(BS_)*col_idx[i] +j]);
            }
            x[buffer_offset + Index(BS_)*row + j] = sum;
          }
        }
      }

      /**
       * \brief Performs a gather-axpy-operation on a primal vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A primal vector whose entries are to be gathered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_ >
      void gather_axpy_prim(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::DenseVectorBlocked<Mem::Main, Ty_, Iy_, BS_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        Tx_* x(buffer.elements());
        const Ty_* y(vector.template elements<Perspective::pod>());
        const IndexType* col_idx(_mirror_gather.col_ind());
        const DataType* val(_mirror_gather.val());
        const IndexType* row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

        ASSERT_(Index(BS_)*num_rows + buffer_offset <= buffer.size());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          for (Index j(0); j < Index(BS_); ++j)
          {
            Tx_ sum(Tx_(0));
            for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
            {
              sum += Tx_(val[i]) * Tx_(y[Index(BS_)*col_idx[i] + j]);
            }
            x[buffer_offset + Index(BS_)*row + j] += alpha*sum;
          }
        }
      }

      /**
       * \brief Performs a scatter-operation on a primal vector.
       *
       * \param[in,out] vector
       * A reference to a primal vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_ >
      void scatter_prim(
        LAFEM::DenseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        Tx_* x(vector.template elements<Perspective::pod>());
        const Ty_* y(buffer.elements());
        const IndexType* col_idx(_mirror_scatter.col_ind());
        const DataType* val(_mirror_scatter.val());
        const IndexType* row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(Index(BS_)*num_cols + buffer_offset <= buffer.size());
#ifndef DEBUG
        (void)num_cols;
#endif

        // loop over all scatter-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          for(Index j(0); j < Index(BS_); ++j)
          {
            Tx_ sum(Tx_(0));
            for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
            {
              sum += Tx_(val[i]) * Tx_(y[buffer_offset + Index(BS_)*col_idx[i] + j]);
            }
            x[Index(BS_)*row + j] = sum;
          }
        }
      }

      /**
       * \brief Performs a scatter-axpy-operation on a primal vector.
       *
       * \param[in,out] vector
       * A reference to a primal vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template< typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_ >
      void scatter_axpy_prim(
        LAFEM::DenseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        Tx_* x(vector.template elements<Perspective::pod>());
        const Ty_* y(buffer.elements());
        const IndexType* col_idx(_mirror_scatter.col_ind());
        const DataType* val(_mirror_scatter.val());
        const IndexType* row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(Index(BS_)*num_cols + buffer_offset <= buffer.size());
#ifndef DEBUG
        (void)num_cols;
#endif

        // loop over all scatter-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          for(Index j(0); j < Index(BS_); ++j)
          {
            Tx_ sum(Tx_(0));
            for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
              sum += Tx_(val[i]) * Tx_(y[buffer_offset + Index(BS_)*col_idx[i] + j]);

            x[Index(BS_)*row + j] += alpha*sum;
          }
        }
      }

      /**
       * \brief Performs a gather-operation on a dual vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A dual vector whose entries are to be gathered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void gather_dual(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::DenseVectorBlocked<Mem::Main, Ty_, Iy_, BS_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        this->gather_prim<Tx_, Ix_, Ty_, Iy_, BS_>(buffer, vector, buffer_offset);
      }

      /**
       * \brief Performs a gather-axpy-operation on a dual vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A dual vector whose entries are to be gathered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void gather_axpy_dual(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::DenseVectorBlocked<Mem::Main, Ty_, Iy_, BS_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        this->gather_axpy_prim<Tx_, Ix_, Ty_, Iy_, BS_>(buffer, vector, alpha, buffer_offset);
      }

      /**
       * \brief Performs a scatter-operation on a dual vector.
       *
       * \param[in,out] vector
       * A reference to a dual vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void scatter_dual(
        LAFEM::DenseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        this->scatter_prim<Tx_, Ix_, Ty_, Iy_, BS_>(vector, buffer, buffer_offset);
      }

      /**
       * \brief Performs a scatter-axpy-operation on a dual vector.
       *
       * \param[in,out] vector
       * A reference to a dual vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template< typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void scatter_axpy_dual(
        LAFEM::DenseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        this->scatter_axpy_prim<Tx_, Ix_, Ty_, Iy_, BS_>(vector, buffer, alpha, buffer_offset);
      }

      /// \copydoc gather_prim()
      template< typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void gather_prim(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::SparseVectorBlocked<Mem::Main, Ty_, Iy_, BS_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        typedef typename LAFEM::SparseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>::ValueType ValueType;
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        Tx_* x(buffer.elements());
        const Ty_* y(vector.template elements<Perspective::pod>());
        const IndexType* col_idx(_mirror_gather.col_ind());
        const DataType* val(_mirror_gather.val());
        const IndexType* row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

        ASSERT_(Index(BS_)*num_rows + buffer_offset <= buffer.size());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          // Because the index array of the SparseVector AND the row_ptr are sorted, we remember the last position
          // in the SparseVector which contained a required entry OR told us that the required entry is not there.
          Index istart(0);

          ValueType sum(Tx_(0));
          // Loop over all columns. For every column, we need to find the corresponding entry in the SparseVector
          // (if any).
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            for(Index isparse(istart); isparse < vector.used_elements(); ++isparse)
            {
              // This is only safe because vector is in Mem::Main
              Index iv(vector.indices()[isparse]);

              if(iv < col_idx[i])
                continue;

              if(iv == col_idx[i])
              {
                // Address sum explicitly here so we can cast the elements of y to Tx_
                for(int j(0); j < BS_; ++j)
                  sum(j) += Tx_(val[i]) * Tx_(y[Index(BS_)*isparse + Index(j)]);
              }

              // If we come to here, we know we found something, so set istart and break
              istart = isparse;
              break;
            }
          }

          for(int j(0); j < BS_; ++j)
            x[buffer_offset + Ix_(BS_)*row + Ix_(j)] = Tx_(sum(j));
        }
      }

      /// \copydoc gather_axpy_prim()
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_ >
      void gather_axpy_prim(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::SparseVectorBlocked<Mem::Main, Ty_, Iy_, BS_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        typedef typename LAFEM::SparseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>::ValueType ValueType;

        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        Tx_* x(buffer.elements());
        const Ty_* y(vector.template elements<Perspective::pod>());
        const IndexType* col_idx(_mirror_gather.col_ind());
        const DataType* val(_mirror_gather.val());
        const IndexType* row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

        ASSERT_(Index(BS_)*num_rows + buffer_offset <= buffer.size());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          // Because the index array of the SparseVector AND the row_ptr are sorted, we remember the last position
          // in the SparseVector which contained a required entry OR told us that the required entry is not there.
          Index istart(0);

          ValueType sum(Tx_(0));
          // Loop over all columns. For every column, we need to find the corresponding entry in the SparseVector
          // (if any).
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            for (Index isparse(istart); isparse < vector.used_elements(); ++isparse)
            {
              // This is only safe because vector is in Mem::Main
              Index iv(vector.indices()[isparse]);

              if(iv < col_idx[i])
                continue;

              if(iv == col_idx[i])
              {
                // Address sum explicitly here so we can cast the elements of y to Tx_
                for(int j(0); j < BS_; ++j)
                  sum(j) += Tx_(val[i]) * Tx_(y[Index(BS_)*isparse + Index(j)]);
              }

              // If we come to here, we know we found something, so set istart and break
              istart = isparse;
              break;
            }

            for (int j(0); j < BS_; ++j)
              x[buffer_offset + Index(BS_)*row + Index(j)] += alpha*Tx_(sum(j));
          }
        }
      }

      /// \copydoc scatter_prim
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_ >
      void scatter_prim(
        LAFEM::SparseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        typedef typename LAFEM::SparseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>::ValueType ValueType;
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        auto* x = vector.template elements<Perspective::native>();
        const Ty_* y(buffer.elements());
        const IndexType* col_idx(_mirror_scatter.col_ind());
        const DataType* val(_mirror_scatter.val());
        const IndexType* row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(Index(BS_)*num_cols + buffer_offset <= buffer.size());
        ASSERT_(num_rows >= vector.indices()[vector.used_elements()-1]);
#ifndef DEBUG
        (void)num_cols;
        (void)num_rows;
#endif

        for (Index isparse(0) ; isparse < vector.used_elements(); ++isparse)
        {
          Index row(vector.indices()[isparse]);

          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          ValueType sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            for(int j(0); j < BS_; ++j)
              sum(j) += Tx_(val[i]) * Tx_(y[buffer_offset + Iy_(BS_)*col_idx[i] + Iy_(j)]);
          }

          x[isparse] = sum;
        }
      }

      /// \copydoc scatter_axpy_prim
      template< typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_ >
      void scatter_axpy_prim(
        LAFEM::SparseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        typedef typename LAFEM::SparseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>::ValueType ValueType;
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        ValueType* x = vector.template elements<Perspective::native>();
        const Ty_* y(buffer.elements());
        const IndexType* col_idx(_mirror_scatter.col_ind());
        const DataType* val(_mirror_scatter.val());
        const IndexType* row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(Index(BS_)*num_cols + buffer_offset <= buffer.size());
        ASSERT_(num_rows >= vector.indices()[vector.used_elements()-1]);
#ifndef DEBUG
        (void)num_cols;
        (void)num_rows;
#endif

        for(Index isparse(0) ; isparse < vector.used_elements() ; ++isparse)
        {
          Index row(vector.indices()[isparse]);

          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          ValueType sum(Tx_(0));
          for(int j(0); j < BS_; ++j)
          {
            for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
              sum(j) += Tx_(val[i]) * Tx_(y[buffer_offset + Index(BS_)*col_idx[i] + Index(j)]);
          }

          x[isparse] += alpha*sum;
        }
      }

      /// \copydoc gather_dual()
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void gather_dual(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::SparseVectorBlocked<Mem::Main, Ty_, Iy_, BS_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        this->gather_prim<Tx_, Ix_, Ty_, Iy_, BS_>(buffer, vector, buffer_offset);
      }

      /// \copydoc gather_axpy_dual()
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void gather_axpy_dual(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::SparseVectorBlocked<Mem::Main, Ty_, Iy_, BS_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        this->gather_axpy_prim<Tx_, Ix_, Ty_, Iy_, BS_>(buffer, vector, alpha, buffer_offset);
      }

      /// \copydoc scatter_dual()
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void scatter_dual(
        LAFEM::SparseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        this->scatter_prim<Tx_, Ix_, Ty_, Iy_, BS_>(vector, buffer, buffer_offset);
      }

      /// \copydoc scatter_axpy_dual()
      template< typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void scatter_axpy_dual(
        LAFEM::SparseVectorBlocked<Mem::Main, Tx_, Ix_, BS_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        this->scatter_axpy_prim<Tx_, Ix_, Ty_, Iy_, BS_>(vector, buffer, alpha, buffer_offset);
      }
    }; // class VectorMirrorBlocked<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_VECTOR_MIRROR_HPP
