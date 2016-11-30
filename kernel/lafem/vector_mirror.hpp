#pragma once
#ifndef KERNEL_LAFEM_VECTOR_MIRROR_HPP
#define KERNEL_LAFEM_VECTOR_MIRROR_HPP 1

// includes, FEAT
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_vector.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/arch/scatter_axpy_prim.hpp>

namespace FEAT
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
     * \author Jordi Paul
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

      /// gather-mirror matrix
      MirrorMatrixType _mirror_gather;
      /// scatter-mirror matrix
      MirrorMatrixType _mirror_scatter;

      /// Our 'base' class type
      template <typename Mem2_, typename DataType2_ = DataType_, typename IndexType2_ = IndexType_>
      using MirrorType = class VectorMirror<Mem2_, DataType2_, IndexType2_>;

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
      VectorMirror clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        return VectorMirror(_mirror_gather.clone(clone_mode), _mirror_scatter.clone(clone_mode));
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _mirror_gather.bytes() + _mirror_scatter.bytes();
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
      const MirrorMatrixType& get_gather() const
      {
        return _mirror_gather;
      }

      const MirrorMatrixType& get_scatter() const
      {
        return _mirror_scatter;
      }

      MirrorMatrixType& get_gather()
      {
        return _mirror_gather;
      }

      MirrorMatrixType& get_scatter()
      {
        return _mirror_scatter;
      }
      /// \endcond

      /**
       * \brief Computes the required buffer size for a DenseVector.
       *
       * \tparam[in] vector
       * The vector whose buffer size is to be computed.
       */
      template<typename Mem2_, typename DT2_, typename IT2_>
      Index buffer_size(const DenseVector<Mem2_, DT2_, IT2_>& DOXY(vector)) const
      {
        return _mirror_gather.rows();
      }

      /**
       * \brief Computes the required buffer size for a DenseVectorBlocked.
       *
       * \tparam[in] vector
       * The vector whose buffer size is to be computed.
       */
      template<typename Mem2_, typename DT2_, typename IT2_, int block_size_>
      Index buffer_size(const DenseVectorBlocked<Mem2_, DT2_, IT2_, block_size_>& DOXY(vector)) const
      {
        return _mirror_gather.rows()*Index(block_size_);
      }

      /**
       * \brief Computes the required buffer size for a SparseVector.
       *
       * \tparam[in] vector
       * The vector whose buffer size is to be computed.
       */
      template<typename Mem2_, typename DT2_, typename IT2_>
      Index buffer_size(const SparseVector<Mem2_, DT2_, IT2_>& DOXY(vector)) const
      {
        return _mirror_gather.rows();
      }

      /**
       * \brief Computes the required buffer size for a SparseVectorBlocked.
       *
       * \tparam[in] vector
       * The vector whose buffer size is to be computed.
       */
      template<typename Mem2_, typename DT2_, typename IT2_, int block_size_>
      Index buffer_size(const SparseVectorBlocked<Mem2_, DT2_, IT2_, block_size_>& DOXY(vector)) const
      {
        return _mirror_gather.rows()*Index(block_size_);
      }

      /**
       * \brief Creates a new buffer vector for a vector.
       *
       * \tparam[in] vector
       * The vector for which the buffer is to be created.
       */
      template<typename Vector_>
      DenseVector<Mem::Main, DataType, IndexType> create_buffer(const Vector_& vector) const
      {
        return DenseVector<Mem::Main, DataType, IndexType>(buffer_size(vector), Pinning::disabled);
      }

      /**
       * \brief Gathers the buffer entries from a DenseVector.
       *
       * \param[in,out] buffer
       * A reference to the buffer.
       *
       * \param[in] vector
       * A vector whose entries are to be gathered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer.
       */
      void gather(
        LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>& buffer,
        const LAFEM::DenseVector<Mem_, DataType_, IndexType_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        Index num_rows(_mirror_gather.rows());
        XASSERTM(num_rows + buffer_offset <= buffer.size(), "buffer vector size mismatch");
        LAFEM::DenseVector<Mem_, DataType_, IndexType_> mem_buffer(num_rows);

        Arch::Apply<Mem_>::csr(mem_buffer.elements(), DataType_(1), vector.elements(), DataType_(0), mem_buffer.elements(), _mirror_gather.val(),
            _mirror_gather.col_ind(), _mirror_gather.row_ptr(), _mirror_gather.rows(), _mirror_gather.columns(), _mirror_gather.used_elements(), false);

        //download
        DenseVector<Mem::Main, DataType_, IndexType_> buffer_range(buffer, num_rows, buffer_offset);
        buffer_range.copy(mem_buffer);
      }

      /**
       * \brief Scatters the buffer entries onto a DenseVector.
       *
       * \param[in,out] vector
       * A reference to the vector.
       *
       * \param[in] buffer
       * A reference to the buffer whose entries are to be scattered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer.
       */
      void scatter_axpy(
        LAFEM::DenseVector<Mem_, DataType_, IndexType_>& vector,
        const LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>& buffer,
        const DataType_ alpha = DataType_(1),
        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        const Index num_cols(_mirror_scatter.columns());
        XASSERTM(num_cols + buffer_offset <= buffer.size(), "buffer vector size mismatch");

        DenseVector<Mem_, DataType_, IndexType_> mem_buffer(num_cols);
        //download
        DenseVector<Mem::Main, DataType_, IndexType_> buffer_range(buffer, num_cols, buffer_offset);
        mem_buffer.copy(buffer_range);

        DataType_ * x(vector.elements());
        const DataType_ * y(mem_buffer.elements());
        const IndexType_ * col_idx(_mirror_scatter.col_ind());
        const DataType_* val(_mirror_scatter.val());
        const IndexType_ * row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());

        Arch::ScatterAxpyPrim<Mem_>::dv_csr(x, y, col_idx, val, row_ptr, alpha, num_rows);
      }

      /**
       * \brief Gathers the buffer entries from a DenseVectorBlocked.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A vector whose entries are to be gathered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<int block_size_>
      void gather(
        LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>& buffer,
        const LAFEM::DenseVectorBlocked<Mem_, DataType_, IndexType_, block_size_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_gather.empty())
          return;

        const Index data_size = Index(block_size_) * _mirror_gather.rows();
        XASSERTM(data_size + buffer_offset <= buffer.size(), "buffer vector size mismatch");

        LAFEM::DenseVector<Mem_, DataType_, IndexType_> mem_buffer(data_size);

        Arch::Apply<Mem_>::template csrsb<DataType_, IndexType_, block_size_>(
          mem_buffer.elements(), DataType_(1), vector.template elements<Perspective::pod>(),
          DataType_(0), mem_buffer.elements(), _mirror_gather.val(), _mirror_gather.col_ind(),
          _mirror_gather.row_ptr(), _mirror_gather.rows(), _mirror_gather.columns(),
          _mirror_gather.used_elements());

        DenseVector<Mem::Main, DataType_, IndexType_> buffer_range(buffer, data_size, buffer_offset);

        //download
        buffer_range.copy(mem_buffer);
      }

      /**
       * \brief Scatters the buffer entries onto a DenseVectorBlocked.
       *
       * \param[in,out] vector
       * A reference to the vector.
       *
       * \param[in] buffer
       * A reference to the buffer whose entries are to be scattered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer.
       */
      template<int block_size_>
      void scatter_axpy(
        LAFEM::DenseVectorBlocked<Mem_, DataType_, IndexType_, block_size_>& vector,
        const LAFEM::DenseVector<Mem::Main, DataType_, IndexType_>& buffer,
        const DataType_ alpha = DataType_(1),
        const Index buffer_offset = Index(0)) const
      {
        // skip on empty mirror
        if(_mirror_scatter.empty())
          return;

        const Index data_size = Index(block_size_) * _mirror_scatter.columns();
        XASSERTM(data_size + buffer_offset <= buffer.size(), "buffer vector size mismatch");

        DenseVector<Mem_, DataType_, IndexType_> mem_buffer(data_size);
        //download
        DenseVector<Mem::Main, DataType_, IndexType_> buffer_range(buffer, data_size, buffer_offset);
        mem_buffer.copy(buffer_range);

        Arch::ScatterAxpyPrim<Mem_>::template dvb_csr<DataType_, IndexType_, block_size_>(
          vector.template elements<Perspective::pod>(), mem_buffer.elements(),
          _mirror_scatter.col_ind(), _mirror_scatter.val(), _mirror_scatter.row_ptr(),
          alpha, _mirror_scatter.rows());
      }

      /**
       * \brief Gathers the buffer entries from a SparseVector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A vector whose entries are to be gathered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_>
      void gather(
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

        XASSERTM(num_rows + buffer_offset <= buffer.size(), "buffer vector size mismatch");

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

      /**
       * \brief Scatters the buffer entries onto a SparseVector.
       *
       * \param[in,out] vector
       * A reference to the vector.
       *
       * \param[in] buffer
       * A reference to the buffer whose entries are to be scattered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer.
       */
      template<typename Tx_, typename Ix_, typename Ty_, typename Iy_>
      void scatter_axpy(
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

        XASSERTM(num_cols + buffer_offset <= buffer.size(), "buffer vector size mismatch");
        XASSERTM(num_rows >= vector.indices()[vector.used_elements()-1], "vector size mismatch");

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

      /**
       * \brief Gathers the buffer entries from a SparseVector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A vector whose entries are to be gathered.
       *
       * \param[in] buffer_offset
       * The offset within the buffer vector.
       */
      template< typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_>
      void gather(
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
        const Index num_rows(_mirror_gather.rows());

        XASSERTM(Index(BS_)*num_rows + buffer_offset <= buffer.size(), "buffer vector size mismatch");

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

      /**
       * \brief Scatters the buffer entries onto a SparseVectorBlocked.
       *
       * \param[in,out] vector
       * A reference to the vector.
       *
       * \param[in] buffer
       * A reference to the buffer whose entries are to be scattered.
       *
       * \param[in] alpha
       * A scaling factor for the operation.
       *
       * \param[in] buffer_offset
       * The offset within the buffer.
       */
      template< typename Tx_, typename Ix_, typename Ty_, typename Iy_, int BS_ >
      void scatter_axpy(
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

        XASSERTM(Index(BS_)*num_cols + buffer_offset <= buffer.size(), "buffer vector size mismatch");
        XASSERTM(num_rows >= vector.indices()[vector.used_elements()-1], "vector size mismatch");

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
    }; // class VectorMirror<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_VECTOR_MIRROR_HPP
