#pragma once
#ifndef KERNEL_LAFEM_VECTOR_MIRROR_HPP
#define KERNEL_LAFEM_VECTOR_MIRROR_HPP 1

// includes, FEAST
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Vector-Mirror class template
     *
     * \author Peter Zajac
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

    protected:
      /// gather-mirror matrix
      MirrorMatrixType _mirror_gather;
      /// scatter-mirror matrix
      MirrorMatrixType _mirror_scatter;

    public:
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
        typename Algo_,
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_prim(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        Tx_ * x(buffer.elements());
        const Ty_ * y(vector.elements());
        const Index * col_idx(_mirror_gather.col_ind());
        const DataType_* val(_mirror_gather.val());
        const Index * row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

        ASSERT_(num_rows + buffer_offset <= buffer.size());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          Tx_ sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            sum += Tx_(val[i]) * Tx_(y[col_idx[i]]);
          }
          x[buffer_offset + row] = sum;
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
      template<
        typename Algo_,
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
        Tx_ * x(buffer.elements());
        const Ty_ * y(vector.elements());
        const Index * col_idx(_mirror_gather.col_ind());
        const DataType_* val(_mirror_gather.val());
        const Index * row_ptr(_mirror_gather.row_ptr());
        Index num_rows(_mirror_gather.rows());

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
        typename Algo_,
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_prim(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        Tx_ * x(vector.elements());
        const Ty_ * y(buffer.elements());
        const Index * col_idx(_mirror_scatter.col_ind());
        const DataType_* val(_mirror_scatter.val());
        const Index * row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(num_cols + buffer_offset <= buffer.size());

        // loop over all scatter-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
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
      template<
        typename Algo_,
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
        Tx_ * x(vector.elements());
        const Ty_ * y(buffer.elements());
        const Index * col_idx(_mirror_scatter.col_ind());
        const DataType_* val(_mirror_scatter.val());
        const Index * row_ptr(_mirror_scatter.row_ptr());
        const Index num_rows(_mirror_scatter.rows());
        const Index num_cols(_mirror_scatter.columns());

        ASSERT_(num_cols + buffer_offset <= buffer.size());

        // loop over all scatter-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          // skip empty rows
          if(row_ptr[row] >= row_ptr[row + 1])
            continue;

          Tx_ sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_ptr[row + 1] ; ++i)
          {
            sum += Tx_(val[i]) * Tx_(y[buffer_offset + col_idx[i]]);
          }
          x[row] += alpha*sum;
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
      template<
        typename Algo_,
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_dual(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        this->gather_prim<Algo_, Tx_, Ix_, Ty_, Iy_>(buffer, vector, buffer_offset);
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
        typename Algo_,
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void gather_axpy_dual(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        this->gather_axpy_prim<Algo_, Tx_, Ix_, Ty_, Iy_>(buffer, vector, alpha, buffer_offset);
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
        typename Algo_,
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_dual(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        this->scatter_prim<Algo_, Tx_, Ix_, Ty_, Iy_>(vector, buffer, buffer_offset);
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
        typename Algo_,
        typename Tx_,
        typename Ix_,
        typename Ty_,
        typename Iy_>
      void scatter_axpy_dual(
        LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_, Iy_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        this->scatter_axpy_prim<Algo_, Tx_, Ix_, Ty_, Iy_>(vector, buffer, alpha, buffer_offset);
      }
    }; // class VectorMirror<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_VECTOR_MIRROR_HPP
