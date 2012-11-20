#pragma once
#ifndef KERNEL_LAFEM_VECTOR_MIRROR_HPP
#define KERNEL_LAFEM_VECTOR_MIRROR_HPP 1

// includes, FEAST
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/transposition.hpp>

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
      typename Arch_,
      typename DataType_>
    class VectorMirror
    {
    public:
      /// arch typedef
      typedef Arch_ MemType;
      /// data-type typedef
      typedef DataType_ DataType;

      /// mirror matrix typedef
      typedef SparseMatrixCSR<Arch_, DataType_> MirrorMatrixType;

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

      /**
       * \brief Graph constructor
       *
       * This constructor creates a vector mirror based on a dof-mirror adjacency graph.
       */
      explicit VectorMirror(const Graph& gather_graph) :
        _mirror_gather(gather_graph),
        _mirror_scatter(Graph(Graph::rt_transpose, gather_graph))
      {
        _mirror_gather.clear(DataType(1));
        _mirror_scatter.clear(DataType(1));
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
      /// \endcond

      /**
       * \brief Returns the number of entries in the mirror.
       */
      Index get_num_mirror_entries() const
      {
        return _mirror_gather.rows();
      }

      /**
       * \brief Creates a buffer vector based on a template vector.
       *
       * \param[in] tmpl_vec
       * A reference to the template vector.
       *
       * \returns
       * A new vector of the same data-type and arch as the template vector.
       */
      template<
        typename DataType2_>
      LAFEM::DenseVector<Arch_, DataType2_> create_buffer(const LAFEM::DenseVector<Arch_, DataType2_>& tmpl_vec) const
      {
        return LAFEM::DenseVector<Arch_, DataType2_>(get_num_mirror_entries());
      }

      /**
       * \brief Performs a gather-operation on a primal vector.
       *
       * \param[in,out] buffer
       * A reference to a buffer vector.
       *
       * \param[in] vector
       * A primal vector whose entries are to be gathered.
       */
      template<
        typename Tx_,
        typename Ty_>
      void gather_prim(
        LAFEM::DenseVector<Mem::Main, Tx_>& buffer,
        const LAFEM::DenseVector<Mem::Main, Ty_>& vector) const
      {
        Tx_ * x(buffer.elements());
        const Ty_ * y(vector.elements());
        const Index * col_idx(_mirror_gather.col_ind());
        const DataType_* val(_mirror_gather.val());
        const Index * row_ptr(_mirror_gather.row_ptr());
        const Index * row_end(_mirror_gather.row_ptr_end());
        Index num_rows(_mirror_gather.rows());

        // loop over all gather-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          Tx_ sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_end[row] ; ++i)
          {
            sum += Tx_(val[i]) * Tx_(y[col_idx[i]]);
          }
          x[row] = sum;
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
       */
      template<
        typename Tx_,
        typename Ty_>
      void scatter_prim(
        LAFEM::DenseVector<Mem::Main, Tx_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_>& buffer) const
      {
        Tx_ * x(vector.elements());
        const Ty_ * y(buffer.elements());
        const Index * col_idx(_mirror_scatter.col_ind());
        const DataType_* val(_mirror_scatter.val());
        const Index * row_ptr(_mirror_scatter.row_ptr());
        const Index * row_end(_mirror_scatter.row_ptr_end());
        const Index num_rows(_mirror_scatter.rows());

        // loop over all scatter-matrix rows
        for (Index row(0) ; row < num_rows ; ++row)
        {
          // skip empty rows
          if(row_ptr[row] >= row_end[row])
            continue;

          Tx_ sum(Tx_(0));
          for (Index i(row_ptr[row]) ; i < row_end[row] ; ++i)
          {
            sum += Tx_(val[i]) * Tx_(y[col_idx[i]]);
          }
          x[row] = sum;
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
       */
      template<
        typename Tx_,
        typename Ty_>
      void gather_dual(
        LAFEM::DenseVector<Mem::Main, Tx_>& buffer,
        const LAFEM::DenseVector<Mem::Main, Ty_>& vector) const
      {
        gather_prim(buffer, vector);
      }

      /**
       * \brief Performs a scatter-operation on a dual vector.
       *
       * \param[in,out] vector
       * A reference to a dual vector.
       *
       * \param[in] buffer
       * A reference to the buffer vector whose entries are to be scattered.
       */
      template<
        typename Tx_,
        typename Ty_>
      void scatter_dual(
        LAFEM::DenseVector<Mem::Main, Tx_>& vector,
        const LAFEM::DenseVector<Mem::Main, Ty_>& buffer) const
      {
        scatter_prim(vector, buffer);
      }
    }; // class VectorMirror<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_VECTOR_MIRROR_HPP
