#pragma once
#ifndef KERNEL_LAFEM_NONE_FILTER_HPP
#define KERNEL_LAFEM_NONE_FILTER_HPP 1

// includes, FEAST
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief None Filter class template.
     *
     * This class implements a filter that does nothing.
     *
     * \note Not to be confused with <em>do-nothing boundary conditions</em>.
     *
     * \author Peter Zajac
     */
    template<
      typename MemType_,
      typename DataType_,
      typename IndexType_>
    class NoneFilter
    {
    public:
      /// mem-type typedef
      typedef MemType_ MemType;
      /// data-type typedef
      typedef DataType_ DataType;
      /// index-type typedef
      typedef IndexType_ IndexType;

      /// our supported vector type
      typedef DenseVector<MemType, DataType, IndexType> VectorType;

      /// Creates and returns a (empty) deep copy of this filter.
      NoneFilter clone() const
      {
        return NoneFilter();
      }

      template<typename MT2_, typename DT2_, typename IT2_>
      void convert(const NoneFilter<MT2_, DT2_, IT2_>&)
      {
      }

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType&) const
      {
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      void filter_sol(VectorType&) const
      {
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      void filter_def(VectorType&) const
      {
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      void filter_cor(VectorType&) const
      {
      }
    }; // class NoneFilter<...>

    /**
     * \brief Blocked None Filter class template.
     *
     * This class implements a filter that does nothing.
     *
     * \note Not to be confused with <em>do-nothing boundary conditions</em>.
     *
     * \author Peter Zajac
     */
    template<
      typename MemType_,
      typename DataType_,
      typename IndexType_,
      int BlockSize_
    >
    class NoneFilterBlocked
    {
    public:
      /// mem-type typedef
      typedef MemType_ MemType;
      /// data-type typedef
      typedef DataType_ DataType;
      /// index-type typedef
      typedef IndexType_ IndexType;

      static constexpr int BlockSize = BlockSize_;

      /// our supported vector type
      typedef DenseVectorBlocked<MemType, DataType, IndexType, BlockSize> VectorType;

      /// Creates and returns a (empty) deep copy of this filter.
      NoneFilterBlocked clone() const
      {
        return NoneFilterBlocked();
      }

      template<typename MT2_, typename DT2_, typename IT2_>
      void convert(const NoneFilterBlocked<MT2_, DT2_, IT2_, BlockSize>&)
      {
      }

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType&) const
      {
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      void filter_sol(VectorType&) const
      {
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      void filter_def(VectorType&) const
      {
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      void filter_cor(VectorType&) const
      {
      }
    }; // class NoneFilterBlocked<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_NONE_FILTER_HPP
