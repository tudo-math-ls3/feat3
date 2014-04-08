#pragma once
#ifndef KERNEL_LAFEM_NONE_FILTER_HPP
#define KERNEL_LAFEM_NONE_FILTER_HPP 1

// includes, FEAST
#include <kernel/lafem/dense_vector.hpp>

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

      /// Creates and returns a (empty) deep copy of this filter.
      NoneFilter clone() const
      {
        return NoneFilter();
      }

      template<typename Algo_, typename Matrix_>
      void filter_mat(Matrix_&) const
      {
      }

      template<typename Algo_, typename Matrix_>
      void filter_offdiag_row_mat(Matrix_&) const
      {
      }

      template<typename Algo_, typename Matrix_>
      void filter_offdiag_col_mat(Matrix_&) const
      {
      }

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      template<typename Algo_>
      void filter_rhs(DenseVector<MemType,DataType,IndexType>&) const
      {
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      template<typename Algo_>
      void filter_sol(DenseVector<MemType,DataType,IndexType>&) const
      {
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      template<typename Algo_>
      void filter_def(DenseVector<MemType,DataType,IndexType>&) const
      {
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      template<typename Algo_>
      void filter_cor(DenseVector<MemType,DataType,IndexType>&) const
      {
      }
    }; // class NoneFilter<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_NONE_FILTER_HPP
