#pragma once
#ifndef KERNEL_LAFEM_PRODUCT_HPP
#define KERNEL_LAFEM_PRODUCT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/component_product.hpp>
#include <kernel/lafem/product_matvec.hpp>
#include <kernel/lafem/dot_product.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Product calculations.
     *
     * Convenience class for all multiplication related operations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Algo_>
    struct Product
    {
      template <typename Mem_, typename DT_>
      static void value(DenseVector<Mem_, DT_> & r, const DenseVector<Mem_, DT_> & x, const DenseVector<Mem_, DT_> & y)
      {
        ComponentProduct<Algo_>::value(r, x, y);
      }

      template <typename Mem_, typename DT_, typename SM_>
      static void value(DenseVector<Mem_, DT_> & r, const SM_ & a, const DenseVector<Mem_, DT_> & b)
      {
        ProductMatVec<Algo_>::value(r, a, b);
      }

      template <typename Mem_, typename DT_>
      static DT_ value(const DenseVector<Mem_, DT_> & x, const DenseVector<Mem_, DT_> & y)
      {
        return DotProduct<Algo_>::value(x, y);
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_PRODUCT_HPP
