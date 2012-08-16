#pragma once
#ifndef KERNEL_LAFEM_ELEMENT_PRODUCT_HPP
#define KERNEL_LAFEM_ELEMENT_PRODUCT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Arch_, typename BType_>
    struct ElementProduct
    {
    };

    template <>
    struct ElementProduct<Archs::CPU, Archs::Generic>
    {
      template <typename DT_>
      static void value(DenseVector<Archs::CPU, DT_> & r, const DenseVector<Archs::CPU, DT_> & x, const DenseVector<Archs::CPU, DT_> & y)
      {
        if (x.size() != y.size())
          throw InternalError("Vector size does not match!");
        if (x.size() != r.size())
          throw InternalError("Vector size does not match!");

        const DT_ * xp(x.elements());
        const DT_ * yp(y.elements());
        DT_ * rp(r.elements());
        const Index size(r.size());

        for (Index i(0) ; i < size ; ++i)
        {
          rp[i] = xp[i] * yp[i];
        }
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ELEMENT_PRODUCT_HPP
