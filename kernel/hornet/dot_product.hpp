#pragma once
#ifndef KERNEL_HORNET_DOT_PRODUCT_HPP
#define KERNEL_HORNET_DOT_PRODUCT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/hornet/dense_vector.hpp>



namespace FEAST
{
  template <typename Arch_, typename BType_>
  struct DotProduct
  {
  };

  template <>
  struct DotProduct<Archs::CPU, Archs::Generic>
  {
    template <typename DT_>
    static DT_ value(const DenseVector<Archs::CPU, DT_> & x, const DenseVector<Archs::CPU, DT_> & y)
    {
      if (x.size() != y.size())
        throw InternalError("Vector size does not match!");

      const DT_ * xp(x.elements());
      const DT_ * yp(y.elements());
      DT_ r(0);
      const Index size(x.size());

      for (Index i(0) ; i < size ; ++i)
      {
        r += xp[i] * yp[i];
      }

      return r;
    }
  };

} // namespace FEAST

#endif // KERNEL_HORNET_DOT_PRODUCT_HPP
