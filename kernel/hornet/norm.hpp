#pragma once
#ifndef KERNEL_HORNET_NORM_HPP
#define KERNEL_HORNET_NORM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/hornet/dense_vector.hpp>

#include <cmath>


namespace FEAST
{
  template <typename Arch_, typename BType_>
  struct Norm2
  {
    template <typename DT_>
    static DT_ value(const DenseVector<Arch_, DT_> & x)
    {
      const DT_ * xp(x.elements());
      DT_ r(0);
      const Index size(x.size());

      DT_ xpv;
      for (Index i(0) ; i < size ; ++i)
      {
        xpv = xp[i];
        r += xpv * xpv;
      }

      return sqrt(r);
    }
  };

} // namespace FEAST

#endif // KERNEL_HORNET_NORM2_HPP
