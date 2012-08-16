#pragma once
#ifndef KERNEL_LAFEM_NORM_HPP
#define KERNEL_LAFEM_NORM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>

#include <cmath>


namespace FEAST
{
  namespace LAFEM
  {
    template <typename Arch_, typename BType_>
    struct Norm2
    {
    };

    template <>
    struct Norm2<Archs::CPU, Archs::Generic>
    {
      template <typename DT_>
      static DT_ value(const DenseVector<Archs::CPU, DT_> & x)
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

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_NORM2_HPP
