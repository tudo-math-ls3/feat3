// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_NORM_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_NORM_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_NORM_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <cmath>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      DT_ Norm2::value_generic(const DT_ * const x, const Index size)
      {
        DT_ r(0);
        FEAT_PRAGMA_OMP(parallel for reduction(+:r))
        for (Index i = 0 ; i < size ; ++i)
        {
          r += x[i] * x[i];
        }

        return (DT_)Math::sqrt(r);
      }

      template <typename ValueType_>
      ValueType_ Norm2::value_blocked_generic(const ValueType_ * const x, const Index size)
      {
        ValueType_ r(0);
        FEAT_PRAGMA_OMP(parallel for)
        for(int j = 0; j<ValueType_::n; ++j)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[j] += x[i][j] * x[i][j];
          }
          r[j] = Math::sqrt(r[j]);
        }
        return r;
      }

      template <typename ValueType_>
      ValueType_ Norm2Sqr::value_blocked_generic(const ValueType_ * const x, const Index size)
      {
        ValueType_ r(0);
        FEAT_PRAGMA_OMP(parallel for)
        for(int j = 0; j<ValueType_::n; ++j)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[j] += x[i][j] * x[i][j];
          }
        }
        return r;
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_NORM_GENERIC_HPP
