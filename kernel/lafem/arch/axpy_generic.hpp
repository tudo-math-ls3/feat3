// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_AXPY_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      void Axpy::value_generic(DT_ * r, const DT_ a, const DT_ * const x, const Index size)
      {
        if (r == x)
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < size ; ++i)
          {
            r[i] *= DT_(1) + a;
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for (Index i = 0 ; i < size ; ++i)
          {
            r[i] += a * x[i];
          }
        }
      }

      template <typename ValueType_>
      void Axpy::value_blocked_generic(ValueType_ * r, const ValueType_ a, const ValueType_ * const x, const Index size)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index i = 0 ; i < size ; ++i)
        {
          for(int j(0); j < ValueType_::n; ++j)
            r[i][j] += a[j] * x[i][j];
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
