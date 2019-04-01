// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
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
      void Axpy<Mem::Main>::dv_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
      {
        if (r == y)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] += a * x[i];
          }
        }
        else if (r == x)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] *= a;
            r[i]+= y[i];
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] = (a * x[i]) + y[i];
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_AXPY_GENERIC_HPP
