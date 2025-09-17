// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_MAX_REL_DIFF_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_MAX_REL_DIFF_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_MAX_REL_DIFF_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      DT_ MaxRelDiff::value_generic(const DT_ * const x, const DT_ * const y, const Index size)
      {
        DT_ max(0);
        static DT_ eps = Math::eps<DT_>();
        DT_ tmp;

        for (Index i(0) ; i < size ; ++i)
        {
          // |x_i-y_i|/max(|x_i|+|y_i|, eps)
          tmp = Math::abs(x[i] - y[i]) / Math::max(Math::abs(x[i]) + Math::abs(y[i]), eps);
          if (tmp > max)
          {
            max = tmp;
          }
        }
        return max;
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_MAX_REL_DIFF_GENERIC_HPP
