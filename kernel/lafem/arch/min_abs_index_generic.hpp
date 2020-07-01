// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_MIN_ABS_INDEX_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_MIN_ABS_INDEX_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_MIN_ABS_INDEX_HPP
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
      Index MinAbsIndex<Mem::Main>::value_generic(const DT_ * const x, const Index size)
      {
        DT_ min(Math::abs(x[0]));
        Index min_i(0);

        for (Index i(0) ; i < size ; ++i)
        {
          if (Math::abs(x[i]) < min)
          {
            min = Math::abs(x[i]);
            min_i = i;
          }
        }

        return min_i;
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_MIN_ABS_INDEX_GENERIC_HPP
