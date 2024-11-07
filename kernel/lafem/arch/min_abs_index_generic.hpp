// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
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
      Index MinAbsIndex::value_generic(const DT_ * const x, const Index size)
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

      template <typename ValueType_>
      ValueType_ MinAbsIndex::value_blocked_generic(const ValueType_ * const x, const Index size)
      {
        ValueType_ min(0);

        if(size > 0)
        {
          for(int j(0); j < ValueType_::n; ++j)
          {
            min[j] = Math::abs(x[0][j]);
          }
        }

        for (Index i(0) ; i < size ; ++i)
        {
          for(int j(0); j < ValueType_::n; ++j)
          {
            if(Math::abs(x[i][j]) < Math::abs(min[j]))
            {
              min[j] = Math::abs(x[i][j]);
            }
          }
        }
        return min;
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_MIN_ABS_INDEX_GENERIC_HPP
