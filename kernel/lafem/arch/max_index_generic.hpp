// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_MAX_INDEX_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_MAX_INDEX_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_MAX_INDEX_HPP
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
      Index MaxIndex::value_generic(const DT_ * const x, const Index size)
      {
        DT_ max(x[0]);
        Index max_i(0);

        for (Index i(0) ; i < size ; ++i)
        {
          if (x[i] > max)
          {
            max = x[i];
            max_i = i;
          }
        }

        return max_i;
      }

      template <typename ValueType_>
      ValueType_ MaxIndex::value_blocked_generic(const ValueType_ * const x, const Index size)
      {
        ValueType_ max(0);

        if(size > 0)
        {
          max = x[0];
        }

        for (Index i(0) ; i < size ; ++i)
        {
          for(int j(0); j < ValueType_::n; ++j)
          {
            if(x[i][j] > max[j])
            {
              max[j] = x[i][j];
            }
          }
        }
        return max;
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_MAX_INDEX_GENERIC_HPP
