// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MIN_MAX_VALUE_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  DT_ MinMaxValueDense::exec_generic_impl(const DT_ * const x, const Index size, const bool min_, const bool abs_)
  {
    DT_ result(abs_ ? Math::abs(x[0]) : x[0]);
    if(min_)
    {
      if(abs_)
      {
        FEAT_PRAGMA_OMP(parallel for reduction(min:result))
        for (Index i = 0 ; i < size ; ++i)
        {
          result = Math::min(result, Math::abs(x[i]));
        }
      }
      else
      {
        FEAT_PRAGMA_OMP(parallel for reduction(min:result))
        for (Index i = 0 ; i < size ; ++i)
        {
          result = Math::min(result, x[i]);
        }
      }
    }
    else // max
    {
      if(abs_)
      {
        FEAT_PRAGMA_OMP(parallel for reduction(max:result))
        for (Index i = 0 ; i < size ; ++i)
        {
          result = Math::max(result, Math::abs(x[i]));
        }
      }
      else
      {
        FEAT_PRAGMA_OMP(parallel for reduction(max:result))
        for (Index i = 0 ; i < size ; ++i)
        {
          result = Math::max(result, x[i]);
        }
      }
    }

    return result;
  }
} // namespace FEAT::LAFEM::Arch
