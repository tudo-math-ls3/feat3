// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MIN_MAX_VALUE_BLOCK_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename ValueType_>
  ValueType_ MinMaxValueBlock::exec_generic_impl(const ValueType_ * const x, const Index size, const bool min_, const bool abs_)
  {
    ValueType_ result(x[0]);

    if(abs_)
    {
      for(int j(0); j < ValueType_::n; ++j)
      {
        result[j] = Math::abs(x[0][j]);
      }
    }

    if(min_)
    {
      if(abs_)
      {
        for (Index i(0) ; i < size ; ++i)
        {
          for(int j(0); j < ValueType_::n; ++j)
          {
            result[j] = Math::min(result[j], Math::abs(x[i][j]));
          }
        }
      }
      else
      {
        for (Index i(0) ; i < size ; ++i)
        {
          for(int j(0); j < ValueType_::n; ++j)
          {
            result[j] = Math::min(result[j], x[i][j]);
          }
        }
      }
    }
    else // max
    {
      if(abs_)
      {
        for (Index i(0) ; i < size ; ++i)
        {
          for(int j(0); j < ValueType_::n; ++j)
          {
            result[j] = Math::max(result[j], Math::abs(x[i][j]));
          }
        }
      }
      else
      {
        for (Index i(0) ; i < size ; ++i)
        {
          for(int j(0); j < ValueType_::n; ++j)
          {
            result[j] = Math::max(result[j], x[i][j]);
          }
        }
      }
    }

    return result;
  }
} // namespace FEAT::LAFEM::Arch
