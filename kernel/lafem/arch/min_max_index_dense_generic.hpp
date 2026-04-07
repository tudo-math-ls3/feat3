// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MIN_MAX_INDEX_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  Index MinMaxIndexDense::exec_generic_impl(const DT_ * const x, const Index size, const bool min_, const bool abs_)
  {
    DT_ val(abs_ ? Math::abs(x[0]) : x[0]);
    Index idx(0);
    if(min_)
    {
      if(abs_)
      {
        for (Index i(0) ; i < size ; ++i)
        {
          DT_ ax(Math::abs(x[i]));
          if (ax < val)
          {
            val = ax;
            idx = i;
          }
        }
      }
      else
      {
        for (Index i(0) ; i < size ; ++i)
        {
          if (x[i] < val)
          {
            val = x[i];
            idx = i;
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
          DT_ ax(Math::abs(x[i]));
          if (ax < val)
          {
            val = ax;
            idx = i;
          }
        }
      }
      else
      {
        for (Index i(0) ; i < size ; ++i)
        {
          if (x[i] < val)
          {
            val = x[i];
            idx = i;
          }
        }
      }
    }

    return idx;
  }
} // namespace FEAT::LAFEM::Arch
