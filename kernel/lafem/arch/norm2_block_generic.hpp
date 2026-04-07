// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_NORM2_BLOCK_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename ValueType_>
  ValueType_ Norm2Block::exec_generic_impl(const ValueType_ * const x, const Index size)
  {
    ValueType_ r(0);

    for (Index i(0) ; i < size ; ++i)
    {
      for(int j = 0; j < ValueType_::n; ++j)
      {
        r[j] += x[i][j] * x[i][j];
      }
    }

    for(int j = 0; j < ValueType_::n; ++j)
      r[j] = Math::sqrt(r[j]);

    return r;
  }
} // namespace FEAT::LAFEM::Arch
