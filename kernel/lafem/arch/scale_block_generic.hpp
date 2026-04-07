// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_SCALE_BLOCK_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename ValueType_>
  void ScaleBlock::exec_generic_impl(ValueType_ * r, const ValueType_ * const x, const ValueType_ alpha, const Index size)
  {
    if (x == r)
    {
      FEAT_PRAGMA_OMP(parallel for)
      for (Index i = 0 ; i < size ; ++i)
      {
        for(int j = 0; j < ValueType_::n; ++j)
        {
          r[i][j] *= alpha[j];
        }
      }
    }
    else
    {
      FEAT_PRAGMA_OMP(parallel for)
      for (Index i = 0 ; i < size ; ++i)
      {
        for(int j = 0; j < ValueType_::n; ++j)
        {
          r[i][j] = x[i][j] * alpha[j];
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
