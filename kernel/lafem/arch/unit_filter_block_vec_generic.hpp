// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCK_VEC_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void UnitFilterBlockVec::exec_generic_impl(DT_ * v, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool zero, const bool ign_nans, const int block_size)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for(Index i = 0; i < f_nzes; ++i)
    {
      for(int j = 0 ; j < block_size ; ++j)
      {
        const Index k = Index(block_size) * i + Index(j);

        // skip if filter value is NaN if desired
        if(ign_nans && Math::isnan(f_val[k]))
          continue;

        // write filter value or zero
        v[Index(block_size) * f_idx[i] + Index(j)] = zero ? DT_(0) : f_val[k];
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
