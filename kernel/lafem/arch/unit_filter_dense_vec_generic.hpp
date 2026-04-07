// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_DENSE_VEC_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void UnitFilterDenseVec::exec_generic_impl(DT_ * v, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool zero)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for(Index i = 0; i < f_nzes ; ++i)
    {
      v[f_idx[i]] = zero ? DT_(0) : f_val[i];
    }
  }
} // namespace FEAT::LAFEM::Arch
