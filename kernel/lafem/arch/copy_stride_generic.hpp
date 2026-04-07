// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_COPY_STRIDE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void CopyStride::exec_generic_impl(DT_ * r, const DT_ * const x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for(Index i = 0; i < count; ++i)
    {
      r[i*Index(stride_r) + Index(comp_r)] = x[i*Index(stride_x) + Index(comp_x)];
    }
  }
} // namespace FEAT::LAFEM::Arch
