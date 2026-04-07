// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MIRROR_DENSE_GATHER_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template<typename DT_, typename IT_>
  void MirrorDenseGather::exec_generic_impl(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for(Index i = 0; i < nidx; ++i)
    {
      buf[boff+i] = vec[idx[i]];
    }
  }
} // namespace FEAT::LAFEM::Arch
