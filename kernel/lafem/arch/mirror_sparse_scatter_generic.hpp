// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MIRROR_SPARSE_SCATTER_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template<typename DT_, typename IT_>
  void MirrorSparseScatter::exec_generic_impl(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, const Index nvec, DT_* vval, const IT_* vidx, const DT_ alpha)
  {
    // loop over all mirror indices
    FEAT_PRAGMA_OMP(parallel for)
    for(Index i = 0; i < nidx; ++i)
    {
      // loop over all vector indices and try to find the match
      for(Index j(0); j < nvec; ++j)
      {
        if(idx[i] == vidx[j])
        {
          vval[j] += alpha*buf[boff+i];
          break;
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
