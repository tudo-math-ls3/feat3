// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_SLIP_FILTER_BLOCK_HPP
#error "Do not include this implementation-only header file directly!"
#endif

/// \cond internal
namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void SlipFilterBlock::exec_generic_impl(DT_ * v, const DT_ * const f_nu, const IT_ * const f_idx, const Index f_nzes, const int bs)
  {
    const Index block_size = Index(bs);

    FEAT_PRAGMA_OMP(parallel for)
    for(Index i = 0; i < f_nzes; ++i)
    {
      DT_ sp(DT_(0));
      DT_ scal(DT_(0));
      for(Index j = 0 ; j < block_size ; ++j)
      {
        sp += v[block_size * f_idx[i] + j] *f_nu[block_size * i + j];
        scal += f_nu[block_size * i + j]*f_nu[block_size * i + j];
      }

      sp /= scal;

      for(Index j = 0 ; j < block_size ; ++j)
        v[block_size * f_idx[i] + j] -= sp*f_nu[block_size * i + j];
    }
  }
} // namespace FEAT::LAFEM::Arch
/// \endcond
