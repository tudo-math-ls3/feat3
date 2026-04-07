// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_LUMPING_BCSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void LumpingBCSR::exec_generic_impl(DT_ * lump, const DT_ * const val, const IT_ * const /*col_idx*/,
    const IT_ * const row_ptr, const Index rows, const int block_height, const int block_width)
  {
    const Index bh = Index(block_height);
    const Index bw = Index(block_width);

    FEAT_PRAGMA_OMP(parallel for)
    for (Index row = 0; row < rows; row++)
    {
      Index end = row_ptr[row + 1];

      for(Index i(0); i < bh; ++i)
      {
          lump[bh*row + i] = DT_(0);
      }

      for (Index col = row_ptr[row]; col < end; col++)
      {
        for(Index i(0); i < bh; ++i)
        {
          for(Index j(0); j < bw; ++j)
          {
            lump[bh*row + i] += val[bh*bw*col + i*bw + j];
          }
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
