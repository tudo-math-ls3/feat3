// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_DIAGONAL_BCSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void DiagonalBCSR::exec_generic_impl(DT_ * diag, const DT_ * const val, const IT_ * const col_idx,
    const IT_ * const row_ptr, const Index rows, const int block_size)
  {
    const Index bs = Index(block_size);

    FEAT_PRAGMA_OMP(parallel for)
    for (Index row = 0; row < rows; row++)
    {
      Index end = row_ptr[row + 1];

      for(Index i(0); i < bs; ++i)
      {
        diag[bs*row + i] = DT_(0);
      }

      for (Index col = row_ptr[row]; col < end; col++)
      {
        if(col_idx[col] == row)
        {
          for(Index i(0); i < bs; ++i)
          {
            diag[bs*row + i] = val[bs*bs*col + i*bs + i];
          }
          break;
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
