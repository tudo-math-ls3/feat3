// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_ROW_NORM2_BCSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void RowNorm2BCSR::exec_generic_impl(DT_* row_norms, const DT_* const val, const IT_* const /*col_idx*/,
    const IT_* const row_ptr, const Index rows, const int bh, const int bw, const bool squared)
  {
    Index block_height = Index(bh);
    Index block_width = Index(bw);

    FEAT_PRAGMA_OMP(parallel for)
    for (Index row = 0; row < rows; row++)
    {
      Index end = row_ptr[row + 1];

      for(Index i = 0; i < block_height; ++i)
      {
        row_norms[block_height*row + i] = DT_(0);
      }

      for (Index k = row_ptr[row]; k < end; ++k)
      {
        // Manually compute norm2 of row
        for(Index i = 0; i < block_height; ++i)
        {
          for(Index j = 0; j < block_width; ++j)
          {
            row_norms[block_height*row + i] += Math::sqr(val[block_height*block_width*k + i*block_width + j]);
          }
          // Take the square root
          if(!squared)
            row_norms[block_height*row + i] = Math::sqrt(row_norms[block_height*row + i]);
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
