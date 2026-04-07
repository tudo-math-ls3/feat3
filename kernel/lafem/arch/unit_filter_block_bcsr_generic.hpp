// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCK_BCSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template<typename DT_, typename IT_>
  void UnitFilterBlockBCSR::exec_generic_impl(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for(Index i = 0; i < f_nzes; ++i)
    {
      const IT_ ix(f_idx[i]);
      const DT_* const vx(&f_val[i*Index(block_height)]);

      // replace by unit row
      for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
      {
        // loop over rows in the block
        for(int k(0); k < block_height; ++k)
        {
          // possibly skip row if filter value is NaN
          if(ign_nans && Math::isnan(vx[k]))
            continue;

          // format block row to zero
          for(int l(0); l < block_width; ++l)
            mat[j*IT_(block_height*block_width) + IT_(k*block_width) + IT_(l)] = DT_(0);

          // replace diagonal entry by 1 if requested
          if(unit && (col_idx[j] == ix) && (k < block_width))
            mat[j*IT_(block_height*block_width) + IT_(k*block_width) + IT_(k)] = DT_(1);
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
