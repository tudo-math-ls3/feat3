// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_DENSE_CSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template<typename DT_, typename IT_>
  void UnitFilterDenseCSR::exec_generic_impl(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const f_idx, const Index f_nzes, const bool unit, const int block_width)
  {
    if(unit)
    {
      FEAT_PRAGMA_OMP(parallel for)
      for(Index i = 0; i < f_nzes; ++i)
      {
        const IT_ ix(f_idx[i]);
        for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
        {
          mat[j] = (col_idx[j] == ix) ? DT_(1) : DT_(0);
        }
      }
    }
    else
    {
      FEAT_PRAGMA_OMP(parallel for)
      for(Index i = 0; i < f_nzes; ++i)
      {
        const IT_ ix(f_idx[i]);
        for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
        {
          for(Index k = 0; k < Index(block_width); ++k)
            mat[j*Index(block_width) + k] = DT_(0);
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
