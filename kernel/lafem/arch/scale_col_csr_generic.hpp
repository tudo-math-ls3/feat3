// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_SCALE_COL_CSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void ScaleColCSR::exec_generic_impl(DT_ * r, const DT_ * const a, const IT_ * const col_idx,
    const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for (Index row = 0 ; row < rows ; ++row)
    {
      const IT_ end(row_ptr[row + 1]);
      for (IT_ i = row_ptr[row] ; i < end ; ++i)
      {
        r[i] = a[i] * x[col_idx[i]];
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
