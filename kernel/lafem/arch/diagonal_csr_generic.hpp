// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_DIAGONAL_CSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void DiagonalCSR::exec_generic_impl(DT_ * diag, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for (Index row = 0; row < rows; row++)
    {
      diag[row] = DT_(0);
      const Index end = row_ptr[row + 1];
      for (Index j = row_ptr[row]; j < end; j++)
      {
        if (row == col_idx[j])
        {
          diag[row] = val[j];
          break;
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
