// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_ROW_NORM2_CSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void RowNorm2CSR::exec_generic_impl(DT_* row_norms, const DT_* const a, const IT_* const /*col_idx*/, const IT_* const row_ptr, const Index rows, const bool squared)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for (Index row = 0; row < rows; row++)
    {
      Index end = row_ptr[row + 1];
      DT_ norm(0);

      for (Index j = row_ptr[row]; j < end; j++)
      {
        norm += Math::sqr(a[j]);
      }

      row_norms[row] = (squared ? norm : Math::sqrt(norm));
    }
  }
} // namespace FEAT::LAFEM::Arch
