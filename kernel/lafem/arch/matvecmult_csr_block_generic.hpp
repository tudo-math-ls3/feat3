// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MATVECMULT_CSR_BLOCK_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_, int block_size_>
  void MatVecMultCSRBlock::exec_generic_impl(Tiny::Vector<DT_,block_size_> * r, const DT_ alpha, const Tiny::Vector<DT_,block_size_> * const x, const Tiny::Vector<DT_,block_size_> * const y, const DT_ * const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index /*cols*/, const Index /*nonzeros*/, const bool transposed)
  {
    XASSERTM(!transposed, "transposed mat-vec-mult not implemented yet");

    FEAT_PRAGMA_OMP(parallel for schedule(static, 2000))
    for (Index row = 0 ; row < rows ; ++row)
    {
      Tiny::Vector<DT_, block_size_> bsum(0);
      const IT_ end(row_ptr[row + 1]);
      for (IT_ i(row_ptr[row]) ; i < end ; ++i)
      {
        bsum.axpy(val[i], x[col_idx[i]]);
      }
      if(y)
        r[row] = y[row];
      else
        r[row] = DT_(0);
      r[row].axpy(alpha, bsum);
    }
  }
} // namespace FEAT::LAFEM::Arch
