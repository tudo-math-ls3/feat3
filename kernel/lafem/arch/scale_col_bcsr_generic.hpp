// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_SCALE_COL_BCSR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void ScaleColBCSR::exec_generic_impl(DT_ * r, const DT_ * const a, const IT_ * const col_idx,
    const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index, const int bh, const int bw)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for (Index i = 0 ; i < rows ; ++i)
    {
      const IT_ end(row_ptr[i + 1]);
      for (IT_ j(row_ptr[i]) ; j < end ; ++j)
      {
        for (int h(0) ; h < bh ; ++h)
        {
          for (int w(0) ; w < bw ; ++w)
          {
            r[j * IT_(bh * bw) + IT_(h * bw + w)] = a[j * IT_(bh * bw) + IT_(h * bw + w)] * x[col_idx[j]*IT_(bw) + IT_(w)];
          }
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
