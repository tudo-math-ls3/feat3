// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_AXPY_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void AxpyDense::exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const Index size)
  {
    if (r == x)
    {
      DT_ alpha1 = DT_(1) + alpha;
      FEAT_PRAGMA_OMP(parallel for)
      for (Index i = 0 ; i < size ; ++i)
      {
        r[i] *= alpha1;
      }
    }
    else
    {
      FEAT_PRAGMA_OMP(parallel for)
      for (Index i = 0 ; i < size ; ++i)
      {
        r[i] += alpha * x[i];
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
