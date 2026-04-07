// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_SCALE_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void ScaleDense::exec_generic_impl(DT_ * r, const DT_ * const x, const DT_ alpha, const Index size)
  {
    if (x == r)
    {
      FEAT_PRAGMA_OMP(parallel for)
      for (Index i = 0 ; i < size ; ++i)
      {
        r[i] *= alpha;
      }
    }
    else
    {
      FEAT_PRAGMA_OMP(parallel for)
      for (Index i = 0 ; i < size ; ++i)
      {
        r[i] = x[i] * alpha;
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
