// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_TRIPLE_DOT_PRODUCT_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  DT_ TripleDotProductDense::exec_generic_impl(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size)
  {
    DT_ r(0);

    if (x == y)
    {
      FEAT_PRAGMA_OMP(parallel for reduction(+:r))
      for (Index i = 0 ; i < size ; ++i)
        r += x[i] * x[i] * z[i];
    }
    else if (x == z)
    {
      FEAT_PRAGMA_OMP(parallel for reduction(+:r))
      for (Index i = 0 ; i < size ; ++i)
        r += x[i] * x[i] * y[i];
    }
    else if (y == z)
    {
      FEAT_PRAGMA_OMP(parallel for reduction(+:r))
      for (Index i = 0 ; i < size ; ++i)
        r += x[i] * y[i] * y[i];
    }
    else
    {
      FEAT_PRAGMA_OMP(parallel for reduction(+:r))
      for (Index i = 0 ; i < size ; ++i)
        r += x[i] * y[i] * z[i];
    }

    return r;
  }
} // namespace FEAT::LAFEM::Arch
