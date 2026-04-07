// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MATMATMULT_DENSE_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void MatMatMultDenseDense::exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index cols, const Index inner)
  {
    FEAT_PRAGMA_OMP(parallel for)
    for (Index i = 0 ; i < rows ; ++i)
    {
      for (Index j = 0 ; j < cols ; ++j)
      {
        DT_ sum(0.);
        Index xindex(i * inner);
        Index yindex(j);
        for (Index xcol(0) ; xcol < inner ; ++xcol)
        {
          sum  = sum + x[xindex + xcol] * y[yindex + xcol * cols];
        }
        if(z != nullptr)
          r[i * cols + j] = z[i * cols + j] + alpha * sum;
        else
          r[i * cols + j] = alpha * sum;
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
