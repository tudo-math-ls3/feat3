// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MATVECMULT_DENSE_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void MatVecMultDenseDense::exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index cols, bool transposed)
  {
    if (transposed)
    {
      if(y == nullptr)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < cols ; ++col)
        {
          DT_ sum(0);
          for (Index row(0); row < rows; ++row)
          {
            sum += val[row * cols + col] * x[row];
          }
          r[col] = alpha * sum;
        }
      }
      else
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < cols ; ++col)
        {
          DT_ sum(0);
          for (Index row(0); row < rows; ++row)
          {
            sum += val[row * cols + col] * x[row];
          }
          r[col] = y[col] + alpha * sum;
        }
      }
    }
    else
    {
      if(y == nullptr)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          DT_ sum(0);
          for (Index col(0); col < cols; ++col)
          {
            sum += val[row * cols + col] * x[col];
          }
          r[row] = alpha * sum;
        }
      }
      else
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          DT_ sum(0);
          for (Index col(0); col < cols; ++col)
          {
            sum += val[row * cols + col] * x[col];
          }
          r[row] = y[row] + alpha * sum;
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
