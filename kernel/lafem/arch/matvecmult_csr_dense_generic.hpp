// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MATVECMULT_CSR_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void MatVecMultCSRDense::exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index, const bool transposed)
  {
    if (transposed)
    {
      if(y == nullptr)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < cols ; ++col)
        {
          r[col] = DT_(0);
        }
      }
      else if((r != y) || (alpha != DT_(1)))
      {
        const DT_ ba = DT_(1) / alpha;
        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < cols ; ++col)
        {
          r[col] = ba * y[col];
        }
      }

      for (Index row(0) ; row < rows ; ++row)
      {
        for (Index i(row_ptr[row]) ; i < row_ptr[row+1] ; ++i)
        {
          r[col_idx[i]] += val[i] * x[row];
        }
      }

      if(alpha != DT_(1))
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < cols ; ++col)
        {
          r[col] = alpha * r[col];
        }
      }
    }
    else
    {
      if(y == nullptr)
      {
        FEAT_PRAGMA_OMP(parallel for schedule(static, 2000))
        for (Index row = 0 ; row < rows ; ++row)
        {
          DT_ sum(0);
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i = row_ptr[row] ; i < end ; ++i)
          {
            sum += val[i] * x[col_idx[i]];
          }
          r[row] = sum * alpha;
        }
      }
      else
      {
        FEAT_PRAGMA_OMP(parallel for schedule(static, 2000))
        for (Index row = 0 ; row < rows ; ++row)
        {
          DT_ sum(0);
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i = row_ptr[row] ; i < end ; ++i)
          {
            sum += val[i] * x[col_idx[i]];
          }
          r[row] = y[row] + (sum * alpha);
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
