// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MATVECMULT_CSCR_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void MatVecMultCSCRDense::exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const IT_ * const row_idx, const Index rows, const Index cols, const Index nonzero_rows, const Index, const bool transposed)
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
      else if(r != y)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index col = 0 ; col < cols ; ++col)
        {
          r[col] = y[col];
        }
      }

      for (Index nzrow(0) ; nzrow < nonzero_rows ; ++nzrow)
      {
        const Index row(row_idx[nzrow]);
        for (Index i(row_ptr[nzrow]) ; i < row_ptr[nzrow+1] ; ++i)
        {
          r[col_idx[i]] += alpha * val[i] * x[row];
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
          r[row] = DT_(0);
        }
      }
      else if(r != y)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          r[row] = y[row];
        }
      }

      FEAT_PRAGMA_OMP(parallel for schedule(static, 2000))
      for (Index nzrow = 0 ; nzrow < nonzero_rows ; ++nzrow)
      {
        const Index row(row_idx[nzrow]);
        DT_ sum(0);
        const IT_ end(row_ptr[nzrow + 1]);
        for (IT_ i = row_ptr[nzrow] ; i < end ; ++i)
        {
          sum += val[i] * x[col_idx[i]];
        }
        r[row] += sum * alpha;
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
