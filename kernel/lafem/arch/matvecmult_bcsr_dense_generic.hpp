// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_MATVECMULT_BCSR_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_, int bh_, int bw_>
  void MatVecMultBCSRDense::exec_generic_impl(Tiny::Vector<DT_, bh_> * r, const DT_ alpha, const Tiny::Vector<DT_, bw_> * const x,
    const Tiny::Vector<DT_, bh_> * const y, const Tiny::Matrix<DT_, bh_, bw_>* const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index , const Index)
  {
    if(y != nullptr)
    {
      FEAT_PRAGMA_OMP(parallel for schedule(static, 2000))
      for (Index row = 0 ; row < rows ; ++row)
      {
        Tiny::Vector<DT_, bh_> bsum(0);
        const IT_ end(row_ptr[row + 1]);
        for (IT_ i = row_ptr[row] ; i < end ; ++i)
        {
          bsum.add_mat_vec_mult(val[i], x[col_idx[i]]);
        }
        r[row] = y[row] + (alpha * bsum);
      }
    }
    else
    {
      FEAT_PRAGMA_OMP(parallel for schedule(static, 2000))
      for (Index row = 0 ; row < rows ; ++row)
      {
        Tiny::Vector<DT_, bh_> bsum(0);
        const IT_ end(row_ptr[row + 1]);
        for (IT_ i = row_ptr[row] ; i < end ; ++i)
        {
          bsum.add_mat_vec_mult(val[i], x[col_idx[i]]);
        }
        r[row] = alpha * bsum;
      }
    }
  }

  template <typename DT_, typename IT_, int bh_, int bw_>
  void MatVecMultBCSRDense::exec_generic_transpose_impl(Tiny::Vector<DT_, bw_> * r, const DT_ alpha, const Tiny::Vector<DT_, bh_> * const x,
    const Tiny::Vector<DT_, bw_> * const y, const Tiny::Matrix<DT_, bh_, bw_>* const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index)
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

    for (Index row = 0 ; row < rows ; ++row)
    {
      for (Index i(row_ptr[row]) ; i < row_ptr[row+1] ; ++i)
      {
        r[col_idx[i]].add_vec_mat_mult(x[row], val[i]);
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
} // namespace FEAT::LAFEM::Arch
