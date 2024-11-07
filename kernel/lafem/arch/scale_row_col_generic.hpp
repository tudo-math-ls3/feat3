// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SCALE_ROW_COL_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {

      template <typename DT_, typename IT_>
      void ScaleRows::csr_generic(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/,
                                                    const IT_ * const row_ptr, const DT_ * const x,
                                                    const Index rows, const Index, const Index)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i = row_ptr[row] ; i < end ; ++i)
          {
            r[i] = a[i] * x[row];
          }
        }
      }

      template <int bh_, int bw_, typename DT_, typename IT_>
      void ScaleRows::bcsr_generic(DT_ * r, const DT_ * const a, const IT_ * const /*col_ind*/,
                                                    const IT_ * const row_ptr, const DT_ * const x,
                                                    const Index rows, const Index, const Index)
      {

        Tiny::Matrix<DT_, bh_, bw_> * const br(reinterpret_cast<Tiny::Matrix<DT_, bh_, bw_> *>(r));
        const Tiny::Matrix<DT_, bh_, bw_> * const ba(reinterpret_cast<const Tiny::Matrix<DT_, bh_, bw_> *>(a));
        const Tiny::Vector<DT_, bh_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, bh_> *>(x));
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i = row_ptr[row] ; i < end ; ++i)
          {
            for (int irow(0); irow < bh_; ++ irow )
            {
              for (int icol(0); icol < bw_; ++icol)
              {
                br[i][irow][icol] = ba[i][irow][icol] * bx[row][irow];
              }
            }
          }
        }
      }

      // ***********************************************

      template <typename DT_, typename IT_>
      void ScaleCols::csr_generic(DT_ * r, const DT_ * const a, const IT_ * const col_ind,
                                                    const IT_ * const row_ptr, const DT_ * const x,
                                                    const Index rows, const Index, const Index)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i = row_ptr[row] ; i < end ; ++i)
          {
            r[i] = a[i] * x[col_ind[i]];
          }
        }
      }

      template <int bh_, int bw_, typename DT_, typename IT_>
      void ScaleCols::bcsr_generic(DT_ * r, const DT_ * const a, const IT_ * const col_ind,
                                                    const IT_ * const row_ptr, const DT_ * const x,
                                                    const Index rows, const Index, const Index)
      {

        Tiny::Matrix<DT_, bh_, bw_> * const br(reinterpret_cast<Tiny::Matrix<DT_, bh_, bw_> *>(r));
        const Tiny::Matrix<DT_, bh_, bw_> * const ba(reinterpret_cast<const Tiny::Matrix<DT_, bh_, bw_> *>(a));
        const Tiny::Vector<DT_, bw_> * const bx(reinterpret_cast<const Tiny::Vector<DT_, bw_> *>(x));
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0 ; row < rows ; ++row)
        {
          const IT_ end(row_ptr[row + 1]);
          for (IT_ i = row_ptr[row] ; i < end ; ++i)
          {
            for (int irow = 0; irow < bh_; ++irow)
            {
              for (int icol = 0; icol < bw_; ++icol)
              {
                br[i][irow][icol] = ba[i][irow][icol] * bx[col_ind[i]][icol];
              }
            }
          }
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_SCALE_ROW_COL_GENERIC_HPP
