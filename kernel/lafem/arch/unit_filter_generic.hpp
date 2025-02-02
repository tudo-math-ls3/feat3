// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP
#error "Do not include this implementation-only header file directly!"
#endif

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_>
      void UnitFilter::filter_rhs_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue ; ++i)
        {
          v[sv_indices[i]] = sv_elements[i];
        }
      }

      template <typename DT_, typename IT_>
      void UnitFilter::filter_def_generic(DT_ * v, const IT_ * const sv_indices, const Index ue)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue ; ++i)
        {
          v[sv_indices[i]] = DT_(0);
        }
      }

      template<typename DT_, typename IT_>
      void UnitFilter::filter_unit_mat_generic(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const sv_indices, const Index ue)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue; ++i)
        {
          const IT_ ix(sv_indices[i]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            mat[j] = (col_idx[j] == ix) ? DT_(1) : DT_(0);
          }
        }
      }

      template<typename DT_, typename IT_>
      void UnitFilter::filter_offdiag_row_mat_generic(DT_* mat, const IT_* const row_ptr, int block_width, const IT_ * const sv_indices, const Index ue)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue; ++i)
        {
          const IT_ ix(sv_indices[i]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            for(int l(0); l < block_width; ++l)
              mat[j*block_width + l] = DT_(0);
          }
        }
      }

      template<typename DT_, typename IT_>
      void UnitFilter::filter_weak_matrix_rows_generic(DT_* mat_a, const DT_ * const mat_m, const IT_* const row_ptr,
                                                              const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue; ++i)
        {
          const IT_ ix(sv_indices[i]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            mat_a[j] = sv_elements[i] * mat_m[j];
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
/// \endcond

#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP
