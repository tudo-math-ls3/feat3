// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_>
      void UnitFilterBlocked::filter_rhs_generic(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
      {
        if(ign_nans)
        {
          FEAT_PRAGMA_OMP(parallel for)
          for(Index i = 0; i < ue; ++i)
          {
            for(int j = 0 ; j < block_size ; ++j)
            {
              if(!Math::isnan(sv_elements[block_size * i + Index(j)]))
                v[block_size * sv_indices[i] + j] = sv_elements[block_size * i + Index(j)];
            }
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for(Index i = 0; i < ue; ++i)
          {
            for(int j = 0 ; j < block_size ; ++j)
            {
              v[block_size * sv_indices[i] + j] = sv_elements[block_size * i + Index(j)];
            }
          }
        }
      }

      template <typename DT_, typename IT_>
      void UnitFilterBlocked::filter_def_generic(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
      {
        if(ign_nans)
        {
          FEAT_PRAGMA_OMP(parallel for)
          for(Index i = 0; i < ue; ++i)
          {
            for(int j = 0 ; j < block_size ; ++j)
            {
              if(!Math::isnan(sv_elements[block_size * i + Index(j)]))
                v[block_size * sv_indices[i] + j] = DT_(0);
            }
          }
        }
        else
        {
          FEAT_PRAGMA_OMP(parallel for)
          for(Index i = 0; i < ue; ++i)
          {
            for(int j = 0 ; j < block_size ; ++j)
              v[block_size * sv_indices[i] + j] = DT_(0);
          }
        }
      }

      template<typename DT_, typename IT_>
      void UnitFilterBlocked::filter_unit_mat_generic(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, int block_height, int block_width,
                                                      const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue; ++i)
        {
          const IT_ ix(sv_indices[i]);
          const DT_* const vx(&sv_elements[i*block_height]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            // loop over rows in the block
            for(int k(0); k < block_height; ++k)
            {
              // possibly skip row if filter value is NaN
              if(ign_nans && Math::isnan(vx[k]))
                continue;
              for(int l(0); l < block_width; ++l)
                mat[j*block_height*block_width + k*block_width + l] = DT_(0);
              if((col_idx[j] == ix) && (k < block_width))
                mat[j*block_height*block_width + k*block_width + k] = DT_(1);
            }
          }
        }
      }

      template<typename DT_, typename IT_>
      void UnitFilterBlocked::filter_offdiag_row_mat_generic(DT_* mat, const IT_* const row_ptr, int block_height, int block_width,
                                                              const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue; ++i)
        {
          const IT_ ix(sv_indices[i]);
          const DT_* const vx(&sv_elements[i*block_height]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            // loop over rows in the block
            for(int k(0); k < block_height; ++k)
            {
              // possibly skip row if filter value is NaN
              if(ign_nans && Math::isnan(vx[k]))
                continue;
              for(int l(0); l < block_width; ++l)
                mat[j*block_height*block_width + k*block_width + l] = DT_(0);
            }
          }
        }
      }

      template<typename DT_, typename IT_>
      void UnitFilterBlocked::filter_weak_matrix_rows_generic(DT_* mat_a, const DT_ * const mat_m, const IT_* const row_ptr, int block_height, int block_width,
                                                              const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue; ++i)
        {
          const IT_ ix(sv_indices[i]);
          const DT_* const vx(&sv_elements[i*block_height]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            // loop over rows in the block
            for(int k(0); k < block_height; ++k)
            {
              for(int l(0); l < block_width; ++l)
              {
                mat_a[j*block_height*block_width + k*block_width + l] = vx[k] * mat_m[j*block_height*block_width + k*block_width + l];
              }
            }
          }
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
/// \endcond

#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP
