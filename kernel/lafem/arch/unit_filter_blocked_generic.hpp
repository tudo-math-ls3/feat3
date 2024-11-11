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
      template <int BlockSize_, typename DT_, typename IT_>
      void UnitFilterBlocked::filter_rhs_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
      {
        Index block_size = Index(BlockSize_);
        if(ign_nans)
        {
          FEAT_PRAGMA_OMP(parallel for)
          for(Index i = 0; i < ue; ++i)
          {
            for(Index j = 0 ; j < block_size ; ++j)
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
            for(Index j = 0 ; j < block_size ; ++j)
            {
              v[block_size * sv_indices[i] + j] = sv_elements[block_size * i + Index(j)];
            }
          }
        }
      }

      template <int BlockSize_, typename DT_, typename IT_>
      void UnitFilterBlocked::filter_def_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
      {
        Index block_size = Index(BlockSize_);
        if(ign_nans)
        {
          FEAT_PRAGMA_OMP(parallel for)
          for(Index i = 0; i < ue; ++i)
          {
            for(Index j = 0 ; j < block_size ; ++j)
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
            for(Index j = 0 ; j < BlockSize_ ; ++j)
              v[block_size * sv_indices[i] + j] = DT_(0);
          }
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
/// \endcond

#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP
