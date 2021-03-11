// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
#error "Do not include this implementation-only header file directly!"
#endif

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_, int BlockSize_>
      void UnitFilterBlocked<Mem::Main>::filter_rhs_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
      {
        Index block_size = Index(BlockSize_);
        for(Index i(0); i < ue; ++i)
        {
          for(Index j(0) ; j < block_size ; ++j)
            v[block_size * sv_indices[i] + j] = sv_elements[block_size * i + Index(j)];
        }
      }

      template <typename DT_, typename IT_, int BlockSize_>
      void UnitFilterBlocked<Mem::Main>::filter_def_generic(DT_ * v, const IT_ * const sv_indices, const Index ue)
      {
        Index block_size = Index(BlockSize_);
        for(Index i(0); i < ue; ++i)
        {
          for(Index j(0) ; j < BlockSize_ ; ++j)
            v[block_size * sv_indices[i] + j] = DT_(0);
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
/// \endcond

#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP
