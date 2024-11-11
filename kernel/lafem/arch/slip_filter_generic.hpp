// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_SLIP_FILTER_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SLIP_FILTER_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SLIP_FILTER_HPP
#error "Do not include this implementation-only header file directly!"
#endif

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <int BlockSize_, typename DT_, typename IT_>
      void SlipFilter::filter_rhs_generic(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue)
      {
        Index block_size = Index(BlockSize_);
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue; ++i)
        {
          DT_ sp(DT_(0));
          DT_ scal(DT_(0));
          for(Index j = 0 ; j < block_size ; ++j)
          {
            sp += v[block_size * sv_indices[i] + j] *nu_elements[block_size * i + j];
            scal += nu_elements[block_size * i + j]*nu_elements[block_size * i + j];
          }

          sp /= scal;

          for(Index j = 0 ; j < block_size ; ++j)
            v[block_size * sv_indices[i] + j] -= sp*nu_elements[block_size * i + j];
        }
      }

      template <int BlockSize_, typename DT_, typename IT_>
      void SlipFilter::filter_def_generic(DT_ * v, const DT_* const nu_elements, const IT_ * const sv_indices, const Index ue)
      {
        Index block_size = Index(BlockSize_);
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < ue; ++i)
        {
          DT_ sp(DT_(0));
          DT_ scal(DT_(0));
          for(Index j = 0 ; j < block_size ; ++j)
          {
            sp += v[block_size * sv_indices[i] + j] *nu_elements[block_size * i + j];
            scal += nu_elements[block_size * i + j]*nu_elements[block_size * i + j];
          }

          sp /= scal;

          for(Index j = 0 ; j < block_size ; ++j)
            v[block_size * sv_indices[i] + j] -= sp*nu_elements[block_size * i + j];
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
/// \endcond

#endif // KERNEL_LAFEM_ARCH_SLIP_FILTER_GENERIC_HPP
