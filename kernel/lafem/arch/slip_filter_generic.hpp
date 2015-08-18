#pragma once
#ifndef KERNEL_LAFEM_ARCH_SLIP_FILTER_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SLIP_FILTER_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SLIP_FILTER_HPP
#error "Do not include this implementation-only header file directly!"
#endif

/// \cond internal
namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_, int BlockSize_>
      void SlipFilter<Mem::Main>::filter_rhs_generic(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue)
      {
        Index block_size = Index(BlockSize_);
        for(Index i(0); i < ue; ++i)
        {
          DT_ sp(DT_(0));
          for(Index j(0) ; j < block_size ; ++j)
            sp += v[block_size * sv_indices[i] + j] *nu_elements[block_size * i + j];

          for(Index j(0) ; j < block_size ; ++j)
            v[block_size * sv_indices[i] + j] -= sp*nu_elements[block_size * i + j];
        }
      }

      template <typename DT_, typename IT_, int BlockSize_>
      void SlipFilter<Mem::Main>::filter_def_generic(DT_ * v, const DT_* const nu_elements, const IT_ * const sv_indices, const Index ue)
      {
        Index block_size = Index(BlockSize_);
        for(Index i(0); i < ue; ++i)
        {
          DT_ sp(DT_(0));
          for(Index j(0) ; j < block_size ; ++j)
            sp += v[block_size * sv_indices[i] + j] *nu_elements[block_size * i + j];

          for(Index j(0) ; j < block_size ; ++j)
            v[block_size * sv_indices[i] + j] -= sp*nu_elements[block_size * i + j];
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST
/// \endcond

#endif // KERNEL_LAFEM_ARCH_SLIP_FILTER_GENERIC_HPP
