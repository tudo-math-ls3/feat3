#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP
  #error "Do not include this implementation-only header file directly!"
#endif

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_, Index BlockSize_>
      void UnitFilterBlocked<Mem::Main, Algo::Generic>::filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
      {
        for(Index i(0); i < ue; ++i)
        {
          for(Index j(0); j < BlockSize_; ++j)
          {
             v[ BlockSize_* sv_indices[i] + j ] = sv_elements[BlockSize_*i+j];

          }

        }
      }

      template <typename DT_, typename IT_, Index BlockSize_>
      void UnitFilterBlocked<Mem::Main, Algo::Generic>::filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue)
      {
        for(Index i(0); i < ue; ++i)
        {
          for(Index j(0); j < BlockSize_; ++j)
          v[ BlockSize_*sv_indices[i] + j ] = DT_(0);
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_GENERIC_HPP
