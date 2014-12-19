#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP
#error "Do not include this implementation-only header file directly!"
#endif

/// \cond internal
namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_>
      void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
      {
        for(Index i(0); i < ue ; ++i)
        {
          v[sv_indices[i]] = sv_elements[i];
        }
      }

      template <typename DT_, typename IT_>
      void UnitFilter<Mem::Main, Algo::Generic>::filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue)
      {
        for(Index i(0); i < ue ; ++i)
        {
          v[sv_indices[i]] = DT_(0);
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST
/// \endcond

#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_GENERIC_HPP
