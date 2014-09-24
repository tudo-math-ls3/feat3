#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_, typename Algo_>
      struct UnitFilterBlocked;

      template <>
      struct UnitFilterBlocked<Mem::Main, Algo::Generic>
      {
        template <typename DT_, typename IT_, Index BlockSize_>
        static void filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_, Index BlockSize_>
        static void filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue);
      };

//      extern template void UnitFilterBlocked<Mem::Main, Algo::Generic>::template filter_rhs<float, unsigned long>(float * v, const float * const sv_elements, const unsigned long * const sv_indices, const Index ue);
//      extern template void UnitFilterBlocked<Mem::Main, Algo::Generic>::template filter_rhs<double, unsigned long>(double * v, const double * const sv_elements, const unsigned long * const sv_indices, const Index ue);
//      extern template void UnitFilterBlocked<Mem::Main, Algo::Generic>::template filter_rhs<float, unsigned int>(float * v, const float * const sv_elements, const unsigned int * const sv_indices, const Index ue);
//      extern template void UnitFilterBlocked<Mem::Main, Algo::Generic>::template filter_rhs<double, unsigned int>(double * v, const double * const sv_elements, const unsigned int * const sv_indices, const Index ue);
//
//      extern template void UnitFilterBlocked<Mem::Main, Algo::Generic>::template filter_def<float, unsigned long>(float * v, const unsigned long * const sv_indices, const Index ue);
//      extern template void UnitFilterBlocked<Mem::Main, Algo::Generic>::template filter_def<double, unsigned long>(double * v, const unsigned long * const sv_indices, const Index ue);
//      extern template void UnitFilterBlocked<Mem::Main, Algo::Generic>::template filter_def<float, unsigned int>(float * v, const unsigned int * const sv_indices, const Index ue);
//      extern template void UnitFilterBlocked<Mem::Main, Algo::Generic>::template filter_def<double, unsigned int>(double * v, const unsigned int * const sv_indices, const Index ue);

      template <>
      struct UnitFilterBlocked<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_, typename IT_, Index BlockSize_>
        static void filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, const Index bs);

        template <typename DT_, typename IT_, Index BlockSize_>
        static void filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue, const Index bs);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/unit_filter_blocked_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
