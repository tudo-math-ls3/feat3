#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

/// \cond internal
namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_, typename Algo_>
      struct UnitFilter;

      template <>
      struct UnitFilter<Mem::Main, Algo::Generic>
      {
        template <typename DT_, typename IT_>
        static void filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_>
        static void filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue);
      };

      extern template void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(float * v, const float * const sv_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(double * v, const double * const sv_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(float * v, const float * const sv_elements, const unsigned int * const sv_indices, const Index ue);
      extern template void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(double * v, const double * const sv_elements, const unsigned int * const sv_indices, const Index ue);

      extern template void UnitFilter<Mem::Main, Algo::Generic>::filter_def(float * v, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilter<Mem::Main, Algo::Generic>::filter_def(double * v, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilter<Mem::Main, Algo::Generic>::filter_def(float * v, const unsigned int * const sv_indices, const Index ue);
      extern template void UnitFilter<Mem::Main, Algo::Generic>::filter_def(double * v, const unsigned int * const sv_indices, const Index ue);

      template <>
      struct UnitFilter<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_, typename IT_>
        static void filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_>
        static void filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST
/// \endcond

#ifndef  __CUDACC__
#include <kernel/lafem/arch/unit_filter_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP
