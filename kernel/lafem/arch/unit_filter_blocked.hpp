#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct UnitFilterBlocked;

      template <>
      struct UnitFilterBlocked<Mem::Main>
      {
        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
        {
          filter_rhs_generic<DT_, IT_, BlockSize_>(v, sv_elements, sv_indices, ue);
        }

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue)
        {
          filter_def_generic<DT_, IT_, BlockSize_>(v, sv_indices, ue);
        }

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def_generic(DT_ * v, const IT_ * const sv_indices, const Index ue);
      };

      // Do not instantiate the following templates as this is done in unit_filter_blocked_generic.cpp and then linked
      // into the shared library
      extern template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<float, unsigned long, 2>(float * v, const float * const sv_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<double, unsigned long, 2>(double * v, const double * const sv_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<float, unsigned int, 2>(float * v, const float * const sv_elements, const unsigned int * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<double, unsigned int, 2>(double * v, const double * const sv_elements, const unsigned int * const sv_indices, const Index ue);

      extern template void UnitFilterBlocked<Mem::Main>::filter_def_generic<float, unsigned long, 2>(float * v, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_def_generic<double, unsigned long, 2>(double * v, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_def_generic<float, unsigned int, 2>(float * v, const unsigned int * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_def_generic<double, unsigned int, 2>(double * v, const unsigned int * const sv_indices, const Index ue);

      extern template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<float, unsigned long, 3>(float * v, const float * const sv_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<double, unsigned long, 3>(double * v, const double * const sv_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<float, unsigned int, 3>(float * v, const float * const sv_elements, const unsigned int * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<double, unsigned int, 3>(double * v, const double * const sv_elements, const unsigned int * const sv_indices, const Index ue);

      extern template void UnitFilterBlocked<Mem::Main>::filter_def_generic<float, unsigned long, 3>(float * v, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_def_generic<double, unsigned long, 3>(double * v, const unsigned long * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_def_generic<float, unsigned int, 3>(float * v, const unsigned int * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked<Mem::Main>::filter_def_generic<double, unsigned int, 3>(double * v, const unsigned int * const sv_indices, const Index ue);

      template <>
      struct UnitFilterBlocked<Mem::CUDA>
      {
        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

/// \endcond
#ifndef  __CUDACC__
#include <kernel/lafem/arch/unit_filter_blocked_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
