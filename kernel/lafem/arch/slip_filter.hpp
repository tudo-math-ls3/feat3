#pragma once
#ifndef KERNEL_LAFEM_ARCH_SLIP_FILTER_HPP
#define KERNEL_LAFEM_ARCH_SLIP_FILTER_HPP 1

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
      struct SlipFilter;

      template <>
      struct SlipFilter<Mem::Main>
      {
        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue)
        {
          filter_rhs_generic<DT_, IT_, BlockSize_>(v, nu_elements, sv_indices, ue);
        }

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue)
        {
          filter_def_generic<DT_, IT_, BlockSize_>(v, nu_elements, sv_indices, ue);
        }

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs_generic(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def_generic(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue);
      }; // SlipFilter<Mem::Main>

      // Do not instantiate the following templates as this is done in slip_filter_generic.cpp and then linked
      // into the shared library
      extern template void SlipFilter<Mem::Main>::filter_rhs_generic<float, unsigned long, 2>(float * v, const float * const nu_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_rhs_generic<double, unsigned long, 2>(double * v, const double * const nu_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_rhs_generic<float, unsigned int, 2>(float * v, const float * const nu_elements, const unsigned int * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_rhs_generic<double, unsigned int, 2>(double * v, const double * const nu_elements, const unsigned int * const sv_indices, const Index ue);

      extern template void SlipFilter<Mem::Main>::filter_def_generic<float, unsigned long, 2>(float * v, const float * const nu_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_def_generic<double, unsigned long, 2>(double * v, const double * const nu_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_def_generic<float, unsigned int, 2>(float * v, const float * const nu_elements, const unsigned int * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_def_generic<double, unsigned int, 2>(double * v, const double * const nu_elements, const unsigned int * const sv_indices, const Index ue);

      extern template void SlipFilter<Mem::Main>::filter_rhs_generic<float, unsigned long, 3>(float * v, const float * const nu_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_rhs_generic<double, unsigned long, 3>(double * v, const double * const nu_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_rhs_generic<float, unsigned int, 3>(float * v, const float * const nu_elements, const unsigned int * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_rhs_generic<double, unsigned int, 3>(double * v, const double * const nu_elements, const unsigned int * const sv_indices, const Index ue);

      extern template void SlipFilter<Mem::Main>::filter_def_generic<float, unsigned long, 3>(float * v, const float * const nu_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_def_generic<double, unsigned long, 3>(double * v, const double * const nu_elements, const unsigned long * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_def_generic<float, unsigned int, 3>(float * v, const float * const nu_elements, const unsigned int * const sv_indices, const Index ue);
      extern template void SlipFilter<Mem::Main>::filter_def_generic<double, unsigned int, 3>(double * v, const double * const nu_elements, const unsigned int * const sv_indices, const Index ue);

      template <>
      struct SlipFilter<Mem::CUDA>
      {
        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue);

      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

/// \endcond
#ifndef  __CUDACC__
#include <kernel/lafem/arch/slip_filter_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_SLIP_FILTER_HPP
