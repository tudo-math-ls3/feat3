// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_SLIP_FILTER_HPP
#define KERNEL_LAFEM_ARCH_SLIP_FILTER_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct SlipFilter
      {
        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue)
        {
          filter_rhs_generic<DT_, IT_, BlockSize_>(v, nu_elements, sv_indices, ue);
        }

        template <int BlockSize_>
        static void filter_rhs(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue)
        {
          constexpr auto filter_rhs_generic_float_u64 = &filter_rhs_generic<float, std::uint64_t, BlockSize_>;
          constexpr auto filter_rhs_cuda_float_u64 = &filter_rhs_cuda<float, std::uint64_t, BlockSize_>;
          BACKEND_SKELETON_VOID(filter_rhs_cuda_float_u64, filter_rhs_generic_float_u64, filter_rhs_generic_float_u64, v, nu_elements, sv_indices, ue)
        }

        template <int BlockSize_>
        static void filter_rhs(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue)
        {
          constexpr auto filter_rhs_generic_double_u64 = &filter_rhs_generic<double, std::uint64_t, BlockSize_>;
          constexpr auto filter_rhs_cuda_double_u64 = &filter_rhs_cuda<double, std::uint64_t, BlockSize_>;
          BACKEND_SKELETON_VOID(filter_rhs_cuda_double_u64, filter_rhs_generic_double_u64, filter_rhs_generic_double_u64, v, nu_elements, sv_indices, ue)
        }

        template <int BlockSize_>
        static void filter_rhs(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue)
        {
          constexpr auto filter_rhs_generic_float_u32 = &filter_rhs_generic<float, std::uint32_t, BlockSize_>;
          constexpr auto filter_rhs_cuda_float_u32 = &filter_rhs_cuda<float, std::uint32_t, BlockSize_>;
          BACKEND_SKELETON_VOID(filter_rhs_cuda_float_u32, filter_rhs_generic_float_u32, filter_rhs_generic_float_u32, v, nu_elements, sv_indices, ue)
        }

        template <int BlockSize_>
        static void filter_rhs(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue)
        {
          constexpr auto filter_rhs_generic_double_u32 = &filter_rhs_generic<double, std::uint32_t, BlockSize_>;
          constexpr auto filter_rhs_cuda_double_u32 = &filter_rhs_cuda<double, std::uint32_t, BlockSize_>;
          BACKEND_SKELETON_VOID(filter_rhs_cuda_double_u32, filter_rhs_generic_double_u32, filter_rhs_generic_double_u32, v, nu_elements, sv_indices, ue)
        }

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue)
        {
          filter_def_generic<DT_, IT_, BlockSize_>(v, nu_elements, sv_indices, ue);
        }

        template <int BlockSize_>
        static void filter_def(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue)
        {
          constexpr auto filter_def_generic_float_u64 = &filter_def_generic<float, std::uint64_t, BlockSize_>;
          constexpr auto filter_def_cuda_float_u64 = &filter_def_cuda<float, std::uint64_t, BlockSize_>;
          BACKEND_SKELETON_VOID(filter_def_cuda_float_u64, filter_def_generic_float_u64, filter_def_generic_float_u64, v, nu_elements, sv_indices, ue)
        }

        template <int BlockSize_>
        static void filter_def(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue)
        {
          constexpr auto filter_def_generic_double_u64 = &filter_def_generic<double, std::uint64_t, BlockSize_>;
          constexpr auto filter_def_cuda_double_u64 = &filter_def_cuda<double, std::uint64_t, BlockSize_>;
          BACKEND_SKELETON_VOID(filter_def_cuda_double_u64, filter_def_generic_double_u64, filter_def_generic_double_u64, v, nu_elements, sv_indices, ue)
        }

        template <int BlockSize_>
        static void filter_def(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue)
        {
          constexpr auto filter_def_generic_float_u32 = &filter_def_generic<float, std::uint32_t, BlockSize_>;
          constexpr auto filter_def_cuda_float_u32 = &filter_def_cuda<float, std::uint32_t, BlockSize_>;
          BACKEND_SKELETON_VOID(filter_def_cuda_float_u32, filter_def_generic_float_u32, filter_def_generic_float_u32, v, nu_elements, sv_indices, ue)
        }

        template <int BlockSize_>
        static void filter_def(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue)
        {
          constexpr auto filter_def_generic_double_u32 = &filter_def_generic<double, std::uint32_t, BlockSize_>;
          constexpr auto filter_def_cuda_double_u32 = &filter_def_cuda<double, std::uint32_t, BlockSize_>;
          BACKEND_SKELETON_VOID(filter_def_cuda_double_u32, filter_def_generic_double_u32, filter_def_generic_double_u32, v, nu_elements, sv_indices, ue)
        }

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs_generic(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def_generic(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_rhs_cuda(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_, int BlockSize_>
        static void filter_def_cuda(DT_ * v, const DT_ * const nu_elements, const IT_ * const sv_indices, const Index ue);
      }; // SlipFilter

      // Do not instantiate the following templates as this is done in slip_filter_generic.cpp and then linked
      // into the shared library
#ifdef FEAT_EICKT
      extern template void SlipFilter::filter_rhs_generic<float, std::uint64_t, 2>(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_rhs_generic<double, std::uint64_t, 2>(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_rhs_generic<float, std::uint32_t, 2>(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_rhs_generic<double, std::uint32_t, 2>(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);

      extern template void SlipFilter::filter_def_generic<float, std::uint64_t, 2>(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_def_generic<double, std::uint64_t, 2>(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_def_generic<float, std::uint32_t, 2>(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_def_generic<double, std::uint32_t, 2>(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);

      extern template void SlipFilter::filter_rhs_generic<float, std::uint64_t, 3>(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_rhs_generic<double, std::uint64_t, 3>(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_rhs_generic<float, std::uint32_t, 3>(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_rhs_generic<double, std::uint32_t, 3>(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);

      extern template void SlipFilter::filter_def_generic<float, std::uint64_t, 3>(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_def_generic<double, std::uint64_t, 3>(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_def_generic<float, std::uint32_t, 3>(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
      extern template void SlipFilter::filter_def_generic<double, std::uint32_t, 3>(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
#endif
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

/// \endcond
#ifndef  __CUDACC__
#include <kernel/lafem/arch/slip_filter_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_SLIP_FILTER_HPP
