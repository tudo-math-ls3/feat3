// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct UnitFilterBlocked
      {
        template <int BlockSize_, typename DT_, typename IT_>
        static void filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          filter_rhs_generic<BlockSize_, DT_, IT_>(v, sv_elements, sv_indices, ue, ign_nans);
        }

        template <int BlockSize_>
        static void filter_rhs(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID_T1(BlockSize_, filter_rhs_cuda, filter_rhs_generic, filter_rhs_generic, v, sv_elements, sv_indices, ue, ign_nans)
        }

        template <int BlockSize_>
        static void filter_rhs(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID_T1(BlockSize_, filter_rhs_cuda, filter_rhs_generic, filter_rhs_generic, v, sv_elements, sv_indices, ue, ign_nans)
        }

        template <int BlockSize_>
        static void filter_rhs(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID_T1(BlockSize_, filter_rhs_cuda, filter_rhs_generic, filter_rhs_generic, v, sv_elements, sv_indices, ue, ign_nans)
        }

        template <int BlockSize_>
        static void filter_rhs(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID_T1(BlockSize_, filter_rhs_cuda, filter_rhs_generic, filter_rhs_generic, v, sv_elements, sv_indices, ue, ign_nans)
        }

        template <int BlockSize_, typename DT_, typename IT_>
        static void filter_def(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          filter_def_generic<BlockSize_, DT_, IT_>(v, sv_elements, sv_indices, ue, ign_nans);
        }

        template <int BlockSize_>
        static void filter_def(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID_T1(BlockSize_, filter_def_cuda, filter_def_generic, filter_def_generic, v, sv_elements, sv_indices, ue, ign_nans)
        }

        template <int BlockSize_>
        static void filter_def(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID_T1(BlockSize_, filter_def_cuda, filter_def_generic, filter_def_generic, v, sv_elements, sv_indices, ue, ign_nans)
        }

        template <int BlockSize_>
        static void filter_def(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID_T1(BlockSize_, filter_def_cuda, filter_def_generic, filter_def_generic, v, sv_elements, sv_indices, ue, ign_nans)
        }

        template <int BlockSize_>
        static void filter_def(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID_T1(BlockSize_, filter_def_cuda, filter_def_generic, filter_def_generic, v, sv_elements, sv_indices, ue, ign_nans)
        }

        template <int BlockSize_, typename DT_, typename IT_>
        static void filter_rhs_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template <int BlockSize_, typename DT_, typename IT_>
        static void filter_def_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template <int BlockSize_, typename DT_, typename IT_>
        static void filter_rhs_cuda(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template <int BlockSize_, typename DT_, typename IT_>
        static void filter_def_cuda(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);
      };

      // Do not instantiate the following templates as this is done in unit_filter_blocked_generic.cpp and then linked
      // into the shared library
#ifdef FEAT_EICKT
      extern template void UnitFilterBlocked::filter_rhs_generic<2, float, std::uint64_t>(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<2, double, std::uint64_t>(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<2, float, std::uint32_t>(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<2, double, std::uint32_t>(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

      extern template void UnitFilterBlocked::filter_def_generic<2, float, std::uint64_t>(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<2, double, std::uint64_t>(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<2, float, std::uint32_t>(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<2, double, std::uint32_t>(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

      extern template void UnitFilterBlocked::filter_rhs_generic<3, float, std::uint64_t>(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<3, double, std::uint64_t>(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<3, float, std::uint32_t>(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<3, double, std::uint32_t>(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

      extern template void UnitFilterBlocked::filter_def_generic<3, float, std::uint64_t>(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<3, double, std::uint64_t>(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<3, float, std::uint32_t>(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<3, double, std::uint32_t>(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

/// \endcond
#ifndef  __CUDACC__
#include <kernel/lafem/arch/unit_filter_blocked_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
