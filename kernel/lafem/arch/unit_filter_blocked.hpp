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
        template <typename DT_, typename IT_>
        static void filter_rhs(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          filter_rhs_generic(v, block_size, sv_elements, sv_indices, ue, ign_nans);
        }

        template <typename IT_>
        static void filter_rhs(double * v, int block_size, const double * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID(filter_rhs_cuda, filter_rhs_generic, filter_rhs_generic, v, block_size, sv_elements, sv_indices, ue, ign_nans)
        }

        template <typename IT_>
        static void filter_rhs(float * v, int block_size, const float * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID(filter_rhs_cuda, filter_rhs_generic, filter_rhs_generic, v, block_size, sv_elements, sv_indices, ue, ign_nans)
        }

        template <typename DT_, typename IT_>
        static void filter_def(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          filter_def_generic(v, block_size, sv_elements, sv_indices, ue, ign_nans);
        }

        template <typename IT_>
        static void filter_def(double * v, int block_size, const double * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID(filter_def_cuda, filter_def_generic, filter_def_generic, v, block_size, sv_elements, sv_indices, ue, ign_nans)
        }

        template <typename IT_>
        static void filter_def(float * v, int block_size, const float * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID(filter_def_cuda, filter_def_generic, filter_def_generic, v, block_size, sv_elements, sv_indices, ue, ign_nans)
        }

        template<typename DT_, typename IT_>
        static void filter_unit_mat(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          filter_unit_mat_generic(mat, row_ptr, col_idx, block_height, block_width, sv_elements, sv_indices, ue, ign_nans);
        }

        template<typename IT_>
        static void filter_unit_mat(double* mat, const IT_* const row_ptr, const IT_* const col_idx, int block_height, int block_width, const double * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID(filter_unit_mat_cuda, filter_unit_mat_generic, filter_unit_mat_generic, mat, row_ptr, col_idx, block_height, block_width, sv_elements, sv_indices, ue, ign_nans)
        }

        template<typename IT_>
        static void filter_unit_mat(float* mat, const IT_* const row_ptr, const IT_* const col_idx, int block_height, int block_width, const float * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID(filter_unit_mat_cuda, filter_unit_mat_generic, filter_unit_mat_generic, mat, row_ptr, col_idx, block_height, block_width, sv_elements, sv_indices, ue, ign_nans)
        }

        template<typename DT_, typename IT_>
        static void filter_offdiag_row_mat(DT_* mat, const IT_* const row_ptr, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          filter_offdiag_row_mat_generic(mat, row_ptr, block_height, block_width, sv_elements, sv_indices, ue, ign_nans);
        }

        template<typename IT_>
        static void filter_offdiag_row_mat(double* mat, const IT_* const row_ptr, int block_height, int block_width, const double * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID(filter_offdiag_row_mat_cuda, filter_offdiag_row_mat_generic, filter_offdiag_row_mat_generic, mat, row_ptr, block_height, block_width, sv_elements, sv_indices, ue, ign_nans)
        }

        template<typename IT_>
        static void filter_offdiag_row_mat(float* mat, const IT_* const row_ptr, int block_height, int block_width, const float * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
        {
          BACKEND_SKELETON_VOID(filter_offdiag_row_mat_cuda, filter_offdiag_row_mat_generic, filter_offdiag_row_mat_generic, mat, row_ptr, block_height, block_width, sv_elements, sv_indices, ue, ign_nans)
        }

        template<typename DT_, typename IT_>
        static void filter_weak_matrix_rows(DT_* mat_a, const DT_* const mat_m, const IT_* const row_ptr, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
        {
          filter_weak_matrix_rows_generic(mat_a, mat_m, row_ptr, block_height, block_width, sv_elements, sv_indices, ue);
        }

        template<typename IT_>
        static void filter_weak_matrix_rows(double* mat_a, const double* const mat_m, const IT_* const row_ptr, int block_height, int block_width, const double * const sv_elements, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_weak_matrix_rows_cuda, filter_weak_matrix_rows_generic, filter_weak_matrix_rows_generic, mat_a, mat_m, row_ptr, block_height, block_width, sv_elements, sv_indices, ue)
        }

        template<typename IT_>
        static void filter_weak_matrix_rows(float* mat_a, const float* const mat_m, const IT_* const row_ptr, int block_height, int block_width, const float * const sv_elements, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_weak_matrix_rows_cuda, filter_weak_matrix_rows_generic, filter_weak_matrix_rows_generic, mat_a, mat_m, row_ptr, block_height, block_width, sv_elements, sv_indices, ue)
        }

        template <typename DT_, typename IT_>
        static void filter_rhs_generic(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template <typename DT_, typename IT_>
        static void filter_def_generic(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template <typename DT_, typename IT_>
        static void filter_rhs_cuda(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template <typename DT_, typename IT_>
        static void filter_def_cuda(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template<typename DT_, typename IT_>
        static void filter_unit_mat_generic(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template<typename DT_, typename IT_>
        static void filter_unit_mat_cuda(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template<typename DT_, typename IT_>
        static void filter_offdiag_row_mat_generic(DT_* mat, const IT_* const row_ptr, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template<typename DT_, typename IT_>
        static void filter_offdiag_row_mat_cuda(DT_* mat, const IT_* const row_ptr, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans);

        template<typename DT_, typename IT_>
        static void filter_weak_matrix_rows_generic(DT_* mat_a, const DT_* const mat_m, const IT_* const row_ptr, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template<typename DT_, typename IT_>
        static void filter_weak_matrix_rows_cuda(DT_* mat_a, const DT_* const mat_m, const IT_* const row_ptr, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);
      };

      // Do not instantiate the following templates as this is done in unit_filter_blocked_generic.cpp and then linked
      // into the shared library
#ifdef FEAT_EICKT
      extern template void UnitFilterBlocked::filter_rhs_generic<float, std::uint64_t>(float * v, int block_size, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<double, std::uint64_t>(double * v, int block_size, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<float, std::uint32_t>(float * v, int block_size, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_rhs_generic<double, std::uint32_t>(double * v, int block_size, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

      extern template void UnitFilterBlocked::filter_def_generic<float, std::uint64_t>(float * v, int block_size, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<double, std::uint64_t>(double * v, int block_size, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<float, std::uint32_t>(float * v, int block_size, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_def_generic<double, std::uint32_t>(double * v, int block_size, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

      extern template void UnitFilterBlocked::filter_unit_mat_generic<float, std::uint64_t>(float* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, int block_height, int block_width, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_unit_mat_generic<double, std::uint64_t>(double* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, int block_height, int block_width, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_unit_mat_generic<float, std::uint32_t>(float* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, int block_height, int block_width, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_unit_mat_generic<double, std::uint32_t>(double* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, int block_height, int block_width, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

      extern template void UnitFilterBlocked::filter_offdiag_row_mat_generic<float, std::uint64_t>(float* mat, const std::uint64_t* const row_ptr, int block_height, int block_width, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_offdiag_row_mat_generic<double, std::uint64_t>(double* mat, const std::uint64_t* const row_ptr, int block_height, int block_width, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_offdiag_row_mat_generic<float, std::uint32_t>(float* mat, const std::uint32_t* const row_ptr, int block_height, int block_width, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
      extern template void UnitFilterBlocked::filter_offdiag_row_mat_generic<double, std::uint32_t>(double* mat, const std::uint32_t* const row_ptr, int block_height, int block_width, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

      extern template void UnitFilterBlocked::filter_weak_matrix_rows_generic<float, std::uint64_t>(float* mat_a, const float* const mat_m, const std::uint64_t* const row_ptr, int block_height, int block_width, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked::filter_weak_matrix_rows_generic<double, std::uint64_t>(double* mat_a, const double* const mat_m, const std::uint64_t* const row_ptr, int block_height, int block_width, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked::filter_weak_matrix_rows_generic<float, std::uint32_t>(float* mat_a, const float* const mat_m, const std::uint32_t* const row_ptr, int block_height, int block_width, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue);
      extern template void UnitFilterBlocked::filter_weak_matrix_rows_generic<double, std::uint32_t>(double* mat_a, const double* const mat_m, const std::uint32_t* const row_ptr, int block_height, int block_width, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

/// \endcond
#ifndef  __CUDACC__
#include <kernel/lafem/arch/unit_filter_blocked_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCKED_HPP
