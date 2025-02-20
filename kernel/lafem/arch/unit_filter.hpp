// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP 1

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
      struct UnitFilter
      {
        template <typename DT_, typename IT_>
        static void filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
        {
          filter_rhs_generic(v, sv_elements, sv_indices, ue);
        }

        static void filter_rhs(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_rhs_cuda, filter_rhs_generic, filter_rhs_generic, v, sv_elements, sv_indices, ue)
        }

        static void filter_rhs(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_rhs_cuda, filter_rhs_generic, filter_rhs_generic, v, sv_elements, sv_indices, ue)
        }

        static void filter_rhs(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_rhs_cuda, filter_rhs_mkl, filter_rhs_generic, v, sv_elements, sv_indices, ue)
        }

        static void filter_rhs(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_rhs_cuda, filter_rhs_mkl, filter_rhs_generic, v, sv_elements, sv_indices, ue)
        }

        template <typename DT_, typename IT_>
        static void filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue)
        {
          filter_def_generic(v, sv_indices, ue);
        }

        template <typename IT_>
        static void filter_def(float * v, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_def_cuda, filter_def_generic, filter_def_generic, v, sv_indices, ue)
        }

        template <typename IT_>
        static void filter_def(double * v, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_def_cuda, filter_def_generic, filter_def_generic, v, sv_indices, ue)
        }

        template<typename DT_, typename IT_>
        static void filter_unit_mat(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const sv_indices, const Index ue)
        {
          filter_unit_mat_generic(mat, row_ptr, col_idx, sv_indices, ue);
        }

        template<typename IT_>
        static void filter_unit_mat(double* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_unit_mat_cuda, filter_unit_mat_generic, filter_unit_mat_generic, mat, row_ptr, col_idx, sv_indices, ue)
        }

        template<typename IT_>
        static void filter_unit_mat(float* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_unit_mat_cuda, filter_unit_mat_generic, filter_unit_mat_generic, mat, row_ptr, col_idx, sv_indices, ue)
        }

        template<typename DT_, typename IT_>
        static void filter_offdiag_row_mat(DT_* mat, const IT_* const row_ptr, int block_width, const IT_ * const sv_indices, const Index ue)
        {
          filter_offdiag_row_mat_generic(mat, row_ptr, block_width, sv_indices, ue);
        }

        template<typename IT_>
        static void filter_offdiag_row_mat(double* mat, const IT_* const row_ptr, int block_width, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_offdiag_row_mat_cuda, filter_offdiag_row_mat_generic, filter_offdiag_row_mat_generic, mat, row_ptr, block_width, sv_indices, ue)
        }

        template<typename IT_>
        static void filter_offdiag_row_mat(float* mat, const IT_* const row_ptr, int block_width, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_offdiag_row_mat_cuda, filter_offdiag_row_mat_generic, filter_offdiag_row_mat_generic, mat, row_ptr, block_width, sv_indices, ue)
        }

        template<typename DT_, typename IT_>
        static void filter_weak_matrix_rows(DT_* mat_a, const DT_* const mat_m, const IT_* const row_ptr, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
        {
          filter_weak_matrix_rows_generic(mat_a, mat_m, row_ptr, sv_elements, sv_indices, ue);
        }

        template<typename IT_>
        static void filter_weak_matrix_rows(double* mat_a, const double* const mat_m, const IT_* const row_ptr, const double * const sv_elements, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_weak_matrix_rows_cuda, filter_weak_matrix_rows_generic, filter_weak_matrix_rows_generic, mat_a, mat_m, row_ptr, sv_elements, sv_indices, ue)
        }

        template<typename IT_>
        static void filter_weak_matrix_rows(float* mat_a, const float* const mat_m, const IT_* const row_ptr, const float * const sv_elements, const IT_ * const sv_indices, const Index ue)
        {
          BACKEND_SKELETON_VOID(filter_weak_matrix_rows_cuda, filter_weak_matrix_rows_generic, filter_weak_matrix_rows_generic, mat_a, mat_m, row_ptr, sv_elements, sv_indices, ue)
        }

        template <typename DT_, typename IT_>
        static void filter_rhs_generic(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_>
        static void filter_def_generic(DT_ * v, const IT_ * const sv_indices, const Index ue);

        static void filter_rhs_mkl(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
        static void filter_rhs_mkl(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);

        template <typename DT_, typename IT_>
        static void filter_rhs_cuda(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template <typename DT_, typename IT_>
        static void filter_def_cuda(DT_ * v, const IT_ * const sv_indices, const Index ue);

        template<typename DT_, typename IT_>
        static void filter_unit_mat_generic(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const sv_indices, const Index ue);

        template<typename DT_, typename IT_>
        static void filter_unit_mat_cuda(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const sv_indices, const Index ue);

        template<typename DT_, typename IT_>
        static void filter_offdiag_row_mat_generic(DT_* mat, const IT_* const row_ptr, int block_width, const IT_ * const sv_indices, const Index ue);

        template<typename DT_, typename IT_>
        static void filter_offdiag_row_mat_cuda(DT_* mat, const IT_* const row_ptr, int block_width, const IT_ * const sv_indices, const Index ue);

        template<typename DT_, typename IT_>
        static void filter_weak_matrix_rows_generic(DT_* mat_a, const DT_* const mat_m, const IT_* const row_ptr, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);

        template<typename DT_, typename IT_>
        static void filter_weak_matrix_rows_cuda(DT_* mat_a, const DT_* const mat_m, const IT_* const row_ptr, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue);
      };

#ifdef FEAT_EICKT
      extern template void UnitFilter::filter_rhs_generic(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void UnitFilter::filter_rhs_generic(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
      extern template void UnitFilter::filter_rhs_generic(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue);
      extern template void UnitFilter::filter_rhs_generic(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue);

      extern template void UnitFilter::filter_def_generic(float * v, const std::uint64_t * const sv_indices, const Index ue);
      extern template void UnitFilter::filter_def_generic(double * v, const std::uint64_t * const sv_indices, const Index ue);
      extern template void UnitFilter::filter_def_generic(float * v, const std::uint32_t * const sv_indices, const Index ue);
      extern template void UnitFilter::filter_def_generic(double * v, const std::uint32_t * const sv_indices, const Index ue);
#endif
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
/// \endcond

#ifndef  __CUDACC__
#include <kernel/lafem/arch/unit_filter_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_UNIT_FILTER_HPP
