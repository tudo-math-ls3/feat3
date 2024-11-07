// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_blocked.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

/// \cond internal
template void UnitFilterBlocked::filter_rhs_generic<float, std::uint64_t>(float * v, int block_size, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<double, std::uint64_t>(double * v, int block_size, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<float, std::uint32_t>(float * v, int block_size, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<double, std::uint32_t>(double * v, int block_size, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

template void UnitFilterBlocked::filter_def_generic<float, std::uint64_t>(float * v, int block_size, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<double, std::uint64_t>(double * v, int block_size, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<float, std::uint32_t>(float * v, int block_size, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<double, std::uint32_t>(double * v, int block_size, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

template void UnitFilterBlocked::filter_unit_mat_generic<float, std::uint64_t>(float* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, int block_height, int block_width, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_unit_mat_generic<double, std::uint64_t>(double* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, int block_height, int block_width, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_unit_mat_generic<float, std::uint32_t>(float* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, int block_height, int block_width, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_unit_mat_generic<double, std::uint32_t>(double* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, int block_height, int block_width, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

template void UnitFilterBlocked::filter_offdiag_row_mat_generic<float, std::uint64_t>(float* mat, const std::uint64_t* const row_ptr, int block_height, int block_width, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_offdiag_row_mat_generic<double, std::uint64_t>(double* mat, const std::uint64_t* const row_ptr, int block_height, int block_width, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_offdiag_row_mat_generic<float, std::uint32_t>(float* mat, const std::uint32_t* const row_ptr, int block_height, int block_width, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_offdiag_row_mat_generic<double, std::uint32_t>(double* mat, const std::uint32_t* const row_ptr, int block_height, int block_width, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

template void UnitFilterBlocked::filter_weak_matrix_rows_generic<float, std::uint64_t>(float* mat_a, const float* const mat_m, const std::uint64_t* const row_ptr, int block_height, int block_width, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
template void UnitFilterBlocked::filter_weak_matrix_rows_generic<double, std::uint64_t>(double* mat_a, const double* const mat_m, const std::uint64_t* const row_ptr, int block_height, int block_width, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
template void UnitFilterBlocked::filter_weak_matrix_rows_generic<float, std::uint32_t>(float* mat_a, const float* const mat_m, const std::uint32_t* const row_ptr, int block_height, int block_width, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue);
template void UnitFilterBlocked::filter_weak_matrix_rows_generic<double, std::uint32_t>(double* mat_a, const double* const mat_m, const std::uint32_t* const row_ptr, int block_height, int block_width, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue);
/// \endcond
