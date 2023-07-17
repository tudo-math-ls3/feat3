// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_blocked.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

/// \cond internal
template void UnitFilterBlocked::filter_rhs_generic<2, float, std::uint64_t>(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<2, double, std::uint64_t>(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<2, float, std::uint32_t>(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<2, double, std::uint32_t>(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

template void UnitFilterBlocked::filter_def_generic<2, float, std::uint64_t>(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<2, double, std::uint64_t>(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<2, float, std::uint32_t>(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<2, double, std::uint32_t>(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

template void UnitFilterBlocked::filter_rhs_generic<3, float, std::uint64_t>(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<3, double, std::uint64_t>(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<3, float, std::uint32_t>(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_generic<3, double, std::uint32_t>(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);

template void UnitFilterBlocked::filter_def_generic<3, float, std::uint64_t>(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<3, double, std::uint64_t>(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<3, float, std::uint32_t>(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
template void UnitFilterBlocked::filter_def_generic<3, double, std::uint32_t>(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue, bool ign_nans);
/// \endcond
