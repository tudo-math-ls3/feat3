// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter.hpp>

/// \cond internal
using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void UnitFilter::filter_rhs_generic(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
template void UnitFilter::filter_rhs_generic(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue);
template void UnitFilter::filter_rhs_generic(float * v, const float * const sv_elements, const std::uint32_t * const sv_indices, const Index ue);
template void UnitFilter::filter_rhs_generic(double * v, const double * const sv_elements, const std::uint32_t * const sv_indices, const Index ue);

template void UnitFilter::filter_def_generic(float * v, const std::uint64_t * const sv_indices, const Index ue);
template void UnitFilter::filter_def_generic(double * v, const std::uint64_t * const sv_indices, const Index ue);
template void UnitFilter::filter_def_generic(float * v, const std::uint32_t * const sv_indices, const Index ue);
template void UnitFilter::filter_def_generic(double * v, const std::uint32_t * const sv_indices, const Index ue);
/// \endcond
