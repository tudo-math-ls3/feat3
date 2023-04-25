// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/slip_filter.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

/// \cond internal
template void SlipFilter::filter_rhs_generic<2, float, std::uint64_t>(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
template void SlipFilter::filter_rhs_generic<2, double, std::uint64_t>(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
template void SlipFilter::filter_rhs_generic<2, float, std::uint32_t>(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
template void SlipFilter::filter_rhs_generic<2, double, std::uint32_t>(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);

template void SlipFilter::filter_def_generic<2, float, std::uint64_t>(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
template void SlipFilter::filter_def_generic<2, double, std::uint64_t>(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
template void SlipFilter::filter_def_generic<2, float, std::uint32_t>(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
template void SlipFilter::filter_def_generic<2, double, std::uint32_t>(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);

template void SlipFilter::filter_rhs_generic<3, float, std::uint64_t>(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
template void SlipFilter::filter_rhs_generic<3, double, std::uint64_t>(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
template void SlipFilter::filter_rhs_generic<3, float, std::uint32_t>(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
template void SlipFilter::filter_rhs_generic<3, double, std::uint32_t>(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);

template void SlipFilter::filter_def_generic<3, float, std::uint64_t>(float * v, const float * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
template void SlipFilter::filter_def_generic<3, double, std::uint64_t>(double * v, const double * const nu_elements, const std::uint64_t * const sv_indices, const Index ue);
template void SlipFilter::filter_def_generic<3, float, std::uint32_t>(float * v, const float * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);
template void SlipFilter::filter_def_generic<3, double, std::uint32_t>(double * v, const double * const nu_elements, const std::uint32_t * const sv_indices, const Index ue);

/// \endcond
