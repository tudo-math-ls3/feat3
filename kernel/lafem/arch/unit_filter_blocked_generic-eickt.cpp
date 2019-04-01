// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/unit_filter_blocked.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

/// \cond internal
template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<float, unsigned long, 2>(float * v, const float * const sv_elements, const unsigned long * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<double, unsigned long, 2>(double * v, const double * const sv_elements, const unsigned long * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<float, unsigned int, 2>(float * v, const float * const sv_elements, const unsigned int * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<double, unsigned int, 2>(double * v, const double * const sv_elements, const unsigned int * const sv_indices, const Index ue);

template void UnitFilterBlocked<Mem::Main>::filter_def_generic<float, unsigned long, 2>(float * v, const unsigned long * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_def_generic<double, unsigned long, 2>(double * v, const unsigned long * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_def_generic<float, unsigned int, 2>(float * v, const unsigned int * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_def_generic<double, unsigned int, 2>(double * v, const unsigned int * const sv_indices, const Index ue);

template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<float, unsigned long, 3>(float * v, const float * const sv_elements, const unsigned long * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<double, unsigned long, 3>(double * v, const double * const sv_elements, const unsigned long * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<float, unsigned int, 3>(float * v, const float * const sv_elements, const unsigned int * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_rhs_generic<double, unsigned int, 3>(double * v, const double * const sv_elements, const unsigned int * const sv_indices, const Index ue);

template void UnitFilterBlocked<Mem::Main>::filter_def_generic<float, unsigned long, 3>(float * v, const unsigned long * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_def_generic<double, unsigned long, 3>(double * v, const unsigned long * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_def_generic<float, unsigned int, 3>(float * v, const unsigned int * const sv_indices, const Index ue);
template void UnitFilterBlocked<Mem::Main>::filter_def_generic<double, unsigned int, 3>(double * v, const unsigned int * const sv_indices, const Index ue);
/// \endcond
