// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/unit_filter.hpp>

/// \cond internal
using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void UnitFilter<Mem::Main>::filter_rhs_generic(float * v, const float * const sv_elements, const unsigned long * const sv_indices, const Index ue);
template void UnitFilter<Mem::Main>::filter_rhs_generic(double * v, const double * const sv_elements, const unsigned long * const sv_indices, const Index ue);
template void UnitFilter<Mem::Main>::filter_rhs_generic(float * v, const float * const sv_elements, const unsigned int * const sv_indices, const Index ue);
template void UnitFilter<Mem::Main>::filter_rhs_generic(double * v, const double * const sv_elements, const unsigned int * const sv_indices, const Index ue);

template void UnitFilter<Mem::Main>::filter_def_generic(float * v, const unsigned long * const sv_indices, const Index ue);
template void UnitFilter<Mem::Main>::filter_def_generic(double * v, const unsigned long * const sv_indices, const Index ue);
template void UnitFilter<Mem::Main>::filter_def_generic(float * v, const unsigned int * const sv_indices, const Index ue);
template void UnitFilter<Mem::Main>::filter_def_generic(double * v, const unsigned int * const sv_indices, const Index ue);
/// \endcond
