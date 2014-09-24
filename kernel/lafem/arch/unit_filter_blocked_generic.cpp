// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/unit_filter.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(float * v, const float * const sv_elements, const unsigned long * const sv_indices, const Index ue, const Index bs);
template void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(double * v, const double * const sv_elements, const unsigned long * const sv_indices, const Index ue, const Index bs);
template void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(float * v, const float * const sv_elements, const unsigned int * const sv_indices, const Index ue, const Index bs);
template void UnitFilter<Mem::Main, Algo::Generic>::filter_rhs(double * v, const double * const sv_elements, const unsigned int * const sv_indices, const Index ue, const Index bs);

template void UnitFilter<Mem::Main, Algo::Generic>::filter_def(float * v, const unsigned long * const sv_indices, const Index ue, const Index bs);
template void UnitFilter<Mem::Main, Algo::Generic>::filter_def(double * v, const unsigned long * const sv_indices, const Index ue, const Index bs);
template void UnitFilter<Mem::Main, Algo::Generic>::filter_def(float * v, const unsigned int * const sv_indices, const Index ue, const Index bs);
template void UnitFilter<Mem::Main, Algo::Generic>::filter_def(double * v, const unsigned int * const sv_indices, const Index ue, const Index bs);
