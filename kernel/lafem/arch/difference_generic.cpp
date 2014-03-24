// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/difference.hpp>


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void Difference<Mem::Main, Algo::Generic>::value(float *, const float * const, const float * const, const Index);
template void Difference<Mem::Main, Algo::Generic>::value(double *, const double * const, const double * const, const Index);
