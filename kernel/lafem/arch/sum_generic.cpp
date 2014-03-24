// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/sum.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void Sum<Mem::Main, Algo::Generic>::value(float *, const float * const, const float * const, const Index);
template void Sum<Mem::Main, Algo::Generic>::value(double *, const double * const, const double * const, const Index);
