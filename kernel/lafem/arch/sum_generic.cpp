// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/sum.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void Sum<Mem::Main>::value_generic(float *, const float * const, const float * const, const Index);
template void Sum<Mem::Main>::value_generic(double *, const double * const, const double * const, const Index);
