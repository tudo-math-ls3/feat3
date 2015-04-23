// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/scale.hpp>


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void Scale<Mem::Main>::value_generic(float *, const float * const, const float, const Index);
template void Scale<Mem::Main>::value_generic(double *, const double * const, const double, const Index);
