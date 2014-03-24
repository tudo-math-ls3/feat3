// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/component_product.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void ComponentProduct<Mem::Main, Algo::Generic>::value(float *, const float * const, const float * const, const Index);
template void ComponentProduct<Mem::Main, Algo::Generic>::value(double *, const double * const, const double * const, const Index);
