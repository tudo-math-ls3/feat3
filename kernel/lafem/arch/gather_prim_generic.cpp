// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/gather_prim.hpp>

#include <cstring>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void GatherPrim<Mem::Main>::dv_csr_generic(float*, const float*, const Index*, const float*, const Index*, const Index, const Index);
template void GatherPrim<Mem::Main>::dv_csr_generic(double*, const double*, const Index*, const double*, const Index*, const Index, const Index);
