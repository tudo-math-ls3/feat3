// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/gather_axpy_prim.hpp>

#include <cstring>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void GatherAxpyPrim<Mem::Main>::dv_csr_generic(float*, const float*, const Index*, const float*, const Index*, const float alpha, const Index, const Index);
template void GatherAxpyPrim<Mem::Main>::dv_csr_generic(double*, const double*, const Index*, const double*, const Index*, const double alpha, const Index, const Index);
