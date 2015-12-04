// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/gather_axpy_prim.hpp>

#include <cstring>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void GatherAxpyPrim<Mem::Main>::dv_csr_generic(float*, const float*, const unsigned long*, const float*, const unsigned long*, const float alpha, const Index, const Index);
template void GatherAxpyPrim<Mem::Main>::dv_csr_generic(double*, const double*, const unsigned long*, const double*, const unsigned long*, const double alpha, const Index, const Index);
template void GatherAxpyPrim<Mem::Main>::dv_csr_generic(float*, const float*, const unsigned int*, const float*, const unsigned int*, const float alpha, const Index, const Index);
template void GatherAxpyPrim<Mem::Main>::dv_csr_generic(double*, const double*, const unsigned int*, const double*, const unsigned int*, const double alpha, const Index, const Index);
