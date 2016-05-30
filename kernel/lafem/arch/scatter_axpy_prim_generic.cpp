// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/scatter_axpy_prim.hpp>

#include <cstring>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(float*, const float*, const unsigned long*, const float*, const unsigned long*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(double*, const double*, const unsigned long*, const double*, const unsigned long*, const double alpha, const Index);
template void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(float*, const float*, const unsigned int*, const float*, const unsigned int*, const float alpha, const Index);
template void ScatterAxpyPrim<Mem::Main>::dv_csr_generic(double*, const double*, const unsigned int*, const double*, const unsigned int*, const double alpha, const Index);
