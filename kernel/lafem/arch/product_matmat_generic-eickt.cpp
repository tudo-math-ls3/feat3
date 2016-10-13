// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/product_matmat.hpp>

#include <cstring>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void ProductMatMat<Mem::Main>::dense_generic(float *, const float * const, const float * const, const Index, const Index, const Index);
template void ProductMatMat<Mem::Main>::dense_generic(double *, const double * const, const double * const, const Index, const Index, const Index);
