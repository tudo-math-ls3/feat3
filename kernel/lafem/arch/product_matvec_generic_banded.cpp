// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>

#include <cstring>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void ProductMatVec<Mem::Main, Algo::Generic>::banded(float *, const float * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::banded(double *, const double * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::banded(float *, const float * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::Main, Algo::Generic>::banded(double *, const double * const, const unsigned int * const, const double * const, const Index, const Index, const Index);
