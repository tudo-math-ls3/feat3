// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/component_product.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

void ComponentProduct<Mem::Main, Algo::MKL>::value(float * r, const float * const x, const float * const y, const Index size)
{
  vsMul((MKL_INT)size, x, y, r);
}

void ComponentProduct<Mem::Main, Algo::MKL>::value(double * r, const double * const x, const double * const y, const Index size)
{
  vdMul((MKL_INT)size, x, y, r);
}
