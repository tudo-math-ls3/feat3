// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/component_product.hpp>

#include <mkl.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void ComponentProduct<Mem::Main>::value_mkl(float * r, const float * const x, const float * const y, const Index size)
{
  vsMul((MKL_INT)size, x, y, r);
}

void ComponentProduct<Mem::Main>::value_mkl(double * r, const double * const x, const double * const y, const Index size)
{
  vdMul((MKL_INT)size, x, y, r);
}
