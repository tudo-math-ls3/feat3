// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/sum.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

void Sum<Mem::Main>::value_mkl(float * r, const float * const x, const float * const y, const Index size)
{
  vsAdd((MKL_INT)size, x, y, r);
}

void Sum<Mem::Main>::value_mkl(double * r, const double * const x, const double * const y, const Index size)
{
  vdAdd((MKL_INT)size, x, y, r);
}
