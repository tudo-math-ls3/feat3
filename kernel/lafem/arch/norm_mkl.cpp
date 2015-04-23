// includes, FEAST
#include <kernel/lafem/arch/norm.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

float Norm2<Mem::Main>::value_mkl(const float * const x, const Index size)
{
  return cblas_snrm2((MKL_INT)size, x, 1);
}

double Norm2<Mem::Main>::value_mkl(const double * const x, const Index size)
{
  return cblas_dnrm2((MKL_INT)size, x, 1);
}
