// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/dot_product.hpp>

#include <mkl.h>


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

float DotProduct<Mem::Main>::value_mkl(const float * const x, const float * const y, const Index size)
{
  return cblas_sdot((MKL_INT)size, x, 1, y, 1);
}

double DotProduct<Mem::Main>::value_mkl(const double * const x, const double * const y, const Index size)
{
  return cblas_ddot((MKL_INT)size, x, 1, y, 1);
}

float TripleDotProduct<Mem::Main>::value_mkl(const float * const x, const float * const y, const float * const z, const Index size)
{
  float * temp(new float[size]);
  vsMul((MKL_INT)size, y, z, temp);
  float result = cblas_sdot((MKL_INT)size, x, 1, temp, 1);
  delete[] temp;
  return result;
}

double TripleDotProduct<Mem::Main>::value_mkl(const double * const x, const double * const y, const double * const z, const Index size)
{
  double * temp(new double[size]);
  vdMul((MKL_INT)size, y, z, temp);
  double result = cblas_ddot((MKL_INT)size, x, 1, temp, 1);
  delete[] temp;
  return result;
}
