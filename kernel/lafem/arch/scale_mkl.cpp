// includes, FEAST
#include <kernel/lafem/arch/scale.hpp>
#include <cstring>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

void Scale<Mem::Main, Algo::MKL>::value(float * r, const float * const x, const float s, const Index size)
{
  if (r == x)
  {
    cblas_sscal((MKL_INT)size, s, r, 1);
  }
  else
  {
    std::memcpy(r, x, size * sizeof(float));
    cblas_sscal((MKL_INT)size, s, r, 1);
  }
}

void Scale<Mem::Main, Algo::MKL>::value(double * r, const double * const x, const double s, const Index size)
{
  if (r == x)
  {
    cblas_dscal((MKL_INT)size, s, r, 1);
  }
  else
  {
    std::memcpy(r, x, size * sizeof(double));
    cblas_dscal((MKL_INT)size, s, r, 1);
  }
}
