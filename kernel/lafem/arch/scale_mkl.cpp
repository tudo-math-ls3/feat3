// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/scale.hpp>

#include <cstring>
#include <mkl.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Scale<Mem::Main>::value_mkl(float * r, const float * const x, const float s, const Index size)
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

void Scale<Mem::Main>::value_mkl(double * r, const double * const x, const double s, const Index size)
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
