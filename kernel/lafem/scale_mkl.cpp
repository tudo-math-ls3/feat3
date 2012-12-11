// includes, FEAST
#include <kernel/lafem/scale.hpp>
#include <kernel/lafem/algorithm.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void Scale<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const float s)
{
  if (r.elements() == x.elements())
  {
    cblas_sscal(r.size(), s, r.elements(), 1);
  }
  else
  {
    copy(r, x);
    cblas_sscal(r.size(), s, r.elements(), 1);
  }
}

void Scale<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const double s)
{
  if (r.elements() == x.elements())
  {
    cblas_dscal(r.size(), s, r.elements(), 1);
  }
  else
  {
    copy(r, x);
    cblas_dscal(r.size(), s, r.elements(), 1);
  }
}
