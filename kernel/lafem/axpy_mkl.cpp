// includes, FEAST
#include <kernel/lafem/axpy.hpp>
#include <kernel/lafem/algorithm.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void Axpy<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const float a, const DenseVector<Mem::Main, float> & x, const DenseVector<Mem::Main, float> & y)
{
  if (r.elements() == y.elements())
  {
    cblas_saxpy(x.size(), a, x.elements(), 1, r.elements(), 1);
  }
  else if (r.elements() == x.elements())
  {
    DenseVector<Mem::Main, float> temp(r.size());
    copy(temp, y);
    cblas_saxpy(x.size(), a, x.elements(), 1, temp.elements(), 1);
    r = temp;
  }
  else
  {
    copy(r, y);
    cblas_saxpy(x.size(), a, x.elements(), 1, r.elements(), 1);
  }
}

void Axpy<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const double a, const DenseVector<Mem::Main, double> & x, const DenseVector<Mem::Main, double> & y)
{
  if (r.elements() == y.elements())
  {
    cblas_daxpy(x.size(), a, x.elements(), 1, r.elements(), 1);
  }
  else if (r.elements() == x.elements())
  {
    DenseVector<Mem::Main, double> temp(r.size());
    copy(temp, y);
    cblas_daxpy(x.size(), a, x.elements(), 1, temp.elements(), 1);
    r = temp;
  }
  else
  {
    copy(r, y);
    cblas_daxpy(x.size(), a, x.elements(), 1, r.elements(), 1);
  }
}
