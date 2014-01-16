// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/sum.hpp>

#ifdef FEAST_GMP
#include <gmpxx.h>
#include <mpfr.h>
#endif

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Sum<Mem::Main, Algo::Generic>::value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size)
{
  if (r == x)
  {
  // hack for clang, not working correctly with mpf_class
#ifdef FEAST_GMP
    if (typeid(DT_) == typeid(mpf_class))
    {
      for (Index i(0) ; i < size ; ++i)
      {
        DT_ t(y[i]);
        t += r[i];
        r[i] = DT_(t);
      }
    }
    else
#endif
      for (Index i(0) ; i < size ; ++i)
        r[i] += y[i];
  }
else if (r == y)
{
  for (Index i(0) ; i < size ; ++i)
  {
    r[i] += x[i];
    }
  }
  else if (r == x && r == y)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] += r[i];
    }
  }
  else
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] = x[i] + y[i];
    }
  }
}

template void Sum<Mem::Main, Algo::Generic>::value(float *, const float * const, const float * const, const Index);
template void Sum<Mem::Main, Algo::Generic>::value(double *, const double * const, const double * const, const Index);
#ifdef FEAST_GMP
template void Sum<Mem::Main, Algo::Generic>::value(mpf_class *, const mpf_class * const, const mpf_class * const, const Index);
#endif
