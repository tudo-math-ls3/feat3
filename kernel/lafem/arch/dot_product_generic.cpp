// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/dot_product.hpp>

#ifdef FEAST_GMP
#include <gmpxx.h>
#include <mpfr.h>
#endif

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
DT_ DotProduct<Mem::Main, Algo::Generic>::value(const DT_ * const x, const DT_ * const y, const Index size)
{
  DT_ r(0);

  if(x == y)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r += x[i] * x[i];
    }
  }
  else
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r += x[i] * y[i];
    }
  }

  return r;
}

template float DotProduct<Mem::Main, Algo::Generic>::value(const float * const, const float * const, const Index);
template double DotProduct<Mem::Main, Algo::Generic>::value(const double * const, const double * const, const Index);
#ifdef FEAST_GMP
template mpf_class DotProduct<Mem::Main, Algo::Generic>::value(const mpf_class * const, const mpf_class * const, const Index);
#endif
