// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/component_product.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void ComponentProduct<Mem::Main, Algo::Generic>::value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size)
{
  if (r == x)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] *= y[i];
    }
  }
  else if (r == y)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] *= x[i];
    }
  }
  else if (r == x && r == y)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] *= r[i];
    }
  }
  else
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] = x[i] * y[i];
    }
  }
}

template void ComponentProduct<Mem::Main, Algo::Generic>::value(float *, const float * const, const float * const, const Index);
template void ComponentProduct<Mem::Main, Algo::Generic>::value(double *, const double * const, const double * const, const Index);
