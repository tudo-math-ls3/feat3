// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/scale.hpp>


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Scale<Mem::Main, Algo::Generic>::value(DT_ * r, const DT_ * const x, const DT_ s, const Index size)
{
  if (x == r)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] *= s;
    }
  }
  else
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r[i] = x[i] * s;
    }
  }
}

template void Scale<Mem::Main, Algo::Generic>::value(float *, const float * const, const float, const Index);
template void Scale<Mem::Main, Algo::Generic>::value(double *, const double * const, const double, const Index);
