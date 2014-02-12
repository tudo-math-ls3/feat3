// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
DT_ Norm2<Mem::Main, Algo::Generic>::value(const DT_ * const x, const Index size)
{
  DT_ r(0);
  for (Index i(0) ; i < size ; ++i)
  {
    r += x[i] * x[i];
  }

  return (DT_)Math::sqrt(r);
}

template float Norm2<Mem::Main, Algo::Generic>::value(const float * const, const Index);
template double Norm2<Mem::Main, Algo::Generic>::value(const double * const, const Index);
