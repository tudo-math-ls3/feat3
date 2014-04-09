// includes, FEAST
#include <kernel/lafem/arch/component_invert.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template void ComponentInvert<Mem::Main, Algo::Generic>::value(float *, const float * const, const float, const Index);
      template void ComponentInvert<Mem::Main, Algo::Generic>::value(double *, const double * const, const double, const Index);
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST
