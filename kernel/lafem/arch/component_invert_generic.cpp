// includes, FEAT
#include <kernel/lafem/arch/component_invert.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template void ComponentInvert<Mem::Main>::value_generic(float *, const float * const, const float, const Index);
      template void ComponentInvert<Mem::Main>::value_generic(double *, const double * const, const double, const Index);
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
