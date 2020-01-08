// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
