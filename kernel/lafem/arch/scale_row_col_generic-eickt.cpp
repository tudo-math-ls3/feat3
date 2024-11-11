// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/lafem/arch/scale_row_col.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template void ScaleRows::csr_generic(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index);
      template void ScaleRows::csr_generic(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index);
      template void ScaleRows::csr_generic(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index);
      template void ScaleRows::csr_generic(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index);

      template void ScaleCols::csr_generic(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index);
      template void ScaleCols::csr_generic(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index);
      template void ScaleCols::csr_generic(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index);
      template void ScaleCols::csr_generic(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index);
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
