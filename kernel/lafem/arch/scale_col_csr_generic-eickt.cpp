// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/lafem/arch/scale_col_csr.hpp>

namespace FEAT::LAFEM::Arch
{
  template void ScaleColCSR::exec_generic_impl(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index);
  template void ScaleColCSR::exec_generic_impl(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index);
  template void ScaleColCSR::exec_generic_impl(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index);
  template void ScaleColCSR::exec_generic_impl(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index);
} // namespace FEAT::LAFEM::Arch
