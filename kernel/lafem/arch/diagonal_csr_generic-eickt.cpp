// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/diagonal_csr.hpp>

namespace FEAT::LAFEM::Arch
{
  template void DiagonalCSR::exec_generic_impl(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
  template void DiagonalCSR::exec_generic_impl(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index);

  template void DiagonalCSR::exec_generic_impl(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
  template void DiagonalCSR::exec_generic_impl(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
} // namespace FEAT::LAFEM::Arch
