// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/diagonal_bcsr.hpp>

namespace FEAT::LAFEM::Arch
{
  template void DiagonalBCSR::exec_generic_impl(float *, const float * const, const Index * const, const Index * const, const Index, const int);
  template void DiagonalBCSR::exec_generic_impl(double *, const double * const, const Index * const, const Index * const, const Index, const int);
} // namespace FEAT::LAFEM::Arch
