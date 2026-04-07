// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/row_norm2_csr.hpp>

namespace FEAT::LAFEM::Arch
{
  template void RowNorm2CSR::exec_generic_impl(float*, const float* const, const Index* const, const Index * const, const Index, const bool);
  template void RowNorm2CSR::exec_generic_impl(double*, const double* const, const Index* const, const Index* const, const Index, const bool);
} // namespace FEAT::LAFEM::Arch
