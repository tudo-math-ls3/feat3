// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/diagonal.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void Diagonal<Mem::Main>::csr_generic(float *, const float * const, const Index * const, const Index * const, const Index);
template void Diagonal<Mem::Main>::csr_generic(double *, const double * const, const Index * const, const Index * const, const Index);

template void Diagonal<Mem::Main>::ell_generic(float *, const float * const, const Index * const,
  const Index * const, const Index * const, const Index, const Index);
template void Diagonal<Mem::Main>::ell_generic(double *, const double * const, const Index * const,
  const Index * const, const Index * const, const Index, const Index);
