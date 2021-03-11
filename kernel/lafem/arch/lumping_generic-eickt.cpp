// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/lumping.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void Lumping<Mem::Main>::csr_generic(float *, const float * const, const Index * const, const Index * const, const Index);
template void Lumping<Mem::Main>::csr_generic(double *, const double * const, const Index * const, const Index * const, const Index);

template void Lumping<Mem::Main>::bcsr_generic(float *, const float * const, const Index * const, const Index * const, const Index, const int, const int);
template void Lumping<Mem::Main>::bcsr_generic(double *, const double * const, const Index * const, const Index * const, const Index, const int, const int);

template void Lumping<Mem::Main>::ell_generic(float *, const float * const, const Index * const,
  const Index * const, const Index * const, const Index, const Index);
template void Lumping<Mem::Main>::ell_generic(double *, const double * const, const Index * const,
  const Index * const, const Index * const, const Index, const Index);
