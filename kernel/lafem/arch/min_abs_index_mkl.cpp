// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/min_abs_index.hpp>

#include <mkl.h>


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

Index MinAbsIndex<Mem::Main>::value_mkl(const float * const x, const Index size)
{
  return cblas_isamin((MKL_INT)size, x, 1);
}

Index MinAbsIndex<Mem::Main>::value_mkl(const double * const x, const Index size)
{
  return cblas_idamin((MKL_INT)size, x, 1);
}
