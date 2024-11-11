// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/component_product.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void ComponentProduct::value_mkl(float * r, const float * const x, const float * const y, const Index size)
{
  vsMul((MKL_INT)size, x, y, r);
}

void ComponentProduct::value_mkl(double * r, const double * const x, const double * const y, const Index size)
{
  vdMul((MKL_INT)size, x, y, r);
}
