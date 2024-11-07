// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/component_copy.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void ComponentCopy::value_mkl(float * r, const float * const x, const int stride, const int block, const Index size)
{
  cblas_scopy(static_cast<MKL_INT>(size), x, 1, &r[block], static_cast<MKL_INT>(stride));
}

void ComponentCopy::value_mkl(double * r, const double * const x, const int stride, const int block, const Index size)
{
  cblas_dcopy(static_cast<MKL_INT>(size), x, 1, &r[block], static_cast<MKL_INT>(stride));
}

void ComponentCopy::value_to_mkl(const float * const r, float * x, const int stride, const int block, const Index size)
{
  cblas_scopy(static_cast<MKL_INT>(size), &r[block], static_cast<MKL_INT>(stride), x, 1);
}

void ComponentCopy::value_to_mkl(const double * const r, double * x, const int stride, const int block, const Index size)
{
  cblas_dcopy(static_cast<MKL_INT>(size), &r[block], static_cast<MKL_INT>(stride), x, 1);
}
