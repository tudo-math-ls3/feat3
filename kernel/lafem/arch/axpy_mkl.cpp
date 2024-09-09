// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/axpy.hpp>

#include <cstring>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
#include <mkl_spblas.h>
FEAT_RESTORE_WARNINGS

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Axpy::value_mkl(float * r, const float a, const float * const x, const Index size)
{
  cblas_saxpy((MKL_INT)size, a, x, 1, r, 1);
}

void Axpy::value_mkl(double * r, const double a, const double * const x, const Index size)
{
  cblas_daxpy((MKL_INT)size, a, x, 1, r, 1);
}
