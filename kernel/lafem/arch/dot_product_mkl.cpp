// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/dot_product.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

float DotProduct::value_mkl(const float * const x, const float * const y, const Index size)
{
  return cblas_sdot((MKL_INT)size, x, 1, y, 1);
}

double DotProduct::value_mkl(const double * const x, const double * const y, const Index size)
{
  return cblas_ddot((MKL_INT)size, x, 1, y, 1);
}

float TripleDotProduct::value_mkl(const float * const x, const float * const y, const float * const z, const Index size)
{
  float * temp(new float[size]);
  vsMul((MKL_INT)size, y, z, temp);
  float result = cblas_sdot((MKL_INT)size, x, 1, temp, 1);
  delete[] temp;
  return result;
}

double TripleDotProduct::value_mkl(const double * const x, const double * const y, const double * const z, const Index size)
{
  double * temp(new double[size]);
  vdMul((MKL_INT)size, y, z, temp);
  double result = cblas_ddot((MKL_INT)size, x, 1, temp, 1);
  delete[] temp;
  return result;
}
