// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/axpy.hpp>

#include <cstring>

#include <mkl.h>
#include <mkl_spblas.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Axpy<Mem::Main>::dv_mkl(float * r, const float a, const float * const x, const float * const y, const Index size)
{
  if (r == y)
  {
    cblas_saxpy((MKL_INT)size, a, x, 1, r, 1);
  }
  else if (r == x)
  {
    float * t = new float[size];
    memcpy(t, x, sizeof(float) * size);
    memcpy(r, y, sizeof(float) * size);
    cblas_saxpy((MKL_INT)size, a, t, 1, r, 1);
    delete[] t;
  }
  else
  {
    memcpy(r, y, sizeof(float) * size);
    cblas_saxpy((MKL_INT)size, a, x, 1, r, 1);
  }
}

void Axpy<Mem::Main>::dv_mkl(double * r, const double a, const double * const x, const double * const y, const Index size)
{
  if (r == y)
  {
    cblas_daxpy((MKL_INT)size, a, x, 1, r, 1);
  }
  else if (r == x)
  {
    double * t = new double[size];
    memcpy(t, x, sizeof(double) * size);
    memcpy(r, y, sizeof(double) * size);
    cblas_daxpy((MKL_INT)size, a, t, 1, r, 1);
    delete[] t;
  }
  else
  {
    memcpy(r, y, sizeof(double) * size);
    cblas_daxpy((MKL_INT)size, a, x, 1, r, 1);
  }
}
