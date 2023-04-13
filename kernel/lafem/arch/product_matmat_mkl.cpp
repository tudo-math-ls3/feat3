// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/arch/product_matmat.hpp>

#include <mkl_blas.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void ProductMatMat::dense_mkl(float * r, const float alpha, const float beta, const float * const x, const float * const y, const float * const z, const Index rows, const Index columns, const Index inner)
{
  XASSERT(r == z);
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT minner = (MKL_INT)inner;
  char trans = 'N';
  sgemm(&trans, &trans, &mcolumns, &mrows, &minner, &alpha, y, &mcolumns, x, &minner, &beta, r, &mcolumns);
}

void ProductMatMat::dense_mkl(double * r, const double alpha, const double beta, const double * const x, const double * const y, const double * const z, const Index rows, const Index columns, const Index inner)
{
  XASSERT(r == z);
  MKL_INT mrows = (MKL_INT)rows;
  MKL_INT mcolumns = (MKL_INT)columns;
  MKL_INT minner = (MKL_INT)inner;
  char trans = 'N';
  dgemm(&trans, &trans, &mcolumns, &mrows, &minner, &alpha, y, &mcolumns, x, &minner, &beta, r, &mcolumns);
}
