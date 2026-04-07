// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/arch/matmatmult_dense_dense.hpp>

#include <mkl_blas.h>

namespace FEAT::LAFEM::Arch
{
  void MatMatMultDenseDense::exec_mkl_impl(float * r, const float alpha, const float * const x, const float * const y, const float * const z, const Index rows, const Index cols, const Index inner)
  {
    XASSERT((r == z) || (z == nullptr));
    MKL_INT mrows = (MKL_INT)rows;
    MKL_INT mcolumns = (MKL_INT)cols;
    MKL_INT minner = (MKL_INT)inner;
    char trans = 'N';
    float beta(z == nullptr ? 0.0f :1.0f);
    sgemm(&trans, &trans, &mcolumns, &mrows, &minner, &alpha, y, &mcolumns, x, &minner, &beta, r, &mcolumns);
  }

  void MatMatMultDenseDense::exec_mkl_impl(double * r, const double alpha, const double * const x, const double * const y, const double * const z, const Index rows, const Index cols, const Index inner)
  {
    XASSERT((r == z) || (z == nullptr));
    MKL_INT mrows = (MKL_INT)rows;
    MKL_INT mcolumns = (MKL_INT)cols;
    MKL_INT minner = (MKL_INT)inner;
    char trans = 'N';
    double beta(z == nullptr ? 0.0 :1.0);
    dgemm(&trans, &trans, &mcolumns, &mrows, &minner, &alpha, y, &mcolumns, x, &minner, &beta, r, &mcolumns);
  }
} // namespace FEAT::LAFEM::Arch
