// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_dense_dense.hpp>
#include <kernel/util/exception.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
#include <mkl_spblas.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  void MatVecMultDenseDense::exec_mkl_impl(float * r, const float alpha, const float * const y, const float * const val, const float * const x, const Index rows, const Index cols, bool transposed)
  {
    MKL_INT mrows = static_cast<MKL_INT>(rows);
    MKL_INT mcolumns = static_cast<MKL_INT>(cols);

    if(y == nullptr)
      Memory::memset_main(r, 0, sizeof(float) * (transposed ? cols : rows));
    else if(r != y)
    {
      MKL_INT one = 1;
      MKL_INT n = transposed ? mcolumns : mrows;
      scopy(&n, y, &one, r, &one);
    }

    float beta(1);
    cblas_sgemv(CblasRowMajor, transposed ? CblasTrans : CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
  }

  void MatVecMultDenseDense::exec_mkl_impl(double * r, const double alpha, const double * const y, const double * const val, const double * const x, const Index rows, const Index cols, bool transposed)
  {
    MKL_INT mrows = static_cast<MKL_INT>(rows);
    MKL_INT mcolumns = static_cast<MKL_INT>(cols);

    if(y == nullptr)
      Memory::memset_main(r, 0, sizeof(double) * (transposed ? cols : rows));
    else if(r != y)
    {
      MKL_INT one = 1;
      MKL_INT n = transposed ? mcolumns : mrows;
      dcopy(&n, y, &one, r, &one);
    }

    double beta(1);
    cblas_dgemv(CblasRowMajor, transposed ? CblasTrans : CblasNoTrans, mrows, mcolumns, alpha, val, mcolumns, x, 1, beta, r, 1);
  }
} // namespace FEAT::LAFEM::Arch
