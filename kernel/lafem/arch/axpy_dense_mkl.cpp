// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/axpy_dense.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
#include <mkl_spblas.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  void AxpyDense::exec_mkl_impl(float * r, const float alpha, const float * const x, const Index size)
  {
    cblas_saxpy(static_cast<MKL_INT>(size), alpha, x, 1, r, 1);
  }

  void AxpyDense::exec_mkl_impl(double * r, const double alpha, const double * const x, const Index size)
  {
    cblas_daxpy(static_cast<MKL_INT>(size), alpha, x, 1, r, 1);
  }
} // namespace FEAT::LAFEM::Arch
