// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/transpose_dense.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  void TransposeDense::exec_mkl_impl(float * r, const float * const x, const Index rows_x, const Index cols_x)
  {
    if (r != x)
    {
      mkl_somatcopy('R', 'T', rows_x, cols_x, float(1), x, cols_x, r, rows_x);
    }
    else
    {
      mkl_simatcopy('R', 'T', rows_x, cols_x, float(1), r, cols_x, rows_x);
    }
  }

  void TransposeDense::exec_mkl_impl(double * r, const double * const x, const Index rows_x, const Index cols_x)
  {
    if (r != x)
    {
      mkl_domatcopy('R', 'T', rows_x, cols_x, double(1), x, cols_x, r, rows_x);
    }
    else
    {
      mkl_dimatcopy('R', 'T', rows_x, cols_x, double(1), r, cols_x, rows_x);
    }
  }
} // namespace FEAT::LAFEM::Arch
