// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/copy_stride.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  void CopyStride::exec_mkl_impl(float * r, const float * const x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index size)
  {
    cblas_scopy(static_cast<MKL_INT>(size), &x[comp_x], static_cast<MKL_INT>(stride_x), &r[comp_r], static_cast<MKL_INT>(stride_r));
  }

  void CopyStride::exec_mkl_impl(double * r, const double * const x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index size)
  {
    cblas_dcopy(static_cast<MKL_INT>(size), &x[comp_x], static_cast<MKL_INT>(stride_x), &r[comp_r], static_cast<MKL_INT>(stride_r));
  }
} // namespace FEAT::LAFEM::Arch
