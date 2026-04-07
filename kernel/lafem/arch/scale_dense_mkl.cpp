// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/scale_dense.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  void ScaleDense::exec_mkl_impl(float * r, const float * const x, const float alpha, const Index size)
  {
    if (r == x)
    {
      cblas_sscal(static_cast<MKL_INT>(size), alpha, r, 1);
    }
    else
    {
      Memory::memcopy_main(r, x, size * sizeof(float));
      cblas_sscal(static_cast<MKL_INT>(size), alpha, r, 1);
    }
  }

  void ScaleDense::exec_mkl_impl(double * r, const double * const x, const double alpha, const Index size)
  {
    if (r == x)
    {
      cblas_dscal(static_cast<MKL_INT>(size), alpha, r, 1);
    }
    else
    {
      Memory::memcopy_main(r, x, size * sizeof(double));
      cblas_dscal(static_cast<MKL_INT>(size), alpha, r, 1);
    }
  }
} // namespace FEAT::LAFEM::Arch
