// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/component_product_dense.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  void ComponentProductDense::exec_mkl_impl(float * r, const float * const x, const float * const y, const Index size)
  {
    vsMul(static_cast<MKL_INT>(size), x, y, r);
  }

  void ComponentProductDense::exec_mkl_impl(double * r, const double * const x, const double * const y, const Index size)
  {
    vdMul(static_cast<MKL_INT>(size), x, y, r);
  }
} // namespace FEAT::LAFEM::Arch
