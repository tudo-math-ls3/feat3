// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/min_max_index_dense.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  Index MinMaxIndexDense::exec_mkl_impl(const float * const x, const Index size, const bool min_, const bool abs_)
  {
    // MKL only supports absolute min/max
    if(!abs_)
      return MinMaxIndexDense::exec_generic_impl(x, size, min_, abs_);

    if(min_)
      return cblas_isamin(static_cast<MKL_INT>(size), x, 1);
    else
      return cblas_isamax(static_cast<MKL_INT>(size), x, 1);
  }

  Index MinMaxIndexDense::exec_mkl_impl(const double * const x, const Index size, const bool min_, const bool abs_)
  {
    // MKL only supports absolute min/max
    if(!abs_)
      return MinMaxIndexDense::exec_generic_impl(x, size, min_, abs_);

    if(min_)
      return cblas_idamin(static_cast<MKL_INT>(size), x, 1);
    else
      return cblas_idamax(static_cast<MKL_INT>(size), x, 1);
  }
} // namespace FEAT::LAFEM::Arch
