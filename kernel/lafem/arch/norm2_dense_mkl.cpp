// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/norm2_dense.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{

float Norm2Dense::exec_mkl_impl(const float * const x, const Index size)
{
  return cblas_snrm2(static_cast<MKL_INT>(size), x, 1);
}

double Norm2Dense::exec_mkl_impl(const double * const x, const Index size)
{
  return cblas_dnrm2(static_cast<MKL_INT>(size), x, 1);
}
} // namespace FEAT::LAFEM::Arch
