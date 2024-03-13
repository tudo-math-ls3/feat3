// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter.hpp>

#include <cstring>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
#include <mkl_spblas.h>
FEAT_RESTORE_WARNINGS

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void UnitFilter::filter_rhs_mkl(float * v, const float * const sv_elements, const std::uint64_t * const sv_indices, const Index ue)
{
  cblas_ssctr((MKL_INT)ue, sv_elements, (const MKL_INT*)sv_indices, v);
}

void UnitFilter::filter_rhs_mkl(double * v, const double * const sv_elements, const std::uint64_t * const sv_indices, const Index ue)
{
  cblas_dsctr((MKL_INT)ue, sv_elements, (const MKL_INT*)sv_indices, v);
}
