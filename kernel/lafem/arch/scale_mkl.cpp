// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/scale.hpp>

#include <cstring>
#include <mkl.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Scale<Mem::Main>::value_mkl(float * r, const float * const x, const float s, const Index size)
{
  if (r == x)
  {
    cblas_sscal((MKL_INT)size, s, r, 1);
  }
  else
  {
    std::memcpy(r, x, size * sizeof(float));
    cblas_sscal((MKL_INT)size, s, r, 1);
  }
}

void Scale<Mem::Main>::value_mkl(double * r, const double * const x, const double s, const Index size)
{
  if (r == x)
  {
    cblas_dscal((MKL_INT)size, s, r, 1);
  }
  else
  {
    std::memcpy(r, x, size * sizeof(double));
    cblas_dscal((MKL_INT)size, s, r, 1);
  }
}
