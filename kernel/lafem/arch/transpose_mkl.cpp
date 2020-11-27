// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/transpose.hpp>

#include <cstring>
#include <mkl.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Transpose<Mem::Main>::value_mkl(float * r, const float * const x, const Index rows_x, const Index columns_x)
{
  if (r != x)
  {
    mkl_somatcopy('R', 'T', rows_x, columns_x, float(1), x, columns_x, r, rows_x);
  }
  else
  {
    mkl_simatcopy('R', 'T', rows_x, columns_x, float(1), r, columns_x, rows_x);
  }
}

void Transpose<Mem::Main>::value_mkl(double * r, const double * const x, const Index rows_x, const Index columns_x)
{
  if (r != x)
  {
    mkl_domatcopy('R', 'T', rows_x, columns_x, double(1), x, columns_x, r, rows_x);
  }
  else
  {
    mkl_dimatcopy('R', 'T', rows_x, columns_x, double(1), r, columns_x, rows_x);
  }
}
