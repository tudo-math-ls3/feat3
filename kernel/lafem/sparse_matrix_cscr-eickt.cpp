// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_cscr.hpp>

namespace FEAT
{
  namespace LAFEM
  {
     template class SparseMatrixCSCR<float, unsigned int>;
     template class SparseMatrixCSCR<double, unsigned int>;
     template class SparseMatrixCSCR<float, unsigned long>;
     template class SparseMatrixCSCR<double, unsigned long>;
  }
}
