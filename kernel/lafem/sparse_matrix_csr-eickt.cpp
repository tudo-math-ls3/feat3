// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace FEAT
{
  namespace LAFEM
  {
     template class SparseMatrixCSR<Mem::Main, float, unsigned int>;
     template class SparseMatrixCSR<Mem::Main, double, unsigned int>;
     template class SparseMatrixCSR<Mem::Main, float, unsigned long>;
     template class SparseMatrixCSR<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
     template class SparseMatrixCSR<Mem::CUDA, float, unsigned int>;
     template class SparseMatrixCSR<Mem::CUDA, double, unsigned int>;
     template class SparseMatrixCSR<Mem::CUDA, float, unsigned long>;
     template class SparseMatrixCSR<Mem::CUDA, double, unsigned long>;
#endif
  }
}
