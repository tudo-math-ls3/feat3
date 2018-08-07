// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_cscr.hpp>

namespace FEAT
{
  namespace LAFEM
  {
     template class SparseMatrixCSCR<Mem::Main, float, unsigned int>;
     template class SparseMatrixCSCR<Mem::Main, double, unsigned int>;
     template class SparseMatrixCSCR<Mem::Main, float, unsigned long>;
     template class SparseMatrixCSCR<Mem::Main, double, unsigned long>;
/*#ifdef FEAT_HAVE_CUDA
     template class SparseMatrixCSCR<Mem::CUDA, float, unsigned int>;
     template class SparseMatrixCSCR<Mem::CUDA, double, unsigned int>;
     template class SparseMatrixCSCR<Mem::CUDA, float, unsigned long>;
     template class SparseMatrixCSCR<Mem::CUDA, double, unsigned long>;
#endif*/
  }
}
