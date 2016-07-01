// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>

namespace FEAT
{
  namespace LAFEM
  {
     template class SparseMatrixELL<Mem::Main, float, unsigned int>;
     template class SparseMatrixELL<Mem::Main, double, unsigned int>;
     template class SparseMatrixELL<Mem::Main, float, unsigned long>;
     template class SparseMatrixELL<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
     template class SparseMatrixELL<Mem::CUDA, float, unsigned int>;
     template class SparseMatrixELL<Mem::CUDA, double, unsigned int>;
     template class SparseMatrixELL<Mem::CUDA, float, unsigned long>;
     template class SparseMatrixELL<Mem::CUDA, double, unsigned long>;
#endif
  }
}
