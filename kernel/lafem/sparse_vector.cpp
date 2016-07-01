// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_vector.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    template class SparseVector<Mem::Main, float, unsigned int>;
    template class SparseVector<Mem::Main, double, unsigned int>;
    template class SparseVector<Mem::Main, float, unsigned long>;
    template class SparseVector<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
    template class SparseVector<Mem::CUDA, float, unsigned int>;
    template class SparseVector<Mem::CUDA, double, unsigned int>;
    template class SparseVector<Mem::CUDA, float, unsigned long>;
    template class SparseVector<Mem::CUDA, double, unsigned long>;
#endif
  }
}
