// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    template class DenseVector<Mem::Main, float, unsigned int>;
    template class DenseVector<Mem::Main, double, unsigned int>;
    template class DenseVector<Mem::Main, float, unsigned long>;
    template class DenseVector<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
    template class DenseVector<Mem::CUDA, float, unsigned int>;
    template class DenseVector<Mem::CUDA, double, unsigned int>;
    template class DenseVector<Mem::CUDA, float, unsigned long>;
    template class DenseVector<Mem::CUDA, double, unsigned long>;
#endif
  }
}
