// includes, FEAST
#include <kernel/lafem/dot_product.hpp>

// includes, CUDA
#include <cublas.h>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      float cuda_dot_product(const float * x, const float * y, const Index size)
      {
        return cublasSdot(size, x, 1, y, 1);
      }

      double cuda_dot_product(const double * x, const double * y, const Index size)
      {
        return cublasDdot(size, x, 1, y, 1);
      }
    }
  }
}

using namespace FEAST;
using namespace FEAST::LAFEM;

template <typename DT_>
DT_ DotProduct<Archs::GPU, Archs::CUDA>::value(const DenseVector<Archs::GPU, DT_> & x, const DenseVector<Archs::GPU, DT_> & y)
{
  if (x.size() != y.size())
    throw InternalError("Vector size does not match!");

  const DT_ * x_gpu(x.elements());
  const DT_ * y_gpu(y.elements());
  cublasInit();
  DT_ result = Intern::cuda_dot_product(x_gpu, y_gpu, x.size());
  cublasShutdown();
  return result;
}

template float DotProduct<Archs::GPU, Archs::CUDA>::value(const DenseVector<Archs::GPU, float> &, const DenseVector<Archs::GPU, float>&);
template double DotProduct<Archs::GPU, Archs::CUDA>::value(const DenseVector<Archs::GPU, double> &, const DenseVector<Archs::GPU, double>&);
