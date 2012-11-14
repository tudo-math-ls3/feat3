// includes, FEAST
#include <kernel/lafem/norm.hpp>

// includes, CUDA
#include <cublas.h>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      float cuda_norm2(const float * x, const Index size)
      {
        return cublasSnrm2(size, x, 1);
      }

      double cuda_norm2(const double * x, const Index size)
      {
        return cublasDnrm2(size, x, 1);
      }
    }
  }
}

using namespace FEAST;
using namespace FEAST::LAFEM;

template <typename DT_>
DT_ Norm2<Algo::CUDA>::value(const DenseVector<Mem::CUDA, DT_> & x)
{
  const DT_ * x_gpu(x.elements());
  cublasInit();
  DT_ result = Intern::cuda_norm2(x_gpu, x.size());
  cublasShutdown();
  return result;
}

template float Norm2<Algo::CUDA>::value(const DenseVector<Mem::CUDA, float> &);
template double Norm2<Algo::CUDA>::value(const DenseVector<Mem::CUDA, double> &);
