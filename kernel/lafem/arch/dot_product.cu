// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/dot_product.hpp>

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
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
DT_ DotProduct<Mem::CUDA, Algo::CUDA>::value(const DT_ * const x, const DT_ * const y, const Index size)
{
  cublasInit();
  DT_ result = Intern::cuda_dot_product(x, y, size);
  cublasShutdown();
  return result;
}

template float DotProduct<Mem::CUDA, Algo::CUDA>::value(const float * const, const float * const, const Index);
template double DotProduct<Mem::CUDA, Algo::CUDA>::value(const double * const, const double * const, const Index);
