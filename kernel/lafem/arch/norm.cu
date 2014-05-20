// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/util/exception.hpp>

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
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
DT_ Norm2<Mem::CUDA, Algo::CUDA>::value(const DT_ * const x, const Index size)
{
  cublasInit();
  DT_ result = Intern::cuda_norm2(x, size);
  cublasShutdown();
  return result;
}

template float Norm2<Mem::CUDA, Algo::CUDA>::value(const float * const, const Index);
template double Norm2<Mem::CUDA, Algo::CUDA>::value(const double * const, const Index);
