// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/util/exception.hpp>

// includes, CUDA
#include <cublas_v2.h>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      extern cublasHandle_t cublas_handle;

      float cuda_norm2(const float * x, const Index size)
      {
        float result;
        cublasSnrm2(cublas_handle, size, x, 1, &result);
        return result;
      }

      double cuda_norm2(const double * x, const Index size)
      {
        double result;
        cublasDnrm2(cublas_handle, size, x, 1, &result);
        return result;
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
  DT_ result = Intern::cuda_norm2(x, size);
  return result;
}

template float Norm2<Mem::CUDA, Algo::CUDA>::value(const float * const, const Index);
template double Norm2<Mem::CUDA, Algo::CUDA>::value(const double * const, const Index);
