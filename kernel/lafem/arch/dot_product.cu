// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/dot_product.hpp>
#include <kernel/lafem/arch/component_product.hpp>
#include <kernel/util/exception.hpp>

// includes, CUDA
#include <cublas_v2.h>

namespace FEAST
{
  namespace Util
  {
    namespace Intern
    {
      extern cublasHandle_t cublas_handle;
    }
  }

  namespace LAFEM
  {
    namespace Intern
    {

      float cuda_dot_product(const float * x, const float * y, const Index size)
      {
        float result;
        if (CUBLAS_STATUS_SUCCESS != cublasSdot(Util::Intern::cublas_handle, size, x, 1, y, 1, &result))
          throw InternalError(__func__, __FILE__, __LINE__, "cublasSdot failed!");
        cudaDeviceSynchronize();
        return result;
      }

      double cuda_dot_product(const double * x, const double * y, const Index size)
      {
        double result;
        if (CUBLAS_STATUS_SUCCESS != cublasDdot(Util::Intern::cublas_handle, size, x, 1, y, 1, &result))
          throw InternalError(__func__, __FILE__, __LINE__, "cublasDdot failed!");
        cudaDeviceSynchronize();
        return result;
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
  DT_ result = Intern::cuda_dot_product(x, y, size);
  return result;
}

template float DotProduct<Mem::CUDA, Algo::CUDA>::value(const float * const, const float * const, const Index);
template double DotProduct<Mem::CUDA, Algo::CUDA>::value(const double * const, const double * const, const Index);

template <typename DT_>
DT_ TripleDotProduct<Mem::CUDA, Algo::CUDA>::value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size)
{
  DT_ * temp;
  cudaMalloc((void **) &temp, size * sizeof(DT_));
  ComponentProduct<Mem::CUDA, Algo::CUDA>::value(temp, y, z, size);
  DT_ result = Intern::cuda_dot_product(x, temp, size);
  cudaFree(temp);
  return result;
}

template float TripleDotProduct<Mem::CUDA, Algo::CUDA>::value(const float * const x, const float * const y, const float * const z, const Index size);
template double TripleDotProduct<Mem::CUDA, Algo::CUDA>::value(const double * const x, const double * const y, const double * const z, const Index size);
