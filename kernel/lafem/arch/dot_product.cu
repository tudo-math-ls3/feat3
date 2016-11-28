// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/dot_product.hpp>
#include <kernel/lafem/arch/component_product.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

// includes, CUDA
#include <cublas_v2.h>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {

      float cuda_dot_product(const float * x, const float * y, const Index size)
      {
        float result;
        cublasStatus_t status;
        status = cublasSdot(Util::Intern::cublas_handle, size, x, 1, y, 1, &result);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasdot failed with status code "+ stringify(status));
        return result;
      }

      double cuda_dot_product(const double * x, const double * y, const Index size)
      {
        double result;
        cublasStatus_t status;
        status = cublasDdot(Util::Intern::cublas_handle, size, x, 1, y, 1, &result);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasdot failed with status code "+ stringify(status));
        return result;
      }
    }
  }
}

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
DT_ DotProduct<Mem::CUDA>::value(const DT_ * const x, const DT_ * const y, const Index size)
{
  DT_ result = Intern::cuda_dot_product(x, y, size);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  return result;
}

template float DotProduct<Mem::CUDA>::value(const float * const, const float * const, const Index);
template double DotProduct<Mem::CUDA>::value(const double * const, const double * const, const Index);

template <typename DT_>
DT_ TripleDotProduct<Mem::CUDA>::value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size)
{
  DT_ * temp;
  cudaMalloc((void **) &temp, size * sizeof(DT_));
  ComponentProduct<Mem::CUDA>::value(temp, y, z, size);
  DT_ result = Intern::cuda_dot_product(x, temp, size);
  cudaFree(temp);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  return result;
}

template float TripleDotProduct<Mem::CUDA>::value(const float * const x, const float * const y, const float * const z, const Index size);
template double TripleDotProduct<Mem::CUDA>::value(const double * const x, const double * const y, const double * const z, const Index size);
