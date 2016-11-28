// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/product_matmat.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

#include "cublas_v2.h"

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      void cublas_product_matmat_dense(int m, int n, int k,
                                       const float * x,
                                       const float * y, float * r)
      {
        cublasStatus_t status;
        const float one(1.f);
        const float zero(0.f);
        float * temp;
        cudaMalloc((void**)&temp, m*n*sizeof(float));

        status = cublasSgemm(Util::Intern::cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, n, m, k, &one, y, n, x, k, &zero, r, n);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasSgemm failed with status code: " + stringify(status));

        cudaFree(temp);
      }

      void cublas_product_matmat_dense(int m, int n, int k,
                                       const double * x,
                                       const double * y, double * r)
      {
        cublasStatus_t status;
        const double one(1.);
        const double zero(0.);
        double * temp;
        cudaMalloc((void**)&temp, m*n*sizeof(double));

        status = cublasDgemm(Util::Intern::cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, n, m, k, &one, y, n, x, k, &zero, r, n);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasDgemm failed with status code: " + stringify(status));

        cudaFree(temp);
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
void ProductMatMat<Mem::CUDA>::dense(DT_ * r, const DT_ * const x, const DT_ * const y, const Index rows, const Index columns, const Index inner)
{
  FEAT::LAFEM::Intern::cublas_product_matmat_dense((int)rows, (int)columns, (int)inner, x, y, r);

#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ProductMatMat<Mem::CUDA>::dense(float *, const float * const, const float * const, const Index, const Index, const Index);
template void ProductMatMat<Mem::CUDA>::dense(double *, const double * const, const double * const, const Index, const Index, const Index);
