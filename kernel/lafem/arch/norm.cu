// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/norm.hpp>
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
      float cuda_norm2(const float * x, const Index size)
      {
        float result;
        cublasStatus_t status;
        status = cublasSnrm2(Util::Intern::cublas_handle, size, x, 1, &result);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasnrm2 failed with status code: " + stringify(status));
        return result;
      }

      double cuda_norm2(const double * x, const Index size)
      {
        double result;
        cublasStatus_t status;
        status = cublasDnrm2(Util::Intern::cublas_handle, size, x, 1, &result);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasnrm2 failed with status code: " + stringify(status));
        return result;
      }
    }
  }
}

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
DT_ Norm2<Mem::CUDA>::value(const DT_ * const x, const Index size)
{
  DT_ result = Intern::cuda_norm2(x, size);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  return result;
}

template float Norm2<Mem::CUDA>::value(const float * const, const Index);
template double Norm2<Mem::CUDA>::value(const double * const, const Index);
