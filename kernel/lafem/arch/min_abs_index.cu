// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/min_abs_index.hpp>
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

      Index cuda_min_abs_index(const float * x, const Index size)
      {
        int result;
        cublasStatus_t status;
        status = cublasIsamin(Util::Intern::cublas_handle, size, x, 1, &result);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasmin failed with status code "+ stringify(status));
        return (Index)result - 1;
      }

      Index cuda_min_abs_index(const double * x, const Index size)
      {
        int result;
        cublasStatus_t status;
        status = cublasIdamin(Util::Intern::cublas_handle, size, x, 1, &result);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cublasmin failed with status code "+ stringify(status));
        return (Index)result - 1;
      }
    }
  }
}

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
Index MinAbsIndex<Mem::CUDA>::value(const DT_ * const x, const Index size)
{
  Index result = Intern::cuda_min_abs_index(x, size);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  return result;
}

template Index MinAbsIndex<Mem::CUDA>::value(const float * const, const Index);
template Index MinAbsIndex<Mem::CUDA>::value(const double * const, const Index);
