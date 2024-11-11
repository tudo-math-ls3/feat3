// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/max_abs_index.hpp>
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

      Index cuda_max_abs_index(const float * x, const Index size)
      {
        int result;
        cublasStatus_t status;
        status = cublasIsamax(Util::Intern::cublas_handle, int(size), x, 1, &result);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
        return (Index)result - 1;
      }

      Index cuda_max_abs_index(const double * x, const Index size)
      {
        int result;
        cublasStatus_t status;
        status = cublasIdamax(Util::Intern::cublas_handle, int(size), x, 1, &result);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
        return (Index)result - 1;
      }
    }
  }
}

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
Index MaxAbsIndex::value_cuda(const DT_ * const x, const Index size)
{
  Index result = Intern::cuda_max_abs_index(x, size);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  return result;
}

template Index MaxAbsIndex::value_cuda(const float * const, const Index);
template Index MaxAbsIndex::value_cuda(const double * const, const Index);
