// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/min_max_index_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

// includes, CUDA
#include <cublas_v2.h>

namespace FEAT::LAFEM::Arch
{
  Index cuda_min_max_index(const float * x, const Index size, const bool min_, const bool abs_)
  {
    XASSERTM(abs_, "only absolute min/max implemented for cuda!");
    int result;
    cublasStatus_t status;
    if(min_)
      status = cublasIsamin(Util::Intern::cublas_handle, int(size), x, 1, &result);
    else
      status = cublasIsamax(Util::Intern::cublas_handle, int(size), x, 1, &result);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
    return (Index)result - 1;
  }

  Index cuda_min_max_index(const double * x, const Index size, const bool min_, const bool abs_)
  {
    XASSERTM(abs_, "only absolute min/max implemented for cuda!");
    int result;
    cublasStatus_t status;
    if(min_)
      status = cublasIdamin(Util::Intern::cublas_handle, int(size), x, 1, &result);
    else
      status = cublasIdamax(Util::Intern::cublas_handle, int(size), x, 1, &result);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
    return (Index)result - 1;
  }

  template <typename DT_>
  Index MinMaxIndexDense::exec_cuda_impl(const DT_ * const x, const Index size, const bool min_, const bool abs_)
  {
    Index result = cuda_min_max_index(x, size, min_, abs_);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
    return result;
  }

  template Index MinMaxIndexDense::exec_cuda_impl(const float * const, const Index, const bool, const bool);
  template Index MinMaxIndexDense::exec_cuda_impl(const double * const, const Index, const bool, const bool);
} // namespace FEAT::LAFEM::Arch
