// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/min_max_value_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

// includes, CUDA
#include <cublas_v2.h>

namespace FEAT::LAFEM::Arch
{
  float cuda_min_max_value(const float * x, const Index size, const bool min_, const bool abs_)
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

    float retval(0.0f);
    cudaMemcpy(&retval, &x[result - 1], sizeof(float), cudaMemcpyDeviceToHost);
    return std::abs(retval);
  }

  double cuda_min_max_value(const double * x, const Index size, const bool min_, const bool abs_)
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

    double retval(0.0);
    cudaMemcpy(&retval, &x[result - 1], sizeof(double), cudaMemcpyDeviceToHost);
    return std::abs(retval);
  }

  template <typename DT_>
  DT_ MinMaxValueDense::exec_cuda_impl(const DT_ * const x, const Index size, const bool min_, const bool abs_)
  {
    DT_ result = cuda_min_max_value(x, size, min_, abs_);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
    return result;
  }

  template float MinMaxValueDense::exec_cuda_impl(const float * const, const Index, const bool, const bool);
  template double MinMaxValueDense::exec_cuda_impl(const double * const, const Index, const bool, const bool);
} // namespace FEAT::LAFEM::Arch
