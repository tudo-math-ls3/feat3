// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/scale_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void ScaleDense::exec_cuda_impl(DT_ * r, const DT_ * const x, const DT_ alpha, const Index size)
  {
    cudaDataType dt;
    cudaDataType et;
    if (typeid(DT_) == typeid(double))
    {
      dt = CUDA_R_64F;
      et = CUDA_R_64F;
    }
    else if (typeid(DT_) == typeid(float))
    {
      dt = CUDA_R_32F;
      et = CUDA_R_32F;
    }
  #ifdef FEAT_HAVE_HALFMATH
    else if (typeid(DT_) == typeid(Half))
    {
      dt = CUDA_R_16F;
      et = CUDA_R_32F;
    }
  #endif
    else
      throw InternalError(__func__, __FILE__, __LINE__, "unsupported data type!");

    if (r != x)
      Memory::memcopy_cuda(r, x, size * sizeof(DT_));

    cublasStatus_t status;

    status = cublasScalEx(Util::Intern::cublas_handle, int(size), &alpha, et, r, dt, 1, et);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }
  #ifdef FEAT_HAVE_HALFMATH
  template void ScaleDense::exec_cuda_impl(Half *, const Half * const, const Half, const Index);
  #endif
  template void ScaleDense::exec_cuda_impl(float *, const float * const, const float, const Index);
  template void ScaleDense::exec_cuda_impl(double *, const double * const, const double, const Index);
} // namespace FEAT::LAFEM::Arch
