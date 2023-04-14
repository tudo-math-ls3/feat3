// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
void Axpy::value_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
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

  cublasStatus_t status;
  void* temp_x;

  if (r == x)
  {
    cudaMalloc((void**)&temp_x, sizeof(DT_) * size);
    cudaMemcpy(temp_x, x, size * sizeof(DT_), cudaMemcpyDefault);
  }
  else
  {
    temp_x = (void*)x;
  }

  if (r != y)
  {
    ///\todo cuse cublasCopyEx when available
    cudaMemcpy(r, y, size * sizeof(DT_), cudaMemcpyDefault);
  }

  status = cublasAxpyEx(Util::Intern::cublas_handle, size, &a, et, temp_x, dt, 1, r, dt, 1, et);
  if (r == x)
    cudaFree (temp_x);

  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
#ifdef FEAT_HAVE_HALFMATH
template void Axpy::value_cuda(Half *, const Half, const Half * const, const Half * const, const Index);
#endif
template void Axpy::value_cuda(float *, const float, const float * const, const float * const, const Index);
template void Axpy::value_cuda(double *, const double, const double * const, const double * const, const Index);
