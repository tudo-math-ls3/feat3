// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

// includes, CUDA
#include <cublas_v2.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
DT_ Norm2::value_cuda(const DT_ * const x, const Index size)
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
  DT_ result(42.);

  status = cublasNrm2Ex(Util::Intern::cublas_handle, size, x, dt, 1, &result, dt, et);
  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  return result;
}

#ifdef FEAT_HAVE_HALFMATH
template Half Norm2::value_cuda(const Half * const, const Index);
#endif
template float Norm2::value_cuda(const float * const, const Index);
template double Norm2::value_cuda(const double * const, const Index);
