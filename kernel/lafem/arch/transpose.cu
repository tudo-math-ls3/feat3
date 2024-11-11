// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/transpose.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Transpose::value_cuda(float * r, const float * const x, Index rows_x, Index columns_x)
{
  cublasStatus_t status;
  float one(1);
  float zero(0);
  float* temp = nullptr;

  if (r == x)
  {
    temp = (float*)Util::cuda_malloc(rows_x * columns_x * sizeof(float));
    Util::cuda_copy_device_to_device(temp, x, rows_x * columns_x * sizeof(float));
    status = cublasSgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, int(rows_x), int(columns_x), &one, temp, int(columns_x), &zero, nullptr, int(columns_x), r, int(rows_x));
  }
  else
  {
    status = cublasSgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, int(rows_x), int(columns_x), &one, x, int(columns_x), &zero, nullptr, int(columns_x), r, int(rows_x));
  }
  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));


  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  //free ptr
  Util::cuda_free(temp);
}

void Transpose::value_cuda(double * r, const double * const x, Index rows_x, Index columns_x)
{
  cublasStatus_t status;
  double one(1);
  double zero(0);
  double *temp = nullptr;

  if (r == x)
  {
    temp = (double*)Util::cuda_malloc(rows_x * columns_x * sizeof(double));
    Util::cuda_copy(temp, x, rows_x * columns_x * sizeof(double));
    status = cublasDgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, int(rows_x), int(columns_x), &one, temp, int(columns_x), &zero, nullptr, int(columns_x), r, int(rows_x));
  }
  else
  {
    status = cublasDgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, int(rows_x), int(columns_x), &one, x, int(columns_x), &zero, nullptr, int(columns_x), r, int(rows_x));
  }

  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  Util::cuda_free(temp);
}
