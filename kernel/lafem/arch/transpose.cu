// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/transpose.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void Transpose<Mem::CUDA>::value(float * r, const float * const x, Index rows_x, Index columns_x)
{
  cublasStatus_t status;
  float one(1);
  float zero(0);

  if (r == x)
  {
    float * temp;
    cudaMalloc((void**)&temp, rows_x * columns_x * sizeof(float));
    cudaMemcpy(temp, x, rows_x * columns_x * sizeof(float), cudaMemcpyDeviceToDevice);
    status = cublasSgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, rows_x, columns_x, &one, temp, columns_x, &zero, nullptr, columns_x, r, rows_x);
  }
  else
  {
    status = cublasSgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, rows_x, columns_x, &one, x, columns_x, &zero, nullptr, columns_x, r, rows_x);
  }
  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cublasSgeam failed with status code: " + stringify(status));

#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

void Transpose<Mem::CUDA>::value(double * r, const double * const x, Index rows_x, Index columns_x)
{
  cublasStatus_t status;
  double one(1);
  double zero(0);

  if (r == x)
  {
    double * temp;
    cudaMalloc((void**)&temp, rows_x * columns_x * sizeof(double));
    cudaMemcpy(temp, x, rows_x * columns_x * sizeof(double), cudaMemcpyDefault);
    status = cublasDgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, rows_x, columns_x, &one, temp, columns_x, &zero, nullptr, columns_x, r, rows_x);
  }
  else
  {
    status = cublasDgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, rows_x, columns_x, &one, x, columns_x, &zero, nullptr, columns_x, r, rows_x);
  }

  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cublasDgeam failed with status code: " + stringify(status));

#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
