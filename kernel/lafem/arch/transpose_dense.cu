// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/transpose_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  void TransposeDense::exec_cuda_impl(float * r, const float * const x, Index rows_x, Index cols_x)
  {
    cublasStatus_t status;
    float one(1);
    float zero(0);
    float* temp = nullptr;

    if (r == x)
    {
      temp = (float*)Memory::alloc_cuda(rows_x * cols_x * sizeof(float));
      Memory::memcopy_cuda(temp, x, rows_x * cols_x * sizeof(float));
      status = cublasSgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, int(rows_x), int(cols_x), &one, temp, int(cols_x), &zero, nullptr, int(cols_x), r, int(rows_x));
    }
    else
    {
      status = cublasSgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, int(rows_x), int(cols_x), &one, x, int(cols_x), &zero, nullptr, int(cols_x), r, int(rows_x));
    }
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));


  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
    //free ptr
    if(temp)
      Memory::free_cuda(temp);
  }

  void TransposeDense::exec_cuda_impl(double * r, const double * const x, Index rows_x, Index cols_x)
  {
    cublasStatus_t status;
    double one(1);
    double zero(0);
    double *temp = nullptr;

    if (r == x)
    {
      temp = (double*)Memory::alloc_cuda(rows_x * cols_x * sizeof(double));
      Memory::memcopy_cuda(temp, x, rows_x * cols_x * sizeof(double));
      status = cublasDgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, int(rows_x), int(cols_x), &one, temp, int(cols_x), &zero, nullptr, int(cols_x), r, int(rows_x));
    }
    else
    {
      status = cublasDgeam(Util::Intern::cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, int(rows_x), int(cols_x), &one, x, int(cols_x), &zero, nullptr, int(cols_x), r, int(rows_x));
    }

    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
    if(temp)
      Memory::free_cuda(temp);
  }
} // namespace FEAT::LAFEM::Arch
