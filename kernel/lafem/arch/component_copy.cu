// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/component_copy.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

// includes, CUDA
#include <cublas_v2.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

void ComponentCopy::value_cuda(float * r, const float * const x, const int stride, const int block, const Index size)
{
  cublasCopyEx(Util::Intern::cublas_handle, int(size), x, CUDA_R_32F, 1, &r[block], CUDA_R_32F, stride);
  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

void ComponentCopy::value_to_cuda(const float * const r, float * x, const int stride, const int block, const Index size)
{
  cublasCopyEx(Util::Intern::cublas_handle, int(size), &r[block], CUDA_R_32F, stride, x, CUDA_R_32F, 1);
  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

void ComponentCopy::value_cuda(double * r, const double * const x, const int stride, const int block, const Index size)
{
  cublasCopyEx(Util::Intern::cublas_handle, int(size), x, CUDA_R_64F, 1, &r[block], CUDA_R_64F, stride);
  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

void ComponentCopy::value_to_cuda(const double * const r, double * x, const int stride, const int block, const Index size)
{
  cublasCopyEx(Util::Intern::cublas_handle, int(size), &r[block], CUDA_R_64F, stride, x, CUDA_R_64F, 1);
  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

#ifdef FEAT_HAVE_HALFMATH
void ComponentCopy::value_cuda(Half * r, const Half * const x, const int stride, const int block, const Index size)
{
  cublasCopyEx(Util::Intern::cublas_handle, int(size), x, CUDA_R_16F, 1, &r[block], CUDA_R_16F, stride);
  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
void ComponentCopy::value_to_cuda(const Half * const r, Half * x, const int stride, const int block, const Index size)
{
  cublasCopyEx(Util::Intern::cublas_handle, int(size), &r[block], CUDA_R_16F, stride, x, CUDA_R_16F, 1);
  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
#endif
