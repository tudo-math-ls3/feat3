// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_dense_dense.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/half.hpp>

#include "cusparse_v2.h"

namespace FEAT::LAFEM::Arch
{
  cublasStatus_t cublas_apply_dense(cublasOperation_t trans,
    int m, int n,
    const float * const alpha,
    const float * const val,
    const float * const x, const float * const beta, float * y)
  {
    return cublasSgemv(Util::Intern::cublas_handle, trans, n, m, alpha, val, n, x, 1, beta, y, 1);
  }

  cublasStatus_t cublas_apply_dense(cublasOperation_t trans,
    int m, int n,
    const double * const alpha,
    const double * const val,
    const double * const x, const double * const beta, double * y)
  {
    return cublasDgemv(Util::Intern::cublas_handle, trans, n, m, alpha, val, n, x, 1, beta, y, 1);
  }

  template <typename DT_>
  void MatVecMultDenseDense::exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index cols, bool transposed)
  {
    // cublas expects column major, so swap transpose and normal
    cublasOperation_t trans = (transposed ? CUBLAS_OP_N : CUBLAS_OP_T);

    if(y == nullptr)
      Memory::memset_cuda(r, 0, sizeof(DT_) * (transposed ? cols : rows));
    else if(r != y)
      Memory::memcopy_cuda(r, y, sizeof(DT_) * (transposed ? cols : rows));

    DT_ beta = (y == nullptr) ? DT_(0) : DT_(1);
    cublasStatus_t status = cublas_apply_dense(trans, (int)rows, (int)cols, &alpha, val, x, &beta, r);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void MatVecMultDenseDense::exec_cuda_impl(float *, const float, const float * const, const float * const, const float * const, const Index, const Index, bool);
  template void MatVecMultDenseDense::exec_cuda_impl(double *, const double, const double * const, const double * const, const double * const, const Index, const Index, bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultDenseDense::exec_cuda_impl(Half *, const Half, const Half * const, const Half * const, const Half * const, const Index, const Index, bool);
  #endif
} // namespace FEAT::LAFEM::Arch
