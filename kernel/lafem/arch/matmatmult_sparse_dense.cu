// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matmatmult_sparse_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

#include <cublas_v2.h>
#include <cublasLt.h>
#include <cusparse_v2.h>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void MatMatMultSparseDense::exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index nonzeros,
    const DT_ * y, const DT_ * z, const Index rows, const Index cols, const Index inner)
  {
    XASSERT((r == z) || (z == nullptr));
    cudaDataType dt;
    cudaDataType ct; //compute type
    if (typeid(DT_) == typeid(double))
    {
      dt = CUDA_R_64F;
      ct = CUDA_R_64F;
    }
    else if (typeid(DT_) == typeid(float))
    {
      dt = CUDA_R_32F;
      ct = CUDA_R_32F;
    }
  #ifdef FEAT_HAVE_HALFMATH
    else if (typeid(DT_) == typeid(Half))
    {
      dt = CUDA_R_16F;
      ct = CUDA_R_32F; //cusparseSpMM does not support computation in half, yet
    }
  #endif
    else
    {
      throw InternalError(__func__, __FILE__, __LINE__, "unsupported data type!");
    }

    cusparseIndexType_t it;
    if(sizeof(IT_) == 4u)
      it = CUSPARSE_INDEX_32I;
    else if(sizeof(IT_) == 8u)
      it = CUSPARSE_INDEX_64I;
    else
    {
      throw InternalError(__func__, __FILE__, __LINE__, "unsupported index type!");
    }

    cusparseStatus_t status;

    cusparseDnMatDescr_t descr_r=0;
    status = cusparseCreateDnMat(&descr_r, rows, cols, cols, (void*)r, dt, CUSPARSE_ORDER_ROW);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

    cusparseSpMatDescr_t descr_x=0;
    status = cusparseCreateCsr(&descr_x, rows, inner, nonzeros, (void*)row_ptr, (void*)col_idx, (void*)val, it, it, CUSPARSE_INDEX_BASE_ZERO, dt);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

    cusparseDnMatDescr_t descr_y=0;
    status = cusparseCreateDnMat(&descr_y, inner, cols, cols, (void*)y, dt, CUSPARSE_ORDER_ROW);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

    cusparseOperation_t trans = CUSPARSE_OPERATION_NON_TRANSPOSE;
    size_t buffer_size(0);
    DT_ beta = (z == nullptr ? DT_(0) : DT_(1));
    status = cusparseSpMM_bufferSize(Util::Intern::cusparse_handle, trans, trans, &alpha, descr_x, descr_y, &beta, descr_r, ct, CUSPARSE_SPMM_CSR_ALG2, &buffer_size);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrmvex_buffersize failed with status code: " + stringify(cusparseGetErrorString(status)));

    void* buffer = Util::cuda_get_static_memory(buffer_size);

    status = cusparseSpMM(Util::Intern::cusparse_handle, trans, trans, &alpha, descr_x, descr_y, &beta, descr_r, ct, CUSPARSE_SPMM_CSR_ALG2, buffer);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMM failed with status code: " + stringify(cusparseGetErrorString(status)));

    cusparseDestroyDnMat(descr_r);
    cusparseDestroySpMat(descr_x);
    cusparseDestroyDnMat(descr_y);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  #ifdef FEAT_HAVE_HALFMATH
  template void MatMatMultSparseDense::exec_cuda_impl(Half *, const Half, const Half * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Half *,  const Half *, const Index, const Index, const Index);
  #endif
  template void MatMatMultSparseDense::exec_cuda_impl(float *, const float, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const float *, const float *, const Index, const Index, const Index);
  template void MatMatMultSparseDense::exec_cuda_impl(double *, const double, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const double *, const double *, const Index, const Index, const Index);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatMatMultSparseDense::exec_cuda_impl(Half *, const Half, const Half * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Half *,  const Half *, const Index, const Index, const Index);
  #endif
  template void MatMatMultSparseDense::exec_cuda_impl(float *, const float, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const float *, const float *, const Index, const Index, const Index);
  template void MatMatMultSparseDense::exec_cuda_impl(double *, const double, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const double *, const double *, const Index, const Index, const Index);
} // namespace FEAT::LAFEM::Arch
