// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_csr_dense.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/half.hpp>

#include "cusparse_v2.h"

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  void MatVecMultCSRDense::exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
  {
    cusparseOperation_t trans;
    if (transposed)
      trans = CUSPARSE_OPERATION_TRANSPOSE;
    else
      trans = CUSPARSE_OPERATION_NON_TRANSPOSE;

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
      ct = CUDA_R_32F; //cusparseSpMV does not support computation in half, yet
    }
  #endif
    else
      XABORTM("unsupported data type!");

    cusparseIndexType_t it;
    if(sizeof(IT_) == 4u)
      it = CUSPARSE_INDEX_32I;
    else if(sizeof(IT_) == 8u)
      it = CUSPARSE_INDEX_64I;
    else
      XABORTM("unsupported index type!");

    cusparseStatus_t status;

    cusparseSpMatDescr_t descr=0;
    status = cusparseCreateCsr(&descr, rows, cols, nonzeros, (void*)row_ptr, (void*)col_idx, (void*)val, it, it, CUSPARSE_INDEX_BASE_ZERO, dt);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

    cusparseDnVecDescr_t dx;
    status = cusparseCreateDnVec(&dx, (transposed?rows:cols), (void*)x, dt);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

    cusparseDnVecDescr_t dr;
    status = cusparseCreateDnVec(&dr, (transposed?cols:rows), (void*)r, dt);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

    if(y == nullptr)
      Memory::memset_cuda(r, 0, sizeof(DT_) * (transposed ? cols : rows));
    else if(r != y)
      Memory::memcopy_cuda(r, y, sizeof(DT_) * (transposed ? cols : rows));

    size_t buffer_size(0);

    DT_ beta = DT_(1);

    status = cusparseSpMV_bufferSize(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_CSR_ALG1, &buffer_size);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMV_bufferSize failed with status code: " + stringify(cusparseGetErrorString(status)));

    void* buffer = Util::cuda_get_static_memory(buffer_size);

    status = cusparseSpMV(Util::Intern::cusparse_handle, trans, &alpha, descr, dx, &beta, dr, ct, CUSPARSE_SPMV_CSR_ALG1, buffer);
    if (status != CUSPARSE_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMV failed with status code: " + stringify(cusparseGetErrorString(status)));

    cusparseDestroySpMat(descr);
    cusparseDestroyDnVec(dx);
    cusparseDestroyDnVec(dr);
  }

  template void MatVecMultCSRDense::exec_cuda_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRDense::exec_cuda_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultCSRDense::exec_cuda_impl(Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  #endif

  template void MatVecMultCSRDense::exec_cuda_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRDense::exec_cuda_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultCSRDense::exec_cuda_impl(Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  #endif
} // namespace FEAT::LAFEM::Arch
