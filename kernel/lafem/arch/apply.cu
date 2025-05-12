// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/apply.hpp>
#include <kernel/lafem/arch/component_product.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/half.hpp>

#include "cusparse_v2.h"

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_apply_banded(DT_ * r, const DT_ alpha, const DT_ * x, const DT_ beta, const DT_ * val, const IT_ * offsets, const Index num_of_offsets, const Index rows, const Index columns)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows)
          return;

        const Index k1(rows - 1);
        const Index k2(rows + columns - 1);

        Index start(0);

        while (k1 > offsets[start] + idx)
        {
          ++start;
        }

        Index end(start);

        while (end < num_of_offsets && idx + offsets[end] < k2)
        {
          ++end;
        }

        DT_ sum(DT_(0.0));
        for (Index diag(start); diag < end; ++diag)
        {
          sum += val[rows * diag + idx] * x[idx + offsets[diag] - rows + 1];
        }
        r[idx] = (sum*alpha) + beta * r[idx];
      }

      void cusparse_apply_bcsr(cusparseDirection_t dir, cusparseOperation_t trans,
                                       int m, int n, int nnz,
                                       const float * alpha, const cusparseMatDescr_t descrA,
                                       const float * csrVal, const int * csrRowPtr, const int *csrColInd,
                                       int block_dim,
                                       const float * x, const float * beta, float * y)
      {
        cusparseStatus_t status;
        status = cusparseSbsrmv(Util::Intern::cusparse_handle, dir, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
                       csrColInd, block_dim, x, beta, y);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrmv failed with status code: " + stringify(status));
      }

      void cusparse_apply_bcsr(cusparseDirection_t dir, cusparseOperation_t trans,
                                       int m, int n, int nnz,
                                       const double * alpha, const cusparseMatDescr_t descrA,
                                       const double * csrVal, const int * csrRowPtr, const int *csrColInd,
                                       int block_dim,
                                       const double * x, const double * beta, double * y)
      {
        cusparseStatus_t status;
        status = cusparseDbsrmv(Util::Intern::cusparse_handle, dir, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
                       csrColInd, block_dim, x, beta, y);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrmv failed with status code: " + stringify(status));
      }

      void cublas_apply_dense(cublasOperation_t trans,
                                       int m, int n,
                                       const float * alpha,
                                       const float * val,
                                       const float * x, const float * beta, float * y)
      {
        cublasStatus_t status;
        status = cublasSgemv(Util::Intern::cublas_handle, trans, n, m, alpha, val, n, x, 1, beta, y, 1);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
      }

      void cublas_apply_dense(cublasOperation_t trans,
                                       int m, int n,
                                       const double * alpha,
                                       const double * val,
                                       const double * x, const double * beta, double * y)
      {
        cublasStatus_t status;
        status = cublasDgemv(Util::Intern::cublas_handle, trans, n, m, alpha, val, n, x, 1, beta, y, 1);
        if (status != CUBLAS_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
      }

      template <int BlockSize_, typename DT_, typename IT_>
      __global__ void cuda_apply_csrsb(DT_ * r, const DT_ a, const DT_ * x, const DT_ b, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ bsum[BlockSize_];
        for (int j(0) ; j < BlockSize_ ; ++j)
        {
          bsum[j] = DT_(0);
        }
        const IT_ end(row_ptr[idx + 1]);
        for (IT_ i(row_ptr[idx]) ; i < end ; ++i)
        {
          const DT_ vali(val[i]);
          const IT_ coli(col_ind[i] * BlockSize_);
          for (int j(0) ; j < BlockSize_ ; ++j)
          {
            bsum[j] += vali * x[coli + j];
          }
        }
        for (int j(0) ; j < BlockSize_ ; ++j)
        {
          r[idx * BlockSize_ + j] = (bsum[j] * a) + b * r[idx * BlockSize_ + j];
        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_apply_bcsr(DT_ * r, const DT_ a, const DT_ * x, const DT_ b, const DT_ * val, const IT_ * col_ind,
          const IT_ * row_ptr, const Index count, const int BlockHeight, const int BlockWidth)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        /// \todo remove hardcoded number
        DT_ bsum[10];
        for (int j(0) ; j < BlockHeight ; ++j)
        {
          bsum[j] = DT_(0);
        }
        const IT_ end(row_ptr[idx + 1]);
        for (IT_ i(row_ptr[idx]) ; i < end ; ++i)
        {
          for (int h(0) ; h < BlockHeight ; ++h)
          {
            for (int w(0) ; w < BlockWidth ; ++w)
            {
              bsum[h] += val[i * BlockHeight * BlockWidth + h * BlockWidth + w] * x[col_ind[i] * BlockWidth + w];
            }
          }
        }
        for (int j(0) ; j < BlockHeight ; ++j)
        {
          r[idx * BlockHeight + j] = (bsum[j] * a) + b * r[idx * BlockHeight + j];
        }
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_, typename IT_>
void Apply::csr_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const bool transposed)
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
    throw InternalError(__func__, __FILE__, __LINE__, "unsupported data type!");

  cusparseIndexType_t it;
  if(sizeof(IT_) == 4u)
    it = CUSPARSE_INDEX_32I;
  else if(sizeof(IT_) == 8u)
    it = CUSPARSE_INDEX_64I;
  else
    throw InternalError(__func__, __FILE__, __LINE__, "unsupported index type!");

  cusparseStatus_t status;

  cusparseSpMatDescr_t descr=0;
  status = cusparseCreateCsr(&descr, rows, columns, used_elements, (void*)row_ptr, (void*)col_ind, (void*)val, it, it, CUSPARSE_INDEX_BASE_ZERO, dt);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

  cusparseDnVecDescr_t dx;
  status = cusparseCreateDnVec(&dx, (transposed?rows:columns), (void*)x, dt);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

  cusparseDnVecDescr_t dr;
  status = cusparseCreateDnVec(&dr, (transposed?columns:rows), (void*)r, dt);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));


  if (r != y)
  {
    MemoryPool::copy(r, y, (transposed?columns:rows));
  }

  size_t buffer_size(0);
  status = cusparseSpMV_bufferSize(Util::Intern::cusparse_handle, trans, &a, descr, dx, &b, dr, ct, CUSPARSE_SPMV_CSR_ALG1, &buffer_size);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMV_bufferSize failed with status code: " + stringify(cusparseGetErrorString(status)));

  // void* buffer = Util::cuda_malloc(buffer_size);
  void* buffer = Util::cuda_get_static_memory(buffer_size);

  status = cusparseSpMV(Util::Intern::cusparse_handle, trans, &a, descr, dx, &b, dr, ct, CUSPARSE_SPMV_CSR_ALG1, buffer);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMV failed with status code: " + stringify(cusparseGetErrorString(status)));

  cusparseDestroySpMat(descr);
  cusparseDestroyDnVec(dx);
  cusparseDestroyDnVec(dr);
  // Util::cuda_free(buffer);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Apply::csr_cuda(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
template void Apply::csr_cuda(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
#ifdef FEAT_HAVE_HALFMATH
template void Apply::csr_cuda(Half *, const Half, const Half * const, const Half, const Half * const, const Half * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
#endif

template void Apply::csr_cuda(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
template void Apply::csr_cuda(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
#ifdef FEAT_HAVE_HALFMATH
template void Apply::csr_cuda(Half *, const Half, const Half * const, const Half, const Half * const, const Half * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
#endif

//silence the compiler by pretending to accept any IT_ but hopefully only 'std::uint32_t' calls will be made
//this circumnavigates the missing static_if in bcsr_wrapper
template <typename DT_, typename IT_>
void Apply::bcsr_intern_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const int BlockSize)
{
  if (r != y)
  {
    Util::cuda_copy_device_to_device(r, y, rows * BlockSize * sizeof(DT_));
  }

  cusparseMatDescr_t descr=0;
  cusparseCreateMatDescr(&descr);
  cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

  FEAT::LAFEM::Intern::cusparse_apply_bcsr(CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE, (int)rows, (int)columns, (int)used_elements, &a, descr, val, (int*)row_ptr, (int*)col_ind,
      BlockSize, x, &b, r);

  cusparseDestroyMatDescr(descr);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template <typename DT_, typename IT_>
void Apply::bcsr_intern_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index /*columns*/, const Index /*used_elements*/, const int BlockHeight, const int BlockWidth)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  if (Math::abs(b) < Math::eps<DT_>())
  {
    MemoryPool::set_memory(r, DT_(0), /*(transposed?columns:rows)*/ rows * BlockHeight);
  }
  else if (r != y)
  {
    MemoryPool::copy(r, y, /*(transposed?columns:rows)*/ rows * BlockHeight);
  }

  FEAT::LAFEM::Intern::cuda_apply_bcsr<DT_, IT_><<<grid, block>>>(r, a, x, b, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template <typename DT_, typename IT_>
void Apply::bcsr_wrapper_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const int BlockHeight, const int BlockWidth)
{
  //CUSPARSE
  if ((BlockHeight == BlockWidth) && (sizeof(IT_) == 4u))
  {
    // CUSPARSE BCSR kernel only supports block sizes > 1; call scalar CSR in this case instead
    if(BlockHeight > 1)
      bcsr_intern_cuda<DT_, IT_>(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, BlockHeight);
    else
      csr_cuda<DT_, IT_>(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, false);
  }
  //GENERIC
  else
  {
    bcsr_intern_cuda<DT_, IT_>(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, BlockHeight, BlockWidth);
  }
}
template void Apply::bcsr_wrapper_cuda<float, std::uint32_t>(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const int, const int);
template void Apply::bcsr_wrapper_cuda<double, std::uint32_t>(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const int, const int);
template void Apply::bcsr_wrapper_cuda<float, std::uint64_t>(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const int, const int);
template void Apply::bcsr_wrapper_cuda<double, std::uint64_t>(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const int, const int);

template <int BlockSize_, typename DT_, typename IT_>
void Apply::csrsb_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows,
    const Index /*columns*/, const Index /*used_elements*/)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  if (Math::abs(b) < Math::eps<DT_>())
  {
    MemoryPool::set_memory(r, DT_(0), /*(transposed?columns:rows)*/ rows * BlockSize_);
  }
  else if (r != y)
  {
    MemoryPool::copy(r, y, /*(transposed?columns:rows)*/ rows * BlockSize_);
  }

  FEAT::LAFEM::Intern::cuda_apply_csrsb<BlockSize_, DT_, IT_><<<grid, block>>>(r, a, x, b, val, col_ind, row_ptr, rows);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Apply::csrsb_cuda<1, float, std::uint64_t>
  (float *, const float, const float *, const float, const float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<1, double, std::uint64_t>
  (double *, const double, const double *, const double, const double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<1, float, std::uint32_t>
  (float *, const float, const float *, const float, const float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<1, double, std::uint32_t>
  (double *, const double, const double *, const double, const double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<2, float, std::uint64_t>
  (float *, const float, const float *, const float, const float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<2, double, std::uint64_t>
  (double *, const double, const double *, const double, const double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<2, float, std::uint32_t>
  (float *, const float, const float *, const float, const float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<2, double, std::uint32_t>
  (double *, const double, const double *, const double, const double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<3, float, std::uint64_t>
  (float *, const float, const float *, const float, const float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<3, double, std::uint64_t>
  (double *, const double, const double *, const double, const double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<3, float, std::uint32_t>
  (float *, const float, const float *, const float, const float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
template void Apply::csrsb_cuda<3, double, std::uint32_t>
  (double *, const double, const double *, const double, const double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void Apply::banded_cuda(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y, const DT_ * const val, const IT_ * const offsets, const Index num_of_offsets, const Index rows, const Index columns)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  if (Math::abs(beta) < Math::eps<DT_>())
  {
    MemoryPool::set_memory(r, DT_(0), rows);
  }
  else if (r != y)
  {
    MemoryPool::copy(r, y, rows);
  }

  FEAT::LAFEM::Intern::cuda_apply_banded<<<grid, block>>>(r, alpha, x, beta, val, offsets, num_of_offsets, rows, columns);
}
template void Apply::banded_cuda(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint32_t * const, const Index, const Index, const Index);
template void Apply::banded_cuda(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const Index, const Index, const Index);
template void Apply::banded_cuda(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint64_t * const, const Index, const Index, const Index);
template void Apply::banded_cuda(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint64_t * const, const Index, const Index, const Index);

template <typename DT_>
void Apply::dense_cuda(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
{
  if (r != y)
  {
    Util::cuda_copy_device_to_device(r, y, rows * sizeof(DT_));
  }

  FEAT::LAFEM::Intern::cublas_apply_dense(CUBLAS_OP_T, (int)rows, (int)columns, &alpha, val, x, &beta, r);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Apply::dense_cuda(float * r, const float, const float, const float * const, const float * const, const float * const, const Index, const Index);
template void Apply::dense_cuda(double * r, const double, const double, const double * const, const double * const, const double * const, const Index, const Index);
