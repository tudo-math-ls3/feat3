// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/apply.hpp>
#include <kernel/lafem/arch/component_product.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/math.hpp>

#include "cusparse_v2.h"

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_apply_mv_csr(DT_ * r, const DT_ a, const DT_ * x, const DT_ b, const DT_ * val,
          const unsigned long * col_ind, const unsigned long * row_ptr, const Index count, const bool transposed)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ sum(0);
        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          sum += val[i] * x[col_ind[i]];
        }
        r[idx] = (sum * a) + b * r[idx];
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_apply_mv_ell(DT_ * r, const DT_ a, const DT_ * x, const DT_ b, const DT_ * val, const IT_ * col_ind,
                                       const IT_ * cs, const IT_ * cl, const Index rows, const Index C)
      {
        const Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows)
          return;


        DT_ sum(0);
        const Index chunk(idx / C);
        const Index local_row(idx % C);
        const Index chunk_end(cs[chunk+1]);

        for (Index pcol(cs[chunk] + local_row) ; pcol < chunk_end ; pcol+=C)
        {
          sum += val[pcol] * x[col_ind[pcol]];
        }
        r[idx] = sum * a + b * r[idx];

      }

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

      void cusparse_apply_csr(cusparseOperation_t trans,
                                       int m, int n, int nnz,
                                       const float * alpha, const cusparseMatDescr_t descrA,
                                       const float * csrVal, const int * csrRowPtr, const int *csrColInd,
                                       const float * x, const float * beta, float * y)
      {
        cusparseStatus_t status;
        status = cusparseScsrmv(Util::Intern::cusparse_handle, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
                       csrColInd, x, beta, y);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrmv failed with status code: " + stringify(status));
      }

      void cusparse_apply_csr(cusparseOperation_t trans,
                                       int m, int n, int nnz,
                                       const double * alpha, const cusparseMatDescr_t descrA,
                                       const double * csrVal, const int * csrRowPtr, const int *csrColInd,
                                       const double * x, const double * beta, double * y)
      {
        cusparseStatus_t status;
        status = cusparseDcsrmv(Util::Intern::cusparse_handle, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
                       csrColInd, x, beta, y);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrmv failed with status code: " + stringify(status));
      }

      void cusparse_apply_csrb(cusparseDirection_t dir, cusparseOperation_t trans,
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

      void cusparse_apply_csrb(cusparseDirection_t dir, cusparseOperation_t trans,
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
          throw InternalError(__func__, __FILE__, __LINE__, "cublasSgemv failed with status code: " + stringify(status));
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
          throw InternalError(__func__, __FILE__, __LINE__, "cublasDgemv failed with status code: " + stringify(status));
      }

      template <typename DT_, typename IT_, int BlockSize_>
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
      __global__ void cuda_apply_csrb(DT_ * r, const DT_ a, const DT_ * x, const DT_ b, const DT_ * val, const IT_ * col_ind,
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
        const unsigned long end(row_ptr[idx + 1]);
        for (unsigned long i(row_ptr[idx]) ; i < end ; ++i)
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

template <typename DT_>
void Apply<Mem::CUDA>::csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows, const Index columns, const Index used_elements, const bool transposed)
{
  if (transposed)
    throw InternalError(__func__, __FILE__, __LINE__, "transposedd csr product not supported for IT=unsigned long!");

  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  if (Math::abs(b) < Math::eps<DT_>())
  {
    MemoryPool<Mem::CUDA>::set_memory(r, DT_(0), (transposed?columns:rows));
  }
  else if (r != y)
  {
    MemoryPool<Mem::CUDA>::copy(r, y, (transposed?columns:rows));
  }

  FEAT::LAFEM::Intern::cuda_apply_mv_csr<<<grid, block>>>(r, a, x, b, val, col_ind, row_ptr, rows, transposed);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Apply<Mem::CUDA>::csr(float *, const float, const float * const, const float, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index, const bool);
template void Apply<Mem::CUDA>::csr(double *, const double, const double * const, const double, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index, const bool);

template <typename DT_>
void Apply<Mem::CUDA>::csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows, const Index columns, const Index used_elements, const bool transposed)
{
  cusparseOperation_t trans;
  if (transposed)
    trans = CUSPARSE_OPERATION_TRANSPOSE;
  else
    trans = CUSPARSE_OPERATION_NON_TRANSPOSE;

  cusparseMatDescr_t descr=0;
  cusparseCreateMatDescr(&descr);
  cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

  if (r != y)
  {
    MemoryPool<Mem::CUDA>::copy(r, y, (transposed?columns:rows));
  }

  FEAT::LAFEM::Intern::cusparse_apply_csr(trans, (int)rows, (int)columns, (int)used_elements, &a, descr, val, (int*)row_ptr, (int*)col_ind, x, &b, r);

  cusparseDestroyMatDescr(descr);

#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Apply<Mem::CUDA>::csr(float *, const float, const float * const, const float, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index, const bool);
template void Apply<Mem::CUDA>::csr(double *, const double, const double * const, const double, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index, const bool);

//silence the compiler by pretending to accept any IT_ but hopefully only 'unsigned int' calls will be made
//this circumnavigates the missing static_if in csrb_wrapper
template <typename DT_, typename IT_>
void Apply<Mem::CUDA>::csrb_intern(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const int BlockSize)
{
  if (r != y)
  {
    cudaMemcpy(r, y, rows * BlockSize * sizeof(DT_), cudaMemcpyDeviceToDevice);
  }

  cusparseMatDescr_t descr=0;
  cusparseCreateMatDescr(&descr);
  cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

  FEAT::LAFEM::Intern::cusparse_apply_csrb(CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE, (int)rows, (int)columns, (int)used_elements, &a, descr, val, (int*)row_ptr, (int*)col_ind,
      BlockSize, x, &b, r);

  cusparseDestroyMatDescr(descr);

#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template <typename DT_, typename IT_>
void Apply<Mem::CUDA>::csrb_intern(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const int BlockHeight, const int BlockWidth)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  if (Math::abs(b) < Math::eps<DT_>())
  {
    MemoryPool<Mem::CUDA>::set_memory(r, DT_(0), /*(transposed?columns:rows)*/ rows * BlockHeight);
  }
  else if (r != y)
  {
    MemoryPool<Mem::CUDA>::copy(r, y, /*(transposed?columns:rows)*/ rows * BlockHeight);
  }

  FEAT::LAFEM::Intern::cuda_apply_csrb<DT_, IT_><<<grid, block>>>(r, a, x, b, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template <typename DT_, typename IT_>
void Apply<Mem::CUDA>:: csrb_wrapper(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const Index columns, const Index used_elements, const int BlockHeight, const int BlockWidth)
{
  //CUSPARSE
  if (BlockHeight == BlockWidth && typeid(unsigned int) == typeid(IT_))
  {
    csrb_intern<DT_, IT_>(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, BlockHeight);
  }
  //GENERIC
  else
  {
    csrb_intern<DT_, IT_>(r, a, x, b, y, val, col_ind, row_ptr, rows, columns, used_elements, BlockHeight, BlockWidth);
  }
}
template void Apply<Mem::CUDA>::csrb_wrapper<float, unsigned int>(float *, const float, const float * const, const float, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index, const int, const int);
template void Apply<Mem::CUDA>::csrb_wrapper<double, unsigned int>(double *, const double, const double * const, const double, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index, const int, const int);
template void Apply<Mem::CUDA>::csrb_wrapper<float, unsigned long>(float *, const float, const float * const, const float, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index, const int, const int);
template void Apply<Mem::CUDA>::csrb_wrapper<double, unsigned long>(double *, const double, const double * const, const double, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index, const int, const int);

template <typename DT_, typename IT_, int BlockSize_>
void Apply<Mem::CUDA>::csrsb(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows,
    const Index columns, const Index used_elements)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  if (Math::abs(b) < Math::eps<DT_>())
  {
    MemoryPool<Mem::CUDA>::set_memory(r, DT_(0), /*(transposed?columns:rows)*/ rows * BlockSize_);
  }
  else if (r != y)
  {
    MemoryPool<Mem::CUDA>::copy(r, y, /*(transposed?columns:rows)*/ rows * BlockSize_);
  }

  FEAT::LAFEM::Intern::cuda_apply_csrsb<DT_, IT_, BlockSize_><<<grid, block>>>(r, a, x, b, val, col_ind, row_ptr, rows);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Apply<Mem::CUDA>::csrsb<float, unsigned long, 1>
  (float *, const float, const float *, const float, const float *, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<double, unsigned long, 1>
  (double *, const double, const double *, const double, const double *, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<float, unsigned int, 1>
  (float *, const float, const float *, const float, const float *, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<double, unsigned int, 1>
  (double *, const double, const double *, const double, const double *, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<float, unsigned long, 2>
  (float *, const float, const float *, const float, const float *, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<double, unsigned long, 2>
  (double *, const double, const double *, const double, const double *, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<float, unsigned int, 2>
  (float *, const float, const float *, const float, const float *, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<double, unsigned int, 2>
  (double *, const double, const double *, const double, const double *, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<float, unsigned long, 3>
  (float *, const float, const float *, const float, const float *, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<double, unsigned long, 3>
  (double *, const double, const double *, const double, const double *, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<float, unsigned int, 3>
  (float *, const float, const float *, const float, const float *, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::csrsb<double, unsigned int, 3>
  (double *, const double, const double *, const double, const double *, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void Apply<Mem::CUDA>::ell(DT_ * r, const DT_ a, const DT_ * const x, const DT_ b, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs, const IT_ * const cl, const Index C, const Index rows)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  if (Math::abs(b) < Math::eps<DT_>())
  {
    MemoryPool<Mem::CUDA>::set_memory(r, DT_(0), rows);
  }
  else if (r != y)
  {
    MemoryPool<Mem::CUDA>::copy(r, y, rows);
  }

  FEAT::LAFEM::Intern::cuda_apply_mv_ell<<<grid, block>>>(r, a, x, b, val, col_ind, cs, cl, rows, C);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Apply<Mem::CUDA>::ell(float *, const float, const float * const, const float, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Apply<Mem::CUDA>::ell(double *, const double, const double * const, const double, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Apply<Mem::CUDA>::ell(float *, const float, const float * const, const float, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
template void Apply<Mem::CUDA>::ell(double *, const double, const double * const, const double, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const Index, const Index);

template <typename DT_, typename IT_>
void Apply<Mem::CUDA>::banded(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ beta, const DT_ * const y, const DT_ * const val, const IT_ * const offsets, const Index num_of_offsets, const Index rows, const Index columns)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  if (Math::abs(beta) < Math::eps<DT_>())
  {
    MemoryPool<Mem::CUDA>::set_memory(r, DT_(0), rows);
  }
  else if (r != y)
  {
    MemoryPool<Mem::CUDA>::copy(r, y, rows);
  }

  FEAT::LAFEM::Intern::cuda_apply_banded<<<grid, block>>>(r, alpha, x, beta, val, offsets, num_of_offsets, rows, columns);
}
template void Apply<Mem::CUDA>::banded(float *, const float, const float * const, const float, const float * const, const float * const, const unsigned int * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::banded(double *, const double, const double * const, const double, const double * const, const double * const, const unsigned int * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::banded(float *, const float, const float * const, const float, const float * const, const float * const, const unsigned long * const, const Index, const Index, const Index);
template void Apply<Mem::CUDA>::banded(double *, const double, const double * const, const double, const double * const, const double * const, const unsigned long * const, const Index, const Index, const Index);

template <typename DT_>
void Apply<Mem::CUDA>::dense(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index columns)
{
  if (r != y)
  {
    cudaMemcpy(r, y, rows * sizeof(DT_), cudaMemcpyDeviceToDevice);
  }

  FEAT::LAFEM::Intern::cublas_apply_dense(CUBLAS_OP_T, (int)rows, (int)columns, &alpha, val, x, &beta, r);

#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Apply<Mem::CUDA>::dense(float * r, const float, const float, const float * const, const float * const, const float * const, const Index, const Index);
template void Apply<Mem::CUDA>::dense(double * r, const double, const double, const double * const, const double * const, const double * const, const Index, const Index);
