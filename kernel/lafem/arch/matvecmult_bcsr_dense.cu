// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_bcsr_dense.hpp>
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
  __global__ void cuda_apply_bcsr(DT_ * r, const DT_ alpha, const DT_ * x, const DT_ * val, const IT_ * col_idx,
    const IT_ * row_ptr, const Index count, const int block_height, const int block_width)
  {
    // grid strided for loop
    for(Index idx = threadIdx.x + blockDim.x * blockIdx.x; idx < count; idx += blockDim.x * gridDim.x)
    {
      /// \todo remove hardcoded number
      DT_ bsum[10];
      for (int j(0) ; j < block_height ; ++j)
      {
        bsum[j] = DT_(0);
      }
      const IT_ end(row_ptr[idx + 1]);
      for (IT_ i(row_ptr[idx]) ; i < end ; ++i)
      {
        for (int h(0) ; h < block_height ; ++h)
        {
          for (int w(0) ; w < block_width ; ++w)
          {
            bsum[h] += val[i * block_height * block_width + h * block_width + w] * x[col_idx[i] * block_width + w];
          }
        }
      }
      for (int j(0) ; j < block_height ; ++j)
      {
        r[idx * block_height + j] += alpha * bsum[j];
      }
    }
  }

  template <typename DT_, typename IT_>
  void bcsr_intern_cuda(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
    const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index /*nonzeros*/,
    const int block_height, const int block_width, const bool transposed)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    Index rsize = (transposed ? cols*Index(block_width) : rows*Index(block_height));
    if(y == nullptr)
      Memory::memset_cuda(r, 0, rsize * sizeof(DT_));
    else if (r != y)
      Memory::memcopy_cuda(r, y, rsize * sizeof(DT_));

    XASSERTM(!transposed, "transposed BCSR matvec for CUDA not implemented");
    XASSERT(block_height < 10); // see cuda_apply_bcsr implementation

    cuda_apply_bcsr<DT_, IT_><<<grid, block>>>(r, alpha, x, val, col_idx, row_ptr, rows, block_height, block_width);

#ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
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

  //silence the compiler by pretending to accept any IT_ but hopefully only 'std::uint32_t' calls will be made
  //this circumnavigates the missing static_if in bcsr_wrapper
  template <typename DT_, typename IT_>
  void bcsr_intern_cusparse(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const int block_size, const bool transposed)
  {
    Index rsize = (transposed ? cols : rows) * Index(block_size);
    if(y == nullptr)
      Memory::memset_cuda(r, 0, rsize * sizeof(DT_));
    else if (r != y)
      Memory::memcopy_cuda(r, y, rsize * sizeof(DT_));

    cusparseMatDescr_t descr=0;
    cusparseCreateMatDescr(&descr);
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

    cusparseOperation_t trans = (transposed ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);

    DT_ b = DT_(y == nullptr ? 0 : 1);
    cusparse_apply_bcsr(CUSPARSE_DIRECTION_ROW, trans, (int)rows, (int)cols, (int)nonzeros, &alpha, descr, val, (int*)row_ptr, (int*)col_idx,
      block_size, x, &b, r);

    cusparseDestroyMatDescr(descr);

#ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  }


  template <typename DT_, typename IT_>
  void MatVecMultBCSRDense::exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const int block_height, const int block_width, const bool transposed)
  {
    if((block_height == 1) && (block_width == 1))
    {
      MatVecMultCSRDense::exec_cuda_impl<DT_, IT_>(r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, transposed);
    }
    else if((sizeof(IT_) == 4u) && (block_height == block_width))
    {
      bcsr_intern_cusparse<DT_, IT_>(r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, block_height, transposed);
    }
    else
    {
      bcsr_intern_cuda<DT_, IT_>(r, alpha, x, y, val, col_idx, row_ptr, rows, cols, nonzeros, block_height, block_width, transposed);
    }
  }

  template void MatVecMultBCSRDense::exec_cuda_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const int, const int, const bool);
  template void MatVecMultBCSRDense::exec_cuda_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const int, const int, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultBCSRDense::exec_cuda_impl(Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const int, const int, const bool);
  #endif

  template void MatVecMultBCSRDense::exec_cuda_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const int, const int, const bool);
  template void MatVecMultBCSRDense::exec_cuda_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const int, const int, const bool);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultBCSRDense::exec_cuda_impl(Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const int, const int, const bool);
  #endif
} // namespace FEAT::LAFEM::Arch
