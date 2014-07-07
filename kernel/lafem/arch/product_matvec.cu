// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/util/exception.hpp>

#include "cusparse_v2.h"

namespace FEAST
{
  namespace Util
  {
    namespace Intern
    {
      extern cusparseHandle_t cusparse_handle;
    }
  }

  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_product_matvec_csr(DT_ * r, const DT_ * b, const DT_ * val, const unsigned long * col_ind,
          const unsigned long * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ sum(0);
        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          sum += val[i] * b[col_ind[i]];
        }
        r[idx] = sum;
      }

      void cusparse_product_matvec_csr(cusparseOperation_t trans,
          int m, int n, int nnz,
          const float * alpha, const cusparseMatDescr_t descrA,
          const float * csrVal, const int * csrRowPtr, const int *csrColInd,
          const float * x, const float * beta, float * y)
      {
        cusparseScsrmv(Util::Intern::cusparse_handle, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
            csrColInd, x, beta, y);
      }

      void cusparse_product_matvec_csr(cusparseOperation_t trans,
          int m, int n, int nnz,
          const double * alpha, const cusparseMatDescr_t descrA,
          const double * csrVal, const int * csrRowPtr, const int *csrColInd,
          const double * x, const double * beta, double * y)
      {
        cusparseDcsrmv(Util::Intern::cusparse_handle, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
            csrColInd, x, beta, y);
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_product_matvec_ell(DT_ * r, const DT_ * b, const DT_ * Ax, const IT_ * Aj,
          const IT_ * Arl, const Index stride, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index row(idx);
        const IT_ * tAj(Aj);
        const DT_ * tAx(Ax);
        DT_ sum(0);
        tAj += row;
        tAx += row;

        const Index max(Arl[row]);
        for(Index n(0); n < max ; n++)
        {
          const DT_ A_ij = *tAx;

          const IT_ col = *tAj;
          sum += A_ij * b[col];

          tAj += stride;
          tAx += stride;
        }
        r[row] = sum;
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_product_matvec_banded(DT_ * r, const DT_ * x, const DT_ * val, const IT_ * offsets, const Index num_of_offsets, const Index rows, const Index columns)
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
        r[idx] = sum;
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_product_matvec_q1(DT_ * r, const DT_ * x, const DT_ * val, const IT_ * offsets, const Index num_of_offsets, const Index rows, const Index columns)
      {
        const Index m(offsets[4] - offsets[1]);
        const Index n(rows);

        const DT_ * ll(val);
        const DT_ * ld(val + rows);
        const DT_ * lu(val + 2 * rows);
        const DT_ * dl(val + 3 * rows);
        const DT_ * dd(val + 4 * rows);
        const DT_ * du(val + 5 * rows);
        const DT_ * ul(val + 6 * rows);
        const DT_ * ud(val + 7 * rows);
        const DT_ * uu(val + 8 * rows);

        extern __shared__ char  vsmvf_cache[];
        DT_ * smvf_cache = (DT_*)vsmvf_cache;


        unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x;

        // runs from 0 to blockDim.x-1
        unsigned long lindex = threadIdx.x;

        DT_* Dcache = smvf_cache;
        DT_* Lcache = smvf_cache + blockDim.x + 2;
        DT_* Ucache = smvf_cache + 2 * (blockDim.x + 2);

        // prefetch chunks from iteration vector
        //
        //
        // data needed for DD, DU, DL: each thread loads one element, the first and last one load the border cases
        // x_0 ... x_blockdim-1 into c_1...c_blockdim
        if (idx < n) Dcache[lindex + 1] = x[idx];
        if (idx >= m && idx - m < n) Lcache[lindex + 1] = x[idx - m];
        if (idx + m < n) Ucache[lindex + 1] = x[idx + m];
        if (lindex == 0)
        {
          // x_-1 in c_0
          if (blockDim.x * blockIdx.x - 1 < n) Dcache[0] = x[blockDim.x * blockIdx.x - 1];
          if (blockDim.x * blockIdx.x - m - 1 < n) Lcache[0] = x[blockDim.x * blockIdx.x - m - 1];
          if (blockDim.x * blockIdx.x + m - 1 < n) Ucache[0] = x[blockDim.x * blockIdx.x + m - 1];
        }
        if (lindex == blockDim.x - 1)
        {
          // x_blockdim in c_blockdim+1
          if (blockDim.x * (blockIdx.x + 1) < n) Dcache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1)];
          if (blockDim.x * (blockIdx.x + 1) - m < n) Lcache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1) - m];
          if (blockDim.x * (blockIdx.x + 1) + m  < n) Ucache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1) + m];
        }
        __syncthreads();
        // now, compute
        if (idx < n)
        {
          DT_ ytemp1 = dd[idx] * Dcache[lindex + 1];
          if (idx > 0) ytemp1 += dl[idx] * Dcache[lindex];
          if (idx < n - 1) ytemp1 += du[idx] * Dcache[lindex + 2];

          if (idx > m) ytemp1 += ll[idx] * Lcache[lindex];
          if (idx > m - 1) ytemp1 += ld[idx] * Lcache[lindex + 1];
          if (idx > m - 2) ytemp1 += lu[idx] * Lcache[lindex + 2];

          if (idx < n - m + 1) ytemp1 += ul[idx] * Ucache[lindex];
          if (idx < n - m) ytemp1 += ud[idx] * Ucache[lindex + 1];
          if (idx < n - m - 1) ytemp1 += uu[idx] * Ucache[lindex + 2];
          r[idx] = ytemp1;
        }
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_product_matvec_csr<<<grid, block>>>(r, x, val, col_ind, row_ptr, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);

template <typename DT_>
void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
{
  cusparseMatDescr_t descr=0;
  cusparseCreateMatDescr(&descr);
  cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

  DT_ one(1);
  DT_ zero(0);
  FEAST::LAFEM::Intern::cusparse_product_matvec_csr(CUSPARSE_OPERATION_NON_TRANSPOSE, (int)rows, (int)columns, (int)used_elements, &one, descr, val, (int*)row_ptr, (int*)col_ind, x, &zero, r);

  cusparseDestroyMatDescr(descr);
}
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);


template <typename DT_, typename IT_>
void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_product_matvec_ell<<<grid, block>>>(r, x, Ax, Aj, Arl, stride, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

template <typename DT_, typename IT_>
void ProductMatVec<Mem::CUDA, Algo::CUDA>::banded(DT_ * r, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  //if (num_of_offsets == 9)
  //  FEAST::LAFEM::Intern::cuda_product_matvec_q1<<<grid, block, 3 * (block.x + 2) * sizeof(DT_)>>>(r, x, val, offsets, num_of_offsets, rows, columns);
  //else
    FEAST::LAFEM::Intern::cuda_product_matvec_banded<<<grid, block>>>(r, x, val, offsets, num_of_offsets, rows, columns);
}
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::banded(float *, const float * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::banded(double *, const double * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::banded(float *, const float * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::banded(double *, const double * const, const unsigned int * const, const double * const, const Index, const Index, const Index);
