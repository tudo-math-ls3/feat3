// includes, FEAST
#include <kernel/lafem/arch/product_matvec.hpp>
#include "cusparse_v2.h"

namespace FEAST
{
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

      void cusparse_product_matvec_csr(cusparseHandle_t handle, cusparseOperation_t trans,
          int m, int n, int nnz,
          const float * alpha, const cusparseMatDescr_t descrA,
          const float * csrVal, const int * csrRowPtr, const int *csrColInd,
          const float * x, const float * beta, float * y)
      {
        cusparseScsrmv(handle, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
            csrColInd, x, beta, y);
      }

      void cusparse_product_matvec_csr(cusparseHandle_t handle, cusparseOperation_t trans,
          int m, int n, int nnz,
          const double * alpha, const cusparseMatDescr_t descrA,
          const double * csrVal, const int * csrRowPtr, const int *csrColInd,
          const double * x, const double * beta, double * y)
      {
        cusparseDcsrmv(handle, trans, m, n, nnz, alpha, descrA, csrVal, csrRowPtr,
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
}
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);

template <typename DT_>
void ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
{
  cusparseHandle_t handle=0;
  cusparseMatDescr_t descr=0;
  cusparseCreate(&handle);
  cusparseCreateMatDescr(&descr);
  cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

  DT_ one(1);
  DT_ zero(0);
  FEAST::LAFEM::Intern::cusparse_product_matvec_csr(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, (int)rows, (int)columns, (int)used_elements, &one, descr, val, (int*)row_ptr, (int*)col_ind, x, &zero, r);

  cusparseDestroyMatDescr(descr);
  cusparseDestroy(handle);
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
}
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
template void ProductMatVec<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);
