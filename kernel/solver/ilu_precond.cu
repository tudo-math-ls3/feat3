// includes, FEAST
#include <kernel/base_header.hpp>

#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

#include "cusparse_v2.h"

// http://docs.nvidia.com/cuda/cusparse/#cusparse-lt-t-gt-csrilu02_solve


using namespace FEAST;

namespace FEAST
{
  namespace Util
  {
    namespace Intern
    {
      extern cusparseHandle_t cusparse_handle;
    }
  }

  namespace Solver
  {
    namespace Intern
    {

      struct CudaIluSolveInfo
      {
        cusparseMatDescr_t descr_M;
        cusparseMatDescr_t descr_L;
        cusparseMatDescr_t descr_U;
        csrilu02Info_t info_M;
        csrsv2Info_t info_L;
        csrsv2Info_t info_U;
        cusparseOperation_t trans_L;
        cusparseOperation_t trans_U;
        cusparseSolvePolicy_t policy_M;
        cusparseSolvePolicy_t policy_L;
        cusparseSolvePolicy_t policy_U;
        void * pBuffer;
        double * z;
        int m;
        int nnz;
      };

      void * cuda_ilu_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd)
      {
        double * z;
        cudaMalloc((void**)&z, m * sizeof(double));

        cusparseMatDescr_t descr_M = 0;
        cusparseMatDescr_t descr_L = 0;
        cusparseMatDescr_t descr_U = 0;
        csrilu02Info_t info_M  = 0;
        csrsv2Info_t  info_L  = 0;
        csrsv2Info_t  info_U  = 0;
        int pBufferSize_M;
        int pBufferSize_L;
        int pBufferSize_U;
        int pBufferSize;
        void *pBuffer = 0;
        int structural_zero;
        const cusparseSolvePolicy_t policy_M = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        const cusparseSolvePolicy_t policy_L = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        const cusparseSolvePolicy_t policy_U = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
        const cusparseOperation_t trans_U  = CUSPARSE_OPERATION_NON_TRANSPOSE;

        cusparseCreateMatDescr(&descr_M);
        cusparseSetMatIndexBase(descr_M, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);

        cusparseCreateMatDescr(&descr_L);
        cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
        cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_UNIT);

        cusparseCreateMatDescr(&descr_U);
        cusparseSetMatIndexBase(descr_U, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(descr_U, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatFillMode(descr_U, CUSPARSE_FILL_MODE_UPPER);
        cusparseSetMatDiagType(descr_U, CUSPARSE_DIAG_TYPE_NON_UNIT);

        cusparseCreateCsrilu02Info(&info_M);
        cusparseCreateCsrsv2Info(&info_L);
        cusparseCreateCsrsv2Info(&info_U);

        cusparseDcsrilu02_bufferSize(Util::Intern::cusparse_handle, m, nnz,
                descr_M, csrVal, csrRowPtr, csrColInd, info_M, &pBufferSize_M);
        cusparseDcsrsv2_bufferSize(Util::Intern::cusparse_handle, trans_L, m, nnz,
                descr_L, csrVal, csrRowPtr, csrColInd, info_L, &pBufferSize_L);
        cusparseDcsrsv2_bufferSize(Util::Intern::cusparse_handle, trans_U, m, nnz,
                descr_U, csrVal, csrRowPtr, csrColInd, info_U, &pBufferSize_U);

        pBufferSize = max(pBufferSize_M, max(pBufferSize_L, pBufferSize_U));

        cudaMalloc((void**)&pBuffer, pBufferSize);

        cusparseDcsrilu02_analysis(Util::Intern::cusparse_handle, m, nnz, descr_M,
                csrVal, csrRowPtr, csrColInd, info_M,
                    policy_M, pBuffer);
        cusparseStatus_t status = cusparseXcsrilu02_zeroPivot(Util::Intern::cusparse_handle, info_M, &structural_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "CUSPARSE ZERO PIVOT ERROR!");
        }

        cusparseDcsrsv2_analysis(Util::Intern::cusparse_handle, trans_L, m, nnz, descr_L,
                csrVal, csrRowPtr, csrColInd, info_L, policy_L, pBuffer);

        cusparseDcsrsv2_analysis(Util::Intern::cusparse_handle, trans_U, m, nnz, descr_U,
                csrVal, csrRowPtr, csrColInd, info_U, policy_U, pBuffer);


#ifdef FEAST_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        CudaIluSolveInfo * info = new CudaIluSolveInfo;
        info->descr_M = descr_M;
        info->descr_L = descr_L;
        info->descr_U = descr_U;
        info->info_M  = info_M;
        info->info_L  = info_L;
        info->info_U  = info_U;
        info->trans_L = trans_L;
        info->trans_U = trans_U;
        info->policy_M = policy_M;
        info->policy_L = policy_L;
        info->policy_U = policy_U;
        info->pBuffer = pBuffer;
        info->z = z;
        info->m = m;
        info->nnz = nnz;

        return (void*)info;
      }

      void cuda_ilu_init_numeric(double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo)
      {
        CudaIluSolveInfo * info = (CudaIluSolveInfo *) vinfo;

        cusparseDcsrilu02(Util::Intern::cusparse_handle, info->m, info->nnz, info->descr_M,
                csrVal, csrRowPtr, csrColInd, info->info_M, info->policy_M, info->pBuffer);
        int numerical_zero;
        cusparseStatus_t status = cusparseXcsrilu02_zeroPivot(Util::Intern::cusparse_handle, info->info_M, &numerical_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "CUSPARSE ZERO PIVOT ERROR!");
        }
      }

      int cuda_ilu_apply(double * y, const double * x, double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo)
      {
        CudaIluSolveInfo * info = (CudaIluSolveInfo *) vinfo;
        const double alpha = 1.;

        cusparseDcsrsv2_solve(Util::Intern::cusparse_handle, info->trans_L, info->m, info->nnz, &alpha, info->descr_L,
               csrVal, csrRowPtr, csrColInd, info->info_L,
                  x, info->z, info->policy_L, info->pBuffer);

        cusparseDcsrsv2_solve(Util::Intern::cusparse_handle, info->trans_U, info->m, info->nnz, &alpha, info->descr_U,
               csrVal, csrRowPtr, csrColInd, info->info_U,
                  info->z, y, info->policy_U, info->pBuffer);

#ifdef FEAST_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }

      void cuda_ilu_done(void * vinfo)
      {
        CudaIluSolveInfo * info = (CudaIluSolveInfo *) vinfo;

        cudaFree(info->z);
        cudaFree(info->pBuffer);
        cusparseDestroyMatDescr(info->descr_M);
        cusparseDestroyMatDescr(info->descr_L);
        cusparseDestroyMatDescr(info->descr_U);
        cusparseDestroyCsrilu02Info(info->info_M);
        cusparseDestroyCsrsv2Info(info->info_L);
        cusparseDestroyCsrsv2Info(info->info_U);

        delete info;
      }
    }
  }
}
