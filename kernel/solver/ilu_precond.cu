// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>

#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

#include "cusparse_v2.h"

// http://docs.nvidia.com/cuda/cusparse/#cusparse-lt-t-gt-csrilu02_solve


using namespace FEAT;

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    namespace Intern
    {
      // CSR
      struct CudaIluSolveInfo
      {
        cusparseMatDescr_t descr_M;
        cusparseMatDescr_t descr_L;
        cusparseMatDescr_t descr_U;
        csrilu02Info_t info_M;
        cusparseSolveAnalysisInfo_t info_L;
        cusparseSolveAnalysisInfo_t info_U;
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
        cusparseSolveAnalysisInfo_t info_L  = 0;
        cusparseSolveAnalysisInfo_t info_U  = 0;
        int pBufferSize_M;
        int pBufferSize;
        void *pBuffer = 0;
        int structural_zero;
        const cusparseSolvePolicy_t policy_M = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        const cusparseSolvePolicy_t policy_L = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        const cusparseSolvePolicy_t policy_U = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
        const cusparseOperation_t trans_U  = CUSPARSE_OPERATION_NON_TRANSPOSE;

        cusparseStatus_t status;

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
        cusparseCreateSolveAnalysisInfo(&info_L);
        cusparseCreateSolveAnalysisInfo(&info_U);

        status = cusparseDcsrilu02_bufferSize(Util::Intern::cusparse_handle, m, nnz,
                descr_M, csrVal, csrRowPtr, csrColInd, info_M, &pBufferSize_M);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrilu02_bufferSize failed with status code: " + stringify(status));

        pBufferSize = pBufferSize_M;

        cudaMalloc((void**)&pBuffer, pBufferSize);

        status = cusparseDcsrilu02_analysis(Util::Intern::cusparse_handle, m, nnz, descr_M,
                csrVal, csrRowPtr, csrColInd, info_M,
                    policy_M, pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrilu02_analysis failed with status code: " + stringify(status));
        status = cusparseXcsrilu02_zeroPivot(Util::Intern::cusparse_handle, info_M, &structural_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "CUSPARSE ZERO PIVOT ERROR!");
        }

        status = cusparseDcsrsv_analysis(Util::Intern::cusparse_handle, trans_L, m, nnz, descr_L,
                csrVal, csrRowPtr, csrColInd, info_L);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparse_csrv_analysis failed with status code: " + stringify(status));

        status = cusparseDcsrsv_analysis(Util::Intern::cusparse_handle, trans_U, m, nnz, descr_U,
                csrVal, csrRowPtr, csrColInd, info_U);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparse_csrv_analysis failed with status code: " + stringify(status));


#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
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

        cusparseStatus_t status = cusparseDcsrilu02(Util::Intern::cusparse_handle, info->m, info->nnz, info->descr_M,
                csrVal, csrRowPtr, csrColInd, info->info_M, info->policy_M, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrilu02 failed with status code: " + stringify(status));
        int numerical_zero;
        status = cusparseXcsrilu02_zeroPivot(Util::Intern::cusparse_handle, info->info_M, &numerical_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "CUSPARSE ZERO PIVOT ERROR!");
        }
      }

      int cuda_ilu_apply(double * y, const double * x, double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo)
      {
        CudaIluSolveInfo * info = (CudaIluSolveInfo *) vinfo;
        const double alpha = 1.;

        cusparseStatus_t status = cusparseDcsrsv_solve(Util::Intern::cusparse_handle, info->trans_L, info->m, &alpha, info->descr_L,
               csrVal, csrRowPtr, csrColInd, info->info_L,
                  x, info->z);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsr2_solve failed with status code: " + stringify(status));

        status = cusparseDcsrsv_solve(Util::Intern::cusparse_handle, info->trans_U, info->m, &alpha, info->descr_U,
               csrVal, csrRowPtr, csrColInd, info->info_U,
                  info->z, y);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsr2_solve failed with status code: " + stringify(status));

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }

      void cuda_ilu_done_symbolic(void * vinfo)
      {
        CudaIluSolveInfo * info = (CudaIluSolveInfo *) vinfo;

        cudaFree(info->z);
        cudaFree(info->pBuffer);
        cusparseDestroyMatDescr(info->descr_M);
        cusparseDestroyMatDescr(info->descr_L);
        cusparseDestroyMatDescr(info->descr_U);
        cusparseDestroyCsrilu02Info(info->info_M);
        cusparseDestroySolveAnalysisInfo(info->info_L);
        cusparseDestroySolveAnalysisInfo(info->info_U);

        delete info;
      }

      // BCSR
      struct CudaIluBSolveInfo
      {
        cusparseMatDescr_t descr_M;
        cusparseMatDescr_t descr_L;
        cusparseMatDescr_t descr_U;
        bsrilu02Info_t info_M;
        bsrsv2Info_t info_L;
        bsrsv2Info_t info_U;
        cusparseOperation_t trans_L;
        cusparseOperation_t trans_U;
        cusparseDirection_t dir;
        cusparseSolvePolicy_t policy_M;
        cusparseSolvePolicy_t policy_L;
        cusparseSolvePolicy_t policy_U;
        void * pBuffer;
        double * z;
        int m;
        int nnz;
        int blocksize;
      };

      void * cuda_ilub_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd, const int blocksize)
      {
        double * z;
        cudaMalloc((void**)&z, m * blocksize * sizeof(double));

        cusparseMatDescr_t descr_M = 0;
        cusparseMatDescr_t descr_L = 0;
        cusparseMatDescr_t descr_U = 0;
        bsrilu02Info_t info_M  = 0;
        bsrsv2Info_t  info_L  = 0;
        bsrsv2Info_t  info_U  = 0;
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
        const cusparseDirection_t dir = CUSPARSE_DIRECTION_ROW;

        cusparseStatus_t status;

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

        cusparseCreateBsrilu02Info(&info_M);
        cusparseCreateBsrsv2Info(&info_L);
        cusparseCreateBsrsv2Info(&info_U);

        cusparseDbsrilu02_bufferSize(Util::Intern::cusparse_handle, dir, m, nnz,
                descr_M, csrVal, csrRowPtr, csrColInd, blocksize, info_M, &pBufferSize_M);
        cusparseDbsrsv2_bufferSize(Util::Intern::cusparse_handle, dir, trans_L, m, nnz,
                descr_L, csrVal, csrRowPtr, csrColInd, blocksize, info_L, &pBufferSize_L);
        cusparseDbsrsv2_bufferSize(Util::Intern::cusparse_handle, dir, trans_U, m, nnz,
                descr_U, csrVal, csrRowPtr, csrColInd, blocksize, info_U, &pBufferSize_U);

        pBufferSize = max(pBufferSize_M, max(pBufferSize_L, pBufferSize_U));

        cudaMalloc((void**)&pBuffer, pBufferSize);

        status = cusparseDbsrilu02_analysis(Util::Intern::cusparse_handle, dir, m, nnz, descr_M,
                csrVal, csrRowPtr, csrColInd, blocksize, info_M,
                    policy_M, pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrilu02 failed with status code: " + stringify(status));
        status = cusparseXbsrilu02_zeroPivot(Util::Intern::cusparse_handle, info_M, &structural_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "CUSPARSE ZERO PIVOT ERROR!");
        }

        status = cusparseDbsrsv2_analysis(Util::Intern::cusparse_handle, dir, trans_L, m, nnz, descr_L,
                csrVal, csrRowPtr, csrColInd, blocksize, info_L, policy_L, pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrv2_analysis failed with status code: " + stringify(status));

        status = cusparseDbsrsv2_analysis(Util::Intern::cusparse_handle, dir, trans_U, m, nnz, descr_U,
                csrVal, csrRowPtr, csrColInd, blocksize, info_U, policy_U, pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrv2_analysis failed with status code: " + stringify(status));


#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        CudaIluBSolveInfo * info = new CudaIluBSolveInfo;
        info->descr_M = descr_M;
        info->descr_L = descr_L;
        info->descr_U = descr_U;
        info->info_M  = info_M;
        info->info_L  = info_L;
        info->info_U  = info_U;
        info->trans_L = trans_L;
        info->trans_U = trans_U;
        info->dir = dir;
        info->policy_M = policy_M;
        info->policy_L = policy_L;
        info->policy_U = policy_U;
        info->pBuffer = pBuffer;
        info->z = z;
        info->m = m;
        info->nnz = nnz;
        info->blocksize = blocksize;

        return (void*)info;
      }

      void cuda_ilub_init_numeric(double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo)
      {
        CudaIluBSolveInfo * info = (CudaIluBSolveInfo *) vinfo;

        cusparseStatus_t status = cusparseDbsrilu02(Util::Intern::cusparse_handle, info->dir, info->m, info->nnz, info->descr_M,
                csrVal, csrRowPtr, csrColInd, info->blocksize, info->info_M, info->policy_M, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrilu02 failed with status code: " + stringify(status));
        int numerical_zero;
        status = cusparseXbsrilu02_zeroPivot(Util::Intern::cusparse_handle, info->info_M, &numerical_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "CUSPARSE ZERO PIVOT ERROR!");
        }
      }

      int cuda_ilub_apply(double * y, const double * x, double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo)
      {
        CudaIluBSolveInfo * info = (CudaIluBSolveInfo *) vinfo;
        const double alpha = 1.;

        cusparseStatus_t status = cusparseDbsrsv2_solve(Util::Intern::cusparse_handle, info->dir, info->trans_L, info->m, info->nnz, &alpha, info->descr_L,
               csrVal, csrRowPtr, csrColInd, info->blocksize, info->info_L,
                  x, info->z, info->policy_L, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrsv2_solve failed with status code: " + stringify(status));

        status = cusparseDbsrsv2_solve(Util::Intern::cusparse_handle, info->dir, info->trans_U, info->m, info->nnz, &alpha, info->descr_U,
               csrVal, csrRowPtr, csrColInd, info->blocksize, info->info_U,
                  info->z, y, info->policy_U, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsebsrsv2_solve failed with status code: " + stringify(status));

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }

      void cuda_ilub_done_symbolic(void * vinfo)
      {
        CudaIluBSolveInfo * info = (CudaIluBSolveInfo *) vinfo;

        cudaFree(info->z);
        cudaFree(info->pBuffer);
        cusparseDestroyMatDescr(info->descr_M);
        cusparseDestroyMatDescr(info->descr_L);
        cusparseDestroyMatDescr(info->descr_U);
        cusparseDestroyBsrilu02Info(info->info_M);
        cusparseDestroyBsrsv2Info(info->info_L);
        cusparseDestroyBsrsv2Info(info->info_U);

        delete info;
      }
    } // namespace Intern
    /// \endcond
  } // namespace Solver
} // namespace FEAT
