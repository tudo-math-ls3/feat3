// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>

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
#if CUSPARSE_VER_MAJOR < 12
        cusparseMatDescr_t descr_L;
        cusparseMatDescr_t descr_U;
#else
        cusparseSpMatDescr_t descr_L;
        cusparseSpMatDescr_t descr_U;
        //in this case, we also need handler for input and output vectors
        //cusparseDnVecDescr_t descr_X;  //should not be needed, since buffersize and anaylsis accept NULL as vectors...
        //cusparseDnVecDescr_t descr_Y;

#endif
        csrilu02Info_t info_M;
#if CUSPARSE_VER_MAJOR < 12
        csrsv2Info_t  info_L;
        csrsv2Info_t  info_U;
#else
        cusparseSpSVDescr_t  info_L;
        cusparseSpSVDescr_t  info_U;
#endif

        int pBufferSize_M;
#if CUSPARSE_VER_MAJOR < 12
        int pBufferSize_L;
        int pBufferSize_U;
#else
        size_t pBufferSize_L;
        size_t pBufferSize_U;
#endif
        int pBufferSize;
        void *pBuffer;
        int structural_zero;
        int numerical_zero;
        const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
        const cusparseOperation_t trans_U  = CUSPARSE_OPERATION_NON_TRANSPOSE;
        const cusparseSolvePolicy_t policy_M = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
        const cusparseSolvePolicy_t policy_L = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
        const cusparseSolvePolicy_t policy_U = CUSPARSE_SOLVE_POLICY_USE_LEVEL; //why?
        double * z;
        int m;
        int nnz;
      };

      void * cuda_ilu_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd)
      {
        CudaIluSolveInfo * info = new CudaIluSolveInfo;
        info->m = m;
        info->nnz = nnz;

        info->z = (double*)Util::cuda_malloc(m * sizeof(double));


        cusparseStatus_t status;

        cusparseCreateMatDescr(&(info->descr_M));
        cusparseSetMatIndexBase(info->descr_M, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(info->descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);

#if CUSPARSE_VER_MAJOR < 12
        cusparseCreateMatDescr(&(info->descr_L));
        cusparseSetMatIndexBase(info->descr_L, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(info->descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatFillMode(info->descr_L, CUSPARSE_FILL_MODE_LOWER);
        cusparseSetMatDiagType(info->descr_L, CUSPARSE_DIAG_TYPE_UNIT);
#else
        //assertion if int is 32 bits
        static_assert(sizeof(int) == 4u, "ERROR: Size of int is not 32 bits");
        cusparseCreateCsr(&(info->descr_L), info->m, info->m, info->nnz,
                                      csrRowPtr, csrColInd, csrVal,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,   //use typedef somewhere? Since this goes wrong, if int is something other than 32 bits...
                                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);     //Also variable size in theroy...
        //set attributes
        {
          cusparseFillMode_t fillmode = CUSPARSE_FILL_MODE_LOWER;
          cusparseDiagType_t diagtype = CUSPARSE_DIAG_TYPE_UNIT;
          cusparseSpMatSetAttribute(info->descr_L, CUSPARSE_SPMAT_FILL_MODE, &fillmode, sizeof(cusparseFillMode_t)); //set relevant data, rest should be set by default due to new CSR implementation...
          cusparseSpMatSetAttribute(info->descr_L, CUSPARSE_SPMAT_DIAG_TYPE, &diagtype, sizeof(cusparseSpMatAttribute_t));
        }
#endif

#if CUSPARSE_VER_MAJOR < 12
        cusparseCreateMatDescr(&(info->descr_U));
        cusparseSetMatIndexBase(info->descr_U, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(info->descr_U, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatFillMode(info->descr_U, CUSPARSE_FILL_MODE_UPPER);
        cusparseSetMatDiagType(info->descr_U, CUSPARSE_DIAG_TYPE_NON_UNIT);
#else
        cusparseCreateCsr(&(info->descr_U), info->m, info->m, info->nnz,
                                      csrRowPtr, csrColInd, csrVal,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
        {
          cusparseFillMode_t fillmode = CUSPARSE_FILL_MODE_UPPER;
          cusparseDiagType_t diagtype = CUSPARSE_DIAG_TYPE_NON_UNIT;
          cusparseSpMatSetAttribute(info->descr_U, CUSPARSE_SPMAT_FILL_MODE, &fillmode, sizeof(cusparseFillMode_t)); //set relevant data, rest should be set by default due to new CSR implementation...
          cusparseSpMatSetAttribute(info->descr_U, CUSPARSE_SPMAT_DIAG_TYPE, &diagtype, sizeof(cusparseSpMatAttribute_t));
        }
#endif

        cusparseCreateCsrilu02Info(&(info->info_M));
#if CUSPARSE_VER_MAJOR < 12
        cusparseCreateCsrsv2Info(&(info->info_L));
        cusparseCreateCsrsv2Info(&(info->info_U));
#else
        // create information handler for cuSparseSolver
        cusparseSpSV_createDescr(&(info->info_L));
        cusparseSpSV_createDescr(&(info->info_U));
#endif

// #if CUSPARSE_VER_MAJOR >= 12
//         //for now, we need to create pseudo vector arrays... we will set these later to the real vector by transfering the data pointer
//         cusparseCreateDnVec(&(info->descr_X), info->m, nullptr, CUDA_R_64F);
//         cusparseCreateDnVec(&(info->descr_Y), info->m, nullptr, CUDA_R_64F);
// #endif

        status = cusparseDcsrilu02_bufferSize(Util::Intern::cusparse_handle, m, nnz,
                info->descr_M, csrVal, csrRowPtr, csrColInd, info->info_M, &(info->pBufferSize_M));
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrilu02_bufferSize failed with status code: " + stringify(status));
#if CUSPARSE_VER_MAJOR < 12
        status = cusparseDcsrsv2_bufferSize(Util::Intern::cusparse_handle, info->trans_L, m, nnz,
            info->descr_L, csrVal, csrRowPtr, csrColInd, info->info_L, &(info->pBufferSize_L));
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseDcsrsv2_bufferSize failed with status code: " + stringify(status));

        status = cusparseDcsrsv2_bufferSize(Util::Intern::cusparse_handle, info->trans_U, m, nnz,
            info->descr_U, csrVal, csrRowPtr, csrColInd, info->info_L, &(info->pBufferSize_U)); //TODO: Error using info_L here?
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseDcsrsv2_bufferSize failed with status code: " + stringify(status));
#else
        const double alpha = 1.;
        status = cusparseSpSV_bufferSize(Util::Intern::cusparse_handle, info->trans_L, &alpha,
            info->descr_L, NULL /*info->descr_X*/, NULL /*info->descr_Y*/, CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, info->info_L, &(info->pBufferSize_L));
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpSV_bufferSize failed with status code: " + stringify(status));

        status = cusparseSpSV_bufferSize(Util::Intern::cusparse_handle, info->trans_U, &alpha,
            info->descr_U, NULL /*info->descr_X*/, NULL /*info->descr_Y*/, CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, info->info_U, &(info->pBufferSize_U));
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpSV_bufferSize failed with status code: " + stringify(status));
#endif
        info->pBufferSize = max(info->pBufferSize_M, int(max(info->pBufferSize_L, info->pBufferSize_U)));
        info->pBuffer = Util::cuda_malloc(info->pBufferSize_M);

        status = cusparseDcsrilu02_analysis(Util::Intern::cusparse_handle, m, nnz, info->descr_M,
                csrVal, csrRowPtr, csrColInd, info->info_M,
                    info->policy_M, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrilu02_analysis failed with status code: " + stringify(status));
        status = cusparseXcsrilu02_zeroPivot(Util::Intern::cusparse_handle, info->info_M, &(info->structural_zero));
        if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "CUSPARSE ZERO PIVOT ERROR!");
        }
#if CUSPARSE_VER_MAJOR< 12
        status = cusparseDcsrsv2_analysis(Util::Intern::cusparse_handle, info->trans_L, m, nnz, info->descr_L,
                csrVal, csrRowPtr, csrColInd, info->info_L, info->policy_L, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparse_csrv_analysis failed with status code: " + stringify(status));

        status = cusparseDcsrsv2_analysis(Util::Intern::cusparse_handle, info->trans_U, m, nnz, info->descr_U,
                csrVal, csrRowPtr, csrColInd, info->info_U, info->policy_U, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparse_csrv_analysis failed with status code: " + stringify(status));
#else
        status = cusparseSpSV_analysis(Util::Intern::cusparse_handle, info->trans_L, &alpha,
                              info->descr_L, NULL /*info->descr_X*/, NULL /*info->descr_Y*/, CUDA_R_64F,
                              CUSPARSE_SPSV_ALG_DEFAULT, info->info_L, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpSV_analysis failed with status code: " + stringify(status));

        status = cusparseSpSV_analysis(Util::Intern::cusparse_handle, info->trans_U, &alpha,
                              info->descr_U, NULL /*info->descr_X*/, NULL /*info->descr_Y*/, CUDA_R_64F,
                              CUSPARSE_SPSV_ALG_DEFAULT, info->info_U, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpSV_analysis failed with status code: " + stringify(status));
#endif

        cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return (void*)info;
      }

      void cuda_ilu_init_numeric(double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo)
      {
        CudaIluSolveInfo * info = (CudaIluSolveInfo *) vinfo;

        cusparseStatus_t status = cusparseDcsrilu02(Util::Intern::cusparse_handle, info->m, info->nnz, info->descr_M,
                csrVal, csrRowPtr, csrColInd, info->info_M, info->policy_M, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrilu02 failed with status code: " + stringify(status));
        status = cusparseXcsrilu02_zeroPivot(Util::Intern::cusparse_handle, info->info_M, &(info->numerical_zero));
        if (CUSPARSE_STATUS_ZERO_PIVOT == status)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "CUSPARSE ZERO PIVOT ERROR!");
        }
      }

      int cuda_ilu_apply(double * y, const double * x, double * csrVal, int * csrRowPtr, int * csrColInd, void * vinfo)
      {
        CudaIluSolveInfo * info = (CudaIluSolveInfo *) vinfo;
        const double alpha = 1.;
#if CUSPARSE_VER_MAJOR < 12
        cusparseStatus_t status = cusparseDcsrsv2_solve(Util::Intern::cusparse_handle, info->trans_L, info->m, info->nnz, &alpha, info->descr_L,
               csrVal, csrRowPtr, csrColInd, info->info_L,
                  x, info->z, info->policy_L, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseDcsrsv2_solve failed with status code: " + stringify(status));

        status = cusparseDcsrsv2_solve(Util::Intern::cusparse_handle, info->trans_U, info->m, info->nnz, &alpha, info->descr_U,
               csrVal, csrRowPtr, csrColInd, info->info_U,
                  info->z, y, info->policy_U, info->pBuffer);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsr2_solve failed with status code: " + stringify(status));
#else
        //we have to create vector handlers to use the vector data in the new api (or shift them into our vector handlers... lets see whats necessary)
        cusparseConstDnVecDescr_t descr_X;
        cusparseDnVecDescr_t descr_Y, descr_Z;
        cusparseCreateConstDnVec(&descr_X, info->m, x, CUDA_R_64F);
        cusparseCreateDnVec(&descr_Z, info->m, info->z, CUDA_R_64F); //first write into z...
        cusparseCreateDnVec(&descr_Y, info->m, y, CUDA_R_64F);
        //now solve first triang system
        cusparseStatus_t status = cusparseSpSV_solve(Util::Intern::cusparse_handle, info->trans_L, &alpha,
                           info->descr_L, descr_X, descr_Z, CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, info->info_L);
        if (status != CUSPARSE_STATUS_SUCCESS)
        {
          //delete vecs descr if someone catches the error
          cusparseDestroyDnVec(descr_Y);
          cusparseDestroyDnVec(descr_Z);
          cusparseDestroyDnVec(descr_X);
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpSV_solve failed with status code: " + stringify(status));
        }

        status = cusparseSpSV_solve(Util::Intern::cusparse_handle, info->trans_U, &alpha,
                           info->descr_U, descr_Z, descr_Y, CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, info->info_U);
        if (status != CUSPARSE_STATUS_SUCCESS)
        {
          cusparseDestroyDnVec(descr_Y);
          cusparseDestroyDnVec(descr_Z);
          cusparseDestroyDnVec(descr_X);
          throw InternalError(__func__, __FILE__, __LINE__, "SECONDcusparseSpSV_solve failed with status code: " + stringify(status));
        }
        //destroy descr
        cusparseDestroyDnVec(descr_Y);
        cusparseDestroyDnVec(descr_Z);
        cusparseDestroyDnVec(descr_X);
#endif

        cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }

      void cuda_ilu_done_symbolic(void * vinfo)
      {
        CudaIluSolveInfo * info = (CudaIluSolveInfo *) vinfo;

        Util::cuda_free(info->z);
        Util::cuda_free(info->pBuffer);
        cusparseDestroyMatDescr(info->descr_M);
#if CUSPARSE_VER_MAJOR < 12
        cusparseDestroyMatDescr(info->descr_L);
        cusparseDestroyMatDescr(info->descr_U);
#else
        cusparseDestroySpMat(info->descr_L);
        cusparseDestroySpMat(info->descr_U);
#endif
        cusparseDestroyCsrilu02Info(info->info_M);
#if CUSPARSE_VER_MAJOR < 12
        cusparseDestroyCsrsv2Info(info->info_L);
        cusparseDestroyCsrsv2Info(info->info_U);
#else
        cusparseSpSV_destroyDescr(info->info_L);
        cusparseSpSV_destroyDescr(info->info_L);
#endif

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
        double * z = (double*)Util::cuda_malloc(m * blocksize * sizeof(double));

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

        pBuffer = Util::cuda_malloc(pBufferSize);

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

        cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
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

        cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }

      void cuda_ilub_done_symbolic(void * vinfo)
      {
        CudaIluBSolveInfo * info = (CudaIluBSolveInfo *) vinfo;

        Util::cuda_free(info->z);
        Util::cuda_free(info->pBuffer);
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
