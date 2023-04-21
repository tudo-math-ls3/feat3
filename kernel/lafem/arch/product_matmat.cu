// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/product_matmat.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

#include <cublas_v2.h>
#include <cublasLt.h>
#include <cusparse_v2.h>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
void ProductMatMat::dense_cuda(DT_ * r, const DT_ alpha, const DT_ beta,  const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index columns, const Index inner)
{
  if (r==y || r==x || x==y || z==x || z==y)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda ProductMatMat does not allow r==y or r==x or x==y or z==x or z==y!");

  cublasStatus_t status;

  // inspired by https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuBLASLt/LtSgemm/sample_cublasLt_LtSgemm.cu

  cublasLtMatmulDesc_t operationDesc = NULL;
  cublasLtMatrixLayout_t Rdesc = NULL, Adesc = NULL, Bdesc = NULL, Cdesc = NULL;
  cublasLtMatmulPreference_t preference = NULL;

  int algo_selector = -1;

  cudaDataType dt;
  cublasComputeType_t ct;
  if (typeid(DT_) == typeid(double))
  {
      dt = CUDA_R_64F;
      ct = CUBLAS_COMPUTE_64F;
      algo_selector = (rows > 1 && columns > 1 && inner > 1) ? 0 : 1;
  }
  else if (typeid(DT_) == typeid(float))
  {
      dt = CUDA_R_32F;
#if __CUDA_ARCH__ < 800
      ct = CUBLAS_COMPUTE_32F;
#else
      ct = CUBLAS_COMPUTE_32F_FAST_TF32;
#endif
      algo_selector = (rows > 1 && columns > 1 && inner > 1) ? 2 : 3;
  }
#ifdef FEAT_HAVE_HALFMATH
  else if (typeid(DT_) == typeid(Half))
  {
      dt = CUDA_R_16F;
      ct = CUBLAS_COMPUTE_16F;
      algo_selector = (rows > 1 && columns > 1 && inner > 1) ? 4 : 5;
  }
#endif
  else
    throw InternalError(__func__, __FILE__, __LINE__, "unsupported data type!");

  status = cublasLtMatmulDescCreate(&operationDesc, ct, dt);
  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  cublasLtOrder_t matrix_order = CUBLASLT_ORDER_ROW;
  status = cublasLtMatrixLayoutCreate(&Rdesc, dt, rows, columns, columns);
  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
  status = cublasLtMatrixLayoutSetAttribute(Rdesc, CUBLASLT_MATRIX_LAYOUT_ORDER, &matrix_order, sizeof(cublasLtOrder_t));
  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
  cublasLtMatrixLayoutCreate(&Adesc, dt, rows, inner, inner);
  cublasLtMatrixLayoutSetAttribute(Adesc, CUBLASLT_MATRIX_LAYOUT_ORDER, &matrix_order, sizeof(cublasLtOrder_t));
  cublasLtMatrixLayoutCreate(&Bdesc, dt, inner, columns, columns);
  cublasLtMatrixLayoutSetAttribute(Bdesc, CUBLASLT_MATRIX_LAYOUT_ORDER, &matrix_order, sizeof(cublasLtOrder_t));
  if (r!=z)
  {
    cublasLtMatrixLayoutCreate(&Cdesc, dt, rows, columns, columns);
    cublasLtMatrixLayoutSetAttribute(Cdesc, CUBLASLT_MATRIX_LAYOUT_ORDER, &matrix_order, sizeof(cublasLtOrder_t));
  }
  else // r==z -> in-place multiplication
  {
    Cdesc = Rdesc;
  }

  cublasLtMatmulAlgo_t * algo = NULL;
  if (! FEAT::Util::Intern::cublas_lt_algo_matmat_initialized[algo_selector])
  {
    int num_algos = 0;
    cublasLtMatmulHeuristicResult_t heuristic_algos = {};

    status = cublasLtMatmulPreferenceCreate(&preference);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
    //status = cublasLtMatmulPreferenceSetAttribute(preference, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES, &(FEAT::Util::Intern::cuda_workspace_size), sizeof(FEAT::Util::Intern::cuda_workspace_size));
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

    status = cublasLtMatmulAlgoGetHeuristic((cublasLtHandle_t)Util::Intern::cublas_handle, operationDesc, Adesc, Bdesc, Cdesc, Rdesc, preference, 1, &heuristic_algos, &num_algos);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

    if (num_algos == 0)
      throw InternalError(__func__, __FILE__, __LINE__, "no algo supports our matrices!");

    FEAT::Util::Intern::cublas_lt_algo_matmat[algo_selector] = heuristic_algos.algo;
    FEAT::Util::Intern::cublas_lt_algo_matmat_initialized[algo_selector] = true;
  }

  algo = &(FEAT::Util::Intern::cublas_lt_algo_matmat[algo_selector]);

  //status = cublasLtMatmul((cublasLtHandle_t)Util::Intern::cublas_handle, operationDesc, &alpha, x, Adesc, y, Bdesc, &beta, z, Cdesc, r, Rdesc, algo, FEAT::Util::Intern::cuda_workspace, FEAT::Util::Intern::cuda_workspace_size, 0);
  status = cublasLtMatmul((cublasLtHandle_t)Util::Intern::cublas_handle, operationDesc, &alpha, x, Adesc, y, Bdesc, &beta, z, Cdesc, r, Rdesc, algo, NULL, 0, 0);
  if (status != CUBLAS_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
#ifdef FEAT_HAVE_HALFMATH
template void ProductMatMat::dense_cuda(Half *, const Half, const Half, const Half * const, const Half * const, const Half * const, const Index, const Index, const Index);
#endif
template void ProductMatMat::dense_cuda(float *, const float, const float, const float * const, const float * const, const float * const, const Index, const Index, const Index);
template void ProductMatMat::dense_cuda(double *, const double, const double, const double * const, const double * const, const double * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void ProductMatMat::dsd_cuda(DT_ * r, const DT_ alpha, const DT_ beta, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index used_elements,
    const DT_ * y, const Index rows, const Index columns, const Index inner)
{
  if (r==y)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda ProductMatMat does not allow r==y!");

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
  status = cusparseCreateDnMat(&descr_r, rows, columns, columns, (void*)r, dt, CUSPARSE_ORDER_ROW);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

  cusparseSpMatDescr_t descr_x=0;
  status = cusparseCreateCsr(&descr_x, rows, inner, used_elements, (void*)row_ptr, (void*)col_ind, (void*)val, it, it, CUSPARSE_INDEX_BASE_ZERO, dt);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

  cusparseDnMatDescr_t descr_y=0;
  status = cusparseCreateDnMat(&descr_y, inner, columns, columns, (void*)y, dt, CUSPARSE_ORDER_ROW);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cusparseGetErrorString(status)));

  cusparseOperation_t trans = CUSPARSE_OPERATION_NON_TRANSPOSE;
  size_t buffer_size(0);
  status = cusparseSpMM_bufferSize(Util::Intern::cusparse_handle, trans, trans, &alpha, descr_x, descr_y, &beta, descr_r, ct, CUSPARSE_SPMM_CSR_ALG2, &buffer_size);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cusparsecsrmvex_buffersize failed with status code: " + stringify(cusparseGetErrorString(status)));

  void* buffer;
  cudaMalloc(&buffer, buffer_size);

  status = cusparseSpMM(Util::Intern::cusparse_handle, trans, trans, &alpha, descr_x, descr_y, &beta, descr_r, ct, CUSPARSE_SPMM_CSR_ALG2, buffer);
  if (status != CUSPARSE_STATUS_SUCCESS)
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseSpMM failed with status code: " + stringify(cusparseGetErrorString(status)));

  cusparseDestroyDnMat(descr_r);
  cusparseDestroySpMat(descr_x);
  cusparseDestroyDnMat(descr_y);
  cudaFree(buffer);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
#ifdef FEAT_HAVE_HALFMATH
template void ProductMatMat::dsd_cuda(Half *, const Half, const Half, const Half * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Half *, const Index, const Index, const Index);
#endif
template void ProductMatMat::dsd_cuda(float *, const float, const float, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const float *, const Index, const Index, const Index);
template void ProductMatMat::dsd_cuda(double *, const double, const double, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const double *, const Index, const Index, const Index);
#ifdef FEAT_HAVE_HALFMATH
template void ProductMatMat::dsd_cuda(Half *, const Half, const Half, const Half * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Half *, const Index, const Index, const Index);
#endif
template void ProductMatMat::dsd_cuda(float *, const float, const float, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const float *, const Index, const Index, const Index);
template void ProductMatMat::dsd_cuda(double *, const double, const double, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const double *, const Index, const Index, const Index);
