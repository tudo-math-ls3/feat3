// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matmatmult_dense_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

#include <cublas_v2.h>
#include <cublasLt.h>
#include <cusparse_v2.h>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  void MatMatMultDenseDense::exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index cols, const Index inner)
  {
    if (r==y || r==x || x==y || z==x || z==y)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda MatMatMultDenseDense does not allow r==y or r==x or x==y or z==x or z==y!");

    cublasStatus_t status;

    //XASSERTM(z != nullptr, "case z = nullptr not implemented yet");
    if(z == nullptr)
      Memory::memset_cuda(r, 0, sizeof(DT_)*rows*cols);

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
        algo_selector = (rows > 1 && cols > 1 && inner > 1) ? 0 : 1;
    }
    else if (typeid(DT_) == typeid(float))
    {
        dt = CUDA_R_32F;
  #if __CUDA_ARCH__ < 800
        ct = CUBLAS_COMPUTE_32F;
  #else
        ct = CUBLAS_COMPUTE_32F_FAST_TF32;
  #endif
        algo_selector = (rows > 1 && cols > 1 && inner > 1) ? 2 : 3;
    }
  #ifdef FEAT_HAVE_HALFMATH
    else if (typeid(DT_) == typeid(Half))
    {
        dt = CUDA_R_16F;
        ct = CUBLAS_COMPUTE_16F;
        algo_selector = (rows > 1 && cols > 1 && inner > 1) ? 4 : 5;
    }
  #endif
    else
      throw InternalError(__func__, __FILE__, __LINE__, "unsupported data type!");

    status = cublasLtMatmulDescCreate(&operationDesc, ct, dt);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

    cublasLtOrder_t matrix_order = CUBLASLT_ORDER_ROW;
    status = cublasLtMatrixLayoutCreate(&Rdesc, dt, rows, cols, cols);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
    status = cublasLtMatrixLayoutSetAttribute(Rdesc, CUBLASLT_MATRIX_LAYOUT_ORDER, &matrix_order, sizeof(cublasLtOrder_t));
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));
    cublasLtMatrixLayoutCreate(&Adesc, dt, rows, inner, inner);
    cublasLtMatrixLayoutSetAttribute(Adesc, CUBLASLT_MATRIX_LAYOUT_ORDER, &matrix_order, sizeof(cublasLtOrder_t));
    cublasLtMatrixLayoutCreate(&Bdesc, dt, inner, cols, cols);
    cublasLtMatrixLayoutSetAttribute(Bdesc, CUBLASLT_MATRIX_LAYOUT_ORDER, &matrix_order, sizeof(cublasLtOrder_t));
    if (r!=z)
    {
      cublasLtMatrixLayoutCreate(&Cdesc, dt, rows, cols, cols);
      cublasLtMatrixLayoutSetAttribute(Cdesc, CUBLASLT_MATRIX_LAYOUT_ORDER, &matrix_order, sizeof(cublasLtOrder_t));
    }
    else // r==z -> in-place multiplication
    {
      Cdesc = Rdesc;
    }

    cublasLtMatmulHeuristicResult_t algo_check_result;
    if (! FEAT::Util::Intern::cublas_lt_algo_matmat_initialized[algo_selector] ||
        CUBLAS_STATUS_SUCCESS != cublasLtMatmulAlgoCheck((cublasLtHandle_t)Util::Intern::cublas_handle, operationDesc, Adesc, Bdesc, Cdesc, Rdesc, &(FEAT::Util::Intern::cublas_lt_algo_matmat[algo_selector]), &algo_check_result))
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

    cublasLtMatmulAlgo_t * algo = &(FEAT::Util::Intern::cublas_lt_algo_matmat[algo_selector]);

    DT_ beta = DT_(z != nullptr ? 1 : 0);
    const DT_* zr = (z != nullptr ? z : r);
    //status = cublasLtMatmul((cublasLtHandle_t)Util::Intern::cublas_handle, operationDesc, &alpha, x, Adesc, y, Bdesc, &beta, z, Cdesc, r, Rdesc, algo, FEAT::Util::Intern::cuda_workspace, FEAT::Util::Intern::cuda_workspace_size, 0);
    status = cublasLtMatmul((cublasLtHandle_t)Util::Intern::cublas_handle, operationDesc, &alpha, x, Adesc, y, Bdesc, &beta, zr, Cdesc, r, Rdesc, algo, NULL, 0, 0);
    if (status != CUBLAS_STATUS_SUCCESS)
      throw InternalError(__func__, __FILE__, __LINE__, "cuda error: " + stringify(cublasGetStatusString(status)));

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  #ifdef FEAT_HAVE_HALFMATH
  template void MatMatMultDenseDense::exec_cuda_impl(Half *, const Half, const Half * const, const Half * const, const Half * const, const Index, const Index, const Index);
  #endif
  template void MatMatMultDenseDense::exec_cuda_impl(float *, const float, const float * const, const float * const, const float * const, const Index, const Index, const Index);
  template void MatMatMultDenseDense::exec_cuda_impl(double *, const double, const double * const, const double * const, const double * const, const Index, const Index, const Index);
} // namespace FEAT::LAFEM::Arch
