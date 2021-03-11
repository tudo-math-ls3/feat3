// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>

#ifdef FEAT_HAVE_CUSOLVER
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>

#include <cusolverSp.h>
#include "cusparse_v2.h"

using namespace FEAT;

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      int cuda_lu(int n, int nnzA, const double * csrValA, const int * csrRowPtrA, const int * csrColIndA,
          const double * b, double * x)
      {
        cusolverSpHandle_t handle;
        cusolverSpCreate(&handle);

        cusparseMatDescr_t descr;
        cusparseCreateMatDescr(&descr);
        cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

        int singularity;
        cusolverStatus_t status = cusolverSpDcsrlsvluHost(handle, n, nnzA, descr, csrValA, csrRowPtrA, csrColIndA, b, 0.0, 1, x, &singularity);
        if (status != CUSOLVER_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusolverSpDcsrlsvluHost failed with status code: " + stringify(status));

        cusparseDestroyMatDescr(descr);
        cusolverSpDestroy(handle);

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return (status != CUSOLVER_STATUS_SUCCESS);
      }


      int cuda_qr(int m, int nnz, const double * csrValA, const int * csrRowPtrA, const int * csrColIndA,
          const double * b, double * x)
      {
        cusolverSpHandle_t handle;
        cusolverSpCreate(&handle);

        cusparseMatDescr_t descr;
        cusparseCreateMatDescr(&descr);
        cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

        int singularity;
        cusolverStatus_t status = cusolverSpDcsrlsvqr(handle, m, nnz, descr, csrValA, csrRowPtrA, csrColIndA, b, 0.0, 1, x, &singularity);
        if (status != CUSOLVER_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusolverSPDcsrlvsqr failed with status code: " + stringify(status));

        cusparseDestroyMatDescr(descr);
        cusolverSpDestroy(handle);
#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return (status != CUSOLVER_STATUS_SUCCESS);
      }
    }
  }
}
#endif // FEAT_HAVE_CUSOLVER
