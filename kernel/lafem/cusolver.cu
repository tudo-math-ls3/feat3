// includes, FEAST
#include <kernel/base_header.hpp>

#ifdef FEAST_HAVE_CUSOLVER
#include <kernel/archs.hpp>
#include <kernel/lafem/cusolver.hpp>
#include <kernel/util/exception.hpp>

#include <cusolverSp.h>
#include "cusparse_v2.h"

using namespace FEAST;
using namespace FEAST::LAFEM;

void CuSolverLU::solve_intern(int n, int nnzA, const double * csrValA, const int * csrRowPtrA, const int * csrColIndA,
    const double * b, double * x)
{
  cusolverSpHandle_t handle;
  cusolverSpCreate(&handle);

  cusparseMatDescr_t descr;
  cusparseCreateMatDescr(&descr);
  cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

  int singularity;
  cusolverSpDcsrlsvluHost(handle, n, nnzA, descr, csrValA, csrRowPtrA, csrColIndA, b, 0.0, 1, x, &singularity);

  cusparseDestroyMatDescr(descr);
  cusolverSpDestroy(handle);
}

void CuSolverQR::solve_intern(int m, int nnz, const double * csrValA, const int * csrRowPtrA, const int * csrColIndA,
    const double * b, double * x)
{
  cusolverSpHandle_t handle;
  cusolverSpCreate(&handle);

  cusparseMatDescr_t descr;
  cusparseCreateMatDescr(&descr);
  cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

  int singularity;
  cusolverSpDcsrlsvqr(handle, m, nnz, descr, csrValA, csrRowPtrA, csrColIndA, b, 0.0, 1, x, &singularity);

  cusparseDestroyMatDescr(descr);
  cusolverSpDestroy(handle);
}
#endif // FEAST_HAVE_CUSOLVER
