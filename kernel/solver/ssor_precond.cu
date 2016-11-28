// includes, FEAT
#include <kernel/base_header.hpp>

#ifdef FEAT_HAVE_CUSOLVER
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

#include <vector>
#include <algorithm>

#include "cusparse_v2.h"

using namespace FEAT;

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      __global__ void cuda_ssor_forward_apply_kernel(int m, double * y, const double * x,
          double * csrVal, int * csrRowPtr, int * csrColInd, double omega, int * inverseRowPtr)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= m)
          return;

        int row = inverseRowPtr[idx];
        double d(0.);
        int col;
        for (col = csrRowPtr[idx*2] ; csrColInd[col] < row ; ++col)
        {
          d += csrVal[col] * y[csrColInd[col]];
        }
        y[row] = (x[row] - omega * d) / csrVal[col];
      }

      __global__ void cuda_ssor_backward_apply_kernel(int m, double * y, const double * x,
          double * csrVal, int * csrRowPtr, int * csrColInd, double omega, int * inverseRowPtr)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= m)
          return;

        int row = inverseRowPtr[idx];
        double d(0.);
        int col;
        for (col = csrRowPtr[idx*2+1] - 1 ; csrColInd[col] > row ; --col)
        {
          d += csrVal[col] * y[csrColInd[col]];
        }
        y[row] -= omega * d / csrVal[col];
      }

      int cuda_ssor_forward_apply(int m, double * y, const double * x, double * csrVal, int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        cudaMemset(y, 0, m * sizeof(double));

        int row_offset(0);
        for (int i(0) ; i < ncolors ; ++i)
        {
          Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
          dim3 grid;
          dim3 block;
          block.x = blocksize;
          grid.x = (unsigned)ceil((rows_per_color[i])/(double)(block.x));

          cuda_ssor_forward_apply_kernel<<<grid, block>>>(rows_per_color[i], y, x, csrVal, colored_row_ptr + row_offset * 2, csrColInd, omega, inverse_row_ptr + row_offset);
          row_offset += rows_per_color[i];
        }

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }

      int cuda_ssor_backward_apply(int m, double * y, const double * x, double * csrVal, int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        int row_offset(0);
        for (int i(0) ; i < ncolors ; ++i)
        {
          Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
          dim3 grid;
          dim3 block;
          block.x = blocksize;
          grid.x = (unsigned)ceil((rows_per_color[i])/(double)(block.x));

          cuda_ssor_backward_apply_kernel<<<grid, block>>>(rows_per_color[i], y, x, csrVal, colored_row_ptr + row_offset * 2, csrColInd, omega, inverse_row_ptr + row_offset);
          row_offset += rows_per_color[i];
        }

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }
    }
  }
}
#endif // FEAT_HAVE_CUSOLVER
