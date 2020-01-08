// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>

#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

#include <vector>
#include <algorithm>

using namespace FEAT;

#ifdef FEAT_HAVE_CUSOLVER

#include "cusparse_v2.h"

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      __global__ void cuda_ssor_forward_apply_kernel(int m, double * y, const double * x,
          const double * csrVal, int * csrRowPtr, const int * csrColInd, double omega, int * inverseRowPtr)
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

      __global__ void cuda_ssor_backward_apply_kernel(int m, double * y, const double * /*x*/,
          const double * csrVal, int * csrRowPtr, const int * csrColInd, double omega, int * inverseRowPtr)
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

      template<int n_>
      struct InverseHelper;

      template<>
      struct InverseHelper<2>
      {
        __device__ static void compute(double (&inv)[2][2], const double * a)
        {
          double d = double(1) / (a[0]*a[3] - a[1]*a[2]);
          inv[0][0] =  d*a[3];
          inv[0][1] = -d*a[1];
          inv[1][0] = -d*a[2];
          inv[1][1] =  d*a[0];
        }
      };

      template<>
      struct InverseHelper<3>
      {
        __device__ static void compute(double (&inv)[3][3], const double * a)
        {
          inv[0][0] = a[4]*a[8] - a[5]*a[7];
          inv[1][0] = a[5]*a[6] - a[3]*a[8];
          inv[2][0] = a[3]*a[7] - a[4]*a[6];
          double d = double(1) / (a[0]*inv[0][0] + a[1]*inv[1][0] + a[2]*inv[2][0]);
          inv[0][0] *= d;
          inv[1][0] *= d;
          inv[2][0] *= d;
          inv[0][1] = d*(a[2]*a[7] - a[1]*a[8]);
          inv[1][1] = d*(a[0]*a[8] - a[2]*a[6]);
          inv[2][1] = d*(a[1]*a[6] - a[0]*a[7]);
          inv[0][2] = d*(a[1]*a[5] - a[2]*a[4]);
          inv[1][2] = d*(a[2]*a[3] - a[0]*a[5]);
          inv[2][2] = d*(a[0]*a[4] - a[1]*a[3]);
        }
      };

      template <int BlockSize_>
      __global__ void cuda_ssor_forward_bcsr_apply_kernel(int m, double * y, const double * x,
          const double * csrVal, const int * csrRowPtr, const int * csrColInd, double omega, int * inverseRowPtr)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= m)
          return;

        int row = inverseRowPtr[idx];
        double d[BlockSize_];
        for (int j(0) ; j < BlockSize_ ; ++j)
        {
          d[j] = double(0);
        }
        int col;
        for (col = csrRowPtr[idx*2] ; csrColInd[col] < row ; ++col)
        {
          for (int h(0) ; h < BlockSize_ ; ++h)
          {
            for (int w(0) ; w < BlockSize_ ; ++w)
            {
              d[h] += csrVal[col * BlockSize_ * BlockSize_ + h * BlockSize_ + w] * y[csrColInd[col] * BlockSize_ + w];
            }
          }
        }

        double inv[BlockSize_][BlockSize_];
        InverseHelper<BlockSize_>::compute(inv, csrVal + col * BlockSize_ * BlockSize_);

        double temp[BlockSize_];
        for (int j(0) ; j < BlockSize_ ; ++j)
        {
          temp[j] = x[row * BlockSize_ + j] - (omega * d[j]);
        }

        for (int h(0) ; h < BlockSize_ ; ++h)
        {
          double r(0);
          for (int w(0) ; w < BlockSize_ ; ++w)
          {
            r += inv[h][w] * temp[w];
          }
          y[row * BlockSize_ + h] = r;
        }
      }

      template <int BlockSize_>
      __global__ void cuda_ssor_backward_bcsr_apply_kernel(int m, double * y, const double * x,
          const double * csrVal, const int * csrRowPtr, const int * csrColInd, double omega, int * inverseRowPtr)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= m)
          return;

        int row = inverseRowPtr[idx];
        double d[BlockSize_];
        for (int j(0) ; j < BlockSize_ ; ++j)
        {
          d[j] = double(0);
        }
        int col;
        for (col = csrRowPtr[idx*2+1] - 1 ; csrColInd[col] > row ; --col)
        {
          for (int h(0) ; h < BlockSize_ ; ++h)
          {
            for (int w(0) ; w < BlockSize_ ; ++w)
            {
              d[h] += csrVal[col * BlockSize_ * BlockSize_ + h * BlockSize_ + w] * y[csrColInd[col] * BlockSize_ + w];
            }
          }
        }

        double inv[BlockSize_][BlockSize_];
        InverseHelper<BlockSize_>::compute(inv, csrVal + col * BlockSize_ * BlockSize_);

        for (int h(0) ; h < BlockSize_ ; ++h)
        {
          double r(0);
          for (int w(0) ; w < BlockSize_ ; ++w)
          {
            r += inv[h][w] * d[w];
          }
          y[row * BlockSize_ + h] -= omega * r;
        }
      }

      int cuda_ssor_forward_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
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

      int cuda_ssor_backward_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
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

      template <int BlockSize_>
      int cuda_ssor_forward_bcsr_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
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

          cuda_ssor_forward_bcsr_apply_kernel<BlockSize_><<<grid, block>>>(rows_per_color[i], y, x, csrVal, colored_row_ptr + row_offset * 2, csrColInd, omega, inverse_row_ptr + row_offset);
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
      // cuda_sor_bcsr_apply_kernel is hardcoded for BS_ == 2 && BS_ == 3
      template int cuda_ssor_forward_bcsr_apply<2>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
      template int cuda_ssor_forward_bcsr_apply<3>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);

      template <int BlockSize_>
      int cuda_ssor_backward_bcsr_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
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

          cuda_ssor_backward_bcsr_apply_kernel<BlockSize_><<<grid, block>>>(rows_per_color[i], y, x, csrVal, colored_row_ptr + row_offset * 2, csrColInd, omega, inverse_row_ptr + row_offset);
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
      // cuda_sor_bcsr_apply_kernel is hardcoded for BS_ == 2 && BS_ == 3
      template int cuda_ssor_backward_bcsr_apply<2>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
      template int cuda_ssor_backward_bcsr_apply<3>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
    }
  }
}

#else // FEAT_HAVE_CUSOLVER

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      int cuda_ssor_forward_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA SSOR not supported in cuda before version 7!\n");
        //return 0;
      }

      int cuda_ssor_backward_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA SSOR not supported in cuda before version 7!\n");
        //return 0;
      }

      template<int BlockSize_>
      int cuda_ssor_forward_bcsr_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA SSOR not supported in cuda before version 7!\n");
        //return 0;
      }
      // cuda_sor_bcsr_apply_kernel is hardcoded for BS_ == 2 && BS_ == 3
      template int cuda_ssor_forward_bcsr_apply<2>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
      template int cuda_ssor_forward_bcsr_apply<3>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);

      template<int BlockSize_>
      int cuda_ssor_backward_bcsr_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA SSOR not supported in cuda before version 7!\n");
        //return 0;
      }
      // cuda_sor_bcsr_apply_kernel is hardcoded for BS_ == 2 && BS_ == 3
      template int cuda_ssor_backward_bcsr_apply<2>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
      template int cuda_ssor_backward_bcsr_apply<3>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
    }
  }
}

#endif // FEAT_HAVE_CUSOLVER
