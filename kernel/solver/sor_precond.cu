// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>

#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/cuda_util.hpp>

#include <vector>
#include <algorithm>

using namespace FEAT;


#include "cusparse_v2.h"

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      __global__ void cuda_sor_apply_kernel(int m, double * y, const double * x,
          const double * csrVal, const int * csrRowPtr, const int * csrColInd, double omega, int * inverseRowPtr)
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
        y[row] = omega * (x[row] - d) / csrVal[col];
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
      __global__ void cuda_sor_bcsr_apply_kernel(int m, double * y, const double * x,
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

        //y[row * BlockSize_ + j] = omega * (x[row * BlockSize_ + j] - d[j]) / csrVal[col];
        double temp[BlockSize_];
        for (int j(0) ; j < BlockSize_ ; ++j)
        {
          temp[j] = x[row * BlockSize_ + j] - d[j];
        }

        for (int h(0) ; h < BlockSize_ ; ++h)
        {
          double r(0);
          for (int w(0) ; w < BlockSize_ ; ++w)
          {
            r += inv[h][w] * temp[w];
          }
          y[row * BlockSize_ + h] = omega * r;
        }
      }

      void cuda_sor_init_symbolic(int m, int nnz, const double * csrVal, const int * csrRowPtr, const int * csrColInd, int & ncolors, int* & colored_row_ptr, int* & rows_per_color, int* & inverse_row_ptr)
      {
        cusparseColorInfo_t cinfo;
        cusparseStatus_t status = cusparseCreateColorInfo(&cinfo);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseCreateColorInfo failed with status code: " + stringify(status));

        int * d_coloring = (int*)Util::cuda_malloc(m * sizeof(int));
        double one(1.0);

        cusparseMatDescr_t descr_M = 0;
        cusparseCreateMatDescr(&descr_M);
        cusparseSetMatIndexBase(descr_M, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);

        status = cusparseDcsrcolor(Util::Intern::cusparse_handle, m, nnz, descr_M,
            csrVal, csrRowPtr, csrColInd, &one, &ncolors, d_coloring, NULL, cinfo);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseDcsrcolor failed with status code: " + stringify(status));

        status = cusparseDestroyColorInfo(cinfo);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseDestroyColorInfo failed with status code: " + stringify(status));

        //std::cout<<"pre colors: "<<ncolors<<" rows: "<<m<<std::endl;
        int * coloring = new int[m];
        Util::cuda_copy_device_to_host(coloring, d_coloring, m * sizeof(int));
        Util::cuda_free(d_coloring);

        //decrement from non existing colors
        for (int color(0) ; color < ncolors ; ++color)
        {
          //search for each color
          bool found(false);
          for (int i(0) ; i < m ; ++i)
          {
            if (coloring[i] == color)
            {
              found = true;
              break;
            }
          }
          // if color was not found, decrement all remaining colors and retry search with the same color again
          if (!found)
          {
            for (int i(0) ; i < m ; ++i)
            {
              if (coloring[i] > color)
              {
                --coloring[i];
              }
            }
            --ncolors;
            --color;
          }
        }

        /*std::cout<<"row to color:"<<std::endl;
        for (int i(0) ; i < m ; ++i)
        {
          std::cout<<coloring[i]<<" ";
        }
        std::cout<<std::endl;*/
        //std::cout<<"colors: "<<ncolors<<" rows: "<<m<<std::endl;

        // count rows per color
        //rows_per_color = MemoryPool::template allocate_memory<int>(ncolors);
        rows_per_color = new int[ncolors];
        for (int i(0) ; i < ncolors ; ++i)
        {
          rows_per_color[i] = 0;
        }
        for (int i(0) ; i < m ; ++i)
        {
          rows_per_color[coloring[i]] += 1;
        }

        int * colors_ascending = MemoryPool::template allocate_memory<int>(ncolors);
        //vector of pair<rows, color>
        std::vector<std::pair<int, int>> temp;
        for (int i(0) ; i < ncolors ; ++i)
        {
          temp.push_back(std::make_pair(rows_per_color[i], i));
        }
        std::sort(temp.begin(), temp.end());

        for (int i(0) ; i < ncolors ; ++i)
        {
          colors_ascending[i] = temp[i].second;
        }

        //resort rows per color by ascending row count
        for (int i(0) ; i < ncolors ; ++i)
        {
          rows_per_color[i] = temp[i].first;
        }
        /*std::cout<<"colors ascending: "<<std::endl;
        for (int i(0) ; i < ncolors ; ++i)
          std::cout<<i<<" "<<colors_ascending[i]<<std::endl;
        std::cout<<std::endl;

        std::cout<<"rows per color:"<<std::endl;
        for (int i(0) ; i < ncolors ; ++i)
        {
          std::cout<<i<<" "<<rows_per_color[i]<<std::endl;
        }
        std::cout<<std::endl;*/

        int * host_irp = MemoryPool::template allocate_memory<int>(m);
        int * host_crp = MemoryPool::template allocate_memory<int>(2*m);
        int * host_row_ptr = MemoryPool::template allocate_memory<int>(m+1);
        Util::cuda_copy_device_to_host(host_row_ptr, csrRowPtr, (m+1) * sizeof(int));

        //iterate over all colors, by ascending row count
        int crp_i(0); //index into host_crp
        for (int color_i(0) ; color_i < ncolors ; ++color_i)
        {
          int color = colors_ascending[color_i];
          //search all rows with matching color
          for (int row(0) ; row < m ; ++row)
          {
            if (coloring[row] == color)
            {
              host_crp[crp_i*2] = host_row_ptr[row];
              host_crp[crp_i*2+1] = host_row_ptr[row+1];
              host_irp[crp_i] = row;
              ++crp_i;
            }
          }
        }

        MemoryPool::release_memory(host_row_ptr);

        inverse_row_ptr = (int*)Util::cuda_malloc(m * sizeof(int));
        Util::cuda_copy_host_to_device(inverse_row_ptr, host_irp, m * sizeof(int));

        colored_row_ptr = (int*)Util::cuda_malloc(2 * m * sizeof(int));
        Util::cuda_copy_host_to_device(colored_row_ptr, host_crp, 2 * m * sizeof(int));

        delete[] coloring;
        MemoryPool::release_memory(colors_ascending);
        MemoryPool::release_memory(host_irp);
        MemoryPool::release_memory(host_crp);

        cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
      }

      void cuda_sor_done_symbolic(int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        Util::cuda_free(colored_row_ptr);
        Util::cuda_free(inverse_row_ptr);
        //MemoryPool::release_memory(rows_per_color);
        delete[] rows_per_color;

        cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
      }

      int cuda_sor_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        cudaMemset(y, 0, m * sizeof(double));

        int row_offset(0);
        for (int i(0) ; i < ncolors ; ++i)
        {
          Index blocksize = Util::cuda_blocksize_spmv;
          dim3 grid;
          dim3 block;
          block.x = (unsigned)blocksize;
          grid.x = (unsigned)ceil((rows_per_color[i])/(double)(block.x));

          cuda_sor_apply_kernel<<<grid, block>>>(rows_per_color[i], y, x, csrVal, colored_row_ptr + row_offset * 2, csrColInd, omega, inverse_row_ptr + row_offset);
          row_offset += rows_per_color[i];
        }

        cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }

      template<int BlockSize_>
      int cuda_sor_bcsr_apply(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        cudaMemset(y, 0, m * sizeof(double));

        int row_offset(0);
        for (int i(0) ; i < ncolors ; ++i)
        {
          Index blocksize = Util::cuda_blocksize_spmv;
          dim3 grid;
          dim3 block;
          block.x = (unsigned)blocksize;
          grid.x = (unsigned)ceil((rows_per_color[i])/(double)(block.x));

          cuda_sor_bcsr_apply_kernel<BlockSize_><<<grid, block>>>(rows_per_color[i], y, x, csrVal, colored_row_ptr + row_offset * 2, csrColInd, omega, inverse_row_ptr + row_offset);
          row_offset += rows_per_color[i];
        }

        cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif

        return 0;
      }
      // cuda_sor_bcsr_apply_kernel is hardcoded for BS_ == 2 && BS_ == 3
      template int cuda_sor_bcsr_apply<2>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
      template int cuda_sor_bcsr_apply<3>(int m, double * y, const double * x, const double * csrVal, const int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
    }
  }
}
