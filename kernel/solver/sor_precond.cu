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
      __global__ void cuda_sor_apply_kernel(int m, double * y, const double * x,
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
        y[row] = omega * (x[row] - d) / csrVal[col];
      }

      void cuda_sor_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd, int & ncolors, int* & colored_row_ptr, int* & rows_per_color, int* & inverse_row_ptr)
      {
        cusparseColorInfo_t cinfo;
        cusparseStatus_t status = cusparseCreateColorInfo(&cinfo);
        if (status != CUSPARSE_STATUS_SUCCESS)
          throw InternalError(__func__, __FILE__, __LINE__, "cusparseCreateColorInfo failed with status code: " + stringify(status));

        int * d_coloring;
        cudaMalloc(&d_coloring, m * sizeof(int));
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
        int * coloring = MemoryPool<Mem::Main>::template allocate_pinned_memory<int>(m);
        cudaMemcpy(coloring, d_coloring, m * sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(d_coloring);

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
        rows_per_color = MemoryPool<Mem::Main>::template allocate_memory<int>(ncolors);
        for (int i(0) ; i < ncolors ; ++i)
        {
          rows_per_color[i] = 0;
        }
        for (int i(0) ; i < m ; ++i)
        {
          rows_per_color[coloring[i]] += 1;
        }

        int * colors_ascending = MemoryPool<Mem::Main>::template allocate_memory<int>(ncolors);
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

        int * host_irp = MemoryPool<Mem::Main>::template allocate_pinned_memory<int>(m);
        int * host_crp = MemoryPool<Mem::Main>::template allocate_pinned_memory<int>(2*m);
        int * host_row_ptr = MemoryPool<Mem::Main>::template allocate_pinned_memory<int>(m+1);
        cudaMemcpy(host_row_ptr, csrRowPtr, (m+1) * sizeof(int), cudaMemcpyDeviceToHost);

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

        MemoryPool<Mem::Main>::release_memory(host_row_ptr);

        cudaMalloc(&inverse_row_ptr, m * sizeof(unsigned int));
        cudaMemcpy(inverse_row_ptr, host_irp, m * sizeof(int), cudaMemcpyHostToDevice);

        cudaMalloc(&colored_row_ptr, 2 * m * sizeof(int));
        cudaMemcpy(colored_row_ptr, host_crp, 2 * m * sizeof(int), cudaMemcpyHostToDevice);

        MemoryPool<Mem::Main>::release_memory(coloring);
        MemoryPool<Mem::Main>::release_memory(colors_ascending);
        MemoryPool<Mem::Main>::release_memory(host_irp);
        MemoryPool<Mem::Main>::release_memory(host_crp);

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
      }

      void cuda_sor_done_symbolic(int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        cudaFree(colored_row_ptr);
        cudaFree(inverse_row_ptr);
        MemoryPool<Mem::Main>::release_memory(rows_per_color);

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
      }

      int cuda_sor_apply(int m, double * y, const double * x, double * csrVal, int * csrColInd, int ncolors, double omega,
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

          cuda_sor_apply_kernel<<<grid, block>>>(rows_per_color[i], y, x, csrVal, colored_row_ptr + row_offset * 2, csrColInd, omega, inverse_row_ptr + row_offset);
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
#else
namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      void cuda_sor_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd, int & ncolors, int* & colored_row_ptr, int* & rows_per_color, int* & inverse_row_ptr)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA SOR not supported in cuda before version 7!\n");
      }

      void cuda_sor_done_symbolic(int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA SOR not supported in cuda before version 7!\n");
      }

      int cuda_sor_apply(int m, double * y, const double * x, double * csrVal, int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "CUDA SOR not supported in cuda before version 7!\n");
        return 0;
      }
    }
  }
}
#endif // FEAT_HAVE_CUSOLVER
