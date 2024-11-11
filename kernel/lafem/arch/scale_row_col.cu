// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/scale_row_col.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_scale_rows_csr(DT_ * r, const DT_ * b, const DT_ * val, const IT_ * col_ind,
                                          const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          r[i] = val[i] * b[idx];
        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_scale_rows_bcsr(DT_ * r, const DT_ * b, const DT_ * val, const IT_ * col_ind,
                                          const IT_ * row_ptr, const Index count, const int bh_, const int bw_)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const IT_ end(row_ptr[idx + 1]);
        for (IT_ i(row_ptr[idx]) ; i < end ; ++i)
        {

          for (int h(0) ; h < bh_ ; ++h)
          {
            for (int w(0) ; w < bw_ ; ++w)
            {
              r[i * bh_ * bw_ + h * bw_ + w] = val[i * bh_ * bw_ + h * bw_ + w] * b[idx*bh_ + h];
            }
          }

        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_scale_cols_csr(DT_ * r, const DT_ * b, const DT_ * val, const IT_ * col_ind,
                                          const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          r[i] = val[i] * b[col_ind[i]];
        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_scale_cols_bcsr(DT_ * r, const DT_ * b, const DT_ * val, const IT_ * col_ind,
                                          const IT_ * row_ptr, const Index count, const int bh_, const int bw_)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const IT_ end(row_ptr[idx + 1]);
        for (IT_ i(row_ptr[idx]) ; i < end ; ++i)
        {

          for (int h(0) ; h < bh_ ; ++h)
          {
            for (int w(0) ; w < bw_ ; ++w)
            {
              r[i * bh_ * bw_ + h * bw_ + w] = val[i * bh_ * bw_ + h * bw_ + w] * b[col_ind[i]*bw_ + w];
            }
          }

        }
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_, typename IT_>
void ScaleRows::csr_cuda(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index /*columns*/, const Index /*used_elements*/)
{
  Index blocksize = Util::cuda_blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_scale_rows_csr<<<grid, block>>>(r, x, val, col_ind, row_ptr, rows);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ScaleRows::csr_cuda(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index);
template void ScaleRows::csr_cuda(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index);
template void ScaleRows::csr_cuda(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index);
template void ScaleRows::csr_cuda(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void ScaleRows::bcsr_cuda_intern(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index /*columns*/, const Index /*used_elements*/, const int bh_, const int bw_)
{
  Index blocksize = Util::cuda_blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_scale_rows_bcsr<<<grid, block>>>(r, x, val, col_ind, row_ptr, rows, bh_, bw_);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ScaleRows::bcsr_cuda_intern(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index, const int, const int);
template void ScaleRows::bcsr_cuda_intern(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index, const int, const int);
template void ScaleRows::bcsr_cuda_intern(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index, const int, const int);
template void ScaleRows::bcsr_cuda_intern(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index, const int, const int);



template <typename DT_, typename IT_>
void ScaleCols::csr_cuda(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index /*columns*/, const Index /*used_elements*/)
{
  Index blocksize = Util::cuda_blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_scale_cols_csr<<<grid, block>>>(r, x, val, col_ind, row_ptr, rows);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ScaleCols::csr_cuda(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index);
template void ScaleCols::csr_cuda(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index);
template void ScaleCols::csr_cuda(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index);
template void ScaleCols::csr_cuda(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void ScaleCols::bcsr_cuda_intern(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index /*columns*/, const Index /*used_elements*/, const int bh_, const int bw_)
{
  Index blocksize = Util::cuda_blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_scale_cols_bcsr<<<grid, block>>>(r, x, val, col_ind, row_ptr, rows, bh_, bw_);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ScaleCols::bcsr_cuda_intern(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index, const int, const int);
template void ScaleCols::bcsr_cuda_intern(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index, const int, const int);
template void ScaleCols::bcsr_cuda_intern(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index, const int, const int);
template void ScaleCols::bcsr_cuda_intern(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index, const int, const int);
