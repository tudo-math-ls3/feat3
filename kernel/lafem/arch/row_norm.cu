// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/row_norm.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_norm2_csr(DT_ * row_norms, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ norm(0);
        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          norm += val[i] * val[i];
        }
        row_norms[idx] = sqrt(norm);
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_norm2sqr_csr(DT_ * row_norms, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ norm(0);
        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          norm += val[i] * val[i];
        }
        row_norms[idx] = norm;
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_norm2sqr_scaled_csr(DT_ * row_norms, const DT_ * scal, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ norm(0);
        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          norm += val[i] * val[i] * scal[idx];
        }
        row_norms[idx] = norm;
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_norm2_bcsr(DT_ * row_norms, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index rows,
                                              const int BlockHeight, const int BlockWidth)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows * BlockHeight)
          return;

        Index csr_row = idx / BlockHeight;
        Index block_row = idx % BlockHeight;

        DT_ norm(0);
        const Index end(row_ptr[csr_row + 1]);
        for (Index i(row_ptr[csr_row]) ; i < end ; ++i)
        {
          for (Index w(0) ; w < BlockWidth ; ++w)
          {
            DT_ ival = val[BlockHeight*BlockWidth*i + block_row*BlockWidth + w];
            norm += ival * ival;
          }
        }
        row_norms[idx] = sqrt(norm);
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_norm2sqr_bcsr(DT_ * row_norms, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index rows,
                                              const int BlockHeight, const int BlockWidth)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows * BlockHeight)
          return;

        Index csr_row = idx / BlockHeight;
        Index block_row = idx % BlockHeight;

        DT_ norm(0);
        const Index end(row_ptr[csr_row + 1]);
        for (Index i(row_ptr[csr_row]) ; i < end ; ++i)
        {
          for (Index w(0) ; w < BlockWidth ; ++w)
          {
            DT_ ival = val[BlockHeight*BlockWidth*i + block_row*BlockWidth + w];
            norm += ival * ival;
          }
        }
        row_norms[idx] = norm;
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_norm2sqr_scaled_bcsr(DT_ * row_norms, const DT_ * scal, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index rows,
                                              const int BlockHeight, const int BlockWidth)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows * BlockHeight)
          return;

        Index csr_row = idx / BlockHeight;
        Index block_row = idx % BlockHeight;

        DT_ norm(0);
        const Index end(row_ptr[csr_row + 1]);
        for (Index i(row_ptr[csr_row]) ; i < end ; ++i)
        {
          for (Index w(0) ; w < BlockWidth ; ++w)
          {
            DT_ ival = val[BlockHeight*BlockWidth*i + block_row*BlockWidth + w];
            norm += ival * ival * scal[idx];
          }
        }
        row_norms[idx] = norm;
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_, typename IT_>
void RowNorm::csr_cuda_norm2(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_norm2_csr<<<grid, block>>>(row_norms, val, col_ind, row_ptr, rows);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void RowNorm::csr_cuda_norm2(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void RowNorm::csr_cuda_norm2(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void RowNorm::csr_cuda_norm2(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
template void RowNorm::csr_cuda_norm2(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index);

template <typename DT_, typename IT_>
void RowNorm::csr_cuda_norm2sqr(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_norm2sqr_csr<<<grid, block>>>(row_norms, val, col_ind, row_ptr, rows);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void RowNorm::csr_cuda_norm2sqr(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void RowNorm::csr_cuda_norm2sqr(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void RowNorm::csr_cuda_norm2sqr(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
template void RowNorm::csr_cuda_norm2sqr(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index);

template <typename DT_, typename IT_>
void RowNorm::csr_cuda_scaled_norm2sqr(DT_ * row_norms, const DT_ * const scal, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_norm2sqr_scaled_csr<<<grid, block>>>(row_norms, scal, val, col_ind, row_ptr, rows);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void RowNorm::csr_cuda_scaled_norm2sqr(float *, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void RowNorm::csr_cuda_scaled_norm2sqr(double *, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void RowNorm::csr_cuda_scaled_norm2sqr(float *, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
template void RowNorm::csr_cuda_scaled_norm2sqr(double *, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index);

template <typename DT_, typename IT_>
void RowNorm::bcsr_cuda_norm2(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows * BlockHeight)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_norm2_bcsr<<<grid, block>>>(row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void RowNorm::bcsr_cuda_norm2(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_norm2(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_norm2(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_norm2(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const int, const int);

template <typename DT_, typename IT_>
void RowNorm::bcsr_cuda_norm2sqr(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows * BlockHeight)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_norm2sqr_bcsr<<<grid, block>>>(row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void RowNorm::bcsr_cuda_norm2sqr(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_norm2sqr(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_norm2sqr(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_norm2sqr(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const int, const int);

template <typename DT_, typename IT_>
void RowNorm::bcsr_cuda_scaled_norm2sqr(DT_ * row_norms, const DT_ * const scal, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
{
  Index blocksize = Util::cuda_blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows * BlockHeight)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_norm2sqr_scaled_bcsr<<<grid, block>>>(row_norms, scal, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void RowNorm::bcsr_cuda_scaled_norm2sqr(float *, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_scaled_norm2sqr(double *, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_scaled_norm2sqr(float *, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const int, const int);
template void RowNorm::bcsr_cuda_scaled_norm2sqr(double *, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const int, const int);
