// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/row_norm2_bcsr.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_norm2_bcsr(DT_ * row_norms, const DT_ * val, const IT_ * col_idx,
    const IT_ * row_ptr, const Index rows, const int block_height, bool block_width, const bool squared)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < rows * block_height; i += blockDim.x * gridDim.x)
    {
      Index csr_row = i / block_height;
      Index block_row = i % block_height;

      DT_ norm(0);
      const Index end(row_ptr[csr_row + 1]);
      for (Index j(row_ptr[csr_row]) ; j < end ; ++j)
      {
        for (Index w(0) ; w < block_width ; ++w)
        {
          DT_ ival = val[block_height*block_width*j + block_row*block_width + w];
          norm += ival * ival;
        }
      }
      if(!squared)
        row_norms[i] = sqrt(norm);
    }
  }

  template <typename DT_, typename IT_>
  void RowNorm2BCSR::exec_cuda_impl(DT_ * row_norms, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const int block_height, const int block_width, const bool squared)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    cuda_norm2_bcsr<<<grid, block>>>(row_norms, val, col_idx, row_ptr, rows, block_height, block_width, squared);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void RowNorm2BCSR::exec_cuda_impl(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const int, const int, const bool);
  template void RowNorm2BCSR::exec_cuda_impl(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const int, const int, const bool);
  template void RowNorm2BCSR::exec_cuda_impl(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const int, const int, const bool);
  template void RowNorm2BCSR::exec_cuda_impl(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const int, const int, const bool);
} // namespace FEAT::LAFEM::Arch
