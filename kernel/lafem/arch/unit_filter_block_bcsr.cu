// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_block_bcsr.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/cuda_math.cuh>

namespace FEAT::LAFEM::Arch
{
  template<typename DT_, typename IT_>
  __global__ void cuda_unit_filter_blocked_mat(
    DT_* __restrict__ mat, const IT_* __restrict__ const row_ptr, const IT_* __restrict__ const col_idx,
    const DT_ * __restrict__ const f_val, const IT_ * __restrict__ const f_idx, const Index f_nzes,
    const bool unit, const bool ign_nans, const int block_height, const int block_width)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < f_nzes; i += blockDim.x * gridDim.x)
    {
      const IT_ ix(f_idx[i]);
      const DT_* vx(&f_val[i*block_height]);

      // replace by unit or zero row
      for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
      {
        // loop over rows in the block
        for(int k(0); k < block_height; ++k)
        {
          // possibly skip row if filter value is NaN
          if(ign_nans && CudaMath::cuda_isnan(vx[k]))
            continue;

          // format block row to zero
          for(int l(0); l < block_width; ++l)
            mat[j*block_height*block_width + k*block_width +l] = DT_(0);

          // replace diagonal entry by 1 if requested
          if(unit && (col_idx[j] == ix) && (k < block_width))
            mat[j*block_height*block_width + k*block_width +k] = DT_(1);
        }
      }
    }
  }

  template<typename DT_, typename IT_>
  void UnitFilterBlockBCSR::exec_cuda_impl(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx,
    const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool unit,
    const bool ign_nans, const int block_height, const int block_width)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
    dim3 grid((static_cast<unsigned int>(f_nzes) + block.x - 1u) / block.x); // round up to next integer

    cuda_unit_filter_blocked_mat<<<grid, block >>>(mat, row_ptr, col_idx, f_val, f_idx, f_nzes, unit, ign_nans, block_height, block_width);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void UnitFilterBlockBCSR::exec_cuda_impl(float* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, const float * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  template void UnitFilterBlockBCSR::exec_cuda_impl(float* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, const float * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  template void UnitFilterBlockBCSR::exec_cuda_impl(double* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, const double * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  template void UnitFilterBlockBCSR::exec_cuda_impl(double* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, const double * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
} // namespace FEAT::LAFEM::Arch
