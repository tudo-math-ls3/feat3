// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_block_vec.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/cuda_math.cuh>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_unit_filter_block_vec(DT_ * v, const DT_ * f_val, const IT_ * f_idx, const Index f_nzes, const bool zero, const bool ign_nans, const Index block_size)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < f_nzes; i += blockDim.x * gridDim.x)
    {
      for(Index j(0) ; j < block_size; ++j)
      {
        const Index k = block_size * i + j;

        // skip if filter value is NaN if desired
        if(ign_nans && CudaMath::cuda_isnan(f_val[k]))
          continue;

        // write filter value or zero
        v[Index(block_size) * f_idx[i] + Index(j)] = zero ? DT_(0) : f_val[k];
      }
    }
  }

  template <typename DT_, typename IT_>
  void UnitFilterBlockVec::exec_cuda_impl(DT_ * v, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool zero, const bool ign_nans, const int block_size)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
    dim3 grid((static_cast<unsigned int>(f_nzes) + block.x - 1u) / block.x); // round up to next integer

    cuda_unit_filter_block_vec<<<grid, block>>>(v, f_val, f_idx, f_nzes, zero, ign_nans, Index(block_size));

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void UnitFilterBlockVec::exec_cuda_impl(float *, const float * const, const std::uint64_t * const, const Index, const bool, const bool, const int);
  template void UnitFilterBlockVec::exec_cuda_impl(double *, const double * const, const std::uint64_t * const, const Index, const bool, const bool, const int);
  template void UnitFilterBlockVec::exec_cuda_impl(float *, const float * const, const std::uint32_t * const, const Index, const bool, const bool, const int);
  template void UnitFilterBlockVec::exec_cuda_impl(double *, const double * const, const std::uint32_t * const, const Index, const bool, const bool, const int);
} // namespace FEAT::LAFEM::Arch
