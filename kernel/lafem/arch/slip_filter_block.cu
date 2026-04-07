// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/slip_filter_block.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_slip_filter(DT_ * v, const DT_ * f_nu, const IT_ * f_idx, const Index f_nzes, const Index block_size)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < f_nzes; i += blockDim.x * gridDim.x)
    {
      DT_ sp(DT_(0));
      DT_ scal(DT_(0));

      for(Index j(0) ; j < block_size; ++j)
      {
        sp += v[block_size* f_idx[i] + j]*f_nu[block_size * i + j];
        scal += f_nu[block_size * i + j]*f_nu[block_size * i + j];
      }

      sp /= scal;

      for(Index j(0) ; j < block_size; ++j)
        v[block_size* f_idx[i] + j] -= sp*f_nu[block_size * i + j];
    }
  }

  template <typename DT_, typename IT_>
  void SlipFilterBlock::exec_cuda_impl(DT_ * v, const DT_ * const f_nu, const IT_ * const f_idx, const Index f_nzes, const int bs)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
    dim3 grid((static_cast<unsigned int>(f_nzes) + block.x - 1u) / block.x); // round up to next integer

    cuda_slip_filter<DT_, IT_><<<grid, block>>>(v, f_nu, f_idx, f_nzes,  Index(bs));

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void SlipFilterBlock::exec_cuda_impl<float, std::uint32_t>(float *, const float * const, const std::uint32_t * const, const Index, const int);
  template void SlipFilterBlock::exec_cuda_impl<float, std::uint64_t>(float *, const float * const, const std::uint64_t * const, const Index, const int);
  template void SlipFilterBlock::exec_cuda_impl<double, std::uint32_t>(double *, const double * const, const std::uint32_t * const, const Index, const int);
  template void SlipFilterBlock::exec_cuda_impl<double, std::uint64_t>(double *, const double * const, const std::uint64_t * const, const Index, const int);

  /// \endcond
} // namespace FEAT::LAFEM::Arch
