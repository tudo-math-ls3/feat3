// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/mirror_block_gather.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template<typename DT_, typename IT_>
  __global__ void cuda_mirror_gather_dvb(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < nidx; i += blockDim.x * gridDim.x)
    {
      for(Index k(0); k < bs; ++k)
      {
        buf[boff + i*bs + k] = vec[idx[i]*bs + k];
      }
    }
  }

  template<typename DT_, typename IT_>
  void MirrorBlockGather::exec_cuda_impl(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
    dim3 grid((static_cast<unsigned int>(nidx) + block.x - 1u) / block.x); // round up to next integer

    cuda_mirror_gather_dvb<<<grid, block>>>(bs, boff, nidx, idx, buf, vec);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void MirrorBlockGather::exec_cuda_impl(const Index, const Index, const Index, const std::uint64_t*, float*, const float*);
  template void MirrorBlockGather::exec_cuda_impl(const Index, const Index, const Index, const std::uint64_t*, double*, const double*);
  template void MirrorBlockGather::exec_cuda_impl(const Index, const Index, const Index, const std::uint32_t*, float*, const float*);
  template void MirrorBlockGather::exec_cuda_impl(const Index, const Index, const Index, const std::uint32_t*, double*, const double*);
} // namespace FEAT::LAFEM::Arch
