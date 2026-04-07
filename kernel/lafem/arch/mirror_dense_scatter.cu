// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/mirror_dense_scatter.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template<typename DT_, typename IT_>
  __global__ void cuda_mirror_scatter_dv(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < nidx; i += blockDim.x * gridDim.x)
      vec[idx[i]] += alpha*buf[boff + i];
  }

  template<typename DT_, typename IT_>
  void MirrorDenseScatter::exec_cuda_impl(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
    dim3 grid((static_cast<unsigned int>(nidx) + block.x - 1u) / block.x); // round up to next integer

    cuda_mirror_scatter_dv<<<grid, block>>>(boff, nidx, idx, buf, vec, alpha);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void MirrorDenseScatter::exec_cuda_impl(const Index, const Index, const std::uint64_t*, const float*, float*, float);
  template void MirrorDenseScatter::exec_cuda_impl(const Index, const Index, const std::uint64_t*, const double*, double*, double);
  template void MirrorDenseScatter::exec_cuda_impl(const Index, const Index, const std::uint32_t*, const float*, float*, float);
  template void MirrorDenseScatter::exec_cuda_impl(const Index, const Index, const std::uint32_t*, const double*, double*, double);
} // namespace FEAT::LAFEM::Arch
