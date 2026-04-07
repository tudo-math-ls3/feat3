// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_dense_vec.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_unit_filter_dense_vec(DT_ * v, const DT_ * f_val, const IT_ * idx, const Index f_nzes, const bool zero)
  {
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < f_nzes; i += blockDim.x * gridDim.x)
    {
      v[idx[i]] = zero ? DT_(0) : f_val[i];
    }
  }

  template <typename DT_, typename IT_>
  void UnitFilterDenseVec::exec_cuda_impl(DT_ * v, const DT_ * const f_val, const IT_ * const idx, const Index f_nzes, const bool zero)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
    dim3 grid((static_cast<unsigned int>(f_nzes) + block.x - 1u) / block.x); // round up to next integer

    cuda_unit_filter_dense_vec<<<grid, block>>>(v, f_val, idx, f_nzes, zero);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void UnitFilterDenseVec::exec_cuda_impl(float *, const float * const, const std::uint64_t * const, const Index, const bool);
  template void UnitFilterDenseVec::exec_cuda_impl(double *, const double * const, const std::uint64_t * const, const Index, const bool);
  template void UnitFilterDenseVec::exec_cuda_impl(float *, const float * const, const std::uint32_t * const, const Index, const bool);
  template void UnitFilterDenseVec::exec_cuda_impl(double *, const double * const, const std::uint32_t * const, const Index, const bool);
} // namespace FEAT::LAFEM::Arch
