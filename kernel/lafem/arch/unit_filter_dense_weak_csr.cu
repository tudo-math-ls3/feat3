// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_dense_weak_csr.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template<typename DT_, typename IT_>
  __global__ void cuda_unit_filter_dense_weak_csr(
    DT_ * __restrict__ mat_a, const DT_ * __restrict__ const mat_m, const IT_* __restrict__ const row_ptr,
    const DT_ * __restrict__ const f_val, const IT_ * __restrict__ const f_idx, const Index f_nzes)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < f_nzes; i += blockDim.x * gridDim.x)
    {
      const IT_ ix(f_idx[i]);
      const DT_ vx(f_val[i]);
      for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
      {
        mat_a[j] = vx * mat_m[j];
      }
    }
  }

  template<typename DT_, typename IT_>
  void UnitFilterDenseWeakCSR::exec_cuda_impl(DT_* mat_a, const DT_ * const mat_m, const IT_* const row_ptr,
    const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_misc));
    dim3 grid((static_cast<unsigned int>(f_nzes) + block.x - 1u) / block.x); // round up to next integer

    cuda_unit_filter_dense_weak_csr<<<grid, block>>>(mat_a, mat_m, row_ptr, f_val, f_idx, f_nzes);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void UnitFilterDenseWeakCSR::exec_cuda_impl<float, std::uint64_t>(float *, const float* const, const std::uint64_t * const, const float * const, const std::uint64_t * const, const Index);
  template void UnitFilterDenseWeakCSR::exec_cuda_impl<float, std::uint32_t>(float *, const float* const, const std::uint32_t * const, const float * const, const std::uint32_t * const, const Index);
  template void UnitFilterDenseWeakCSR::exec_cuda_impl<double, std::uint64_t>(double *, const double* const, const std::uint64_t * const, const double * const, const std::uint64_t * const, const Index);
  template void UnitFilterDenseWeakCSR::exec_cuda_impl<double, std::uint32_t>(double *, const double* const, const std::uint32_t * const, const double * const, const std::uint32_t * const, const Index);
} // namespace FEAT::LAFEM::Arch
