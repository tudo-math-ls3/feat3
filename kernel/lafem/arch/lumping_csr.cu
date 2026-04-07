// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/lumping_csr.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_lumping_csr(DT_ * lump, const DT_ * val, const IT_ * col_idx, const IT_ * row_ptr, const Index rows)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < rows; i += blockDim.x * gridDim.x)
    {
      DT_ sum(0);
      const Index end(row_ptr[i + 1]);
      for (Index j(row_ptr[i]) ; j < end ; ++j)
      {
        sum += val[j];
      }
      lump[i] = sum;
    }
  }

  template <typename DT_, typename IT_>
  void LumpingCSR::exec_cuda_impl(DT_ * lump, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    cuda_lumping_csr<<<grid, block>>>(lump, val, col_idx, row_ptr, rows);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void LumpingCSR::exec_cuda_impl(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
  template void LumpingCSR::exec_cuda_impl(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index);

  template void LumpingCSR::exec_cuda_impl(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
  template void LumpingCSR::exec_cuda_impl(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
} // namespace FEAT::LAFEM::Arch
