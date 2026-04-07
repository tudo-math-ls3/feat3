// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/diagonal_csr.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_diagonal_csr(DT_ * diag, const DT_ * val, const IT_ * col_idx, const IT_ * row_ptr, const Index rows)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < rows; i += blockDim.x * gridDim.x)
    {
      diag[i] = DT_(0);
      const Index end(row_ptr[i + 1]);
      for (Index j(row_ptr[i]); j < end; ++j)
      {
        if(i == col_idx[j])
        {
          diag[i] = val[j];
          return;
        }
      }
    }
  }

  template <typename DT_, typename IT_>
  void DiagonalCSR::exec_cuda_impl(DT_ * diag, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    cuda_diagonal_csr<<<grid, block>>>(diag, val, col_idx, row_ptr, rows);

  #ifdef DEBUG
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void DiagonalCSR::exec_cuda_impl(float*, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
  template void DiagonalCSR::exec_cuda_impl(float*, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
  template void DiagonalCSR::exec_cuda_impl(double*, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
  template void DiagonalCSR::exec_cuda_impl(double*, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
} // namespace FEAT::LAFEM::Arch
