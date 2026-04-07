// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/scale_row_bcsr.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_scale_rows_bcsr(DT_ * r, const DT_ * b, const DT_ * val, const IT_ * col_idx,
    const IT_ * row_ptr, const Index rows, const int bh_, const int bw_)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < rows; i += blockDim.x * gridDim.x)
    {
      const IT_ end(row_ptr[i + 1]);
      for (IT_ j(row_ptr[i]) ; j < end ; ++j)
      {
        for (int h(0) ; h < bh_ ; ++h)
        {
          for (int w(0) ; w < bw_ ; ++w)
          {
            r[j * bh_ * bw_ + h * bw_ + w] = val[j * bh_ * bw_ + h * bw_ + w] * b[i*bh_ + h];
          }
        }
      }
    }
  }

  template <typename DT_, typename IT_>
  void ScaleRowBCSR::exec_cuda_impl(DT_ * r, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index /*cols*/, const Index /*nonzeros*/, const int bh, const int bw)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    cuda_scale_rows_bcsr<<<grid, block>>>(r, x, val, col_idx, row_ptr, rows, bh, bw);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
  }

  template void ScaleRowBCSR::exec_cuda_impl(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index, const int, const int);
  template void ScaleRowBCSR::exec_cuda_impl(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index, const int, const int);
  template void ScaleRowBCSR::exec_cuda_impl(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index, const int, const int);
  template void ScaleRowBCSR::exec_cuda_impl(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index, const int, const int);
} // namespace FEAT::LAFEM::Arch
