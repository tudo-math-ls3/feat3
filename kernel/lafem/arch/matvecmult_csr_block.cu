// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_csr_block.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/half.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_, int block_size_>
  __global__ void cuda_apply_csrsb(DT_ * r, const DT_ alpha, const DT_ * x, const DT_ * y, const DT_ * val, const IT_ * col_idx, const IT_ * row_ptr, const Index count)
  {
    Index idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx >= count)
      return;

    DT_ bsum[block_size_];
    for (int j(0) ; j < block_size_ ; ++j)
    {
      bsum[j] = DT_(0);
    }
    const IT_ end(row_ptr[idx + 1]);
    for (IT_ i(row_ptr[idx]) ; i < end ; ++i)
    {
      const DT_ vali(val[i]);
      const IT_ coli(col_idx[i] * block_size_);
      for (int j(0) ; j < block_size_ ; ++j)
      {
        bsum[j] += vali * x[coli + j];
      }
    }

    if(y != nullptr)
    {
      for (int j(0) ; j < block_size_ ; ++j)
      {
        r[idx * block_size_ + j] = y[idx * block_size_ + j] + alpha * bsum[j];
      }
    }
    else
    {
      for (int j(0) ; j < block_size_ ; ++j)
      {
        r[idx * block_size_ + j] = alpha * bsum[j];
      }
    }
  }

  template <typename DT_, typename IT_, int bs_>
  void MatVecMultCSRBlock::exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
  {
    XASSERTM(!transposed, "transposed mat-vec-mult not implemented yet");
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    cuda_apply_csrsb<DT_, IT_, bs_><<<grid, block>>>(r, alpha, x, y, val, col_idx, row_ptr, rows);

#ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  }

  template void MatVecMultCSRBlock::exec_cuda_impl<float, std::uint64_t,1>(float *, const float, const float *, const float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<double, std::uint64_t,1>(double *, const double, const double *, const double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<float, std::uint32_t,1>(float *, const float, const float *, const float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<double, std::uint32_t,1>(double *, const double, const double *, const double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<float, std::uint64_t,2>(float *, const float, const float *, const float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<double, std::uint64_t,2>(double *, const double, const double *, const double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<float, std::uint32_t,2>(float *, const float, const float *, const float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<double, std::uint32_t,2>(double *, const double, const double *, const double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<float, std::uint64_t,3>(float *, const float, const float *, const float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<double, std::uint64_t,3>(double *, const double, const double *, const double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<float, std::uint32_t,3>(float *, const float, const float *, const float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  template void MatVecMultCSRBlock::exec_cuda_impl<double, std::uint32_t,3>(double *, const double, const double *, const double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
} // namespace FEAT::LAFEM::Arch
