// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_banded_dense.hpp>
#include <kernel/util/memory_aux.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/half.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_, typename IT_>
  __global__ void cuda_apply_banded(DT_ * r, const DT_ alpha, const DT_ * x, const DT_ * val, const IT_ * offsets, const Index bands, const Index rows, const Index cols)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < rows; i += blockDim.x * gridDim.x)
    {
      const Index k1(rows - 1);
      const Index k2(rows + cols - 1);

      Index start(0);

      while ((start < bands) && (k1 > offsets[start] + i))
      {
        ++start;
      }

      // this may happen if the first rows of the matrix are empty
      if(bands <= start)
        continue;

      Index end(start);

      while (end < bands && i + offsets[end] < k2)
      {
        ++end;
      }

      DT_ sum(DT_(0.0));
      for (Index diag(start); diag < end; ++diag)
      {
        sum += val[rows * diag + i] * x[i + offsets[diag] - rows + 1];
      }
      r[i] += (sum*alpha);
    }
  }

  template <typename DT_, typename IT_>
  void MatVecMultBandedDense::exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y,
    const DT_ * const val, const IT_ * const offsets, const Index bands, const Index rows, const Index cols)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_spmv));
    dim3 grid((static_cast<unsigned int>(rows) + block.x - 1u) / block.x); // round up to next integer

    if(y == nullptr)
      Memory::memset_cuda(r, 0, rows * sizeof(DT_));
    else if (r != y)
      Memory::memcopy_cuda(r, y, rows * sizeof(DT_));

    cuda_apply_banded<<<grid, block>>>(r, alpha, x, val, offsets, bands, rows, cols);
  }

  template void MatVecMultBandedDense::exec_cuda_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint32_t * const, const Index, const Index, const Index);
  template void MatVecMultBandedDense::exec_cuda_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint32_t * const, const Index, const Index, const Index);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultBandedDense::exec_cuda_impl(Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint32_t * const, const Index, const Index, const Index);
  #endif

  template void MatVecMultBandedDense::exec_cuda_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint64_t * const, const Index, const Index, const Index);
  template void MatVecMultBandedDense::exec_cuda_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint64_t * const, const Index, const Index, const Index);
  #ifdef FEAT_HAVE_HALFMATH
  template void MatVecMultBandedDense::exec_cuda_impl(Half *, const Half, const Half * const, const Half * const, const Half * const, const std::uint64_t * const,const Index, const Index, const Index);
  #endif
} // namespace FEAT::LAFEM::Arch
