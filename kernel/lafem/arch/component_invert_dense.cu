// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/component_invert_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  __global__ void cuda_component_invert_dense(DT_ * r, const DT_ * x, const DT_ alpha, const Index size)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < size; i += blockDim.x * gridDim.x)
    {
      r[i] = alpha / x[i];
    }
  }

#ifdef FEAT_HAVE_HALFMATH
  __global__ void cuda_component_invert_dense(Half * r, const Half * x, const Half alpha, const Index size)
  {
    // grid strided for loop
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < size; i += blockDim.x * gridDim.x)
    {
      r[i] = __half2float(alpha / x[i]);
    }
  }
#endif

  template <typename DT_>
  void ComponentInvertDense::exec_cuda_impl(DT_ * r, const DT_ * const x, const DT_ alpha, const Index size)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_axpy));
    dim3 grid((static_cast<unsigned int>(size) + block.x - 1u) / block.x); // round up to next integer

    cuda_component_invert_dense<<<grid, block>>>(r, x, alpha, size);

#ifdef DEBUG
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  }

#ifdef FEAT_HAVE_HALFMATH
  template void ComponentInvertDense::exec_cuda_impl(Half *, const Half * const, const Half, const Index);
#endif
  template void ComponentInvertDense::exec_cuda_impl(float *, const float * const, const float, const Index);
  template void ComponentInvertDense::exec_cuda_impl(double *, const double * const, const double, const Index);

} // namespace FEAT::LAFEM::Arch
