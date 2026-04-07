// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/component_product_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  __global__ void cuda_component_product_dense(DT_ * r, const DT_ * x, const DT_ * y, const Index size)
  {
    for(Index i = threadIdx.x + blockDim.x * blockIdx.x; i < size; i += blockDim.x * gridDim.x)
    {
      r[i] = x[i] * y[i];
    }
  }

  template <typename DT_>
  void ComponentProductDense::exec_cuda_impl(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size)
  {
    dim3 block(static_cast<unsigned int>(Util::cuda_blocksize_axpy));
    dim3 grid((static_cast<unsigned int>(size) + block.x - 1u) / block.x); // round up to next integer

    cuda_component_product_dense<<<grid, block>>>(r, x, y, size);

#ifdef DEBUG
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
  }

#ifdef FEAT_HAVE_HALFMATH
  template void ComponentProductDense::exec_cuda_impl(Half *, const Half * const, const Half * const, const Index);
#endif
  template void ComponentProductDense::exec_cuda_impl(float *, const float * const, const float * const, const Index);
  template void ComponentProductDense::exec_cuda_impl(double *, const double * const, const double * const, const Index);
} // namespace FEAT::LAFEM::Arch
