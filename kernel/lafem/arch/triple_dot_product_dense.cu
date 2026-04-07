// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/triple_dot_product_dense.hpp>
#include <kernel/lafem/arch/dot_product_dense.hpp>
#include <kernel/lafem/arch/component_product_dense.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

namespace FEAT::LAFEM::Arch
{
  template <typename DT_>
  DT_ TripleDotProductDense::exec_cuda_impl(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size)
  {
    DT_ * temp = (DT_*)Util::cuda_get_static_memory(size * sizeof(DT_));
    ComponentProductDense::exec_cuda_impl(temp, y, z, size);
    DT_ result = DotProductDense::exec_cuda_impl(x, temp, size);

  #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t last_error(cudaGetLastError());
    if (cudaSuccess != last_error)
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
  #endif
    return result;
  }

  #ifdef FEAT_HAVE_HALFMATH
  template Half TripleDotProductDense::exec_cuda_impl(const Half * const x, const Half * const y, const Half * const z, const Index size);
  #endif
  template float TripleDotProductDense::exec_cuda_impl(const float * const x, const float * const y, const float * const z, const Index size);
  template double TripleDotProductDense::exec_cuda_impl(const double * const x, const double * const y, const double * const z, const Index size);
} // namespace FEAT::LAFEM::Arch
