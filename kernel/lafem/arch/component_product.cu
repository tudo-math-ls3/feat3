// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/component_product.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/half.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_ComponentProduct(DT_ * r, const DT_ * x, const DT_ * y, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        r[idx] = x[idx] * y[idx];
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_>
void ComponentProduct::value_cuda(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size)
{
  Index blocksize = Util::cuda_blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_ComponentProduct<<<grid, block>>>(r, x, y, size);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
#ifdef FEAT_HAVE_HALFMATH
template void ComponentProduct::value_cuda(Half *, const Half * const, const Half * const, const Index);
#endif
template void ComponentProduct::value_cuda(float *, const float * const, const float * const, const Index);
template void ComponentProduct::value_cuda(double *, const double * const, const double * const, const Index);
