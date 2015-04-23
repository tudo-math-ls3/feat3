// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/util/exception.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_scale(DT_ * r, const DT_ * x, const DT_ s, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        r[idx] = x[idx] * s;
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Scale<Mem::CUDA>::value(DT_ * r, const DT_ * const x, const DT_ s, const Index size)
{
  Index blocksize(256);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_scale<<<grid, block>>>(r, x, s, size);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Scale<Mem::CUDA>::value(float *, const float * const, const float, const Index);
template void Scale<Mem::CUDA>::value(double *, const double * const, const double, const Index);
