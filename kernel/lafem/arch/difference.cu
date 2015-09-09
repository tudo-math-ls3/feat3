// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_difference(DT_ * r, const DT_ * x, const DT_ * y, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        r[idx] = x[idx] - y[idx];
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Difference<Mem::CUDA>::value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size)
{
  Index blocksize = Util::MemoryPool<Mem::CUDA>::blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_difference<<<grid, block>>>(r, x, y, size);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void Difference<Mem::CUDA>::value(float *, const float * const, const float * const, const Index);
template void Difference<Mem::CUDA>::value(double *, const double * const, const double * const, const Index);
