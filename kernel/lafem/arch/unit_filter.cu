// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/unit_filter.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

/// \cond internal
namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_unit_filter_rhs(DT_ * v, const DT_ * sv_elements, const IT_ * sv_indices, const Index ue)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= ue)
          return;
        v[sv_indices[idx]] = sv_elements[idx];
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_unit_filter_def(DT_ * v, const IT_ * sv_indices, const Index ue)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= ue)
          return;
        v[sv_indices[idx]] = DT_(0);
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_, typename IT_>
void UnitFilter<Mem::CUDA>::filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::MemoryPool<Mem::CUDA>::instance()->blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_unit_filter_rhs<<<grid, block>>>(v, sv_elements, sv_indices, ue);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilter<Mem::CUDA>::filter_rhs(float *, const float * const, const unsigned long * const, const Index);
template void UnitFilter<Mem::CUDA>::filter_rhs(double *, const double * const, const unsigned long * const, const Index);
template void UnitFilter<Mem::CUDA>::filter_rhs(float *, const float * const, const unsigned int * const, const Index);
template void UnitFilter<Mem::CUDA>::filter_rhs(double *, const double * const, const unsigned int * const, const Index);

template <typename DT_, typename IT_>
void UnitFilter<Mem::CUDA>::filter_def(DT_ * v, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::MemoryPool<Mem::CUDA>::instance()->blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_unit_filter_def<<<grid, block>>>(v, sv_indices, ue);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilter<Mem::CUDA>::filter_def(float *, const unsigned long * const, const Index);
template void UnitFilter<Mem::CUDA>::filter_def(double *,  const unsigned long * const, const Index);
template void UnitFilter<Mem::CUDA>::filter_def(float *,  const unsigned int * const, const Index);
template void UnitFilter<Mem::CUDA>::filter_def(double *,  const unsigned int * const, const Index);

/// \endcond
