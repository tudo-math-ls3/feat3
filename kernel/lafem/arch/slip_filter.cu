// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/slip_filter.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_, int BlockSize_>
      __global__ void cuda_slip_filter_rhs(DT_ * v, const DT_ * sv_elements, const IT_ * sv_indices, const Index ue)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= ue)
          return;

        Index block_size = Index(BlockSize_);
        DT_ sp(DT_(0));
        DT_ scal(DT_(0));

        for(Index j(0) ; j < block_size; ++j)
        {
          sp += v[block_size* sv_indices[idx] + j]*sv_elements[block_size * idx + j];
          scal += sv_elements[block_size * idx + j]*sv_elements[block_size * idx + j];
        }

        sp /= scal;

        for(Index j(0) ; j < block_size; ++j)
          v[block_size* sv_indices[idx] + j] -= sp*sv_elements[block_size * idx + j];
      }

      template <typename DT_, typename IT_, int BlockSize_>
      __global__ void cuda_slip_filter_def(DT_ * v, const DT_ * sv_elements, const IT_ * sv_indices, const Index ue)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= ue)
          return;

        Index block_size = Index(BlockSize_);
        DT_ sp(DT_(0));
        DT_ scal(DT_(0));

        for(Index j(0) ; j < block_size; ++j)
        {
          sp += v[block_size* sv_indices[idx] + j]*sv_elements[block_size * idx + j];
          scal += sv_elements[block_size * idx + j]*sv_elements[block_size * idx + j];
        }

        sp /= scal;

        for(Index j(0) ; j < block_size; ++j)
          v[block_size* sv_indices[idx] + j] -= sp*sv_elements[block_size * idx + j];
      }

    } // namespace Intern
  } // namespace LAFEM
} // namespace FEAT

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_, typename IT_, int BlockSize_>
void SlipFilter<Mem::CUDA>::filter_rhs(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_slip_filter_rhs<DT_, IT_, BlockSize_><<<grid, block>>>(v, sv_elements, sv_indices, ue);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void SlipFilter<Mem::CUDA>::filter_rhs<float, unsigned long, 2>(float *, const float * const, const unsigned long * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_rhs<double, unsigned long, 2>(double *, const double * const, const unsigned long * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_rhs<float, unsigned int, 2>(float *, const float * const, const unsigned int * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_rhs<double, unsigned int, 2>(double *, const double * const, const unsigned int * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_rhs<float, unsigned long, 3>(float *, const float * const, const unsigned long * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_rhs<double, unsigned long, 3>(double *, const double * const, const unsigned long * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_rhs<float, unsigned int, 3>(float *, const float * const, const unsigned int * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_rhs<double, unsigned int, 3>(double *, const double * const, const unsigned int * const, const Index);

template <typename DT_, typename IT_, int BlockSize_>
void SlipFilter<Mem::CUDA>::filter_def(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_slip_filter_def<DT_, IT_, BlockSize_><<<grid, block>>>(v, sv_elements, sv_indices, ue);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void SlipFilter<Mem::CUDA>::filter_def<float, unsigned long, 2>(float *, const float * const, const unsigned long * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_def<double, unsigned long, 2>(double *, const double * const, const unsigned long * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_def<float, unsigned int, 2>(float *, const float * const, const unsigned int * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_def<double, unsigned int, 2>(double *, const double * const, const unsigned int * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_def<float, unsigned long, 3>(float *, const float * const, const unsigned long * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_def<double, unsigned long, 3>(double *, const double * const, const unsigned long * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_def<float, unsigned int, 3>(float *, const float * const, const unsigned int * const, const Index);
template void SlipFilter<Mem::CUDA>::filter_def<double, unsigned int, 3>(double *, const double * const, const unsigned int * const, const Index);

/// \endcond
