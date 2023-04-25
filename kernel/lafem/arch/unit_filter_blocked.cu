// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_blocked.hpp>
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
      __global__ void cuda_unit_filter_blocked_rhs(DT_ * v, const DT_ * sv_elements, const IT_ * sv_indices, const Index ue)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= ue)
          return;

        Index block_size = Index(BlockSize_);
        for(Index j(0) ; j < block_size; ++j)
          v[block_size* sv_indices[idx] + j] = sv_elements[block_size * idx + j];
      }

      template <typename DT_, typename IT_, int BlockSize_>
      __global__ void cuda_unit_filter_blocked_def(DT_ * v, const IT_ * sv_indices, const Index ue)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= ue)
          return;

        Index block_size = Index(BlockSize_);
        for(Index j(0) ; j < block_size; ++j)
          v[block_size * sv_indices[idx] + j] = DT_(0);
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <int BlockSize_, typename DT_, typename IT_>
void UnitFilterBlocked::filter_rhs_cuda(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_blocked_rhs<DT_, IT_, BlockSize_><<<grid, block>>>(v, sv_elements, sv_indices, ue);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilterBlocked::filter_rhs_cuda<1, float, std::uint64_t>(float *, const float * const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<1, double, std::uint64_t>(double *, const double * const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<1, float, std::uint32_t>(float *, const float * const, const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<1, double, std::uint32_t>(double *, const double * const, const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<2, float, std::uint64_t>(float *, const float * const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<2, double, std::uint64_t>(double *, const double * const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<2, float, std::uint32_t>(float *, const float * const, const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<2, double, std::uint32_t>(double *, const double * const, const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<3, float, std::uint64_t>(float *, const float * const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<3, double, std::uint64_t>(double *, const double * const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<3, float, std::uint32_t>(float *, const float * const, const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<3, double, std::uint32_t>(double *, const double * const, const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<4, float, std::uint64_t>(float *, const float * const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<4, double, std::uint64_t>(double *, const double * const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<4, float, std::uint32_t>(float *, const float * const, const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_rhs_cuda<4, double, std::uint32_t>(double *, const double * const, const std::uint32_t * const, const Index);

template <int BlockSize_, typename DT_, typename IT_>
void UnitFilterBlocked::filter_def_cuda(DT_ * v, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_blocked_def<DT_, IT_, BlockSize_><<<grid, block>>>(v, sv_indices, ue);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilterBlocked::filter_def_cuda<1, float, std::uint64_t>(float *, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<1, double, std::uint64_t>(double *,  const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<1, float, std::uint32_t>(float *,  const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<1, double, std::uint32_t>(double *,  const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<2, float, std::uint64_t>(float *, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<2, double, std::uint64_t>(double *,  const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<2, float, std::uint32_t>(float *,  const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<2, double, std::uint32_t>(double *,  const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<3, float, std::uint64_t>(float *, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<3, double, std::uint64_t>(double *,  const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<3, float, std::uint32_t>(float *,  const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<3, double, std::uint32_t>(double *,  const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<4, float, std::uint64_t>(float *, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<4, double, std::uint64_t>(double *,  const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<4, float, std::uint32_t>(float *,  const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_def_cuda<4, double, std::uint32_t>(double *,  const std::uint32_t * const, const Index);

/// \endcond
