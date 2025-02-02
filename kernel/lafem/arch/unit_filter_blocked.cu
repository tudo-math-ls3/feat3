// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_blocked.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/cuda_math.cuh>

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_unit_filter_blocked_rhs(DT_ * __restrict__ v, int block_size, const DT_ * __restrict__ sv_elements, const IT_ * __restrict__ sv_indices,
                                                  const Index ue, bool ign_nans)
      {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          if(ign_nans)
          {
            for(Index j(0) ; j < block_size; ++j)
            {
              if(!isnan(sv_elements[block_size * idx + j]))
                v[block_size* sv_indices[i] + j] = sv_elements[block_size * i + j];
            }
          }
          else
          {
            for(Index j(0) ; j < block_size; ++j)
              v[block_size* sv_indices[i] + j] = sv_elements[block_size * i + j];
          }
        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_unit_filter_blocked_def(DT_ * __restrict__ v, int block_size, const DT_ * __restrict__ sv_elements, const IT_ * __restrict__ sv_indices,
                                                    const Index ue, bool ign_nans)
      {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          if(ign_nans)
          {
            for(Index j(0) ; j < block_size; ++j)
            {
              if(!isnan(sv_elements[block_size * i + j]))
                v[block_size* sv_indices[i] + j] = DT_(0);
            }
          }
          else
          {
            for(Index j(0) ; j < block_size; ++j)
              v[block_size * sv_indices[i] + j] = DT_(0);
          }
        }
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_unit_filter_blocked_mat(DT_* __restrict__ mat, const IT_* __restrict__ const row_ptr, const IT_* __restrict__ const col_idx, int block_height,
                                          int block_width, const DT_ * __restrict__ const sv_elements, const IT_ * __restrict__ const sv_indices, const Index ue, bool ign_nans)
      {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          const IT_ ix(sv_indices[i]);
          const DT_* vx(&sv_elements[i*block_height]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            // loop over rows in the block
            for(int k(0); k < block_height; ++k)
            {
              // possibly skip row if filter value is NaN
              if(ign_nans && CudaMath::cuda_isnan(vx[k]))
                continue;
              for(int l(0); l < block_width; ++l)
                mat[j*block_height*block_width + k*block_width +l] = DT_(0);
              if((col_idx[j] == ix) && (k < block_width))
                mat[j*block_height*block_width + k*block_width +k] = DT_(1);
            }
          }
        }
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_unit_filter_blocked_offdiag_row_mat(DT_* __restrict__ mat, const IT_* __restrict__ const row_ptr, int block_height, int block_width,
                                                  const DT_ * __restrict__ const sv_elements, const IT_ * __restrict__ const sv_indices, const Index ue, bool ign_nans)
      {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          const IT_ ix(sv_indices[i]);
          const DT_* vx(&sv_elements[i*block_height]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            // loop over rows in the block
            for(int k(0); k < block_height; ++k)
            {
              // possibly skip row if filter value is NaN
              if(ign_nans && CudaMath::cuda_isnan(vx[k]))
                continue;
              for(int l(0); l < block_width; ++l)
                mat[j*block_height*block_width + k*block_width +l] = DT_(0);
            }
          }
        }
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_unit_filter_blocked_weak_matrix_rows(DT_ * __restrict__ mat_a, const DT_ * __restrict__ const mat_m, const IT_* __restrict__ const row_ptr, int block_height,
                                                   int block_width, const DT_ * __restrict__ const sv_elements, const IT_ * __restrict__ const sv_indices, const Index ue)
      {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          const IT_ ix(sv_indices[i]);
          const DT_* vx(&sv_elements[i*block_height]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            // loop over rows in the block
            for(int k(0); k < block_height; ++k)
            {
              for(int l(0); l < block_width; ++l)
              {
                mat_a[j*block_height*block_width + k*block_width + l] = vx[k] * mat_m[j*block_height*block_width + k*block_width + l];
              }
            }
          }
        }
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_, typename IT_>
void UnitFilterBlocked::filter_rhs_cuda(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_blocked_rhs<DT_, IT_><<<grid, block>>>(v, block_size, sv_elements, sv_indices, ue, ign_nans);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilterBlocked::filter_rhs_cuda<float, std::uint64_t>(float *, int, const float * const, const std::uint64_t * const, const Index, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_cuda<double, std::uint64_t>(double *, int, const double * const, const std::uint64_t * const, const Index, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_cuda<float, std::uint32_t>(float *, int, const float * const, const std::uint32_t * const, const Index, bool ign_nans);
template void UnitFilterBlocked::filter_rhs_cuda<double, std::uint32_t>(double *, int, const double * const, const std::uint32_t * const, const Index, bool ign_nans);

template <typename DT_, typename IT_>
void UnitFilterBlocked::filter_def_cuda(DT_ * v, int block_size, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_blocked_def<DT_, IT_><<<grid, block>>>(v, block_size, sv_elements, sv_indices, ue, ign_nans);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilterBlocked::filter_def_cuda<float, std::uint64_t>(float *, int, const float* const, const std::uint64_t * const, const Index, bool ign_nans);
template void UnitFilterBlocked::filter_def_cuda<double, std::uint64_t>(double *, int, const double* const, const std::uint64_t * const, const Index, bool ign_nans);
template void UnitFilterBlocked::filter_def_cuda<float, std::uint32_t>(float *, int, const float* const, const std::uint32_t * const, const Index, bool ign_nans);
template void UnitFilterBlocked::filter_def_cuda<double, std::uint32_t>(double *, int, const double* const, const std::uint32_t * const, const Index, bool ign_nans);

template<typename DT_, typename IT_>
void UnitFilterBlocked::filter_unit_mat_cuda(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_blocked_mat<DT_, IT_><<<grid, block >>>(mat, row_ptr, col_idx, block_height, block_width, sv_elements, sv_indices, ue, ign_nans);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilterBlocked::filter_unit_mat_cuda<float, std::uint64_t>(float *, const std::uint64_t * const, const std::uint64_t * const, int, int, const float* const, const std::uint64_t * const, const Index, bool);
template void UnitFilterBlocked::filter_unit_mat_cuda<float, std::uint32_t>(float *, const std::uint32_t * const, const std::uint32_t * const, int, int, const float* const, const std::uint32_t * const, const Index, bool);
template void UnitFilterBlocked::filter_unit_mat_cuda<double, std::uint64_t>(double *, const std::uint64_t * const, const std::uint64_t * const, int, int, const double* const, const std::uint64_t * const, const Index, bool);
template void UnitFilterBlocked::filter_unit_mat_cuda<double, std::uint32_t>(double *, const std::uint32_t * const, const std::uint32_t * const, int, int, const double* const, const std::uint32_t * const, const Index, bool);

template<typename DT_, typename IT_>
void UnitFilterBlocked::filter_offdiag_row_mat_cuda(DT_* mat, const IT_* const row_ptr, int block_height, int block_width, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue, bool ign_nans)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_blocked_offdiag_row_mat<DT_, IT_><<<grid, block >>>(mat, row_ptr, block_height, block_width, sv_elements, sv_indices, ue, ign_nans);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilterBlocked::filter_offdiag_row_mat_cuda<float, std::uint64_t>(float *, const std::uint64_t * const, int, int, const float* const, const std::uint64_t * const, const Index, bool);
template void UnitFilterBlocked::filter_offdiag_row_mat_cuda<float, std::uint32_t>(float *, const std::uint32_t * const, int, int, const float* const, const std::uint32_t * const, const Index, bool);
template void UnitFilterBlocked::filter_offdiag_row_mat_cuda<double, std::uint64_t>(double *, const std::uint64_t * const, int, int, const double* const, const std::uint64_t * const, const Index, bool);
template void UnitFilterBlocked::filter_offdiag_row_mat_cuda<double, std::uint32_t>(double *, const std::uint32_t * const, int, int, const double* const, const std::uint32_t * const, const Index, bool);

template<typename DT_, typename IT_>
void UnitFilterBlocked::filter_weak_matrix_rows_cuda(DT_* mat_a, const DT_ * const mat_m, const IT_* const row_ptr, int block_height, int block_width,
                                                        const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_blocked_weak_matrix_rows<DT_, IT_><<<grid, block >>>(mat_a, mat_m, row_ptr, block_height, block_width, sv_elements, sv_indices, ue);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilterBlocked::filter_weak_matrix_rows_cuda<float, std::uint64_t>(float *, const float* const, const std::uint64_t * const, int, int, const float* const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_weak_matrix_rows_cuda<float, std::uint32_t>(float *, const float* const, const std::uint32_t * const, int, int, const float* const, const std::uint32_t * const, const Index);
template void UnitFilterBlocked::filter_weak_matrix_rows_cuda<double, std::uint64_t>(double *, const double* const, const std::uint64_t * const, int, int, const double* const, const std::uint64_t * const, const Index);
template void UnitFilterBlocked::filter_weak_matrix_rows_cuda<double, std::uint32_t>(double *, const double* const, const std::uint32_t * const, int, int, const double* const, const std::uint32_t * const, const Index);
/// \endcond
