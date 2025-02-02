// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_unit_filter_rhs(DT_ * v, const DT_ * sv_elements, const IT_ * sv_indices, const Index ue)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          v[sv_indices[i]] = sv_elements[i];
        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_unit_filter_def(DT_ * v, const IT_ * sv_indices, const Index ue)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          v[sv_indices[i]] = DT_(0);
        }
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_unit_filter_mat(DT_* __restrict__ mat, const IT_* __restrict__ const row_ptr, const IT_* __restrict__ const col_idx,
                                          const IT_ * __restrict__ const sv_indices, const Index ue)
      {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          const IT_ ix(sv_indices[i]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            mat[j] = (col_idx[j] == ix) ? DT_(1) : DT_(0);
          }
        }
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_unit_filter_offdiag_row_mat(DT_* __restrict__ mat, const IT_* __restrict__ const row_ptr, int block_width,
                                                  const IT_ * __restrict__ const sv_indices, const Index ue)
      {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          const IT_ ix(sv_indices[i]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
              for(int l(0); l < block_width; ++l)
                mat[j*block_width +l] = DT_(0);
          }
        }
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_unit_filter_weak_matrix_rows(DT_ * __restrict__ mat_a, const DT_ * __restrict__ const mat_m, const IT_* __restrict__ const row_ptr,
                                                   const DT_ * __restrict__ const sv_elements, const IT_ * __restrict__ const sv_indices, const Index ue)
      {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        // grid strided for loop
        for(int i = idx; idx < ue; idx += blockDim.x * gridDim.x)
        {
          const IT_ ix(sv_indices[i]);
          const DT_ vx(sv_elements[i]);

          // replace by unit row
          for(IT_ j(row_ptr[ix]); j < row_ptr[ix + 1]; ++j)
          {
            mat_a[j] = vx * mat_m[j];
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
void UnitFilter::filter_rhs_cuda(DT_ * v, const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_rhs<<<grid, block>>>(v, sv_elements, sv_indices, ue);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilter::filter_rhs_cuda(float *, const float * const, const std::uint64_t * const, const Index);
template void UnitFilter::filter_rhs_cuda(double *, const double * const, const std::uint64_t * const, const Index);
template void UnitFilter::filter_rhs_cuda(float *, const float * const, const std::uint32_t * const, const Index);
template void UnitFilter::filter_rhs_cuda(double *, const double * const, const std::uint32_t * const, const Index);

template <typename DT_, typename IT_>
void UnitFilter::filter_def_cuda(DT_ * v, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_def<<<grid, block>>>(v, sv_indices, ue);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilter::filter_def_cuda(float *, const std::uint64_t * const, const Index);
template void UnitFilter::filter_def_cuda(double *, const std::uint64_t * const, const Index);
template void UnitFilter::filter_def_cuda(float *, const std::uint32_t * const, const Index);
template void UnitFilter::filter_def_cuda(double *, const std::uint32_t * const, const Index);

template<typename DT_, typename IT_>
void UnitFilter::filter_unit_mat_cuda(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_mat<DT_, IT_><<<grid, block >>>(mat, row_ptr, col_idx, sv_indices, ue);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilter::filter_unit_mat_cuda<float, std::uint64_t>(float *, const std::uint64_t * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void UnitFilter::filter_unit_mat_cuda<float, std::uint32_t>(float *, const std::uint32_t * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
template void UnitFilter::filter_unit_mat_cuda<double, std::uint64_t>(double *, const std::uint64_t * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void UnitFilter::filter_unit_mat_cuda<double, std::uint32_t>(double *, const std::uint32_t * const, const std::uint32_t * const, const std::uint32_t * const, const Index);

template<typename DT_, typename IT_>
void UnitFilter::filter_offdiag_row_mat_cuda(DT_* mat, const IT_* const row_ptr, int block_width, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_offdiag_row_mat<DT_, IT_><<<grid, block >>>(mat, row_ptr, block_width, sv_indices, ue);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilter::filter_offdiag_row_mat_cuda<float, std::uint64_t>(float *, const std::uint64_t * const, int, const std::uint64_t * const, const Index);
template void UnitFilter::filter_offdiag_row_mat_cuda<float, std::uint32_t>(float *, const std::uint32_t * const, int, const std::uint32_t * const, const Index);
template void UnitFilter::filter_offdiag_row_mat_cuda<double, std::uint64_t>(double *, const std::uint64_t * const, int, const std::uint64_t * const, const Index);
template void UnitFilter::filter_offdiag_row_mat_cuda<double, std::uint32_t>(double *, const std::uint32_t * const, int, const std::uint32_t * const, const Index);

template<typename DT_, typename IT_>
void UnitFilter::filter_weak_matrix_rows_cuda(DT_* mat_a, const DT_ * const mat_m, const IT_* const row_ptr,
                                                        const DT_ * const sv_elements, const IT_ * const sv_indices, const Index ue)
{
  Index blocksize = Util::cuda_blocksize_misc;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((ue)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_unit_filter_weak_matrix_rows<DT_, IT_><<<grid, block >>>(mat_a, mat_m, row_ptr, sv_elements, sv_indices, ue);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void UnitFilter::filter_weak_matrix_rows_cuda<float, std::uint64_t>(float *, const float* const, const std::uint64_t * const, const float * const, const std::uint64_t * const, const Index);
template void UnitFilter::filter_weak_matrix_rows_cuda<float, std::uint32_t>(float *, const float* const, const std::uint32_t * const, const float * const, const std::uint32_t * const, const Index);
template void UnitFilter::filter_weak_matrix_rows_cuda<double, std::uint64_t>(double *, const double* const, const std::uint64_t * const, const double * const, const std::uint64_t * const, const Index);
template void UnitFilter::filter_weak_matrix_rows_cuda<double, std::uint32_t>(double *, const double* const, const std::uint32_t * const, const double * const, const std::uint32_t * const, const Index);
/// \endcond
