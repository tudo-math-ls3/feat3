// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/diagonal.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_diagonal_csr(DT_ * diag, const DT_ * val, const IT_ * col_ind, const IT_ * row_ptr, const Index rows)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows)
          return;

        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]); i < end; ++i)
        {
          if (idx == col_ind[i])
          {
            diag[idx] = val[i];
            return;
          }
        }
        diag[idx] = 0;
      }

      template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      __global__ void cuda_diagonal_csrb(DT_ * diag, const DT_ * val, const IT_ * col_ind, const IT_ * row_ptr, const Index rows)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows)
          return;

        const Index row(idx);

        for (Index i(0) ; i < BlockHeight_ ; ++i)
        {
          diag[row * BlockHeight_ + i] = DT_(0);
        }
        const Index end = row_ptr[row + 1];

        for (Index col = row_ptr[row]; col < end; col++)
        {
          if (row == col_ind[col])
          {
            for (Index i(0) ; i < BlockHeight_ ; ++i)
            {
              diag[row * BlockHeight_ + i] = val[(col * BlockHeight_ * BlockWidth_) + i + i * BlockWidth_];
            }
          }
        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_diagonal_ell(DT_ * diag, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * cs, const IT_ * cl, const Index C, const Index rows)
      {
        const Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows)
          return;


        const Index chunk(idx / C);
        const Index local_row(idx % C);
        const Index chunk_end(cs[chunk+1]);

        for (Index pcol(cs[chunk] + local_row) ; pcol < chunk_end ; pcol+=C)
        {
          if (col_ind[pcol] == idx)
          {
            diag[idx] = val[pcol];
            return;
          }
        }
        diag[idx] = DT_(0);

      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_, typename IT_>
void Diagonal<Mem::CUDA>::csr(DT_ * diag, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_diagonal_csr<<<grid, block>>>(diag, val, col_ind, row_ptr, rows);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void Diagonal<Mem::CUDA>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const Index);
template void Diagonal<Mem::CUDA>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const Index);

template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
void Diagonal<Mem::CUDA>::csrb(DT_ * diag, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_diagonal_csrb<DT_, IT_, BlockHeight_, BlockWidth_><<<grid, block>>>(diag, val, col_ind, row_ptr, rows);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Diagonal<Mem::CUDA>::csrb<float, unsigned long, 1, 1>(float *, const float * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<double, unsigned long, 1, 1>(double *, const double * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<float, unsigned int, 1, 1>(float *, const float * const, const unsigned int * const, const unsigned int * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<double, unsigned int, 1, 1>(double *, const double * const, const unsigned int * const, const unsigned int * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<float, unsigned long, 2, 2>(float *, const float * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<double, unsigned long, 2, 2>(double *, const double * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<float, unsigned int, 2, 2>(float *, const float * const, const unsigned int * const, const unsigned int * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<double, unsigned int, 2, 2>(double *, const double * const, const unsigned int * const, const unsigned int * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<float, unsigned long, 3, 3>(float *, const float * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<double, unsigned long, 3, 3>(double *, const double * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<float, unsigned int, 3, 3>(float *, const float * const, const unsigned int * const, const unsigned int * const, const Index);
template void Diagonal<Mem::CUDA>::csrb<double, unsigned int, 3, 3>(double *, const double * const, const unsigned int * const, const unsigned int * const, const Index);

template <typename DT_, typename IT_>
void Diagonal<Mem::CUDA>::ell(DT_ * diag, const DT_ * const val, const IT_ * const col_ind,
    const IT_ * const cs, const IT_ * const cl, const Index C, const Index rows)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_diagonal_ell<<<grid, block>>>(diag, val, col_ind, cs, cl, C, rows);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void Diagonal<Mem::CUDA>::ell(float *, const float * const, const unsigned long * const,
    const unsigned long * const, const unsigned long * const, const Index, const Index);
template void Diagonal<Mem::CUDA>::ell(double *, const double * const, const unsigned long * const,
    const unsigned long * const, const unsigned long * const, const Index, const Index);
template void Diagonal<Mem::CUDA>::ell(float *, const float * const, const unsigned int * const,
    const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Diagonal<Mem::CUDA>::ell(double *, const double * const, const unsigned int * const,
    const unsigned int * const, const unsigned int * const, const Index, const Index);
