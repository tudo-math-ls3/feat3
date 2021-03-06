// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/lumping.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_lumping_csr(DT_ * lump, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ sum(0);
        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          sum += val[i];
        }
        lump[idx] = sum;
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_lumping_ell(DT_ * lump, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * cs, const IT_ * cl, const Index C, const Index rows)
      {
        const Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows)
          return;


        DT_ sum(0);
        const Index chunk(idx / C);
        const Index local_row(idx % C);
        const Index chunk_end(cs[chunk+1]);

        for (Index pcol(cs[chunk] + local_row) ; pcol < chunk_end ; pcol+=C)
        {
          sum += val[pcol];
        }
        lump[idx] = sum;
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_lumping_bcsr(DT_ * lump, const DT_ * val, const IT_ * col_ind,
                                              const IT_ * row_ptr, const Index rows,
                                              const int BlockHeight, const int BlockWidth)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows * BlockHeight)
          return;

        Index csr_row = idx / BlockHeight;
        Index block_row = idx % BlockHeight;

        DT_ sum(0);
        const Index end(row_ptr[csr_row + 1]);
        for (Index i(row_ptr[csr_row]) ; i < end ; ++i)
        {
          for (Index w(0) ; w < BlockWidth ; ++w)
          {
            sum += val[BlockHeight*BlockWidth*i + block_row*BlockWidth + w];
          }
        }
        lump[idx] = sum;
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename DT_, typename IT_>
void Lumping<Mem::CUDA>::csr(DT_ * lump, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_lumping_csr<<<grid, block>>>(lump, val, col_ind, row_ptr, rows);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Lumping<Mem::CUDA>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const Index);
template void Lumping<Mem::CUDA>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const Index);
template void Lumping<Mem::CUDA>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const Index);
template void Lumping<Mem::CUDA>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const Index);


template <typename DT_, typename IT_>
void Lumping<Mem::CUDA>::ell(DT_ * lump, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs, const IT_ * const cl, const Index C, const Index rows)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_lumping_ell<<<grid, block>>>(lump, val, col_ind, cs, cl, C, rows);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Lumping<Mem::CUDA>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
template void Lumping<Mem::CUDA>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
template void Lumping<Mem::CUDA>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Lumping<Mem::CUDA>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const Index, const Index);

template <typename DT_, typename IT_>
void Lumping<Mem::CUDA>::bcsr(DT_ * lump, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows * BlockHeight)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_lumping_bcsr<<<grid, block>>>(lump, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);
#ifdef FEAT_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Lumping<Mem::CUDA>::bcsr(float *, const float * const, const unsigned long * const, const unsigned long * const, const Index, const int, const int);
template void Lumping<Mem::CUDA>::bcsr(double *, const double * const, const unsigned long * const, const unsigned long * const, const Index, const int, const int);
template void Lumping<Mem::CUDA>::bcsr(float *, const float * const, const unsigned int * const, const unsigned int * const, const Index, const int, const int);
template void Lumping<Mem::CUDA>::bcsr(double *, const double * const, const unsigned int * const, const unsigned int * const, const Index, const int, const int);
