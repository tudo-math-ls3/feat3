// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/mirror.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template<typename DT_, typename IT_>
      __global__ void cuda_mirror_gather_dv(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
      {
        Index ii = threadIdx.x + blockDim.x * blockIdx.x;
        if (ii >= nidx)
          return;

        buf[boff+ii] = vec[idx[ii]];
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_mirror_scatter_dv(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
      {
        Index ii = threadIdx.x + blockDim.x * blockIdx.x;
        if (ii >= nidx)
          return;

        vec[idx[ii]] += alpha*buf[boff+ii];
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_mirror_gather_dvb(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
      {
        Index ii = threadIdx.x + blockDim.x * blockIdx.x;
        if (ii >= nidx)
          return;

        for(Index k(0); k < bs; ++k)
        {
          buf[boff+ii*bs+k] = vec[idx[ii]*bs+k];
        }
      }

      template<typename DT_, typename IT_>
      __global__ void cuda_mirror_scatter_dvb(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
      {
        Index ii = threadIdx.x + blockDim.x * blockIdx.x;
        if (ii >= nidx)
          return;

        for(Index k(0); k < bs; ++k)
        {
          vec[idx[ii]*bs+k] += alpha*buf[boff+ii*bs+k];
        }
      }
    } // namespace Intern

    namespace Arch
    {
      template<typename DT_, typename IT_>
      void Mirror<Mem::CUDA>::gather_dv(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
      {
        Index blocksize = MemoryPool<Mem::CUDA>::blocksize_misc;
        dim3 grid;
        dim3 block;
        block.x = blocksize;
        grid.x = (unsigned)ceil((nidx)/(double)(block.x));

        FEAT::LAFEM::Intern::cuda_mirror_gather_dv<<<grid, block>>>(boff, nidx, idx, buf, vec);

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
      }

      template void Mirror<Mem::CUDA>::gather_dv(const Index, const Index, const unsigned long*, float*, const float*);
      template void Mirror<Mem::CUDA>::gather_dv(const Index, const Index, const unsigned long*, double*, const double*);
      template void Mirror<Mem::CUDA>::gather_dv(const Index, const Index, const unsigned int*, float*, const float*);
      template void Mirror<Mem::CUDA>::gather_dv(const Index, const Index, const unsigned int*, double*, const double*);

      template<typename DT_, typename IT_>
      void Mirror<Mem::CUDA>::scatter_dv(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
      {
        Index blocksize = MemoryPool<Mem::CUDA>::blocksize_misc;
        dim3 grid;
        dim3 block;
        block.x = blocksize;
        grid.x = (unsigned)ceil((nidx)/(double)(block.x));

        FEAT::LAFEM::Intern::cuda_mirror_scatter_dv<<<grid, block>>>(boff, nidx, idx, buf, vec, alpha);

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
      }

      template void Mirror<Mem::CUDA>::scatter_dv(const Index, const Index, const unsigned long*, const float*, float*, float);
      template void Mirror<Mem::CUDA>::scatter_dv(const Index, const Index, const unsigned long*, const double*, double*, double);
      template void Mirror<Mem::CUDA>::scatter_dv(const Index, const Index, const unsigned int*, const float*, float*, float);
      template void Mirror<Mem::CUDA>::scatter_dv(const Index, const Index, const unsigned int*, const double*, double*, double);

      template<typename DT_, typename IT_>
      void Mirror<Mem::CUDA>::gather_dvb(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
      {
        Index blocksize = MemoryPool<Mem::CUDA>::blocksize_misc;
        dim3 grid;
        dim3 block;
        block.x = blocksize;
        grid.x = (unsigned)ceil((nidx)/(double)(block.x));

        FEAT::LAFEM::Intern::cuda_mirror_gather_dvb<<<grid, block>>>(bs, boff, nidx, idx, buf, vec);

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
      }

      template void Mirror<Mem::CUDA>::gather_dvb(const Index, const Index, const Index, const unsigned long*, float*, const float*);
      template void Mirror<Mem::CUDA>::gather_dvb(const Index, const Index, const Index, const unsigned long*, double*, const double*);
      template void Mirror<Mem::CUDA>::gather_dvb(const Index, const Index, const Index, const unsigned int*, float*, const float*);
      template void Mirror<Mem::CUDA>::gather_dvb(const Index, const Index, const Index, const unsigned int*, double*, const double*);

      template<typename DT_, typename IT_>
      void Mirror<Mem::CUDA>::scatter_dvb(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
      {
        Index blocksize = MemoryPool<Mem::CUDA>::blocksize_misc;
        dim3 grid;
        dim3 block;
        block.x = blocksize;
        grid.x = (unsigned)ceil((nidx)/(double)(block.x));

        FEAT::LAFEM::Intern::cuda_mirror_scatter_dvb<<<grid, block>>>(bs, boff, nidx, idx, buf, vec, alpha);

#ifdef FEAT_DEBUG_MODE
        cudaDeviceSynchronize();
        cudaError_t last_error(cudaGetLastError());
        if (cudaSuccess != last_error)
          throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
      }

      template void Mirror<Mem::CUDA>::scatter_dvb(const Index, const Index, const Index, const unsigned long*, const float*, float*, float);
      template void Mirror<Mem::CUDA>::scatter_dvb(const Index, const Index, const Index, const unsigned long*, const double*, double*, double);
      template void Mirror<Mem::CUDA>::scatter_dvb(const Index, const Index, const Index, const unsigned int*, const float*, float*, float);
      template void Mirror<Mem::CUDA>::scatter_dvb(const Index, const Index, const Index, const unsigned int*, const double*, double*, double);
    }
  } // namespace LAFEM
} // namespace FEAT

/// \endcond
