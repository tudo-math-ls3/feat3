// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/diagonal.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/cuda_util.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename IT_>
      __global__ void cuda_diagonal_csr(IT_ * diag, const IT_ * col_ind, const IT_ * row_ptr, const Index rows)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows)
          return;

        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]); i < end; ++i)
        {
          if (idx == col_ind[i])
          {
            diag[idx] = i;
            return;
          }
        }
        diag[idx] = row_ptr[rows];
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template <typename IT_>
void Diagonal::csr_cuda(IT_ * diag, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
{
  Index blocksize = Util::cuda_blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = (unsigned)blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAT::LAFEM::Intern::cuda_diagonal_csr<<<grid, block>>>(diag, col_ind, row_ptr, rows);

  cudaDeviceSynchronize();
#ifdef FEAT_DEBUG_MODE
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void Diagonal::csr_cuda(std::uint64_t *, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void Diagonal::csr_cuda(std::uint32_t *, const std::uint32_t * const, const std::uint32_t * const, const Index);
