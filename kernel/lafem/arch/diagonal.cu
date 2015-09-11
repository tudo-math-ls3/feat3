// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/diagonal.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAST
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


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_, typename IT_>
void Diagonal<Mem::CUDA>::csr(DT_ * diag, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_diagonal_csr<<<grid, block>>>(diag, val, col_ind, row_ptr, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void Diagonal<Mem::CUDA>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const Index);
template void Diagonal<Mem::CUDA>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const Index);
template void Diagonal<Mem::CUDA>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const Index);

template <typename DT_, typename IT_>
void Diagonal<Mem::CUDA>::ell(DT_ * diag, const DT_ * const val, const IT_ * const col_ind,
    const IT_ * const cs, const IT_ * const cl, const Index C, const Index rows)
{
  Index blocksize = MemoryPool<Mem::CUDA>::blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_diagonal_ell<<<grid, block>>>(diag, val, col_ind, cs, cl, C, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
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
