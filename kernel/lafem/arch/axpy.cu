// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/component_product.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/memory_pool.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_axpy(DT_ * r, const DT_ a, const DT_ * x, const DT_ * y, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        r[idx] = a * x[idx] + y[idx];
      }

      template <typename DT_>
      __global__ void cuda_axpy_mv_csr(DT_ * r, const DT_ a, const DT_ * x, const DT_ * y, const DT_ * val,
          const unsigned long * col_ind, const unsigned long * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ sum(0);
        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          sum += val[i] * x[col_ind[i]];
        }
        r[idx] = (sum * a) + y[idx];
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_axpy_mv_ell(DT_ * r, const DT_ a, const DT_ * x, const DT_ * y, const DT_ * val, const IT_ * col_ind,
                                       const IT_ * cs, const IT_ * cl, const Index rows, const Index C)
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
          sum += val[pcol] * x[col_ind[pcol]];
        }
        r[idx] = sum * a + y[idx];

      }

      template <typename DT_, typename IT_>
      __global__ void cuda_axpy_banded(DT_ * r, const DT_ alpha, const DT_ * x, const DT_ * y, const DT_ * val, const IT_ * offsets, const Index num_of_offsets, const Index rows, const Index columns)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= rows)
          return;

        const Index k1(rows - 1);
        const Index k2(rows + columns - 1);

        Index start(0);

        while (k1 > offsets[start] + idx)
        {
          ++start;
        }

        Index end(start);

        while (end < num_of_offsets && idx + offsets[end] < k2)
        {
          ++end;
        }

        DT_ sum(DT_(0.0));
        for (Index diag(start); diag < end; ++diag)
        {
          sum += val[rows * diag + idx] * x[idx + offsets[diag] - rows + 1];
        }
        r[idx] = (sum*alpha) + y[idx];
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Axpy<Mem::CUDA>::dv(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
{
  Index blocksize = Util::MemoryPool<Mem::CUDA>::instance()->blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy<<<grid, block>>>(r, a, x, y, size);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}

template void Axpy<Mem::CUDA>::dv(float *, const float, const float * const, const float * const, const Index);
template void Axpy<Mem::CUDA>::dv(double *, const double, const double * const, const double * const, const Index);

template <typename DT_>
void Axpy<Mem::CUDA>::csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows, const Index columns, const Index used_elements)
{
  Index blocksize = Util::MemoryPool<Mem::CUDA>::instance()->blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy_mv_csr<<<grid, block>>>(r, a, x, y, val, col_ind, row_ptr, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Axpy<Mem::CUDA>::csr(float *, const float, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA>::csr(double *, const double, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);

template <typename DT_>
void Axpy<Mem::CUDA>::csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows, const Index columns, const Index used_elements)
{
  FEAST::LAFEM::Arch::ProductMatVec<Mem::CUDA>::csr(r, val, col_ind, row_ptr, x, rows, columns, used_elements);
  FEAST::LAFEM::Arch::Scale<Mem::CUDA>::value(r, r, a, rows);
  FEAST::LAFEM::Arch::Sum<Mem::CUDA>::value(r, r, y, rows);
}
template void Axpy<Mem::CUDA>::csr(float *, const float, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA>::csr(double *, const double, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void Axpy<Mem::CUDA>::ell(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_ind, const IT_ * const cs, const IT_ * const cl, const Index C, const Index rows)
{
  Index blocksize = Util::MemoryPool<Mem::CUDA>::instance()->blocksize_axpy;
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy_mv_ell<<<grid, block>>>(r, a, x, y, val, col_ind, cs, cl, rows, C);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Axpy<Mem::CUDA>::ell(float *, const float, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Axpy<Mem::CUDA>::ell(double *, const double, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Axpy<Mem::CUDA>::ell(float *, const float, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
template void Axpy<Mem::CUDA>::ell(double *, const double, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const Index, const Index);

template <typename DT_, typename IT_>
void Axpy<Mem::CUDA>::banded(DT_ * r, const DT_ * const y, const DT_ alpha, const DT_ * const val, const IT_ * const offsets, const DT_ * const x, const Index num_of_offsets, const Index rows, const Index columns)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy_banded<<<grid, block>>>(r, alpha, x, y, val, offsets, num_of_offsets, rows, columns);
}
template void Axpy<Mem::CUDA>::banded(float *, const float * const, const float, const float * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA>::banded(double *, const double * const, const double, const double * const, const unsigned int * const, const double * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA>::banded(float *, const float * const, const float, const float * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA>::banded(double *, const double * const, const double, const double * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
