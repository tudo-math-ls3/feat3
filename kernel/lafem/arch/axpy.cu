// includes, FEAST
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/component_product.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/sum.hpp>

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
      __global__ void cuda_axpyv(DT_ * r, const DT_ * a, const DT_ * x, const DT_ * y, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        r[idx] = a[idx] * x[idx] + y[idx];
      }

      template <typename DT_>
      __global__ void cuda_axpy_mv_csr(DT_ * r, const DT_ a, const DT_ * x, const DT_ * y, const DT_ * val, const Index * col_ind,
          const Index * row_ptr, const Index count)
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

      template <typename DT_>
      __global__ void cuda_axpy_mv_v_csr(DT_ * r, const DT_ * a, const DT_ * x, const DT_ * y, const DT_ * val, const Index * col_ind,
          const Index * row_ptr, const Index count)
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
        r[idx] = (sum * a[idx]) + y[idx];
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_axpy_mv_ell(DT_ * r, const DT_ a, const DT_ * x, const DT_ * y, const DT_ * Ax, const IT_ * Aj,
          const IT_ * Arl, const Index stride, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index row(idx);
        const IT_ * tAj(Aj);
        const DT_ * tAx(Ax);
        DT_ sum(0);
        tAj += row;
        tAx += row;

        const Index max(Arl[row]);
        for(Index n(0); n < max ; n++)
        {
          const DT_ A_ij = *tAx;

          const Index col = *tAj;
          sum += A_ij * x[col];

          tAj += stride;
          tAx += stride;
        }
        r[row] = (sum * a) + y[idx];
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_axpy_mv_v_ell(DT_ * r, const DT_ * a, const DT_ * x, const DT_ * y, const DT_ * Ax, const IT_ * Aj,
          const IT_ * Arl, const Index stride, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index row(idx);
        const IT_ * tAj(Aj);
        const DT_ * tAx(Ax);
        DT_ sum(0);
        tAj += row;
        tAx += row;

        const Index max(Arl[row]);
        for(Index n(0); n < max ; n++)
        {
          const DT_ A_ij = *tAx;

          const Index col = *tAj;
          sum += A_ij * x[col];

          tAj += stride;
          tAx += stride;
        }
        r[row] = (sum * a[idx]) + y[idx];
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Axpy<Mem::CUDA, Algo::CUDA>::dv(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy<<<grid, block>>>(r, a, x, y, size);
}

template void Axpy<Mem::CUDA, Algo::CUDA>::dv(float *, const float, const float * const, const float * const, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::dv(double *, const double, const double * const, const double * const, const Index);

template <typename DT_>
void Axpy<Mem::CUDA, Algo::CUDA>::dv(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const Index size)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((size)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpyv<<<grid, block>>>(r, a, x, y, size);
}
template void Axpy<Mem::CUDA, Algo::CUDA>::dv(float *, const float * const, const float * const, const float * const, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::dv(double *, const double * const, const double * const, const double * const, const Index);

template <typename DT_>
void Axpy<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows, const Index columns, const Index used_elements)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy_mv_csr<<<grid, block>>>(r, a, x, y, val, col_ind, row_ptr, rows);
}
template void Axpy<Mem::CUDA, Algo::CUDA>::csr(float *, const float, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::csr(double *, const double, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);

template <typename DT_>
void Axpy<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows, const Index columns, const Index used_elements)
{
  FEAST::LAFEM::Arch::ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(r, val, col_ind, row_ptr, x, rows, columns, used_elements);
  FEAST::LAFEM::Arch::Scale<Mem::CUDA, Algo::CUDA>::value(r, r, a, rows);
  FEAST::LAFEM::Arch::Sum<Mem::CUDA, Algo::CUDA>::value(r, r, y, rows);
}
template void Axpy<Mem::CUDA, Algo::CUDA>::csr(float *, const float, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::csr(double *, const double, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);

template <typename DT_>
void Axpy<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const Index rows, const Index columns, const Index used_elements)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy_mv_v_csr<<<grid, block>>>(r, a, x, y, val, col_ind, row_ptr, rows);
}
template void Axpy<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index, const Index);

template <typename DT_>
void Axpy<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const Index rows, const Index columns, const Index used_elements)
{
  FEAST::LAFEM::Arch::ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(r, val, col_ind, row_ptr, x, rows, columns, used_elements);
  FEAST::LAFEM::Arch::ComponentProduct<Mem::CUDA, Algo::CUDA>::value(r, r, a, rows);
  FEAST::LAFEM::Arch::Sum<Mem::CUDA, Algo::CUDA>::value(r, r, y, rows);
}
template void Axpy<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void Axpy<Mem::CUDA, Algo::CUDA>::ell(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const Index stride, const Index rows)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy_mv_ell<<<grid, block>>>(r, a, x, y, Ax, Aj, Arl, stride, rows);
}

template void Axpy<Mem::CUDA, Algo::CUDA>::ell(float *, const float, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::ell(double *, const double, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::ell(float *, const float, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::ell(double *, const double, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index);

template <typename DT_, typename IT_>
void Axpy<Mem::CUDA, Algo::CUDA>::ell(DT_ * r, const DT_ * const a, const DT_ * const x, const DT_ * const y, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const Index stride, const Index rows)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_axpy_mv_v_ell<<<grid, block>>>(r, a, x, y, Ax, Aj, Arl, stride, rows);
}

template void Axpy<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const float * const, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const double * const, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const float * const, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
template void Axpy<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const double * const, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const Index, const Index);
