// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/scale_row_col.hpp>
#include <kernel/util/exception.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_, typename IT_>
      __global__ void cuda_scale_rows_csr(DT_ * r, const DT_ * b, const DT_ * val, const IT_ * col_ind,
          const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          r[i] = val[i] * b[idx];
        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_scale_cols_csr(DT_ * r, const DT_ * b, const DT_ * val, const IT_ * col_ind,
          const IT_ * row_ptr, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ sum(0);
        const Index end(row_ptr[idx + 1]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          sum += val[i] * b[col_ind[i]];
        }
      }


      template <typename DT_, typename IT_>
      __global__ void cuda_scale_rows_ell(DT_ * r, const DT_ * b, const DT_ * Ax, const IT_ * Aj,
          const IT_ * Arl, const Index stride, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index row(idx);
        const DT_ * tAx(Ax);
        DT_ * tr(r);
        tAx += row;
        tr += row;

        const Index max(Arl[row]);
        for(Index n(0); n < max ; n++)
        {
          *tr = *Ax * b[row];

          tAx += stride;
          tr += stride;
        }
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_scale_cols_ell(DT_ * r, const DT_ * b, const DT_ * Ax, const IT_ * Aj,
          const IT_ * Arl, const Index stride, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index row(idx);
        const IT_ * tAj(Aj);
        const DT_ * tAx(Ax);
        DT_ * tr(r);
        tAj += row;
        tAx += row;
        tr += row;

        const Index max(Arl[row]);
        for(Index n(0); n < max ; n++)
        {
          *tr = *Ax * b[*tAj];

          tAj += stride;
          tAx += stride;
          tr += stride;
        }
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_, typename IT_>
void ScaleRows<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_scale_rows_csr<<<grid, block>>>(r, x, val, col_ind, row_ptr, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ScaleRows<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
template void ScaleRows<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
template void ScaleRows<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
template void ScaleRows<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);

template <typename DT_, typename IT_>
void ScaleCols<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_scale_cols_csr<<<grid, block>>>(r, x, val, col_ind, row_ptr, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ScaleCols<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
template void ScaleCols<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);
template void ScaleCols<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
template void ScaleCols<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);


template <typename DT_, typename IT_>
void ScaleRows<Mem::CUDA, Algo::CUDA>::ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_scale_rows_ell<<<grid, block>>>(r, x, Ax, Aj, Arl, stride, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ScaleRows<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
template void ScaleRows<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
template void ScaleRows<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
template void ScaleRows<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);

template <typename DT_, typename IT_>
void ScaleCols<Mem::CUDA, Algo::CUDA>::ell(DT_ * r, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_scale_cols_ell<<<grid, block>>>(r, x, Ax, Aj, Arl, stride, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void ScaleCols<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
template void ScaleCols<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
template void ScaleCols<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
template void ScaleCols<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);
