// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/product_matvec.hpp>
#include <kernel/util/exception.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_defect_csr(DT_ * r, const DT_ * rhs, const DT_ * b, const DT_ * val, const unsigned long * col_ind,
          const unsigned long * row_ptr, const Index count)
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
        r[idx] = rhs[idx] - sum;
      }

      template <typename DT_, typename IT_>
      __global__ void cuda_defect_ell(DT_ * r, const DT_ * rhs, const DT_ * b, const DT_ * Ax, const IT_ * Aj,
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
          sum += A_ij * b[col];

          tAj += stride;
          tAx += stride;
        }
        r[row] = rhs[row] - sum;
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
void Defect<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned long * const col_ind, const unsigned long * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_defect_csr<<<grid, block>>>(r, rhs, x, val, col_ind, row_ptr, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Defect<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index, const Index);
template void Defect<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index, const Index);

template <typename DT_>
void Defect<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const unsigned int * const col_ind, const unsigned int * const row_ptr, const DT_ * const x, const Index rows, const Index columns, const Index used_elements)
{
  FEAST::LAFEM::Arch::ProductMatVec<Mem::CUDA, Algo::CUDA>::csr(r, val, col_ind, row_ptr, x, rows, columns, used_elements);
  FEAST::LAFEM::Arch::Difference<Mem::CUDA, Algo::CUDA>::value(r, rhs, r, rows);
}
template void Defect<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index, const Index);
template void Defect<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index, const Index);


template <typename DT_, typename IT_>
void Defect<Mem::CUDA, Algo::CUDA>::ell(DT_ * r, const DT_ * const rhs, const DT_ * const Ax, const IT_ * const Aj, const IT_ * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_defect_ell<<<grid, block>>>(r, rhs, x, Ax, Aj, Arl, stride, rows);
#ifdef FEAST_DEBUG_MODE
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
#endif
}
template void Defect<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
template void Defect<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
template void Defect<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
template void Defect<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);
