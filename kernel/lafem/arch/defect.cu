// includes, FEAST
#include <kernel/lafem/arch/defect.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_product_matvec_csr(DT_ * r, const DT_ * rhs, const DT_ * b, const DT_ * val, const Index * col_ind,
          const Index * row_ptr, const Index count)
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

      template <typename DT_>
      __global__ void cuda_product_matvec_ell(DT_ * r, const DT_ * rhs, const DT_ * b, const DT_ * Ax, const Index * Aj,
          const Index * Arl, const Index stride, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        const Index row(idx);
        const Index * tAj(Aj);
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
void Defect<Mem::CUDA, Algo::CUDA>::csr(DT_ * r, const DT_ * const rhs, const DT_ * const val, const Index * const col_ind, const Index * const row_ptr, const DT_ * const x, const Index rows)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_product_matvec_csr<<<grid, block>>>(r, rhs, x, val, col_ind, row_ptr, rows);
}
template void Defect<Mem::CUDA, Algo::CUDA>::csr(float *, const float * const, const float * const, const Index * const, const Index * const, const float * const, const Index);
template void Defect<Mem::CUDA, Algo::CUDA>::csr(double *, const double * const, const double * const, const Index * const, const Index * const, const double * const, const Index);


template <typename DT_>
void Defect<Mem::CUDA, Algo::CUDA>::ell(DT_ * r, const DT_ * const rhs, const DT_ * const Ax, const Index * const Aj, const Index * const Arl, const DT_ * const x, const Index stride, const Index rows)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((rows)/(double)(block.x));

  FEAST::LAFEM::Intern::cuda_product_matvec_ell<<<grid, block>>>(r, rhs, x, Ax, Aj, Arl, stride, rows);
}
template void Defect<Mem::CUDA, Algo::CUDA>::ell(float *, const float * const, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void Defect<Mem::CUDA, Algo::CUDA>::ell(double *, const double * const, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
