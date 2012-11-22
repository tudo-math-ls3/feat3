// includes, FEAST
#include <kernel/lafem/product.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_product_csr(DT_ * r, const DT_ * b, const DT_ * val, const Index * col_ind,
          const Index * row_ptr, const Index * row_ptr_end, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ sum(0);
        const Index end(row_ptr_end[idx]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          sum += val[i] * b[col_ind[i]];
        }
        r[idx] = sum;
      }

      template <typename DT_>
      __global__ void cuda_product_ell(DT_ * r, const DT_ * b, const DT_ * Ax, const Index * Aj,
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
        r[row] = sum;
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;

template <typename DT_>
void Product<Algo::CUDA>::value(DenseVector<Mem::CUDA, DT_> & r, const SparseMatrixCSR<Mem::CUDA, DT_> & a, const DenseVector<Mem::CUDA, DT_> & b)
{
  if (b.size() != a.columns())
    throw InternalError("Vector size does not match!");
  if (a.rows() != r.size())
    throw InternalError("Vector size does not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.size())/(double)(block.x));

  DT_ * r_gpu(r.elements());
  const DT_ * b_gpu(b.elements());
  const DT_ * val_gpu(a.val());
  const Index * col_ind_gpu(a.col_ind());
  const Index * row_ptr_gpu(a.row_ptr());
  const Index * row_ptr_end_gpu(a.row_ptr_end());

  FEAST::LAFEM::Intern::cuda_product_csr<<<grid, block>>>(r_gpu, b_gpu, val_gpu, col_ind_gpu, row_ptr_gpu, row_ptr_end_gpu, r.size());
}

template void Product<Algo::CUDA>::value(DenseVector<Mem::CUDA, float> &, const SparseMatrixCSR<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &);
template void Product<Algo::CUDA>::value(DenseVector<Mem::CUDA, double> &, const SparseMatrixCSR<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &);

template <typename DT_>
void Product<Algo::CUDA>::value(DenseVector<Mem::CUDA, DT_> & r, const SparseMatrixELL<Mem::CUDA, DT_> & a, const DenseVector<Mem::CUDA, DT_> & b)
{
  if (b.size() != a.columns())
    throw InternalError("Vector size does not match!");
  if (a.rows() != r.size())
    throw InternalError("Vector size does not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.size())/(double)(block.x));

  DT_ * r_gpu(r.elements());
  const DT_ * b_gpu(b.elements());
  const DT_ * Ax_gpu(a.Ax());
  const Index * Aj_gpu(a.Aj());
  const Index * Arl_gpu(a.Arl());

  FEAST::LAFEM::Intern::cuda_product_ell<<<grid, block>>>(r_gpu, b_gpu, Ax_gpu, Aj_gpu, Arl_gpu, a.stride(), r.size());
}

template void Product<Algo::CUDA>::value(DenseVector<Mem::CUDA, float> &, const SparseMatrixELL<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &);
template void Product<Algo::CUDA>::value(DenseVector<Mem::CUDA, double> &, const SparseMatrixELL<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &);
