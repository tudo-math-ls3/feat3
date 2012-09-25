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
          const Index * row_ptr, const Index * row_ptr_end,  const Index count)
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
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;

template <typename DT_>
void Product<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, DT_> & r, const SparseMatrixCSR<Archs::GPU, DT_> & a, const DenseVector<Archs::GPU, DT_> & b)
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

template void Product<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, float> &, const SparseMatrixCSR<Archs::GPU, float> &, const DenseVector<Archs::GPU, float> &);
template void Product<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, double> &, const SparseMatrixCSR<Archs::GPU, double> &, const DenseVector<Archs::GPU, double> &);
