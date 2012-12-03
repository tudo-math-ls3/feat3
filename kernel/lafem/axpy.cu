// includes, FEAST
#include <kernel/lafem/axpy.hpp>

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
          const Index * row_ptr, const Index * row_ptr_end, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;

        DT_ sum(0);
        const Index end(row_ptr_end[idx]);
        for (Index i(row_ptr[idx]) ; i < end ; ++i)
        {
          sum += val[i] * x[col_ind[i]];
        }
        r[idx] = sum * a + y[idx];
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;

template <typename DT_>
void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, DT_> & r, const DT_ a, const DenseVector<Mem::CUDA, DT_> & x, const DenseVector<Mem::CUDA, DT_> & y)
{
  if (x.size() != y.size())
    throw InternalError("Vector size does not match!");
  if (x.size() != r.size())
    throw InternalError("Vector size does not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.size())/(double)(block.x));

  DT_ * r_gpu(r.elements());
  const DT_ * x_gpu(x.elements());
  const DT_ * y_gpu(y.elements());

  FEAST::LAFEM::Intern::cuda_axpy<<<grid, block>>>(r_gpu, a, x_gpu, y_gpu, r.size());
}

template void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, float> &, const float, const DenseVector<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &);
template void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, double>&, const double, const DenseVector<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &);

template <typename DT_>
void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & a, const DenseVector<Mem::CUDA, DT_> & x, const DenseVector<Mem::CUDA, DT_> & y)
{
  if (x.size() != y.size())
    throw InternalError("Vector size does not match!");
  if (x.size() != r.size())
    throw InternalError("Vector size does not match!");
  if (a.size() != r.size())
    throw InternalError("Vector size does not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.size())/(double)(block.x));

  DT_ * r_gpu(r.elements());
  const DT_ * x_gpu(x.elements());
  const DT_ * y_gpu(y.elements());
  const DT_ * a_gpu(a.elements());

  FEAST::LAFEM::Intern::cuda_axpyv<<<grid, block>>>(r_gpu, a_gpu, x_gpu, y_gpu, r.size());
}

template void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &);
template void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &);

template <typename DT_>
void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, DT_> & r, const DT_ a, const SparseMatrixCSR<Mem::CUDA, DT_> & P, const DenseVector<Mem::CUDA, DT_> & x, const DenseVector<Mem::CUDA, DT_> & y)
{
  if (x.size() != P.columns())
    throw InternalError("Vector size does not match!");
  if (P.rows() != r.size())
    throw InternalError("Vector size does not match!");
  if (y.size() != r.size())
    throw InternalError("Vector size does not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.size())/(double)(block.x));

  DT_ * r_gpu(r.elements());
  const DT_ * x_gpu(x.elements());
  const DT_ * y_gpu(y.elements());
  const DT_ * val_gpu(P.val());
  const Index * col_ind_gpu(P.col_ind());
  const Index * row_ptr_gpu(P.row_ptr());
  const Index * row_ptr_end_gpu(P.row_ptr_end());

  FEAST::LAFEM::Intern::cuda_axpy_mv_csr<<<grid, block>>>(r_gpu, a, x_gpu, y_gpu, val_gpu, col_ind_gpu, row_ptr_gpu, row_ptr_end_gpu, r.size());
}

template void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, float> &, const float, const SparseMatrixCSR<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &);
template void Axpy<Algo::CUDA>::value(DenseVector<Mem::CUDA, double> &, const double, const SparseMatrixCSR<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &);
