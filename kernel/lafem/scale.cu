// includes, FEAST
#include <kernel/lafem/scale.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_scale(DT_ * r, const DT_ * x, const DT_ s, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        r[idx] = x[idx] * s;
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;

template <typename DT_>
void Scale<Algo::CUDA>::value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & x, const DT_ s)
{
  if (x.size() != r.size())
    throw InternalError("Vector size does not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.size())/(double)(block.x));

  DT_ * r_gpu(r.elements());
  const DT_ * x_gpu(x.elements());

  FEAST::LAFEM::Intern::cuda_scale<<<grid, block>>>(r_gpu, x_gpu, s, r.size());
}
template void Scale<Algo::CUDA>::value(DenseVector<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &, const float);
template void Scale<Algo::CUDA>::value(DenseVector<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &, const double);

template <typename DT_>
void Scale<Algo::CUDA>::value(SparseMatrixCOO<Mem::CUDA, DT_> & r, const SparseMatrixCOO<Mem::CUDA, DT_> & x, const DT_ s)
{
  if(x.rows() != r.rows())
    throw InternalError("Matrix Rows doe not match!");
  if(x.columns() != r.columns())
    throw InternalError("Matrix Columns doe not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.used_elements())/(double)(block.x));

  DT_ * r_gpu(r.val());
  const DT_ * x_gpu(x.val());

  FEAST::LAFEM::Intern::cuda_scale<<<grid, block>>>(r_gpu, x_gpu, s, r.used_elements());
}
template void Scale<Algo::CUDA>::value(SparseMatrixCOO<Mem::CUDA, float> &, const SparseMatrixCOO<Mem::CUDA, float> &, const float);
template void Scale<Algo::CUDA>::value(SparseMatrixCOO<Mem::CUDA, double> &, const SparseMatrixCOO<Mem::CUDA, double> &, const double);

template <typename DT_>
void Scale<Algo::CUDA>::value(SparseMatrixCSR<Mem::CUDA, DT_> & r, const SparseMatrixCSR<Mem::CUDA, DT_> & x, const DT_ s)
{
  if(x.rows() != r.rows())
    throw InternalError("Matrix Rows doe not match!");
  if(x.columns() != r.columns())
    throw InternalError("Matrix Columns doe not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.used_elements())/(double)(block.x));

  DT_ * r_gpu(r.val());
  const DT_ * x_gpu(x.val());

  FEAST::LAFEM::Intern::cuda_scale<<<grid, block>>>(r_gpu, x_gpu, s, r.used_elements());
}
template void Scale<Algo::CUDA>::value(SparseMatrixCSR<Mem::CUDA, float> &, const SparseMatrixCSR<Mem::CUDA, float> &, const float);
template void Scale<Algo::CUDA>::value(SparseMatrixCSR<Mem::CUDA, double> &, const SparseMatrixCSR<Mem::CUDA, double> &, const double);

template <typename DT_>
void Scale<Algo::CUDA>::value(SparseMatrixELL<Mem::CUDA, DT_> & r, const SparseMatrixELL<Mem::CUDA, DT_> & x, const DT_ s)
{
  if(x.rows() != r.rows())
    throw InternalError("Matrix Rows doe not match!");
  if(x.columns() != r.columns())
    throw InternalError("Matrix Columns doe not match!");

  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((r.stride() * r.num_cols_per_row())/(double)(block.x));

  DT_ * r_gpu(r.Ax());
  const DT_ * x_gpu(x.Ax());

  FEAST::LAFEM::Intern::cuda_scale<<<grid, block>>>(r_gpu, x_gpu, s, r.stride() * r.num_cols_per_row());
}
template void Scale<Algo::CUDA>::value(SparseMatrixELL<Mem::CUDA, float> &, const SparseMatrixELL<Mem::CUDA, float> &, const float);
template void Scale<Algo::CUDA>::value(SparseMatrixELL<Mem::CUDA, double> &, const SparseMatrixELL<Mem::CUDA, double> &, const double);
