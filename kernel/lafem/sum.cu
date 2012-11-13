// includes, FEAST
#include <kernel/lafem/sum.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_sum(DT_ * r, const DT_ * x, const DT_ * y, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        r[idx] = x[idx] + y[idx];
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;

template <typename DT_>
void Sum<Mem::CUDA, Algo::CUDA>::value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & x, const DenseVector<Mem::CUDA, DT_> & y)
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

  FEAST::LAFEM::Intern::cuda_sum<<<grid, block>>>(r_gpu, x_gpu, y_gpu, r.size());
}

template void Sum<Mem::CUDA, Algo::CUDA>::value(DenseVector<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &, const DenseVector<Mem::CUDA, float> &);
template void Sum<Mem::CUDA, Algo::CUDA>::value(DenseVector<Mem::CUDA, double>&, const DenseVector<Mem::CUDA, double> &, const DenseVector<Mem::CUDA, double> &);
