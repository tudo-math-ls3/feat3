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
void Sum<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, DT_> & r, const DenseVector<Archs::GPU, DT_> & x, const DenseVector<Archs::GPU, DT_> & y)
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

template void Sum<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, float> &, const DenseVector<Archs::GPU, float> &, const DenseVector<Archs::GPU, float> &);
template void Sum<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, double>&, const DenseVector<Archs::GPU, double> &, const DenseVector<Archs::GPU, double> &);
