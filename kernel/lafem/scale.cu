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
