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
void Scale<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, DT_> & r, const DenseVector<Archs::GPU, DT_> & x, const DT_ s)
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

template void Scale<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, float> &, const DenseVector<Archs::GPU, float> &, const float);
template void Scale<Archs::GPU, Archs::CUDA>::value(DenseVector<Archs::GPU, double> &, const DenseVector<Archs::GPU, double> &, const double);
