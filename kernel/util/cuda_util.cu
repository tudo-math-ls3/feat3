// includes, FEAST
#include <kernel/util/cuda_util.hpp>


void FEAST::Util::cuda_set_device(int device)
{
  cudaSetDevice(device);
}
