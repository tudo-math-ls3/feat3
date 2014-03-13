#include <test_system/cuda.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

void reset_device()
{
  cudaDeviceReset();
}
