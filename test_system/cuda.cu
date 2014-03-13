#include <test_system/cuda.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

void FEAST::TestSystem::reset_device()
{
  cudaDeviceReset();
}
