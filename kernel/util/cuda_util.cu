// includes, FEAST
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/exception.hpp>


void FEAST::Util::cuda_set_device(const int device)
{
  cudaSetDevice(device);
}

void * FEAST::Util::cuda_malloc_host(const Index bytes)
{
  void * memory(nullptr);
  if (bytes == 0)
    return memory;

  if (cudaErrorMemoryAllocation == cudaMallocHost((void**)&memory, bytes))
    throw InternalError("MemoryPool<CUDA> cuda pinned allocation error (cudaErrorMemoryAllocation)");
  if (memory == nullptr)
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_malloc_host allocation error (null pointer returned)");
  return memory;
}

void FEAST::Util::cuda_free_host(void * address)
{
  if (address == nullptr)
    return;

  if (cudaSuccess != cudaFreeHost(address))
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_free_host: cudaFreeHost failed!");
}
