// includes, FEAST
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/string.hpp>
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

  if (cudaErrorMemoryAllocation == cudaMallocHost((void**)&memory, bytes, cudaHostAllocMapped))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA> cuda pinned allocation error (cudaErrorMemoryAllocation)");
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

void FEAST::Util::cuda_check_last_error()
{
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
}

void * FEAST::Util::cuda_get_device_pointer(void * host)
{
  void * device(nullptr);
  if (cudaSuccess != cudaHostGetDevicePointer((void**)&device, host, 0))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaHostGetDevicePointer failed!");
  return device;
}
