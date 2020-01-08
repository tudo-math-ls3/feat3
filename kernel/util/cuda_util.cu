// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>


void FEAT::Util::cuda_set_device(const int device)
{
  cudaSetDevice(device);
}

void * FEAT::Util::cuda_malloc_host(const Index bytes)
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

void FEAT::Util::cuda_free_host(void * address)
{
  if (address == nullptr)
    return;

  if (cudaSuccess != cudaFreeHost(address))
    throw InternalError(__func__, __FILE__, __LINE__, "Util::cuda_free_host: cudaFreeHost failed!");
}

void FEAT::Util::cuda_check_last_error()
{
  cudaDeviceSynchronize();
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error)));
}

void * FEAT::Util::cuda_get_device_pointer(void * host)
{
  void * device(nullptr);
  if (cudaSuccess != cudaHostGetDevicePointer((void**)&device, host, 0))
    throw InternalError(__func__, __FILE__, __LINE__, "cudaHostGetDevicePointer failed!");
  return device;
}
