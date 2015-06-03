// includes, FEAST
#include <kernel/util/memory_pool.hpp>

#include <cublas_v2.h>
#include "cusparse_v2.h"

#include <cstdio>

namespace FEAST
{
  namespace Util
  {
    namespace Intern
    {
      cublasHandle_t cublas_handle;
      cusparseHandle_t cusparse_handle;

      template <typename DT_>
      __global__ void cuda_set_memory(DT_ * ptr, const DT_ val, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        ptr[idx] = val;
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::Util;

MemoryPool<Mem::CUDA>::MemoryPool()
{
  if (CUBLAS_STATUS_SUCCESS != cublasCreate(&Intern::cublas_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cublasCreate failed!");
  if (CUSPARSE_STATUS_SUCCESS != cusparseCreate(&Intern::cusparse_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseCreate failed!");
}

MemoryPool<Mem::CUDA>::~MemoryPool()
{
  if (_pool.size() > 0)
  {
    std::cout << stderr << " Error: MemoryPool<CUDA> still contains memory chunks on deconstructor call" << std::endl;
    std::exit(1);
  }

  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "Pending cuda errors occured in execution!\n" + stringify(cudaGetErrorString(last_error)));
}

template <typename DT_>
DT_ * MemoryPool<Mem::CUDA>::allocate_memory(const Index count)
{
  DT_ * memory(nullptr);
  if (count == 0)
    return memory;

  if (cudaErrorMemoryAllocation == cudaMalloc((void**)&memory, count * sizeof(DT_)))
    throw InternalError("MemoryPool<CUDA> cuda allocation error (cudaErrorMemoryAllocation)");
  if (memory == nullptr)
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA> allocation error (null pointer returned)");
  Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = count * sizeof(DT_);
  _pool.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));
  return memory;
}

void MemoryPool<Mem::CUDA>::increase_memory(void * address)
{
  if (address == nullptr)
    return;

  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::increase_memory: Memory address not found!");
  else
  {
    it->second.counter = it->second.counter + 1;
  }
}

void MemoryPool<Mem::CUDA>::release_memory(void * address)
{
  if (address == nullptr)
    return;

  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::relase_memory: Memory address not found!");
  else
  {
    if(it->second.counter == 1)
    {
      if (cudaSuccess != cudaFree(address))
        throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::release_memory: cudaFree failed!");
      _pool.erase(it);
    }
    else
    {
      it->second.counter = it->second.counter - 1;
    }
  }
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::download(DT_ * dest, const DT_ * const src, const Index count)
{
  if (cudaSuccess != cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyDeviceToHost))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::download failed!");
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::upload(DT_ * dest, const DT_ * const src, const Index count)
{
  if (cudaSuccess != cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyHostToDevice))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::upload failed!");
}

template <typename DT_>
DT_ MemoryPool<Mem::CUDA>::get_element(const DT_ * data, const Index index)
{
  const void * src(data + index);
  DT_ value;
  if (cudaSuccess != cudaMemcpy(&value, src, sizeof(DT_), cudaMemcpyDeviceToHost))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::get_element failed!");
  return value;
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::set_memory(DT_ * address, const DT_ val, const Index count)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((count)/(double)(block.x));
  FEAST::Util::Intern::cuda_set_memory<<<grid, block>>>(address, val, count);
  cudaError_t last_error(cudaGetLastError());
  if (cudaSuccess != last_error)
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::set_memory failed!\n" + stringify(cudaGetErrorString(last_error)));
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::copy(DT_ * dest, const DT_ * src, const Index count)
{
  if (dest == src)
    return;

  if (cudaSuccess != cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyDeviceToDevice))
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CUDA>::copy failed!");
}

void MemoryPool<Mem::CUDA>::synchronize()
{
  cudaDeviceSynchronize();
}

void MemoryPool<Mem::CUDA>::reset_device()
{
  cudaDeviceReset();
}

void MemoryPool<Mem::CUDA>::shutdown_device()
{
  if (CUBLAS_STATUS_SUCCESS != cublasDestroy(Intern::cublas_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cublasDestroy failed!");
  if (CUSPARSE_STATUS_SUCCESS != cusparseDestroy(Intern::cusparse_handle))
    throw InternalError(__func__, __FILE__, __LINE__, "cusparseDestroy failed!");

  if (cudaSuccess != cudaDeviceSynchronize())
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSynchronize failed!");
  if (cudaSuccess != cudaDeviceReset())
    throw InternalError(__func__, __FILE__, __LINE__, "cudaDeviceSynchronize failed!");
}

template float * MemoryPool<Mem::CUDA>::allocate_memory<float>(const Index);
template double * MemoryPool<Mem::CUDA>::allocate_memory<double>(const Index);
template unsigned int * MemoryPool<Mem::CUDA>::allocate_memory<unsigned int>(const Index);
template unsigned long * MemoryPool<Mem::CUDA>::allocate_memory<unsigned long>(const Index);

template void MemoryPool<Mem::CUDA>::download<float>(float *, const float * const, const Index);
template void MemoryPool<Mem::CUDA>::download<double>(double *, const double * const, const Index);
template void MemoryPool<Mem::CUDA>::download<unsigned int>(unsigned int *, const unsigned int * const, const Index);
template void MemoryPool<Mem::CUDA>::download<unsigned long>(unsigned long *, const unsigned long * const, const Index);

template void MemoryPool<Mem::CUDA>::upload<float>(float *, const float * const, const Index);
template void MemoryPool<Mem::CUDA>::upload<double>(double *, const double * const , const Index);
template void MemoryPool<Mem::CUDA>::upload<unsigned int>(unsigned int *, const unsigned int * const, const Index);
template void MemoryPool<Mem::CUDA>::upload<unsigned long>(unsigned long *, const unsigned long * const, const Index);

template float MemoryPool<Mem::CUDA>::get_element(const float *, const Index);
template double MemoryPool<Mem::CUDA>::get_element(const double *, const Index);
template unsigned int MemoryPool<Mem::CUDA>::get_element(const unsigned int *, const Index);
template unsigned long MemoryPool<Mem::CUDA>::get_element(const unsigned long *, const Index);

template void MemoryPool<Mem::CUDA>::set_memory(float * , const float, const Index);
template void MemoryPool<Mem::CUDA>::set_memory(double * , const double, const Index);
template void MemoryPool<Mem::CUDA>::set_memory(unsigned int * , const unsigned int, const Index);
template void MemoryPool<Mem::CUDA>::set_memory(unsigned long * , const unsigned long, const Index);

template void MemoryPool<Mem::CUDA>::copy<float>(float *, const float *, const Index);
template void MemoryPool<Mem::CUDA>::copy<double>(double *, const double *, const Index);
template void MemoryPool<Mem::CUDA>::copy<unsigned int>(unsigned int *, const unsigned int *, const Index);
template void MemoryPool<Mem::CUDA>::copy<unsigned long>(unsigned long *, const unsigned long *, const Index);
