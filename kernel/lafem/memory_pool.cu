// includes, FEAST
#include <kernel/lafem/memory_pool.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      template <typename DT_>
      __global__ void cuda_set_memory(DT_ * ptr, const DT_ val, const Index count)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx >= count)
          return;
        ptr[idx] = val;
      }

      __global__ void cuda_generate_hash(char * cd, const Index bytes, void * resultv)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if (idx != 0)
          return;

        unsigned long * result = (unsigned long *)resultv;
        unsigned long t(0);
        for (Index i(0) ; i < bytes ; ++i)
        {
          t += (cd[i] * i) % bytes;
        }
        t = t % bytes;
        result[0] = t;
      }
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;

MemoryPool<Mem::CUDA>::MemoryPool()
{
}

MemoryPool<Mem::CUDA>::~MemoryPool()
{
  if (_pool.size() > 0)
    throw InternalError("Memory Pool<GPU> still contains memory chunks!");
}

template <typename DT_>
DT_ * MemoryPool<Mem::CUDA>::allocate_memory(Index bytes)
{
  DT_ * memory(NULL);
  if (cudaErrorMemoryAllocation == cudaMalloc((void**)&memory, bytes))
    throw InternalError("MemoryPool<GPU> cuda allocation error (cudaErrorMemoryAllocation)");
  if (memory == NULL)
    throw InternalError("MemoryPool<GPU> allocation error (null pointer returned)");
  Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = bytes;
  _pool.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));
  return memory;
}

void MemoryPool<Mem::CUDA>::increase_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError("MemoryPool<GPU>::increase_memory: Memory address not found!");
  else
  {
    it->second.counter = it->second.counter + 1;
  }
}

void MemoryPool<Mem::CUDA>::release_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError("MemoryPool<GPU>::relase_memory: Memory address not found!");
  else
  {
    if(it->second.counter == 1)
    {
      cudaFree(address);
      _pool.erase(it);
    }
    else
    {
      it->second.counter = it->second.counter - 1;
    }
  }
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::download(void * dest, void * src, Index bytes)
{
  cudaMemcpy(dest, src, bytes, cudaMemcpyDeviceToHost);
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::upload(void * dest, void * src, Index bytes)
{
  cudaMemcpy(dest, src, bytes, cudaMemcpyHostToDevice);
}

template <typename DT_>
DT_ MemoryPool<Mem::CUDA>::get_element(const DT_ * data, Index index)
{
  const void * src(data + index);
  DT_ value;
  cudaMemcpy(&value, src, sizeof(DT_), cudaMemcpyDeviceToHost);
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
  FEAST::LAFEM::Intern::cuda_set_memory<<<grid, block>>>(address, val, count);
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::copy(void * dest, const void * src, const Index bytes)
{
  if (dest == src)
    return;

  cudaMemcpy(dest, src, bytes, cudaMemcpyDeviceToDevice);
}

unsigned long MemoryPool<Mem::CUDA>::generate_hash(void * data, Index bytes)
{
  dim3 grid(1,1,1);
  dim3 block(128,1,1);
  unsigned long result(0);
  void * result_gpu(0);
  cudaMalloc((void**)&result_gpu, sizeof(unsigned long));
  char * datac((char *)data);
  FEAST::LAFEM::Intern::cuda_generate_hash<<<grid, block>>>(datac, bytes, result_gpu);
  cudaMemcpy(&result, result_gpu, sizeof(unsigned long), cudaMemcpyDeviceToHost);
  cudaFree(result_gpu);
  return result;
}

template float * MemoryPool<Mem::CUDA>::allocate_memory<float>(Index bytes);
template double * MemoryPool<Mem::CUDA>::allocate_memory<double>(Index bytes);
template unsigned long * MemoryPool<Mem::CUDA>::allocate_memory<unsigned long>(Index bytes);

template void MemoryPool<Mem::CUDA>::download<float>(void * dest, void * src, Index bytes);
template void MemoryPool<Mem::CUDA>::download<double>(void * dest, void * src, Index bytes);
template void MemoryPool<Mem::CUDA>::download<unsigned long>(void * dest, void * src, Index bytes);

template void MemoryPool<Mem::CUDA>::upload<float>(void * dest, void * src, Index bytes);
template void MemoryPool<Mem::CUDA>::upload<double>(void * dest, void * src, Index bytes);
template void MemoryPool<Mem::CUDA>::upload<unsigned long>(void * dest, void * src, Index bytes);

template float MemoryPool<Mem::CUDA>::get_element(const float * data, Index index);
template double MemoryPool<Mem::CUDA>::get_element(const double * data, Index index);
template Index MemoryPool<Mem::CUDA>::get_element(const Index * data, Index index);

template void MemoryPool<Mem::CUDA>::set_memory(float * address , const float val, const Index count);
template void MemoryPool<Mem::CUDA>::set_memory(double * address , const double val, const Index count);
template void MemoryPool<Mem::CUDA>::set_memory(Index * address , const Index val, const Index count);

template void MemoryPool<Mem::CUDA>::copy<float>(void * dest, const void * src, const Index Bytes);
template void MemoryPool<Mem::CUDA>::copy<double>(void * dest, const void * src, const Index Bytes);
template void MemoryPool<Mem::CUDA>::copy<unsigned long>(void * dest, const void * src, const Index Bytes);
