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
    }
  }
}


using namespace FEAST;
using namespace FEAST::LAFEM;

MemoryPool<Archs::GPU>::MemoryPool()
{
}

MemoryPool<Archs::GPU>::~MemoryPool()
{
  if (_pool.size() > 0)
    throw InternalError("Memory Pool still contains memory chunks!");
}

void * MemoryPool<Archs::GPU>::allocate_memory(Index bytes)
{
  void * memory(NULL);
  if (cudaErrorMemoryAllocation == cudaMalloc((void**)&memory, bytes))
    throw InternalError("GPU Allocation Error");
  Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = bytes;
  _pool.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));

  return memory;
}

void MemoryPool<Archs::GPU>::increase_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError("Memory address not found!");
  else
  {
    it->second.counter = it->second.counter + 1;
  }
}

void MemoryPool<Archs::GPU>::release_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError("Memory address not found!");
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

void MemoryPool<Archs::GPU>::download(void * dest, void * src, Index bytes)
{
  cudaMemcpy(dest, src, bytes, cudaMemcpyDeviceToHost);
}

void MemoryPool<Archs::GPU>::upload(void * dest, void * src, Index bytes)
{
  cudaMemcpy(dest, src, bytes, cudaMemcpyHostToDevice);
}

template <typename DT_>
DT_ MemoryPool<Archs::GPU>::get_element(const DT_ * data, Index index)
{
  const void * src(data + index);
  DT_ value;
  cudaMemcpy(&value, src, sizeof(DT_), cudaMemcpyDeviceToHost);
  return value;
}

template <typename DT_>
void MemoryPool<Archs::GPU>::modify_element(DT_ * data, Index index, DT_ value)
{
  void * dest(data + index);
  cudaMemcpy(dest, &value, sizeof(DT_), cudaMemcpyHostToDevice);
}

template <typename DT_>
void MemoryPool<Archs::GPU>::set_memory(DT_ * address, const DT_ val, const Index count)
{
  Index blocksize(128);
  dim3 grid;
  dim3 block;
  block.x = blocksize;
  grid.x = (unsigned)ceil((count)/(double)(block.x));
  FEAST::LAFEM::Intern::cuda_set_memory<<<grid, blocksize>>>(address, val, count);
}

template float MemoryPool<Archs::GPU>::get_element(const float * data, Index index);
template double MemoryPool<Archs::GPU>::get_element(const double * data, Index index);
template Index MemoryPool<Archs::GPU>::get_element(const Index * data, Index index);
template void MemoryPool<Archs::GPU>::modify_element(float * data, Index index, float value);
template void MemoryPool<Archs::GPU>::modify_element(double * data, Index index, double value);
template void MemoryPool<Archs::GPU>::modify_element(Index * data, Index index, Index value);
template void MemoryPool<Archs::GPU>::set_memory(float * address , const float val, const Index count);
template void MemoryPool<Archs::GPU>::set_memory(double * address , const double val, const Index count);
template void MemoryPool<Archs::GPU>::set_memory(Index * address , const Index val, const Index count);
