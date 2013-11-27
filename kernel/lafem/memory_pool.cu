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

MemoryPool<Mem::CUDA>::MemoryPool()
{
}

MemoryPool<Mem::CUDA>::~MemoryPool()
{
  if (_pool.size() > 0)
    throw InternalError("Memory Pool<GPU> still contains memory chunks!");
}

template <typename DT_>
DT_ * MemoryPool<Mem::CUDA>::allocate_memory(const Index count)
{
  DT_ * memory(NULL);
  if (cudaErrorMemoryAllocation == cudaMalloc((void**)&memory, count * sizeof(DT_)))
    throw InternalError("MemoryPool<GPU> cuda allocation error (cudaErrorMemoryAllocation)");
  if (memory == NULL)
    throw InternalError("MemoryPool<GPU> allocation error (null pointer returned)");
  Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = count * sizeof(DT_);
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
void MemoryPool<Mem::CUDA>::download(DT_ * dest, DT_ * src, const Index count)
{
  cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyDeviceToHost);
}

template <typename DT_>
void MemoryPool<Mem::CUDA>::upload(DT_ * dest, DT_ * src, const Index count)
{
  cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyHostToDevice);
}

template <typename DT_>
DT_ MemoryPool<Mem::CUDA>::get_element(const DT_ * data, const Index index)
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
void MemoryPool<Mem::CUDA>::copy(DT_ * dest, const DT_ * src, const Index count)
{
  if (dest == src)
    return;

  cudaMemcpy(dest, src, count * sizeof(DT_), cudaMemcpyDeviceToDevice);
}

template float * MemoryPool<Mem::CUDA>::allocate_memory<float>(const Index);
template double * MemoryPool<Mem::CUDA>::allocate_memory<double>(const Index);
template Index * MemoryPool<Mem::CUDA>::allocate_memory<Index>(const Index);

template void MemoryPool<Mem::CUDA>::download<float>(float *, float *, const Index);
template void MemoryPool<Mem::CUDA>::download<double>(double *, double *, const Index);
template void MemoryPool<Mem::CUDA>::download<Index>(Index *, Index *, const Index);

template void MemoryPool<Mem::CUDA>::upload<float>(flaot *, float *, const Index);
template void MemoryPool<Mem::CUDA>::upload<double>(double *, double * const Index);
template void MemoryPool<Mem::CUDA>::upload<Index>(Index *, Index *, const Index);

template float MemoryPool<Mem::CUDA>::get_element(const float *, const Index);
template double MemoryPool<Mem::CUDA>::get_element(const double *, const Index);
template Index MemoryPool<Mem::CUDA>::get_element(const Index *, const Index);

template void MemoryPool<Mem::CUDA>::set_memory(float * , const float, const const Index);
template void MemoryPool<Mem::CUDA>::set_memory(double * , const double, const const Index);
template void MemoryPool<Mem::CUDA>::set_memory(Index * , const Index, const Index);

template void MemoryPool<Mem::CUDA>::copy<float>(float *, const float *, const Index);
template void MemoryPool<Mem::CUDA>::copy<double>(double *, const double *, const Index);
template void MemoryPool<Mem::CUDA>::copy<Index>(Index *, const Index *, const Index);
