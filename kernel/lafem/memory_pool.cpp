// includes, FEAST
#include <kernel/lafem/memory_pool.hpp>

#include <cstdio>


using namespace FEAST;
using namespace FEAST::LAFEM;

MemoryPool<Mem::Main>::MemoryPool()
{
}

MemoryPool<Mem::Main>::~MemoryPool()
{
  if (_pool.size() > 0)
  {
    std::cout << stderr << " Error: MemoryPool<CPU> still contains memory chunks on deconstructor call" << std::endl;
    std::exit(1);
  }
}

template <typename DT_>
DT_ * MemoryPool<Mem::Main>::allocate_memory(const Index count)
{
  DT_ * memory(NULL);
  memory = (DT_*)::malloc(count * sizeof(DT_));
  if (memory == NULL)
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU> allocation error!");

  Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = count * sizeof(DT_);
  _pool.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));

  return memory;
}

void MemoryPool<Mem::Main>::increase_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU>::increase_memory: Memory address not found!");
  else
  {
    it->second.counter = it->second.counter + 1;
  }
}

void MemoryPool<Mem::Main>::release_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU>::release_memory: Memory address not found!");
  else
  {
    if(it->second.counter == 1)
    {
      ::free(address);
      _pool.erase(it);
    }
    else
    {
      it->second.counter = it->second.counter - 1;
    }
  }
}

template <typename DT_>
void MemoryPool<Mem::Main>::download(DT_ * dest, DT_ * src, const Index count)
{
  if (dest == src)
    return;

  ::memcpy(dest, src, count * sizeof(DT_));
}

template <typename DT_>
void MemoryPool<Mem::Main>::upload(DT_ * dest, DT_ * src, const Index count)
{
  if (dest == src)
    return;

  ::memcpy(dest, src, count * sizeof(DT_));
}

template float * MemoryPool<Mem::Main>::allocate_memory<float>(const Index);
template double * MemoryPool<Mem::Main>::allocate_memory<double>(const Index);
template Index * MemoryPool<Mem::Main>::allocate_memory<Index>(const Index);

template void MemoryPool<Mem::Main>::download<float>(float *, float *, const Index);
template void MemoryPool<Mem::Main>::download<double>(double *, double *, const Index);
template void MemoryPool<Mem::Main>::download<Index>(Index *t, Index *, const Index);

template void MemoryPool<Mem::Main>::upload<float>(float *, float *, const Index);
template void MemoryPool<Mem::Main>::upload<double>(double *, double *, const Index);
template void MemoryPool<Mem::Main>::upload<Index>(Index *, Index *, const Index);
