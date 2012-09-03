// includes, FEAST
#include <kernel/lafem/memory_pool.hpp>

#include <cstdio>


using namespace FEAST;
using namespace FEAST::LAFEM;

MemoryPool<Archs::CPU>::MemoryPool()
{
}

MemoryPool<Archs::CPU>::~MemoryPool()
{
  if (_pool.size() > 0)
    throw InternalError("Memory Pool still contains memory chunks!");
}

void * MemoryPool<Archs::CPU>::allocate_memory(Index bytes)
{
  void * memory(::malloc(bytes));
  Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = bytes;
  _pool.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));

  return memory;
}

void MemoryPool<Archs::CPU>::increase_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError("Memory address not found!");
  else
  {
    it->second.counter = it->second.counter + 1;
  }
}

void MemoryPool<Archs::CPU>::release_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError("Memory address not found!");
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

void MemoryPool<Archs::CPU>::download(void * dest, void * src, Index bytes)
{
  ::memcpy(dest, src, bytes);
}

void MemoryPool<Archs::CPU>::upload(void * dest, void * src, Index bytes)
{
  ::memcpy(dest, src, bytes);
}

template <typename DT_>
void MemoryPool<Archs::CPU>::set_memory(DT_ * address, const DT_ val, const Index count)
{
  for (Index i(0) ; i < count ; ++i)
  {
    address[i] = val;
  }
}

template void MemoryPool<Archs::CPU>::set_memory(float * address , const float val, const Index count);
template void MemoryPool<Archs::CPU>::set_memory(double * address , const double val, const Index count);
template void MemoryPool<Archs::CPU>::set_memory(Index * address , const Index val, const Index count);
