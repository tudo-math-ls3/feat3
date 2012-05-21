// includes, FEAST
#include <kernel/hornet/memory_pool.hpp>

#include <stdio.h>


using namespace FEAST;

MemoryPool<Archs::CPU>::MemoryPool()
{
}

MemoryPool<Archs::CPU>::~MemoryPool()
{
  /// \todo throw error if any memory chunks exist
}

//template<>
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
    }
    else
    {
      it->second.counter = it->second.counter - 1;
    }
  }
}
