// includes, FEAST
#include <kernel/hornet/memory_arbiter.hpp>

#include <stdio.h>


using namespace FEAST;

MemoryArbiter::MemoryArbiter()
{
}

MemoryArbiter::~MemoryArbiter()
{
  /// \todo throw error if any memory chunks exist
}

void * MemoryArbiter::allocate_memory(Index bytes)
{
  void * memory(::malloc(bytes));
  Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = bytes;
  _memory_pool.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));

  return memory;
}

void MemoryArbiter::increase_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_memory_pool.find(address));
  if (it == _memory_pool.end())
    throw InternalError("Memory address not found!");
  else
  {
    it->second.counter = it->second.counter + 1;
  }
}

void MemoryArbiter::release_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_memory_pool.find(address));
  if (it == _memory_pool.end())
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
