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
    throw InternalError("MemoryPool<CPU> still contains memory chunks!");
}

void * MemoryPool<Mem::Main>::allocate_memory(Index bytes)
{
  void * memory(NULL);
  memory = ::malloc(bytes);
  if (memory == NULL)
    throw InternalError("MemoryPool<CPU> allocation error!");
  Intern::MemoryInfo mi;
  mi.counter = 1;
  mi.size = bytes;
  _pool.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));

  return memory;
}

void MemoryPool<Mem::Main>::increase_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError("MemoryPool<CPU>::increase_memory: Memory address not found!");
  else
  {
    it->second.counter = it->second.counter + 1;
  }
}

void MemoryPool<Mem::Main>::release_memory(void * address)
{
  std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
  if (it == _pool.end())
    throw InternalError("MemoryPool<CPU>::release_memory: Memory address not found!");
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

void MemoryPool<Mem::Main>::download(void * dest, void * src, Index bytes)
{
  ::memcpy(dest, src, bytes);
}

void MemoryPool<Mem::Main>::upload(void * dest, void * src, Index bytes)
{
  ::memcpy(dest, src, bytes);
}

template <typename DT_>
void MemoryPool<Mem::Main>::set_memory(DT_ * address, const DT_ val, const Index count)
{
  for (Index i(0) ; i < count ; ++i)
  {
    address[i] = val;
  }
}

unsigned long MemoryPool<Mem::Main>::generate_hash(void * data, Index bytes)
{
  char * cd((char * )data);
  Index t(0);
  for (Index i(0) ; i < bytes ; ++i)
  {
    t += ((Index)cd[i] * i) % bytes;
  }
  t = t % bytes;
  return t;
}

template void MemoryPool<Mem::Main>::set_memory(float * address , const float val, const Index count);
template void MemoryPool<Mem::Main>::set_memory(double * address , const double val, const Index count);
template void MemoryPool<Mem::Main>::set_memory(Index * address , const Index val, const Index count);
