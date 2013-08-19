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

template <typename DT_>
void * MemoryPool<Mem::Main>::allocate_memory(Index bytes)
{
  void * memory(NULL);
  memory = ::malloc(bytes);
  if (memory == NULL)
    throw InternalError("MemoryPool<CPU> allocation error!");
#ifdef FEAST_GMP
  if (typeid(DT_) == typeid(mpf_class))
  {
    for (DT_ * location((DT_ *) memory), * end(location + bytes/sizeof(DT_)) ; location != end ; ++location)
      new (location) DT_(0);
  }
#endif

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

template <typename DT_>
void MemoryPool<Mem::Main>::download(void * dest, void * src, Index bytes)
{
  if (dest == src)
    return;

#ifdef FEAST_GMP
  if (typeid(DT_) == typeid(mpf_class))
  {
    const DT_ * s((DT_ *) src);
    for (DT_ * d((DT_ * )dest), * d_end((DT_ *)dest + bytes/sizeof(DT_)) ; d != d_end ; ++d, ++s)
    {
      *d = *s;
    }
  }
  else
#endif
    ::memcpy(dest, src, bytes);
}

template <typename DT_>
void MemoryPool<Mem::Main>::upload(void * dest, void * src, Index bytes)
{
  if (dest == src)
    return;

#ifdef FEAST_GMP
  if (typeid(DT_) == typeid(mpf_class))
  {
    const DT_ * s((DT_ *) src);
    for (DT_ * d((DT_ * )dest), * d_end((DT_ *)dest + bytes/sizeof(DT_)) ; d != d_end ; ++d, ++s)
    {
      *d = *s;
    }
  }
  else
#endif
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

template void * MemoryPool<Mem::Main>::allocate_memory<float>(Index bytes);
template void * MemoryPool<Mem::Main>::allocate_memory<double>(Index bytes);
#ifdef FEAST_GMP
template void * MemoryPool<Mem::Main>::allocate_memory<mpf_class>(Index bytes);
#endif
template void * MemoryPool<Mem::Main>::allocate_memory<unsigned long>(Index bytes);

template void MemoryPool<Mem::Main>::download<float>(void * dest, void * src, Index bytes);
template void MemoryPool<Mem::Main>::download<double>(void * dest, void * src, Index bytes);
#ifdef FEAST_GMP
template void MemoryPool<Mem::Main>::download<mpf_class>(void * dest, void * src, Index bytes);
#endif
template void MemoryPool<Mem::Main>::download<unsigned long>(void * dest, void * src, Index bytes);

template void MemoryPool<Mem::Main>::upload<float>(void * dest, void * src, Index bytes);
template void MemoryPool<Mem::Main>::upload<double>(void * dest, void * src, Index bytes);
#ifdef FEAST_GMP
template void MemoryPool<Mem::Main>::upload<mpf_class>(void * dest, void * src, Index bytes);
#endif
template void MemoryPool<Mem::Main>::upload<unsigned long>(void * dest, void * src, Index bytes);

template void MemoryPool<Mem::Main>::set_memory(float * address , const float val, const Index count);
template void MemoryPool<Mem::Main>::set_memory(double * address , const double val, const Index count);
#ifdef FEAST_GMP
template void MemoryPool<Mem::Main>::set_memory(mpf_class * address , const mpf_class val, const Index count);
#endif
template void MemoryPool<Mem::Main>::set_memory(Index * address , const Index val, const Index count);
