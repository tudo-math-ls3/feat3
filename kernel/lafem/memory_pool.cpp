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
DT_ * MemoryPool<Mem::Main>::allocate_memory(const Index count)
{
  DT_ * memory(NULL);
  memory = (DT_*)::malloc(count * sizeof(DT_));
  if (memory == NULL)
    throw InternalError("MemoryPool<CPU> allocation error!");
#ifdef FEAST_GMP
  if (typeid(DT_) == typeid(mpf_class))
  {
    for (DT_ * location((DT_ *) memory), * end(location + count) ; location != end ; ++location)
      new (location) DT_(0);
  }
#endif

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
void MemoryPool<Mem::Main>::download(DT_ * dest, DT_ * src, const Index count)
{
  if (dest == src)
    return;

#ifdef FEAST_GMP
  if (typeid(DT_) == typeid(mpf_class))
  {
    const DT_ * s( src);
    for (DT_ * d(dest), * d_end(dest + count) ; d != d_end ; ++d, ++s)
    {
      *d = *s;
    }
  }
  else
#endif
    ::memcpy(dest, src, count * sizeof(DT_));
}

template <typename DT_>
void MemoryPool<Mem::Main>::upload(DT_ * dest, DT_ * src, const Index count)
{
  if (dest == src)
    return;

#ifdef FEAST_GMP
  if (typeid(DT_) == typeid(mpf_class))
  {
    const DT_ * s( src);
    for (DT_ * d(dest), * d_end(dest + count) ; d != d_end ; ++d, ++s)
    {
      *d = *s;
    }
  }
  else
#endif
    ::memcpy(dest, src, count * sizeof(DT_));
}

template <typename DT_>
void MemoryPool<Mem::Main>::set_memory(DT_ * address, const DT_ val, const Index count)
{
  for (Index i(0) ; i < count ; ++i)
  {
    address[i] = val;
  }
}

template float * MemoryPool<Mem::Main>::allocate_memory<float>(const Index);
template double * MemoryPool<Mem::Main>::allocate_memory<double>(const Index);
#ifdef FEAST_GMP
template mpf_class * MemoryPool<Mem::Main>::allocate_memory<mpf_class>(const Index);
#endif
template unsigned long * MemoryPool<Mem::Main>::allocate_memory<unsigned long>(const Index);

template void MemoryPool<Mem::Main>::download<float>(float *, float *, const Index);
template void MemoryPool<Mem::Main>::download<double>(double *, double *, const Index);
#ifdef FEAST_GMP
template void MemoryPool<Mem::Main>::download<mpf_class>(mpf_class *, mpf_class *, const Index);
#endif
template void MemoryPool<Mem::Main>::download<unsigned long>(unsigned long *t, unsigned long *, const Index);

template void MemoryPool<Mem::Main>::upload<float>(float *, float *, const Index);
template void MemoryPool<Mem::Main>::upload<double>(double *, double *, const Index);
#ifdef FEAST_GMP
template void MemoryPool<Mem::Main>::upload<mpf_class>(mpf_class *, mpf_class *, const Index);
#endif
template void MemoryPool<Mem::Main>::upload<unsigned long>(unsigned long *, unsigned long *, const Index);

template void MemoryPool<Mem::Main>::set_memory(float * , const float, const Index);
template void MemoryPool<Mem::Main>::set_memory(double * , const double, const Index);
#ifdef FEAST_GMP
template void MemoryPool<Mem::Main>::set_memory(mpf_class * , const mpf_class, const Index);
#endif
template void MemoryPool<Mem::Main>::set_memory(Index * , const Index, const Index);
