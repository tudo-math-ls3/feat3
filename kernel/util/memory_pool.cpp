// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/memory_pool.hpp>


using namespace FEAST;
using namespace FEAST::Util;

template float * MemoryPool<Mem::Main>::allocate_memory<float>(const Index);
template double * MemoryPool<Mem::Main>::allocate_memory<double>(const Index);
template unsigned int * MemoryPool<Mem::Main>::allocate_memory<unsigned int>(const Index);
template unsigned long * MemoryPool<Mem::Main>::allocate_memory<unsigned long>(const Index);
#ifdef FEAST_HAVE_QUADMATH
template __float128 * MemoryPool<Mem::Main>::allocate_memory<__float128>(const Index);
#endif

template void MemoryPool<Mem::Main>::download<float>(float *, const float * const, const Index);
template void MemoryPool<Mem::Main>::download<double>(double *, const double * const, const Index);
template void MemoryPool<Mem::Main>::download<unsigned int>(unsigned int *, const unsigned int * const, const Index);
template void MemoryPool<Mem::Main>::download<unsigned long>(unsigned long *, const unsigned long * const, const Index);
#ifdef FEAST_HAVE_QUADMATH
template void MemoryPool<Mem::Main>::download<__float128>(__float128 *, const __float128 * const, const Index);
#endif

template void MemoryPool<Mem::Main>::upload<float>(float *, const float * const, const Index);
template void MemoryPool<Mem::Main>::upload<double>(double *, const double * const, const Index);
template void MemoryPool<Mem::Main>::upload<unsigned int>(unsigned int *, const unsigned int * const, const Index);
template void MemoryPool<Mem::Main>::upload<unsigned long>(unsigned long *, const unsigned long * const, const Index);
#ifdef FEAST_HAVE_QUADMATH
template void MemoryPool<Mem::Main>::upload<__float128>(__float128 *, const __float128 * const, const Index);
#endif
