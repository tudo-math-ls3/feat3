#include <kernel/util/memory_pool.hpp>

// static member initialisation
std::map<void*, FEAST::Util::Intern::MemoryInfo> FEAST::Util::MemoryPool<FEAST::Mem::Main>::_pool;
