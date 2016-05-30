#include <kernel/util/memory_pool.hpp>

// static member initialisation
std::map<void*, FEAT::Util::Intern::MemoryInfo> FEAT::MemoryPool<FEAT::Mem::Main>::_pool;
std::map<void*, FEAT::Util::Intern::MemoryInfo> FEAT::MemoryPool<FEAT::Mem::Main>::_pinned_pool;
