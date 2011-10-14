#pragma once
#ifndef KERNEL_HORNET_MEMORY_ARBITER_HPP
#define KERNEL_HORNET_MEMORY_ARBITER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/instantiation_policy.hpp>

#include <map>


namespace FEAST
{
  namespace Intern
  {
    struct MemoryInfo
    {
      Index counter;
      /// \todo backend
    };
  }

  class MemoryArbiter
      : public InstantiationPolicy<MemoryArbiter, Singleton>
  {
    private:
      std::map<void*, Intern::MemoryInfo> _memory_pool;

      /// default CTOR
      MemoryArbiter();

    public:
      ~MemoryArbiter();

      /// pointer to MemoryArbiter singleton
      friend MemoryArbiter* InstantiationPolicy<MemoryArbiter, Singleton>::instance();

      /// allocate new memory
      void * allocate_memory(Index bytes);

      /// increase memory counter
      void increase_memory(void * address);

      /// release memory or decrease reference counter
      void release_memory(void * address);
  };
} // namespace FEAST

#endif // KERNEL_HORNET_MEMORY_ARBITER_HPP
