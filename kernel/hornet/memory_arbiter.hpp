#pragma once
#ifndef KERNEL_HORNET_MEMORY_ARBITER_HPP
#define KERNEL_HORNET_MEMORY_ARBITER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/instantiation_policy.hpp>

#include <map>
#include <stdio.h>


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
      std::map<void*, Intern::MemoryInfo> _memory;

      /// default CTOR
      MemoryArbiter()
      {
      }

    public:
      ~MemoryArbiter()
      {
        /// \todo throw error if any memory chunks exist
      }

      /// pointer to MemoryArbiter singleton
      friend MemoryArbiter* InstantiationPolicy<MemoryArbiter, Singleton>::instance();

      /// allocate new memory
      void * allocate_memory(Index bytes)
      {
        void * memory(::malloc(bytes));
        Intern::MemoryInfo mi;
        mi.counter = 1;
        _memory.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));

        return memory;
      }

      /// increase memory counter
      void increase_memory(void * address)
      {
        std::map<void*, Intern::MemoryInfo>::iterator it(_memory.find(address));
        if (it == _memory.end())
          throw InternalError("Memory address not found!");
        else
        {
          it->second.counter = it->second.counter + 1;
        }
      }

      /// release memory or decrease reference counter
      void release_memory(void * address)
      {
        std::map<void*, Intern::MemoryInfo>::iterator it(_memory.find(address));
        if (it == _memory.end())
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
  };
} // namespace FEAST

#endif // KERNEL_HORNET_MEMORY_ARBITER_HPP
