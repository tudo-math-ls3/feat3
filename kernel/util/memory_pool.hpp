#pragma once
#ifndef KERNEL_UTIL_MEMORY_POOL_HPP
#define KERNEL_UTIL_MEMORY_POOL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/archs.hpp>

#include <map>
#include <cstring>
#include <typeinfo>
#include <cstdio>
#include <cstddef>


namespace FEAST
{
  namespace Util
  {
    /// \cond internal
    namespace Intern
    {
      struct MemoryInfo
      {
        Index counter;
        Index size;
      };
    }
    /// \endcond

    template <typename Mem_>
    class MemoryPool
        : public InstantiationPolicy<MemoryPool<Mem_>, Singleton>
    {
    };

    /**
     * \brief Memory managment.
     *
     * This class manages the used memory chunks and releases them, if neccesary.
     *
     * \author Dirk Ribbrock
     */
    template <>
    class MemoryPool<Mem::Main>
        : public InstantiationPolicy<MemoryPool<Mem::Main>, Singleton>
    {
      private:
        /// Map of all memory chunks in use.
        std::map<void*, Intern::MemoryInfo> _pool;

        /// default CTOR
        MemoryPool()
        {
        }

      public:
        ~MemoryPool()
        {
          if (_pool.size() > 0)
          {
            std::cout << stderr << " Error: MemoryPool<CPU> still contains memory chunks on deconstructor call" << std::endl;
            std::exit(1);
          }
        }

        /// pointer to MemoryPool singleton
        friend MemoryPool* InstantiationPolicy<MemoryPool<Mem::Main>, Singleton>::instance();

        /// allocate new memory
        template <typename DT_>
        DT_ * allocate_memory(const Index count)
        {
          DT_ * memory(nullptr);
          if (count == 0)
            return memory;

          memory = (DT_*)::malloc(count * sizeof(DT_));
          if (memory == nullptr)
            throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU> allocation error!");

          Intern::MemoryInfo mi;
          mi.counter = 1;
          mi.size = count * sizeof(DT_);
          _pool.insert(std::pair<void*, Intern::MemoryInfo>(memory, mi));

          return memory;
        }

        /// increase memory counter
        void increase_memory(void * address)
        {
          if (address == nullptr)
            return;

          std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
          if (it == _pool.end())
            throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU>::increase_memory: Memory address not found!");
          else
          {
            it->second.counter = it->second.counter + 1;
          }
        }

        /// release memory or decrease reference counter
        void release_memory(void * address)
        {
          if (address == nullptr)
            return;

          std::map<void*, Intern::MemoryInfo>::iterator it(_pool.find(address));
          if (it == _pool.end())
            throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU>::release_memory: Memory address not found!");
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

        /// download memory chunk to host memory
        template <typename DT_>
        inline static void download(DT_ * dest, const DT_ * const src, const Index count)
        {
          if (dest == src)
            return;

          ::memcpy(dest, src, count * sizeof(DT_));
        }

        /// upload memory chunk from host memory to device memory
        template <typename DT_>
        inline static void upload(DT_ * dest, const DT_ * const src, const Index count)
        {
          if (dest == src)
            return;

          ::memcpy(dest, src, count * sizeof(DT_));
        }

        /// recieve element
        template <typename DT_>
        inline static const DT_ & get_element(const DT_ * data, const Index index)
        {
          return data[index];
        }

        /// set memory to specific value
        template <typename DT_>
        static void set_memory(DT_ * address, const DT_ val, const Index count = 1)
        {
          for (Index i(0) ; i < count ; ++i)
          {
            address[i] = val;
          }
        }

        /// Copy memory area from src to dest
        template <typename DT_>
        static void copy(DT_ * dest, const DT_ * src, const Index count)
        {
          if (dest == src)
            return;

          ::memcpy(dest, src, count * sizeof(DT_));
        }

        static void synchronize()
        {
        }

        /// Convert datatype DT2_ from src into DT1_ in dest
        template <typename DT1_, typename DT2_>
        static void convert(DT1_ * dest, const DT2_ * src, const Index count)
        {
          for (Index i(0) ; i < count ; ++i)
          {
            dest[i] = (DT1_)src[i];
          }
        }
    };

    template <>
    class MemoryPool<Mem::CUDA>
        : public InstantiationPolicy<MemoryPool<Mem::CUDA>, Singleton>
    {
      private:
        std::map<void*, Intern::MemoryInfo> _pool;

        /// default CTOR
        MemoryPool();

      public:
        ~MemoryPool();

        /// pointer to MemoryPool singleton
        friend MemoryPool* InstantiationPolicy<MemoryPool<Mem::CUDA>, Singleton>::instance();

        /// allocate new memory
        template <typename DT_>
        DT_ * allocate_memory(const Index count);

        /// increase memory counter
        void increase_memory(void * address);

        /// release memory or decrease reference counter
        void release_memory(void * address);

        /// download memory chunk to host memory
        template <typename DT_>
        static void download(DT_ * dest, const DT_ * const src, const Index count);

        /// upload memory chunk from host memory to device memory
        template <typename DT_>
        static void upload(DT_ * dest, const DT_ * const src, const Index count);

        /// recieve element
        template <typename DT_>
        static DT_ get_element(const DT_ * data, const Index index);

        /// set memory to specific value
        template <typename DT_>
        static void set_memory(DT_ * address, const DT_ val, const Index count = 1);

        /// Copy memory area from src to dest
        template <typename DT_>
        static void copy(DT_ * dest, const DT_ * src, const Index count);

        /// Convert datatype DT2_ from src into DT1_ in dest
        template <typename DT1_, typename DT2_>
        static void convert(DT1_ * dest, const DT2_ * src, const Index count);

        static void synchronize();

        static void reset_device();

        /**
         * Explicitly shut down cuda device
         *
         * This is necessary as the last user defined line of code, if cuda-memcheck is used with leak checking.
         * This includes a call to cudaResetDevice().
         *
        **/
        static void shutdown_device();
    };
  } // namespace Util
} // namespace FEAST

#endif // KERNEL_UTIL_MEMORY_POOL_HPP
