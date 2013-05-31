#pragma once
#ifndef KERNEL_LAFEM_MEMORY_POOL_HPP
#define KERNEL_LAFEM_MEMORY_POOL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/archs.hpp>

#include <map>
#include <cstring>


namespace FEAST
{
  namespace LAFEM
  {
    namespace Intern
    {
      struct MemoryInfo
      {
        Index counter;
        Index size;
      };
    }

    template <typename Arch_>
    class MemoryPool
        : public InstantiationPolicy<MemoryPool<Arch_>, Singleton>
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
        MemoryPool();

      public:
        ~MemoryPool();

        /// pointer to MemoryPool singleton
        friend MemoryPool* InstantiationPolicy<MemoryPool<Mem::Main>, Singleton>::instance();

        /// allocate new memory
        void * allocate_memory(Index bytes);

        /// increase memory counter
        void increase_memory(void * address);

        /// release memory or decrease reference counter
        void release_memory(void * address);

        /// download memory chunk to host memory
        static void download(void * dest, void * src, Index bytes);

        /// upload memory chunk from host memory to device memory
        static void upload(void * dest, void * src, Index bytes);

        /// recieve element
        template <typename DT_>
        inline static const DT_ & get_element(const DT_ * data, Index index)
        {
          return data[index];
        }

        /// recieve const element
        template <typename DT_>
        inline static DT_ & get_element(DT_ * data, Index index)
        {
          return data[index];
        }

        /// modify element
        template <typename DT_>
        inline static void modify_element(DT_ * data, Index index, DT_ value)
        {
          data[index] = value;
        }

        /// set memory to specific value
        template <typename DT_>
        void set_memory(DT_ * address, const DT_ val, const Index count);

        /// Copy memory area from src to dest
        static void copy(void * dest, const void * src, const Index bytes)
        {
          if (dest == src)
            return;

          ::memcpy(dest, src, bytes);
        }

        /// Generate hash value for given byte sequence
        unsigned long generate_hash(void * address, Index bytes);
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
        void * allocate_memory(Index bytes);

        /// increase memory counter
        void increase_memory(void * address);

        /// release memory or decrease reference counter
        void release_memory(void * address);

        /// download memory chunk to host memory
        static void download(void * dest, void * src, Index bytes);

        /// upload memory chunk from host memory to device memory
        static void upload(void * dest, void * src, Index bytes);

        /// recieve element
        template <typename DT_>
        static DT_ get_element(const DT_ * data, Index index);

        /// modify element
        template <typename DT_>
        static void modify_element(DT_ * data, Index index, DT_ value);

        /// set memory to specific value
        template <typename DT_>
        void set_memory(DT_ * address, const DT_ val, const Index count);

        /// Copy memory area from src to dest
        static void copy(void * dest, const void * src, const Index bytes);

        /// Generate hash value for given byte sequence
        unsigned long generate_hash(void * address, Index bytes);
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_MEMORY_POOL_HPP
