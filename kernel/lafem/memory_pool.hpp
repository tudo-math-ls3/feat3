#pragma once
#ifndef KERNEL_LAFEM_MEMORY_POOL_HPP
#define KERNEL_LAFEM_MEMORY_POOL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/archs.hpp>

#ifdef FEAST_GMP
#include <gmpxx.h>
#include <mpfr.h>
#endif

#include <map>
#include <cstring>
#include <typeinfo>


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
        template <typename DT_>
        DT_ * allocate_memory(const Index count);

        /// increase memory counter
        void increase_memory(void * address);

        /// release memory or decrease reference counter
        void release_memory(void * address);

        /// download memory chunk to host memory
        template <typename DT_>
        static void download(void * dest, void * src, const Index count);

        /// upload memory chunk from host memory to device memory
        template <typename DT_>
        static void upload(void * dest, void * src, const Index count);

        /// recieve element
        template <typename DT_>
        inline static const DT_ & get_element(const DT_ * data, const Index index)
        {
          return data[index];
        }

        /// set memory to specific value
        template <typename DT_>
        static void set_memory(DT_ * address, const DT_ val, const Index count = 1);

        /// Copy memory area from src to dest
        template <typename DT_>
        static void copy(void * dest, const void * src, const Index count)
        {
          if (dest == src)
            return;

#ifdef FEAST_GMP
          if (typeid(DT_) == typeid(mpf_class))
          {
            const DT_ * s((DT_ *) src);
            for (DT_ * d((DT_ * )dest), * d_end((DT_ *)dest + count) ; d != d_end ; ++d, ++s)
            {
              *d = *s;
            }
          }
          else
#endif
            ::memcpy(dest, src, count * sizeof(DT_));
        }

        /// Generate hash value for given byte sequence
        unsigned long generate_hash(void * address, const Index bytes);
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
        void * allocate_memory(const Index count);

        /// increase memory counter
        void increase_memory(void * address);

        /// release memory or decrease reference counter
        void release_memory(void * address);

        /// download memory chunk to host memory
        template <typename DT_>
        static void download(void * dest, void * src, const Index count);

        /// upload memory chunk from host memory to device memory
        template <typename DT_>
        static void upload(void * dest, void * src, const Index count);

        /// recieve element
        template <typename DT_>
        static DT_ get_element(const DT_ * data, const Index index);

        /// set memory to specific value
        template <typename DT_>
        static void set_memory(DT_ * address, const DT_ val, const Index count = 1);

        /// Copy memory area from src to dest
        template <typename DT_>
        static void copy(void * dest, const void * src, const Index count);

        /// Generate hash value for given byte sequence
        unsigned long generate_hash(void * address, const Index bytes);
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_MEMORY_POOL_HPP
