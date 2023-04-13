// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_MEMORY_POOL_HPP
#define KERNEL_UTIL_MEMORY_POOL_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/runtime.hpp>

#include <map>
#include <cstring>
#include <typeinfo>
#include <cstdio>
#include <cstddef>

#ifdef FEAT_HAVE_MKL
#include <mkl.h>
#endif


namespace FEAT
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
  } // namespace Util

    /**
     * \brief Memory management.
     *
     * This class manages the used memory chunks and releases them, if necessary.
     *
     * \author Dirk Ribbrock
     */
    class MemoryPool
    {
      private:
        /// Map of all memory chunks in use.
        static std::map<void*, Util::Intern::MemoryInfo> _pool;

      public:

        /// Setup memory pools
        static void initialize()
        {
        }

        /// Shutdown memory pool and clean up allocated memory pools
        static void finalize()
        {
          if (_pool.size() > 0)
          {
            std::cout << stderr << " Error: MemoryPool still contains memory chunks on deconstructor call" << std::endl;
            std::exit(1);
          }

#ifdef FEAT_HAVE_MKL
          mkl_free_buffers();
#endif
        }

        /// allocate new memory
        template <typename DT_>
        static DT_ * allocate_memory(Index count)
        {
          DT_ * memory(nullptr);

          if (count == 0)
            return memory;

          if (count%4 != 0)
            count = count + (4ul - count%4);

#ifdef FEAT_HAVE_CUDA
          memory = (DT_*)Util::cuda_malloc_managed(count * sizeof(DT_));
#else
          memory = (DT_*)::malloc(count * sizeof(DT_));
#endif
          if (memory == nullptr)
            XABORTM("MemoryPool allocation error!");

          Util::Intern::MemoryInfo mi;
          mi.counter = 1;
          mi.size = count * sizeof(DT_);
          _pool.insert(std::pair<void*, Util::Intern::MemoryInfo>(memory, mi));

          return memory;
        }

        /// increase memory counter
        static void increase_memory(void * address)
        {
          XASSERT(address != nullptr);

          std::map<void*, Util::Intern::MemoryInfo>::iterator it(_pool.find(address));
          if (it != _pool.end())
          {
            it->second.counter = it->second.counter + 1;
            return;
          }

          XABORTM("MemoryPool::increase_memory: Memory address not found!");
        }

        /// release memory or decrease reference counter
        static void release_memory(void * address)
        {
          if (address == nullptr)
            return;

          std::map<void*, Util::Intern::MemoryInfo>::iterator it(_pool.find(address));
          if (it != _pool.end())
          {
            if(it->second.counter == 1)
            {
#ifdef FEAT_HAVE_CUDA
              Util::cuda_free(address);
#else
              ::free(address);
#endif
              _pool.erase(it);
            }
            else
            {
              it->second.counter = it->second.counter - 1;
            }
            return;
          }

          XABORTM("MemoryPool::release_memory: Memory address not found!");
        }

        /// download memory chunk to host memory
        template <typename DT_>
        [[deprecated("no download necessary in unified memory environment.")]]
        inline static void download(DT_ * dest, const DT_ * const src, const Index count)
        {
          if (dest == src)
            return;

          ::memcpy(dest, src, count * sizeof(DT_));
        }

        /// upload memory chunk from host memory to device memory
        template <typename DT_>
        [[deprecated("no upload necessary in unified memory environment.")]]
        inline static void upload(DT_ * dest, const DT_ * const src, const Index count)
        {
          if (dest == src)
            return;

          ::memcpy(dest, src, count * sizeof(DT_));
        }

        /// receive element
        template <typename DT_>
        [[deprecated("no get_element necessary in unified memory environment.")]]
        inline static const DT_ & get_element(const DT_ * data, const Index index)
        {
          return data[index];
        }

        /// set memory to specific value
        template <typename DT_>
        static void set_memory(DT_ * address, const DT_ val, const Index count = 1)
        {
          switch(Runtime::get_preferred_backend())
          {
#ifdef FEAT_HAVE_CUDA
            case PreferredBackend::cuda:
              FEAT::Util::cuda_set_memory(address, val, count);
              return;
#endif
            case PreferredBackend::generic:
            default:
              for (Index i(0) ; i < count ; ++i)
              {
                address[i] = val;
              }
              return;
          }
        }

        /// Copy memory area from src to dest
        template <typename DT_>
        static void copy(DT_ * dest, const DT_ * src, const Index count)
        {
          if (dest == src)
            return;

          switch(Runtime::get_preferred_backend())
          {
#ifdef FEAT_HAVE_CUDA
            case PreferredBackend::cuda:
              FEAT::Util::cuda_copy(dest, src, count * sizeof(DT_));
              return;
#endif
            case PreferredBackend::generic:
            default:
              ::memcpy(dest, src, count * sizeof(DT_));
              return;
          }
        }

        /// Copy memory area from src to dest
        template <typename DT_>
        static void convert(DT_ * dest, const DT_ * src, const Index count)
        {
          if (dest == src)
            return;

          ::memcpy(dest, src, count * sizeof(DT_));
        }

        /// Convert datatype DT2_ from src into DT1_ in dest
        template <typename DT1_, typename DT2_>
        static void convert(DT1_ * dest, const DT2_ * src, const Index count)
        {

          switch(Runtime::get_preferred_backend())
          {
#ifdef FEAT_HAVE_CUDA
            case PreferredBackend::cuda:
              FEAT::Util::cuda_convert(dest, src, count);
              return;
#endif
            case PreferredBackend::generic:
            default:
              for (Index i(0) ; i < count ; ++i)
              {
                dest[i] = DT1_(src[i]);
              }
              return;
          }
        }

        NOINLINE static void synchronize()
        {
#ifdef FEAT_HAVE_CUDA
          FEAT::Util::cuda_synchronize();
#endif
        }

        static Index allocated_memory()
        {
          Index bytes(0);
          for (auto& i : _pool)
          {
            bytes += i.second.size;
          }
          return bytes;
        }

        static Index allocated_size(void * address)
        {
          std::map<void*, Util::Intern::MemoryInfo>::iterator it(_pool.find(address));
          if (it != _pool.end())
          {
            return it->second.size;
          }
          else
            XABORTM("MemoryPool::allocated_size: Memory address not found!");
        }
    };

} // namespace FEAT

#endif // KERNEL_UTIL_MEMORY_POOL_HPP
