#pragma once
#ifndef KERNEL_UTIL_MEMORY_POOL_HPP
#define KERNEL_UTIL_MEMORY_POOL_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/cuda_util.hpp>

#include <map>
#include <cstring>
#include <typeinfo>
#include <cstdio>
#include <cstddef>

#ifdef FEAT_HAVE_MKL
#include <mkl.h>
#endif

#ifdef __CUDACC__
#include "cusparse_v2.h"
#include <cublas_v2.h>
#endif

namespace FEAT
{
  namespace Util
  {
    /// \cond internal
    namespace Intern
    {
#ifdef __CUDACC__
      extern cusparseHandle_t cusparse_handle;
      extern cublasHandle_t cublas_handle;
#endif

      struct MemoryInfo
      {
        Index counter;
        Index size;
      };
    }
    /// \endcond
  } // namespace Util

    template <typename Mem_>
    class MemoryPool;

    template <>
    class MemoryPool<Mem::CUDA>
    {
      private:
        /// Map of allocated device memory patches
        static std::map<void*, Util::Intern::MemoryInfo> _pool;

      public:
        /// Setup memory pools
        static void initialise(int rank, int ranks_per_node, int ranks_per_uma, int gpus_per_node);

        /// Shutdown memory pool and clean up allocated memory pools
        static void finalise();

        /// allocate new memory
        template <typename DT_>
        static DT_ * allocate_memory(const Index count);

        /// increase memory counter
        static void increase_memory(void * address);

        /// release memory or decrease reference counter
        static void release_memory(void * address);

        /// download memory chunk to host memory
        template <typename DT_>
        static void download(DT_ * dest, const DT_ * const src, const Index count);

        /// upload memory chunk from host memory to device memory
        template <typename DT_>
        static void upload(DT_ * dest, const DT_ * const src, const Index count);

        /// receive element
        template <typename DT_>
        static DT_ get_element(const DT_ * data, const Index index);

        /// set memory to specific value
        template <typename DT_>
        static void set_memory(DT_ * address, const DT_ val, const Index count = 1);

        /// Copy memory area from src to dest
        template <typename DT_>
        static void copy(DT_ * dest, const DT_ * src, const Index count);

        /// Copy memory area from src to dest
        template <typename DT_>
        static void convert(DT_ * dest, const DT_ * src, const Index count);

        /// Convert datatype DT2_ from src into DT1_ in dest
        template <typename DT1_, typename DT2_>
        static void convert(DT1_ * dest, const DT2_ * src, const Index count);

        static void synchronise();

        static void reset_device();

        static void set_blocksize(Index misc, Index reduction, Index spmv, Index axpy);

        /// cuda threading grid blocksize for miscellaneous ops
        static Index blocksize_misc;

        /// cuda threading grid blocksize for reduction type ops
        static Index blocksize_reduction;

        /// cuda threading grid blocksize for blas-2 type ops
        static Index blocksize_spmv;

        /// cuda threading grid blocksize for blas-1 type ops
        static Index blocksize_axpy;

        /// return the overall amount of memory allocated in bytes
        static Index allocated_memory();
    };

    /**
     * \brief Memory management.
     *
     * This class manages the used memory chunks and releases them, if necessary.
     *
     * \author Dirk Ribbrock
     */
    template <>
    class MemoryPool<Mem::Main>
    {
      private:
        /// Map of all memory chunks in use.
        static std::map<void*, Util::Intern::MemoryInfo> _pool;

        /// Map of allocated pinned main memory patches
        static std::map<void*, Util::Intern::MemoryInfo> _pinned_pool;

      public:

        /// Setup memory pools
        static void initialise()
        {
        }

        /// Shutdown memory pool and clean up allocated memory pools
        static void finalise()
        {
          if (_pool.size() > 0 || _pinned_pool.size() > 0)
          {
            std::cout << stderr << " Error: MemoryPool<CPU> still contains memory chunks on deconstructor call" << std::endl;
            std::exit(1);
          }

#ifdef FEAT_HAVE_MKL
          mkl_free_buffers();
#endif
        }

        /// allocate new memory
        template <typename DT_>
        static DT_ * allocate_memory(const Index count)
        {
          XASSERT(count != 0);

          DT_ * memory(nullptr);
          memory = (DT_*)::malloc(count * sizeof(DT_));
          if (memory == nullptr)
            throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU> allocation error!");

          Util::Intern::MemoryInfo mi;
          mi.counter = 1;
          mi.size = count * sizeof(DT_);
          _pool.insert(std::pair<void*, Util::Intern::MemoryInfo>(memory, mi));

          return memory;
        }

#ifdef FEAT_HAVE_CUDA
        /// allocate new pinned memory
        template <typename DT_>
        static DT_ * allocate_pinned_memory(const Index count)
        {
          DT_ * memory(nullptr);
          if (count == 0)
            return memory;

          memory = (DT_*)Util::cuda_malloc_host(count * sizeof(DT_));

          Util::Intern::MemoryInfo mi;
          mi.counter = 1;
          mi.size = count * sizeof(DT_);
          _pinned_pool.insert(std::pair<void*, Util::Intern::MemoryInfo>(memory, mi));

          return memory;
        }
#endif

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

#ifdef FEAT_HAVE_CUDA
          it = _pinned_pool.find(address);
          if (it != _pinned_pool.end())
          {
            it->second.counter = it->second.counter + 1;
            return;
          }
#endif

          throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU>::increase_memory: Memory address not found!");
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
              ::free(address);
              _pool.erase(it);
            }
            else
            {
              it->second.counter = it->second.counter - 1;
            }
            return;
          }

#ifdef FEAT_HAVE_CUDA
          it = _pinned_pool.find(address);
          if (it != _pinned_pool.end())
          {
            if(it->second.counter == 1)
            {
              Util::cuda_free_host(address);
              _pinned_pool.erase(it);
            }
            else
            {
              it->second.counter = it->second.counter - 1;
            }
            return;
          }
#endif

          throw InternalError(__func__, __FILE__, __LINE__, "MemoryPool<CPU>::release_memory: Memory address not found!");
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

        /// receive element
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
          for (Index i(0) ; i < count ; ++i)
          {
            dest[i] = DT1_(src[i]);
          }
        }

        static void synchronise()
        {
        }

        static Index allocated_memory()
        {
          Index bytes(0);
          for (auto& i : _pool)
          {
            bytes += i.second.size;
          }
          for (auto& i : _pinned_pool)
          {
            bytes += i.second.size;
          }
          return bytes;
        }
    };

} // namespace FEAT

#endif // KERNEL_UTIL_MEMORY_POOL_HPP
