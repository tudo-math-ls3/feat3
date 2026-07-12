// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifdef FEAT_HAVE_CUDA

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>

#ifdef __CUDACC__
#include <cusparse_v2.h>
#include <cublas_v2.h>
#include <cublasLt.h>
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
      extern cublasLtMatmulAlgo_t * cublas_lt_algo_matmat;
      extern bool * cublas_lt_algo_matmat_initialized;
      extern size_t cuda_workspace_size;
      extern void * cuda_workspace;
      /// cache whether we we support stream ordered operations
      /// by default false, runtime has to query this the first time
      extern bool stream_ordered_supported;
#endif
    }
    /// \endcond

    /// The device number this process uses
    extern int cuda_device_number;

    /// cuda threading grid blocksize for miscellaneous ops
    extern Index cuda_blocksize_misc;

    /// cuda threading grid blocksize for reduction type ops
    extern Index cuda_blocksize_reduction;

    /// cuda threading grid blocksize for blas-2 type ops
    extern Index cuda_blocksize_spmv;

    /// cuda threading grid blocksize for blas-1 type ops
    extern Index cuda_blocksize_axpy;

    /// cuda threading grid blocksize for scalar typed assembly
    extern Index cuda_blocksize_scalar_assembly;

    /// cuda threading grid blocksize for block typed assembly
    extern Index cuda_blocksize_blocked_assembly;

    /// cuda threading grid blocksize for block typed assembly
    extern Index cuda_blocksize_vanka_assembly;

    /// helper for typicall runtime parameter flags
    struct RuntimeParameter
    {
      /// how much shared mem should be used
      int shared_mem;
      /// 1d blocksize that should be used
      int blocksize;
      /// 1d gridsize that should be used, i.e. overall start block_size*grid_size threads
      int gridsize;
      /// print occupancy of the cuda kernel
      bool print_occupancy;
      /// should the function force a synchronization at the end?
      bool synchronize;
    };

    /// Wrapper for stream parameter
    struct StreamWrapper
    {
      /// default stream is zero
      void* stream = nullptr;
      /// default stream is zero
      std::uint32_t stream_idx = 0;
    };

    void cuda_set_blocksize(Index misc, Index reduction, Index spmv, Index axpy, Index scalar_assembly, Index blocked_assembly);
    void cuda_set_device(const int device);
    void cuda_check_last_error();
    //void * cuda_get_device_pointer(void * host);
    //void * cuda_malloc_managed(const Index bytes);
    //void * cuda_malloc(const Index bytes);
    //void * cuda_malloc_host(const Index bytes);
    void * cuda_get_static_memory(const Index bytes);
    //void cuda_free(void * address);
    //void cuda_free_host(void * address);
    void cuda_free_static_memory();
    void cuda_initialize(int rank, int ranks_per_node, int ranks_per_uma, int gpus_per_node);
    void cuda_finalize();
    /**
     * \brief Registers a callback to be invoked during cuda_finalize().
     *
     * This lets higher-level modules (e.g. LAFEM::Arch) hook the CUDA teardown to release their
     * own device resources without the low-level util layer having to depend on them (which would
     * create a circular link dependency). Callbacks are invoked in registration order at the very
     * start of cuda_finalize(), while the CUDA context is still alive.
     */
    void cuda_register_finalize_callback(void (*callback)());
    NOINLINE void cuda_synchronize();
    NOINLINE void cuda_force_synchronize();
    void cuda_reset_device();
    //void cuda_copy(void * dest, const void * src, const Index bytes);
    //void cuda_copy_host_to_device(void * dest, const void * src, const Index bytes);
    //void cuda_copy_device_to_host(void * dest, const void * src, const Index bytes);
    //void cuda_copy_device_to_device(void * dest, const void * src, const Index bytes);
    void cuda_reset_algos();
    //template <typename DT_>
    //void cuda_set_memory(DT_ * address, const DT_ val, const Index count);
    //template <typename DT1_, typename DT2_>
    //void cuda_convert(DT1_ * dest, const DT2_ * src, const Index count);
    int cuda_get_device_count();
    int cuda_get_device_id();
    String cuda_get_visible_devices();
    std::size_t cuda_get_max_cache_thread();
    void cuda_set_max_cache_thread(const std::size_t bytes);
    void cuda_start_profiling();
    void cuda_stop_profiling();
    std::size_t cuda_get_shared_mem_per_sm();
    std::size_t cuda_get_max_blocks_per_sm();
    std::size_t cuda_get_sm_count();

    /// returns true if the current device supports stream-ordered memory allocation (cudaMallocAsync/cudaFreeAsync)
    /// only valid after init call for internal variable (which is done in cuda_init)
    bool cuda_get_stream_ordered_alloc_supported();

    #ifdef __CUDACC__
    /// returns the number of active blocks of size blocksize for a given cuda kernel T
    /// for this to work, you have to be in teh compile unit of the targeted kernel, i.e.
    /// in the .cu file itself
    template<typename T>
    inline int cuda_get_occupancy(T kernel_func, int blocksize, int shared_memory = 0)
    {
      int num_blocks = 0;
      if(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&num_blocks, kernel_func, blocksize, shared_memory) != cudaSuccess)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "cudaOccupancyMaxActiveBlockPerMultiprocessor failed!");
      }
      return num_blocks;
    }
    #endif

    std::uint64_t cuda_start_range(const char* name);

    void cuda_stop_range(std::uint64_t id);

    void cuda_mark_event(const char* name);
  }
}
#endif // FEAT_HAVE_CUDA

// The following NVTX profiling helpers are always compiled (also without CUDA). When FEAT_HAVE_CUDA
// is not defined they degrade to no-ops, so call sites may use them unconditionally without guards.
#include <cstddef>
#include <cstdint>
#include <vector>

namespace FEAT
{
  namespace Util
  {
    /**
     * \brief Stack manager for nested NVTX profiling ranges.
     *
     * Abstracts the paired \c Util::cuda_start_range / \c Util::cuda_stop_range calls behind a
     * push/pop stack: push() opens a named range and stores its id, pop() closes the most recently
     * opened one (LIFO). The destructor closes any ranges still open, so an early return or a thrown
     * exception cannot leak a dangling range.
     *
     * All operations are no-ops unless \c FEAT_HAVE_CUDA is defined.
     *
     * \author Maximilian Esser
     */
    class NvtxRangeStack
    {
    private:
      /// the ids of the currently open ranges (LIFO order)
      std::vector<std::uint64_t> _ids;

    public:
      NvtxRangeStack() = default;

      /// stacks must not be copied (ownership of open range ids is unique)
      NvtxRangeStack(const NvtxRangeStack&) = delete;
      NvtxRangeStack& operator=(const NvtxRangeStack&) = delete;

      /// closes all ranges that are still open
      ~NvtxRangeStack()
      {
        while(!_ids.empty())
          pop();
      }

      /// opens a new named range and pushes it onto the stack
      void push([[maybe_unused]] const char* name)
      {
#ifdef FEAT_HAVE_CUDA
        _ids.push_back(cuda_start_range(name));
#endif
      }

      /// closes the most recently opened range (no-op if the stack is empty)
      void pop()
      {
#ifdef FEAT_HAVE_CUDA
        if(!_ids.empty())
        {
          cuda_stop_range(_ids.back());
          _ids.pop_back();
        }
#endif
      }

      /// returns the number of currently open ranges
      std::size_t size() const
      {
        return _ids.size();
      }
    };

    /// returns the calling thread's default NvtxRangeStack (used by the single-argument ScopedNvtxRange ctor)
    inline NvtxRangeStack& nvtx_thread_stack()
    {
      thread_local NvtxRangeStack stack;
      return stack;
    }

    /**
     * \brief RAII scope guard for an NVTX profiling range.
     *
     * Built on top of NvtxRangeStack: opens a range on construction and closes it on destruction,
     * so wrapping a scope in a range is a single line. Uses the calling thread's default stack
     * unless an explicit stack is supplied.
     *
     * All operations are no-ops unless \c FEAT_HAVE_CUDA is defined.
     *
     * \author Maximilian Esser
     */
    class ScopedNvtxRange
    {
    private:
      /// the stack this range was pushed onto
      NvtxRangeStack& _stack;

    public:
      /// opens a range \p name on the calling thread's default stack
      explicit ScopedNvtxRange(const char* name) :
        _stack(nvtx_thread_stack())
      {
        _stack.push(name);
      }

      /// opens a range \p name on the given \p stack
      ScopedNvtxRange(NvtxRangeStack& stack, const char* name) :
        _stack(stack)
      {
        _stack.push(name);
      }

      /// scope guards must not be copied or moved
      ScopedNvtxRange(const ScopedNvtxRange&) = delete;
      ScopedNvtxRange& operator=(const ScopedNvtxRange&) = delete;

      /// closes the range
      ~ScopedNvtxRange()
      {
        _stack.pop();
      }
    };

    /// emits an instantaneous NVTX marker event (no-op unless FEAT_HAVE_CUDA is defined)
    inline void nvtx_mark([[maybe_unused]] const char* name)
    {
#ifdef FEAT_HAVE_CUDA
      cuda_mark_event(name);
#endif
    }
  } // namespace Util
} // namespace FEAT
