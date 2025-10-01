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

    void cuda_set_blocksize(Index misc, Index reduction, Index spmv, Index axpy, Index scalar_assembly, Index blocked_assembly);
    void cuda_set_device(const int device);
    void cuda_check_last_error();
    void * cuda_get_device_pointer(void * host);
    void * cuda_malloc_managed(const Index bytes);
    void * cuda_malloc(const Index bytes);
    void * cuda_malloc_host(const Index bytes);
    void * cuda_get_static_memory(const Index bytes);
    void cuda_free(void * address);
    void cuda_free_host(void * address);
    void cuda_free_static_memory();
    void cuda_initialize(int rank, int ranks_per_node, int ranks_per_uma, int gpus_per_node);
    void cuda_finalize();
    NOINLINE void cuda_synchronize();
    NOINLINE void cuda_force_synchronize();
    void cuda_reset_device();
    void cuda_copy(void * dest, const void * src, const Index bytes);
    void cuda_copy_host_to_device(void * dest, const void * src, const Index bytes);
    void cuda_copy_device_to_host(void * dest, const void * src, const Index bytes);
    void cuda_copy_device_to_device(void * dest, const void * src, const Index bytes);
    void cuda_reset_algos();
    template <typename DT_>
    void cuda_set_memory(DT_ * address, const DT_ val, const Index count);
    template <typename DT1_, typename DT2_>
    void cuda_convert(DT1_ * dest, const DT2_ * src, const Index count);
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
  }
}
#endif // FEAT_HAVE_CUDA
