#pragma once
#include <kernel/base_header.hpp>

#ifndef __CUDACC__
  static_assert(false, "Header should only be included from cuda context");
#endif

namespace FEAT
{

  #define cudaCheck(ans, funcStr) if((ans) != cudaSuccess) throw InternalError(__func__, __FILE__, __LINE__, (funcStr) + String(" with error: \n") + stringify(cudaGetErrorString((ans))))

  #define cudaCheckLastErr() { \
    cudaError_t last_error(cudaGetLastError()); \
    if (cudaSuccess != last_error) \
      throw InternalError(__func__, __FILE__, __LINE__, "CUDA error occurred in execution!\n" + stringify(cudaGetErrorString(last_error))); \
  }


  template<typename DT_>
  CUDA_HOST_DEVICE DT_* get_aligned_address(void* ptr)
  {
    // alignment is always a power of 2, so can use the following bit magic
    // first add the (alignment-1) of the datatype, which leads to the next aligned
    // address modulo bits that fall into the alignment.
    // The simply cut off everything from the address larger then our aligned address by using
    // AND with the negative alignment, which is simply all bits flipped of (alignment-1)
    constexpr std::uint64_t align = alignof(DT_);
    std::uint64_t ptr_a = reinterpret_cast<std::uint64_t>(ptr);
    return (DT_*)((ptr_a + align - 1u)&(-align));
  }


  /// reduce shared array across block size arrays so the final value is inside shared_array[0]
  template<typename DT_, int blockSize_>
  CUDA_DEVICE void inline static reduce_sum(DT_* shared_array, int thread_idx)
  {
    __syncthreads();
    if constexpr(blockSize_ >= 1024)
    {
      if(thread_idx < 512)
        shared_array[thread_idx] += shared_array[thread_idx + 512];
      __syncthreads();
    }
    if constexpr(blockSize_ >= 512)
    {
      if(thread_idx < 256)
        shared_array[thread_idx] += shared_array[thread_idx + 256];
      __syncthreads();
    }
    if constexpr(blockSize_ >= 256)
    {
      if(thread_idx < 128)
        shared_array[thread_idx] += shared_array[thread_idx + 128];
      __syncthreads();
    }
    if constexpr(blockSize_ >= 128)
    {
      if(thread_idx < 64)
        shared_array[thread_idx] += shared_array[thread_idx + 64];
      __syncthreads();
    }

    // warp internal sync does not require explicit synchonization
    // from here, we always assume that we have at least a block size of 64, so these adds are safe
    if constexpr(blockSize_ >= 64)
    {
      if(thread_idx < 32) shared_array[thread_idx] += shared_array[thread_idx + 32];
    }
    // final warp is done by warp shuffle intrinsics, else we would need to sync after each shuffle
    if(thread_idx < 32)
    {
      shared_array[thread_idx] += __shfl_down_sync(0xFFFFFFFF, shared_array[thread_idx], 16);
      shared_array[thread_idx] += __shfl_down_sync(0xFFFFFFFF, shared_array[thread_idx], 8);
      shared_array[thread_idx] += __shfl_down_sync(0xFFFFFFFF, shared_array[thread_idx], 4);
      shared_array[thread_idx] += __shfl_down_sync(0xFFFFFFFF, shared_array[thread_idx], 2);
      shared_array[thread_idx] += __shfl_down_sync(0xFFFFFFFF, shared_array[thread_idx], 1);
    }
  }
}