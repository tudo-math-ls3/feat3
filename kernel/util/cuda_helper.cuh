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
}