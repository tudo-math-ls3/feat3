#pragma once
#ifndef KERNEL_UTIL_CUDA_UTIL_HPP
#define KERNEL_UTIL_CUDA_UTIL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace Util
  {
    void cuda_set_device(const int device);
    void * cuda_malloc_host(const Index bytes);
    void cuda_free_host(void * address);
  }
}

#endif //KERNEL_UTIL_CUDA_UTIL_HPP
