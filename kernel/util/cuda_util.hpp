#pragma once
#ifndef KERNEL_UTIL_CUDA_UTIL_HPP
#define KERNEL_UTIL_CUDA_UTIL_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>

namespace FEAT
{
  namespace Util
  {
    void cuda_set_device(const int device);
    void * cuda_malloc_host(const Index bytes);
    void cuda_free_host(void * address);
    void cuda_check_last_error();
    void * cuda_get_device_pointer(void * host);
  }
}

#endif //KERNEL_UTIL_CUDA_UTIL_HPP
