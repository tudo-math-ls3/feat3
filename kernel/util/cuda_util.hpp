#pragma once
#ifndef KERNEL_UTIL_CUDA_UTIL_HPP
#define KERNEL_UTIL_CUDA_UTIL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace Util
  {
    void cuda_set_device(int device);
  }
}

#endif //KERNEL_UTIL_CUDA_UTIL_HPP
