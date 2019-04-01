// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
