#pragma once
#ifndef KERNEL_HORNET_ARCHS_HPP
#define KERNEL_HORNET_ARCHS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/instantiation_policy.hpp>

namespace FEAST
{
    namespace Archs
    {
      /**
       * Tag values for all supported destinations.
       */
      enum TagValue
      {
        tv_nil = 0,
        tv_cpu,
        tv_gpu_cuda,
        tv_opencl_cpu,
        tv_opencl_gpu,
        tv_none
      };

      /**
       * This is an empty tag class which may be used for templates with optional parameters.\n
       * Some template implementations might recognise the usage of a \c None parameter as <em>parameter not given</em>.
       */
      struct None :
        public InstantiationPolicy<None, NonCopyable>
      {
        const static TagValue tag_value = tv_nil;
        const static TagValue memory_value = tv_nil;
        const static String name;
      };

      /**
       * Tag-type for generic/C++-based operations.
       */
      struct CPU :
        public InstantiationPolicy<CPU, NonCopyable>
      {
        const static TagValue tag_value = tv_cpu;
        const static TagValue memory_value = tv_cpu;
        const static String name;
      };

    }

    /**
     * Output operator for tags::TagValue.
     */
    std::ostream & operator<< (std::ostream & left, Archs::TagValue value);
}
#endif
