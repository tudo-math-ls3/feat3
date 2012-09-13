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
        tv_generic,
        tv_gpu,
        tv_cuda,
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
        static String name()
        {
          return "none";
        }
      };

      /**
       * Tag-type for operations in main (CPU) memory
       */
      struct CPU :
        public InstantiationPolicy<CPU, NonCopyable>
      {
        const static TagValue tag_value = tv_cpu;
        const static TagValue memory_value = tv_cpu;
        static String name()
        {
          return "cpu";
        }
      };

      /**
       * Tag-type for generic/C++-based operations.
       */
      struct Generic :
        public InstantiationPolicy<Generic, NonCopyable>
      {
        const static TagValue tag_value = tv_generic;
        const static TagValue memory_value = tv_generic;
        static String name()
        {
          return "generic";
        }
      };

      /**
       * Tag-type for operations in GPU memory
       */
      struct GPU :
        public InstantiationPolicy<GPU, NonCopyable>
      {
        const static TagValue tag_value = tv_gpu;
        const static TagValue memory_value = tv_gpu;
        static String name()
        {
          return "gpu";
        }
      };

      /**
       * Tag-type for cuda-based operations.
       */
      struct CUDA :
        public InstantiationPolicy<CUDA, NonCopyable>
      {
        const static TagValue tag_value = tv_cuda;
        const static TagValue memory_value = tv_cuda;
        static String name()
        {
          return "cuda";
        }
      };

    }

    /**
     * Output operator for tags::TagValue.
     */
    std::ostream & operator<< (std::ostream & left, Archs::TagValue value);
}
#endif
