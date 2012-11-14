#pragma once
#ifndef KERNEL_HORNET_ARCHS_HPP
#define KERNEL_HORNET_ARCHS_HPP 1

/**
 * \file
 * \brief FEAST Architecture header.
 *
 * This file contains all architecture definitions, used in FEAST.
 */

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
        tv_none = 0,
        tv_serial,
        tv_parallel
      };

      /**
       * This is an empty tag class which may be used for templates with optional parameters.\n
       * Some template implementations might recognise the usage of a \c None parameter as <em>parameter not given</em>.
       */
      struct None :
        public InstantiationPolicy<None, NonCopyable>
      {
        const static TagValue tag_value = tv_none;
        const static TagValue memory_value = tv_none;
        static String name()
        {
          return "none";
        }
      };

      /**
       * Tag-type for serial operations.
       */
      struct Serial :
        public InstantiationPolicy<Serial, NonCopyable>
      {
        const static TagValue tag_value = tv_serial;
        const static TagValue memory_value = tv_serial;
        static String name()
        {
          return "serial";
        }
      };

      /**
       * Tag-type for Parallel operations.
       */
      struct Parallel :
        public InstantiationPolicy<Parallel, NonCopyable>
      {
        const static TagValue tag_value = tv_parallel;
        const static TagValue memory_value = tv_parallel;
        static String name()
        {
          return "parallel";
        }
      };
    }


    namespace Mem
    {
      /**
       * Tag values for all supported destinations.
       */
      enum TagValue
      {
        tv_main,
        tv_cuda,
      };

      /**
       * Tag-type for operations in main (CPU) memory
       */
      struct Main :
        public InstantiationPolicy<Main, NonCopyable>
      {
        const static TagValue tag_value = tv_main;
        const static TagValue memory_value = tv_main;
        static String name()
        {
          return "main";
        }
      };

      /**
       * Tag-type for operations in cuda (GPU) memory
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

    namespace Algo
    {
      /**
       * Tag values for all supported destinations.
       */
      enum TagValue
      {
        tv_generic,
        tv_cuda,
      };
      /**
       * Tag-type for generic/C++-based operations.
       */
      struct Generic :
        public InstantiationPolicy<Generic, NonCopyable>
      {
        typedef Mem::Main mem_type;
        const static TagValue tag_value = tv_generic;
        const static TagValue memory_value = tv_generic;
        static String name()
        {
          return "generic";
        }
      };

      /**
       * Tag-type for cuda-based operations.
       */
      struct CUDA :
        public InstantiationPolicy<CUDA, NonCopyable>
      {
        typedef Mem::CUDA mem_type;
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
    std::ostream & operator<< (std::ostream & left, Mem::TagValue value);
    std::ostream & operator<< (std::ostream & left, Algo::TagValue value);
}
#endif
