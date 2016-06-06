#pragma once
#ifndef KERNEL_HORNET_ARCHS_HPP
#define KERNEL_HORNET_ARCHS_HPP 1

/**
 * \file
 * \brief FEAT Architecture header.
 *
 * This file contains all architecture definitions, used in FEAT.
 */

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>

namespace FEAT
{
    /**
     * \brief FEAT::Archs namespace
     */
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
      struct None
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
      struct Serial
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
      struct Parallel
      {
        const static TagValue tag_value = tv_parallel;
        const static TagValue memory_value = tv_parallel;
        static String name()
        {
          return "parallel";
        }
      };
    }


    /**
     * \brief FEAT::Mem namespace
     */
    namespace Mem
    {
      /**
       * Tag values for all supported destinations.
       */
      enum TagValue
      {
        tv_main,
        tv_cuda
      };

      /**
       * Tag-type for operations in main (CPU) memory
       */
      struct Main
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
      struct CUDA
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
     * \brief FEAT::Algo namespace
     *
     * \deprecated {The Algo Tag is not used anymore.
     * Only the micro benchmarks are stuck to this to explicitly address the several backends.}
     */
    namespace Algo
    {
      /**
       * Tag values for all supported destinations.
       */
      enum TagValue
      {
        tv_generic,
        tv_mkl,
        tv_cuda
      };

      /**
       * Tag-type for generic/C++-based operations.
       *
       * \deprecated {The Algo Tag is not used anymore.
       * Only the micro benchmarks are stuck to this to explicitly address the several backends.}
       */
      struct Generic
      {
        typedef Mem::Main MemType;
        const static TagValue tag_value = tv_generic;
        const static Mem::TagValue memory_value = Mem::tv_main;
        static String name()
        {
          return "generic";
        }
      };

      /**
       * Tag-type for MKL-based operations.
       *
       * \deprecated {The Algo Tag is not used anymore.
       * Only the micro benchmarks are stuck to this to explicitly address the several backends.}
       */
      struct MKL
      {
        typedef Mem::Main MemType;
        const static TagValue tag_value = tv_mkl;
        const static Mem::TagValue memory_value = Mem::tv_main;
        static String name()
        {
          return "mkl";
        }
      };

      /**
       * Tag-type for cuda-based operations.
       *
       * \deprecated {The Algo Tag is not used anymore.
       * Only the micro benchmarks are stuck to this to explicitly address the several backends.}
       */
      struct CUDA
      {
        typedef Mem::CUDA MemType;
        const static TagValue tag_value = tv_cuda;
        const static Mem::TagValue memory_value = Mem::tv_cuda;
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
