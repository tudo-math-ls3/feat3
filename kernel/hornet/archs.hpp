#pragma once
#ifndef KERNEL_HORNET_ARCHS_HPP
#define KERNEL_HORNET_ARCHS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/instantiation_policy.hpp>

namespace FEAST
{
  namespace HOrNET
  {
    namespace Archs
    {
      /**
       * Tag values for all supported destinations.
       */
      enum TagValue
      {
        tv_cpu = 1,
        tv_gpu_cuda,
        tv_opencl_cpu,
        tv_opencl_gpu,
        tv_none
      };

      /**
       * Tag-type for generic/C++-based operations.
       *
       * \ingroup grptagscpu
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
}
#endif
