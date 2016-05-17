#pragma once
#ifndef KERNEL_UTIL_MEMORY_USAGE_HPP
#define KERNEL_UTIL_MEMORY_USAGE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>

namespace FEAST
{
  namespace Util
  {
    struct MemUsageInfo
    {
      std::size_t current_physical;
      std::size_t peak_physical;
      std::size_t current_virtual;
      std::size_t peak_virtual;
      std::size_t current_swap;

      MemUsageInfo()
      {
        current_physical = 0;
        peak_physical = 0;
        current_virtual = 0;
        peak_virtual = 0;
        current_swap = 0;
      }
    };

    /**
     * \brief Report current and peak memory usage.
     *
     * Reports current and peak memory usage as reported by the operating system.
     *
     * \return a std::tuple containing the following values in bytes:
     *         * current physical memory
     *         * peak physical memory
     *         * current virtual memory
     *         * peak virtual memory
     *         * current swap memory
     *
     * \note This method only works on *ix systems, by parsing /proc/self/status.
     *
     * \note Backport from FEAST 1 feast1/feast/kernel/arch/sysextra.c 3ccb13f633
     *       and http://stackoverflow.com/questions/1558402/memory-usage-of-current-process-in-c
     *
     */
    MemUsageInfo get_memory_usage();

    /// Retrieve formated memory usage string
    String get_formated_memory_usage();
  } // namespace Util
} // namespace FEAST

#endif // KERNEL_UTIL_MEMORY_USAGE_HPP
