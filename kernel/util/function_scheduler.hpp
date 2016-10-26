#pragma once
#ifndef UTIL_GUARD_FUNCTION_SCHEDULER_HPP
#define UTIL_GUARD_FUNCTION_SCHEDULER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

#include <functional>

namespace FEAT
{
  /// Utility collection
  namespace Util
  {
    /**
     * Supported Schedule modes.
     */
    enum class ScheduleMode
    {
      sequential, /**< Execute function on ranks sequential */
      clustered /**< Execute function on a bunch of ranks at once */
    };

    /// Execute a function on all ranks, scheduled by the given mode
    void schedule_function(std::function<void (void)> func, ScheduleMode mode);

  } // namespace Util
} // namespace FEAT


#endif // UTIL_GUARD_FUNCTION_SCHEDULER_HPP
