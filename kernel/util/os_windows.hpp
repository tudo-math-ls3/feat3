#pragma once
#ifndef KERNEL_UTIL_OS_WINDOWS_HPP
#define KERNEL_UTIL_OS_WINDOWS_HPP
#if defined(_WIN32) || defined(DOXYGEN)

#include <kernel/base_header.hpp>

namespace FEAT
{
  /**
   * \brief Windows OS utility namespace
   *
   * This namespace encapsulates a set of functions which are wrapped from the Win32 API,
   * and provides some additional functions for debugging and profiling.
   */
  namespace Windows
  {
    /**
     * \brief Wraps around the QueryPerformanceCounter function.
     *
     * \see https://msdn.microsoft.com/en-us/library/windows/desktop/ms644904.aspx
     *
     * \note This function is used by the TimeStamp class.
     *
     * \author Peter Zajac
     */
    long long query_performance_counter();

    /**
     * \brief Wraps around the QueryPerformanceCounter function.
     *
     * \see https://msdn.microsoft.com/en-us/library/windows/desktop/ms644905.aspx
     *
     * \note This function is used by the TimeStamp class.
     *
     * \author Peter Zajac
     */
    long long query_performance_frequency();

    /**
     * \brief Queries memory usage information
     *
     * \see https://msdn.microsoft.com/en-us/library/windows/desktop/aa965225.aspx
     *
     * \note This function is used by the Util::get_memory_usage() function.
     *
     * \author Peter Zajac
     */
    void query_memory_usage(
      unsigned long long& work_set_size,
      unsigned long long& work_set_size_peak,
      unsigned long long& page_file_usage,
      unsigned long long& page_file_usage_peak);

    /**
     * \brief Dumps the call-stack to stderr
     *
     * \author Peter Zajac
     */
    void dump_call_stack_to_stderr();

    /**
     * \brief Installs custom Windows Structured-Exception-Handler filter
     *
     * This function installs a handler for the Windows SEH system, which
     * prints an error message and a call-stack dump to stderr.
     *
     * \note This function is called by the automated regression test system.
     *
     * \author Peter Zajac
     */
    void install_seh_filter();

    /**
     * \brief Disables Windows error dialog boxes
     *
     * This functions disables the Popup-Dialog-Boxes when an application crashes.
     *
     * \note This function is called by the automated regression test system.
     *
     * \author Peter Zajac
     */
    void disable_error_prompts();
  } // namespace Windows
} // namespace FEAT

#endif // defined(_WIN32) || defined(DOXYGEN)
#endif // KERNEL_UTIL_OS_WINDOWS_HPP
