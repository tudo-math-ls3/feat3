// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_MEMORY_USAGE_HPP
#define KERNEL_UTIL_MEMORY_USAGE_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/os_windows.hpp>

#ifdef __unix__
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#endif

#include <fstream>
#include <vector>
#include <string>

namespace FEAT
{
  /**
   * \brief Memory usage info object
   *
   * This class is used to store and report memory usage statistics.
   *
   * \note The values reported by the various getter methods are all related to the moment when
   * the object was instantiated or the last time, the stamp method was called.
   * The term 'current' in the variables names has historic reasons.
   *
   * \platformswitch linux, windows and bsd all have different memory statistic tools
   *
   * \note Backport from FEAT 1 feat1/feat/kernel/arch/sysextra.c 3ccb13f633
   *       and http://stackoverflow.com/questions/1558402/memory-usage-of-current-process-in-c
   */
  class MemoryUsage
  {
  private:
    std::size_t _current_physical;
    std::size_t _peak_physical;
    std::size_t _current_virtual;
    std::size_t _peak_virtual;
    std::size_t _current_swap;

  public:
    MemoryUsage() :
      _current_physical(0),
      _peak_physical(0),
      _current_virtual(0),
      _peak_virtual(0),
      _current_swap(0)
    {
      stamp();
    }

    /// update the memory usage statistics with current data
    void stamp()
    {
#if defined(__linux)
      String line;
      std::ifstream status_file("/proc/self/status");
      if (!status_file.is_open())
        XABORTM("could not open /proc/self/status!");

      /* man 5 proc:
       *   /proc/[pid]/status
       *   Provides much of the information in /proc/[pid]/stat and
       *   /proc/[pid]/statm in a format that's easier for humans to parse.
       *   The fields are as follows:
       *     Name: Command run by this process.
       *     State: Current state of the process.  One of "R (running)", "S (sleeping)", "D (disk sleep)", "T (stopped)", "T  (trac-
       *     ing stop)", "Z (zombie)", or "X (dead)".
       *     Tgid: Thread group ID (i.e., Process ID).
       *     Pid: Thread ID (see gettid(2)).
       *     PPid: PID of parent process.
       *     TracerPid: PID of process tracing this process (0 if not being traced).
       *     Uid, Gid: Real, effective, saved set, and file system UIDs (GIDs).
       *     FDSize: Number of file descriptor slots currently allocated.
       *     Groups: Supplementary group list.
       *!    VmPeak: Peak virtual memory size.
       *!    VmSize: Virtual memory size.
       *     VmLck: Locked memory size (see mlock(3)).
       *!    VmHWM: Peak resident set size ("high water mark").
       *!    VmRSS: Resident set size.
       *     VmData, VmStk, VmExe: Size of data, stack, and text segments.
       *     VmLib: Shared library code size.
       *     VmPTE: Page table entries size (since Linux 2.6.10).
       *     Threads: Number of threads in process containing this thread.
       *     [...]
       */
      while (std::getline(status_file, line))
      {
        if (line.starts_with("VmRSS"))
        {
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_current_physical);
          _current_physical *= 1024;
          continue;
        }

        if (line.starts_with("VmHWM"))
        {
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_peak_physical);
          _peak_physical *= 1024;
          continue;
        }

        if (line.starts_with("VmSize"))
        {
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_current_virtual);
          _current_virtual *= 1024;
          continue;
        }

        if (line.starts_with("VmPeak"))
        {
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_peak_virtual);
          _peak_virtual *= 1024;
          continue;
        }

        if (line.starts_with("VmSwap"))
        {
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_current_swap);
          _current_swap *= 1024;
          continue;
        }
      }
      status_file.close();

#elif defined(_WIN32)

      unsigned long long work_set_size(0ull), work_set_size_peak(0ull), page_file_usage(0ull), page_file_usage_peak(0ull);

      Windows::query_memory_usage(work_set_size, work_set_size_peak, page_file_usage, page_file_usage_peak);

      _current_physical = std::size_t(work_set_size);
      _current_virtual = std::size_t(page_file_usage);
      _peak_physical = std::size_t(work_set_size_peak);
      _peak_virtual = std::size_t(page_file_usage_peak);
      _current_swap = std::size_t(0);

#elif defined(__unix__)
      // see:
      // https://linux.die.net/man/2/getrusage
      // https://www.freebsd.org/cgi/man.cgi?query=getrusage
      struct rusage r_usage;
      if (0 != getrusage(RUSAGE_SELF, &r_usage))
        XABORTM("Error in getrusage call!");
      _peak_physical = (std::size_t)r_usage.ru_maxrss;
      _peak_virtual = (std::size_t)r_usage.ru_maxrss;
      _current_physical = (std::size_t)r_usage.ru_maxrss;
      _current_virtual = (std::size_t)r_usage.ru_maxrss;
      _current_swap = std::size_t(0);
      // 'ru_maxrss' is given in kilobytes
      _peak_physical *= 1024u;
      _peak_virtual *= 1024u;
      _current_physical *= 1024u;
      _current_virtual *= 1024u;
#endif
    }

    /// Returns current (last stamp call) physical memory
    std::size_t get_current_physical() const
    {
      return _current_physical;
    }

    /// Returns peak physical memory
    std::size_t get_peak_physical() const
    {
      return _peak_physical;
    }

    /// Returns current (last stamp call) virtual memory
    std::size_t get_current_virtual() const
    {
      return _current_virtual;
    }

    /// Returns peak virtual memory
    std::size_t get_peak_virtual() const
    {
      return _peak_virtual;
    }

    /// Returns current (last stamp call) swap memory
    std::size_t get_current_swap() const
    {
      return _current_swap;
    }

    /// Retrieve formatted memory usage string
    String get_formatted_memory_usage() const
    {
      String r;
      r += String("Current real:").pad_back(20) + stringify(_current_physical / 1024 / 1024) + " MByte\n";
      r += String("Peak real:").pad_back(20) + stringify(_peak_physical / 1024 / 1024) + " MByte\n";
      r += String("Current virtual:").pad_back(20) + stringify(_current_virtual / 1024 / 1024) + " MByte\n";
      r += String("Peak virtual:").pad_back(20) + stringify(_peak_virtual / 1024 / 1024) + " MByte\n";
      r += String("Current swap:").pad_back(20) + stringify(_current_swap / 1024 / 1024) + " MByte\n";

      return r;
    }
  }; //class MemoryUsage
} // namespace FEAT

#endif // KERNEL_UTIL_MEMORY_USAGE_HPP
