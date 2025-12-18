// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/os_windows.hpp>
#include <kernel/util/dist.hpp>

#ifdef __unix__
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#endif

#include <fstream>
#include <vector>
#include <string>
#include <array>

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
  public:
    enum memory_type : std::uint8_t
    {
      peak_physical = 0,
      cur_physical = 1,
      peak_virtual = 2,
      cur_virtual = 3,
      cur_swap = 4,
      size = 5
    };

    static const char* format(memory_type type)
    {
      switch (type) {
        case peak_physical:
          return "Peak Physical";
        case cur_physical:
          return "Current Physical";
        case peak_virtual:
          return "Peak Virtual";
        case cur_virtual:
          return "Current Virtual";
        case cur_swap:
          return "Current Swap";
        default:
          return "Unkown";
      }
    }

  private:
    std::array<std::size_t, memory_type::size> _memory;

  public:
    MemoryUsage() :
     _memory()
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
          std::size_t& _current_physical = _memory[memory_type::cur_physical];
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_current_physical);
          _current_physical *= 1024;
          continue;
        }

        if (line.starts_with("VmHWM"))
        {
          std::size_t& _peak_physical = _memory[memory_type::peak_physical];
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_peak_physical);
          _peak_physical *= 1024;
          continue;
        }

        if (line.starts_with("VmSize"))
        {
          std::size_t& _current_virtual = _memory[memory_type::cur_virtual];
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_current_virtual);
          _current_virtual *= 1024;
          continue;
        }

        if (line.starts_with("VmPeak"))
        {
          std::size_t& _peak_virtual = _memory[memory_type::peak_virtual];
          std::deque<String> v = line.split_by_whitespaces();
          XASSERTM(v.back() == "kB", "get_memory_usage: unit mismatch!");
          v.at(v.size() - 2).parse(_peak_virtual);
          _peak_virtual *= 1024;
          continue;
        }

        if (line.starts_with("VmSwap"))
        {
          std::size_t& _current_swap = _memory[memory_type::cur_swap];
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

      _memory[memory_type::cur_physical]  = std::size_t(work_set_size);
      _memory[memory_type::cur_virtual]   = std::size_t(page_file_usage);
      _memory[memory_type::peak_physical] = std::size_t(work_set_size_peak);
      _memory[memory_type::peak_virtual]  = std::size_t(page_file_usage_peak);
      _memory[memory_type::cur_swap]      = std::size_t(0);

#elif defined(__unix__)
      // see:
      // https://linux.die.net/man/2/getrusage
      // https://www.freebsd.org/cgi/man.cgi?query=getrusage
      struct rusage r_usage;
      if (0 != getrusage(RUSAGE_SELF, &r_usage))
        XABORTM("Error in getrusage call!");
      // 'ru_maxrss' is given in kilobytes
      _memory[memory_type::cur_physical]  = std::size_t(r_usage.ru_maxrss) * 1024u;
      _memory[memory_type::cur_virtual]   = std::size_t(r_usage.ru_maxrss) * 1024u;
      _memory[memory_type::peak_physical] = std::size_t(r_usage.ru_maxrss) * 1024u;
      _memory[memory_type::peak_virtual]  = std::size_t(r_usage.ru_maxrss) * 1024u;
      _memory[memory_type::cur_swap]      = std::size_t(0);
#endif
    }

    /// Returns current (last stamp call) physical memory
    std::size_t get_current_physical() const
    {
      return _memory[memory_type::cur_physical];
    }

    /// Returns peak physical memory
    std::size_t get_peak_physical() const
    {
      return _memory[memory_type::peak_physical];
    }

    /// Returns current (last stamp call) virtual memory
    std::size_t get_current_virtual() const
    {
      return _memory[memory_type::cur_virtual];
    }

    /// Returns peak virtual memory
    std::size_t get_peak_virtual() const
    {
      return _memory[memory_type::peak_virtual];
    }

    /// Returns current (last stamp call) swap memory
    std::size_t get_current_swap() const
    {
      return _memory[memory_type::cur_swap];
    }

    /// Retrieve formatted memory usage string
    String get_formatted_memory_usage() const
    {
      String r;
      for(std::size_t k = 0u; k < memory_type::size; ++k)
      {
        memory_type type = memory_type(k);
        r += (String(format(type)) + String(":")).pad_back(20) + stringify(_memory[type] / 1024 / 1024) + " MByte\n";
      }

      return r;
    }

    String get_formatted_memory_usage(const Dist::Comm& comm, memory_type mem_type) const
    {
      std::uint64_t min_p, max_p, sum_p;
      min_p = max_p = sum_p = _memory[mem_type];
      comm.allreduce(&min_p, &min_p, 1u, Dist::op_min);
      comm.allreduce(&max_p, &max_p, 1u, Dist::op_max);
      comm.allreduce(&sum_p, &sum_p, 1u, Dist::op_sum);
      return String(format(mem_type)).pad_back(20u, '.') + ":   " + stringify_bytes(sum_p, 3, 0) + " [ Max: " + stringify_bytes(max_p, 3, 0) + " / Min: " + stringify_bytes(min_p, 3, 0) + " ]";
    }

    /**
     * \brief Returns the formatted peak physical memory usage over an entire communicator
     *
     * This function formats the sum of the peak memory usage over all ranks in the communicator
     * as well as the maximum and the minimum peak memory usage.
     */
    static String format_peak_physical_usage(const Dist::Comm& comm)
    {
      MemoryUsage mu;
      return mu.get_formatted_memory_usage(comm, memory_type::peak_physical);
    }
  }; //class MemoryUsage
} // namespace FEAT
