#include <kernel/util/memory_usage.hpp>
#include <kernel/util/assertion.hpp>

#include <fstream>
#include <vector>

// Note:
// We cannot include the <Windows.h> and <Psapi.h> headers, as these cannot be
// compiled with disabled compiler extensions. Therefore, we define the required
// function prototypes and structures by hand:
#ifdef _WIN32
extern "C" void* __stdcall GetCurrentProcess(void);
extern "C" int __stdcall K32GetProcessMemoryInfo(void*, void*, unsigned long);
struct PROCESS_MEMORY_COUNTERS
{
  unsigned long cb;
  unsigned long PageFaultCount;
  std::size_t PeakWorkingSetSize;
  std::size_t WorkingSetSize;
  std::size_t QuotaPeakPagedPoolUsage;
  std::size_t QuotaPagedPoolUsage;
  std::size_t QuotaPeakNonPagedPoolUsage;
  std::size_t QuotaNonPagedPoolUsage;
  std::size_t PagefileUsage;
  std::size_t PeakPagefileUsage;
};
#endif // _WIN32

namespace FEAST
{
  namespace Util
  {
    MemUsageInfo get_memory_usage()
    {
      MemUsageInfo info;

#if defined(__linux) || defined(__unix__)
      String line;
      std::ifstream status_file("/proc/self/status");
      if (!status_file.is_open())
        throw InternalError(__func__, __FILE__, __LINE__, "could not open /proc/self/status!");

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
          std::vector<String> v;
          line.split_by_charset(v);
          if (v.at(v.size()-1) != "kB")
            throw InternalError(__func__, __FILE__, __LINE__, "get_memory_usage: unit missmatch!");
          info.current_physical = std::stoul(v.at(v.size()-2));
          info.current_physical *= 1024;
          continue;
        }

        if (line.starts_with("VmHWM"))
        {
          std::vector<String> v;
          line.split_by_charset(v);
          if (v.at(v.size()-1) != "kB")
            throw InternalError(__func__, __FILE__, __LINE__, "get_memory_usage: unit missmatch!");
          info.peak_physical = std::stoul(v.at(v.size()-2));
          info.peak_physical *= 1024;
          continue;
        }

        if (line.starts_with("VmSize"))
        {
          std::vector<String> v;
          line.split_by_charset(v);
          if (v.at(v.size()-1) != "kB")
            throw InternalError(__func__, __FILE__, __LINE__, "get_memory_usage: unit missmatch!");
          info.current_virtual = std::stoul(v.at(v.size()-2));
          info.current_virtual *= 1024;
          continue;
        }

        if (line.starts_with("VmPeak"))
        {
          std::vector<String> v;
          line.split_by_charset(v);
          if (v.at(v.size()-1) != "kB")
            throw InternalError(__func__, __FILE__, __LINE__, "get_memory_usage: unit missmatch!");
          info.peak_virtual = std::stoul(v.at(v.size()-2));
          info.peak_virtual *= 1024;
          continue;
        }

        if (line.starts_with("VmSwap"))
        {
          std::vector<String> v;
          line.split_by_charset(v);
          if (v.at(v.size()-1) != "kB")
            throw InternalError(__func__, __FILE__, __LINE__, "get_memory_usage: unit missmatch!");
          info.current_swap = std::stoul(v.at(v.size()-2));
          info.current_swap *= 1024;
          continue;
        }
      }
      status_file.close();

#elif defined(_WIN32)

    // initialise memory counters structure
    PROCESS_MEMORY_COUNTERS pmc;
    memset(&pmc, 0, sizeof(PROCESS_MEMORY_COUNTERS));

    // get memory info
    if(K32GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(PROCESS_MEMORY_COUNTERS)) != 0)
    {
      info.current_physical = pmc.WorkingSetSize;
      info.current_virtual  = pmc.PagefileUsage;
      info.peak_physical    = pmc.PeakWorkingSetSize;
      info.peak_virtual     = pmc.PeakPagefileUsage;
      info.current_swap     = std::size_t(0);
    }
#endif

      return info;
    }

    String get_formated_memory_usage()
    {
      String r;
      auto m = get_memory_usage();
      r += "Current real: " + stringify(m.current_physical / 1024 / 1024) + " MByte\n";
      r += "Peak real: " + stringify(m.peak_physical / 1024 / 1024) + " MByte\n";
      r += "Current virtual: " + stringify(m.current_virtual / 1024 / 1024) + " MByte\n";
      r += "Peak virtual: " + stringify(m.peak_virtual / 1024 / 1024) + " MByte\n";
      r += "Current swap: " + stringify(m.current_swap / 1024 / 1024) + " MByte\n";

      return r;
    }
  } // namespace Util
} // namespace FEAST
