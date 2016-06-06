
#ifndef _WIN32
#error This file can only be compiler under Windows.
#endif

// Warning:
// This source file MUST NOT include any other FEAT header file!
#include <Windows.h>
#include <DbgHelp.h>
#include <Psapi.h>
#include <cstdio>
#include <cstdlib>

namespace FEAT
{
  namespace Windows
  {
    long long query_performance_counter()
    {
      LARGE_INTEGER li;
      if(QueryPerformanceCounter(&li) == FALSE)
        return 0ll;
      return li.QuadPart;
    }

    long long query_performance_frequency()
    {
      LARGE_INTEGER li;
      if(QueryPerformanceFrequency(&li) == FALSE)
        return 0ll;
      return li.QuadPart;
    }

    void query_memory_usage(
      unsigned long long& work_set_size,
      unsigned long long& work_set_size_peak,
      unsigned long long& page_file_usage,
      unsigned long long& page_file_usage_peak)
    {
      work_set_size = work_set_size_peak = page_file_usage = page_file_usage_peak = 0ull;

      // initialise memory counters structure
      PROCESS_MEMORY_COUNTERS pmc;
      memset(&pmc, 0, sizeof(PROCESS_MEMORY_COUNTERS));

      // get memory info
      if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(PROCESS_MEMORY_COUNTERS)) != FALSE)
      {
        work_set_size        = pmc.WorkingSetSize;
        work_set_size_peak   = pmc.PeakWorkingSetSize;
        page_file_usage      = pmc.PagefileUsage;
        page_file_usage_peak = pmc.PeakPagefileUsage;
      }
    }

    void dump_call_stack(CONTEXT& ctx, FILE* stream)
    {
      DWORD line_dsp;
      DWORD64 addr;
      struct
      {
        SYMBOL_INFO si;
        char name[512];
      } sim;

      // retrieve our process and thread id
      HANDLE prc = GetCurrentProcess();
      HANDLE thr = GetCurrentThread();

      // initialise line info
      IMAGEHLP_LINE64 img_line;
      img_line.SizeOfStruct = sizeof(IMAGEHLP_LINE64);

      // create a symbol information structure
      memset(&sim, 0, sizeof(SYMBOL_INFO) + 512);
      sim.si.MaxNameLen = 512;
      sim.si.SizeOfStruct = sizeof(SYMBOL_INFO);

      // initialise symbol table
      SymInitialize(prc, NULL, TRUE);

      // initialise the stack frame stuff
      STACKFRAME64 stack_frame;
      memset(&stack_frame, 0, sizeof(stack_frame));

#if defined(_WIN64)
      int machine_type = IMAGE_FILE_MACHINE_AMD64;
      stack_frame.AddrPC.Offset = ctx.Rip;
      stack_frame.AddrFrame.Offset = ctx.Rbp;
      stack_frame.AddrStack.Offset = ctx.Rsp;
#else
      int machine_type = IMAGE_FILE_MACHINE_I386;
      stack_frame.AddrPC.Offset = ctx.Eip;
      stack_frame.AddrFrame.Offset = ctx.Ebp;
      stack_frame.AddrStack.Offset = ctx.Esp;
#endif
      stack_frame.AddrPC.Mode = AddrModeFlat;
      stack_frame.AddrFrame.Mode = AddrModeFlat;
      stack_frame.AddrStack.Mode = AddrModeFlat;

      // flush the output stream
      fflush(stream);

      // dump call stack
      fprintf(stream, "\nCall-Stack Back-Trace:\n");
      fprintf(stream,   "----------------------\n");

      // walk the stack
      while (StackWalk64(machine_type, prc, thr, &stack_frame, &ctx, NULL,
                         &SymFunctionTableAccess64, &SymGetModuleBase64, NULL) != FALSE)
      {
        // get stack address
        addr = stack_frame.AddrPC.Offset;

        // print stack address
        fprintf(stream, "0x%0*I64X", 2*sizeof(void*), addr);

        // get symbol info to obtain function name
        if(SymFromAddr(prc, addr, 0, &sim.si) != FALSE)
        {
          // print function name
          fprintf(stream, ": '%s'", sim.si.Name);
        }
        else
        {
          // cannot query function name
          fprintf(stream, ": ???");
        }

        // get source file line info
        if(SymGetLineFromAddr64(prc, addr, &line_dsp, &img_line) != FALSE)
        {
          // print source file and line
          fprintf(stream, " ['%s' @ %i]", img_line.FileName, img_line.LineNumber);
        }

        // next stack entry
        fprintf(stream, "\n");
      }

      // flush stream
      fflush(stream);
    }

    void dump_call_stack_to_stderr()
    {
      // query context
      CONTEXT ctx;
      RtlCaptureContext(&ctx);
      // make a copy
      CONTEXT ctx2(ctx);
      dump_call_stack(ctx2, stderr);
    }

    LONG WINAPI FeatWinExceptionFilter(LPEXCEPTION_POINTERS p)
    {
      // print an error header
      fprintf(stderr, "\n>>>>> FATAL ERROR <<<<<\n\n");
      fprintf(stderr, "Exception at 0x%0*I64X: ", 2*sizeof(void*), (DWORD64)p->ExceptionRecord->ExceptionAddress);

      // print exception info
      switch(p->ExceptionRecord->ExceptionCode)
      {
      // note: There are more than 20 exceptions, but we do not handle all of them here explicitly.
      case EXCEPTION_ACCESS_VIOLATION:
        // This one is also known as 'segmentation fault'...
        fprintf(stderr, "Access Violation when ");
        if(p->ExceptionRecord->ExceptionInformation[0] == 0)
          fprintf(stderr, "reading from");
        else
          fprintf(stderr, "writing to");
        fprintf(stderr, " 0x%0*I64X\n", 2*sizeof(void*), (DWORD64)p->ExceptionRecord->ExceptionInformation[1]);
        break;
      case EXCEPTION_ILLEGAL_INSTRUCTION:
        fprintf(stderr, "Illegal Instruction\n");
        break;
      case EXCEPTION_IN_PAGE_ERROR:
        fprintf(stderr, "Page Error\n");
        break;
      case EXCEPTION_STACK_OVERFLOW:
        fprintf(stderr, "Stack Overflow\n");
        break;
      default:
        // We do not handle all error types, so just dump the error code here...
        fprintf(stderr, "Error 0x%08X\n", p->ExceptionRecord->ExceptionCode);
        break;
      }

      // dump call stack
      CONTEXT ctx = *p->ContextRecord;
      dump_call_stack(ctx, stderr);

      // flush
      fprintf(stderr, "\n>>>>> TERMINATING <<<<<\n");
      fflush(stderr);

      // return
      return EXCEPTION_EXECUTE_HANDLER;
    }

    void install_seh_filter()
    {
      // set our custom exception filter
      SetUnhandledExceptionFilter(FeatWinExceptionFilter);
    }

    void disable_error_prompts()
    {
      // disable error prompts
      SetErrorMode(GetErrorMode() | SEM_FAILCRITICALERRORS|SEM_NOGPFAULTERRORBOX|SEM_NOOPENFILEERRORBOX);
      // disable handling of abort function
      _set_abort_behavior(0, _WRITE_ABORT_MSG);
    }
  } // namespace Windows
} // namespace FEAT
