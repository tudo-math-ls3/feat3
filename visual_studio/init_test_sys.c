// *****************************************************************************
// Error Handler Installer for Visual Studio test system
// -----------------------------------------------------
// This file implements a function named 'InitTestSystem', which configures
// the Windows structured exception handlers (SEH) for the FEAT test system.
//
// This function will has 3 effects:
// 1. It disables the popup dialog which appears when an application crashes.
// 2. It disables the default handling of the 'abort' function, which also
//    causes a popup dialog to appear otherwise.
// 3. It installs an 'unhandled exception filter' for SEH which dumps the
//    error reason as well as a call-stack back-trace to stderr.
//
// X86 compile and link:
// cl.exe /TC /Zl /c init_test_sys.c /Fo"..\obj\init_test_sys-x86.obj"
// lib.exe /MACHINE:X86 /NODEFAULTLIB "..\obj\init_test_sys-x86.obj" /OUT:"..\lib\init_test_sys-x86.lib"
//
// X64 compile and link:
// cl.exe /TC /Zl /c init_test_sys.c /Fo"..\obj\init_test_sys-x64.obj"
// lib.exe /MACHINE:X64 /NODEFAULTLIB "..\obj\init_test_sys-x64.obj" /OUT:"..\lib\init_test_sys-x64.lib"
//
// \author Peter Zajac
// *****************************************************************************
#include <Windows.h>
#include <DbgHelp.h>
#include <stdio.h>
#include <memory.h>

LONG WINAPI WinExceptionFilter(LPEXCEPTION_POINTERS p)
{
  // local variables
  CONTEXT ctx;
  HANDLE prc;
  HANDLE thr;
  STACKFRAME64 stack_frame;
  IMAGEHLP_LINE64 img_line;
  DWORD line_dsp;
  DWORD64 addr;
  struct
  {
    SYMBOL_INFO si;
    char name[512];
  } sim;

  // first of all, flush the output buffer
  fflush(stdout);

  // make a copy of the context record
  ctx = *p->ContextRecord;

  // retrieve our process and thread id
  prc = GetCurrentProcess();
  thr = GetCurrentThread();

  // initialise line info
  img_line.SizeOfStruct = sizeof(IMAGEHLP_LINE64);

  // create a symbol information structure
  memset(&sim, 0, sizeof(SYMBOL_INFO) + 512);
  sim.si.MaxNameLen = 512;
  sim.si.SizeOfStruct = sizeof(SYMBOL_INFO);

  // initialise symbol table
  SymInitialize(prc, NULL, TRUE);

  // initialise the stack frame stuff
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
  fprintf(stderr, "\nCall-Stack Back-Trace:\n");
  fprintf(stderr,   "----------------------\n");

  // walk the stack
  while (StackWalk64(machine_type, prc, thr, &stack_frame, &ctx, NULL,
                     &SymFunctionTableAccess64, &SymGetModuleBase64, NULL) != FALSE)
  {
    // get stack address
    addr = stack_frame.AddrPC.Offset;

    // print stack address
    fprintf(stderr, "0x%0*I64X", 2*sizeof(void*), addr);

    // get symbol info to obtain function name
    if(SymFromAddr(prc, addr, 0, &sim.si) != FALSE)
    {
      // print function name
      fprintf(stderr, ": '%s'", sim.si.Name);
    }
    else
    {
      // cannot query function name
      fprintf(stderr, ": ???");
    }

    // get source file line info
    if(SymGetLineFromAddr64(prc, addr, &line_dsp, &img_line) != FALSE)
    {
      // print source file and line
      fprintf(stderr, " ['%s' @ %i]", img_line.FileName, img_line.LineNumber);
    }

    // next stack entry
    fprintf(stderr, "\n");
  }

  // flush
  fprintf(stderr, "\n>>>>> TERMINATING <<<<<\n");
  fflush(stderr);

  // return
  return EXCEPTION_EXECUTE_HANDLER;
}

void __stdcall InitTestSystem()
{
  // disable any error prompts for testing
  SetErrorMode(GetErrorMode() | 0x8003);
  // disable handling of abort function
  _set_abort_behavior(0, 0x1);
  // set our exception filter
  SetUnhandledExceptionFilter(WinExceptionFilter);
}
