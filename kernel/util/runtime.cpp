// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/os_windows.hpp>
#ifdef FEAT_HAVE_DEATH_HANDLER
#include <death_handler.h>
#endif

#include <cstdlib>
#include <fstream>

#if defined(__linux) || defined(__unix__)
#include <execinfo.h>
#endif

#ifdef FEAT_HAVE_MPI
#include <mpi.h>
#endif

using namespace FEAT;

// static member initialization
bool Runtime::_initialized = false;
bool Runtime::_finished = false;

void Runtime::initialize(int& argc, char**& argv)
{
  /// \platformswitch
  /// On Windows, these two function calls MUST come before anything else,
  /// otherwise one the following calls may cause the automated regression
  /// test system to halt with an error prompt awaiting user interaction.
#if defined(_WIN32) && defined(FEAT_TESTING_VC)
  Windows::disable_error_prompts();
  Windows::install_seh_filter();
#endif

  if (_initialized)
  {
    std::cerr << "ERROR: Runtime::initialize called twice!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }
  if (_finished)
  {
    std::cerr << "ERROR: Runtime::initialize called after Runtime::finalize!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }

  // initialize Dist operations
  if(!Dist::initialize(argc, argv))
  {
    std::cerr << "ERROR: Failed to initialize Dist operations!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }

  // initialize memory pool for main memory
  MemoryPool<Mem::Main>::initialize();

#ifdef FEAT_HAVE_CUDA
  // initialize memory pool for CUDA memory
  int rank = Dist::Comm::world().rank();
  MemoryPool<Mem::CUDA>::initialize(rank, 1, 1, 1);
  MemoryPool<Mem::CUDA>::set_blocksize(256, 256, 256, 256);
#endif

  _initialized = true;
}

void Runtime::abort(bool dump_call_stack)
{
  if(dump_call_stack)
  {
#if defined(__linux) || defined(__unix__)
#if defined(FEAT_HAVE_DEATH_HANDLER) and not defined(FEAT_HAVE_MPI)
    Debug::DeathHandler death_handler;
#else
    // https://www.gnu.org/software/libc/manual/html_node/Backtraces.html
    void* buffer[1024];
    auto bt_size = backtrace(buffer, 1024);
    char** bt_symb = backtrace_symbols(buffer, bt_size);
    if((bt_size > 0) && (bt_symb != nullptr))
    {
      fprintf(stderr, "\nCall-Stack Back-Trace:\n");
      fprintf(stderr,   "----------------------\n");
      for(decltype(bt_size) i(0); i < bt_size; ++i)
        fprintf(stderr, "%s\n", bt_symb[i]);
      fflush(stderr);
    }
#endif
#elif defined(_WIN32)
    Windows::dump_call_stack_to_stderr();
#endif
#ifdef FEAT_HAVE_MPI
  ::MPI_Abort(MPI_COMM_WORLD, 1);
#endif

  std::abort();
  }
  else
  {
#ifdef FEAT_HAVE_MPI
    ::MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    std::abort();
  }
}

int Runtime::finalize()
{
  if (!_initialized)
  {
    std::cerr << "ERROR: Runtime::finalize called before Runtime::initialize!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }
  if (_finished)
  {
    std::cerr << "ERROR: Runtime::finalize called twice!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }

  MemoryPool<Mem::Main>::finalize();
#ifdef FEAT_HAVE_CUDA
  MemoryPool<Mem::CUDA>::finalize();
#endif

  // finalize Dist operations
  Dist::finalize();

  _finished = true;

  // return successful exit code
  return EXIT_SUCCESS;
}
