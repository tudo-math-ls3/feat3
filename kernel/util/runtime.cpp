#include <kernel/util/runtime.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/os_windows.hpp>

#include <cstdlib>
#include <fstream>

#if defined(__linux) || defined(__unix__)
#include <execinfo.h>
#endif

#ifdef FEAT_HAVE_MPI
#include <mpi.h>
#endif

using namespace FEAT;

// static member initialisation
bool Runtime::_initialised = false;
bool Runtime::_finished = false;

void Runtime::initialise(int& argc, char**& argv)
{
  /// \platformswitch
  /// On Windows, these two function calls MUST come before anything else,
  /// otherwise one the following calls may cause the automated regression
  /// test system to halt with an error prompt awaiting user interaction.
#if defined(_WIN32) && defined(FEAT_TESTING_VC)
  Windows::disable_error_prompts();
  Windows::install_seh_filter();
#endif

  if (_initialised)
  {
    std::cerr << "ERROR: Runtime::initialise called twice!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }
  if (_finished)
  {
    std::cerr << "ERROR: Runtime::initialise called after Runtime::finalise!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }

  // initialise Dist operations
  if(!Dist::initialise(argc, argv))
  {
    std::cerr << "ERROR: Failed to initialise Dist operations!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }

  // initialise memory pool for main memory
  MemoryPool<Mem::Main>::initialise();

#ifdef FEAT_HAVE_CUDA
  // initialise memory pool for CUDA memory
  int rank = Dist::Comm::world().rank();
  MemoryPool<Mem::CUDA>::initialise(rank, 1, 1, 1);
  MemoryPool<Mem::CUDA>::set_blocksize(256, 256, 256, 256);
#endif

  _initialised = true;
}

void Runtime::abort(bool dump_call_stack)
{
  if(dump_call_stack)
  {
#if defined(__linux) || defined(__unix__)
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
#elif defined(_WIN32)
    Windows::dump_call_stack_to_stderr();
#endif
  }

#ifdef FEAT_HAVE_MPI
  ::MPI_Abort(MPI_COMM_WORLD, 1);
#endif

  std::abort();
}

int Runtime::finalise()
{
  if (!_initialised)
  {
    std::cerr << "ERROR: Runtime::finalise called before Runtime::initialise!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }
  if (_finished)
  {
    std::cerr << "ERROR: Runtime::finalise called twice!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }

  MemoryPool<Mem::Main>::finalise();
#ifdef FEAT_HAVE_CUDA
  MemoryPool<Mem::CUDA>::finalise();
#endif

  // finalise Dist operations
  Dist::finalise();

  _finished = true;

  // return successful exit code
  return EXIT_SUCCESS;
}
