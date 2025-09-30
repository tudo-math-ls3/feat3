// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/os_windows.hpp>
#include <kernel/util/likwid_marker.hpp>
#ifdef FEAT_HAVE_DEATH_HANDLER
#include <death_handler.h>
#endif

#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cstdlib>

#if defined(__linux) || defined(__unix__)
#include <sys/types.h>
#include <unistd.h>
#include <execinfo.h>
#endif

#ifdef FEAT_HAVE_MPI
#include <mpi.h>
#endif

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

using namespace FEAT;

// static member initialization
bool Runtime::_initialized = false;
bool Runtime::_finalized = false;
bool Runtime::SyncGuard::sync_on = true;

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
    std::cerr << "ERROR: Runtime::initialize called twice!\n";
    std::cerr.flush();
    Runtime::abort();
  }
  if (_finalized)
  {
    std::cerr << "ERROR: Runtime::initialize called after Runtime::finalize!\n";
    std::cerr.flush();
    Runtime::abort();
  }

  // initialize Dist operations
  if(!Dist::initialize(argc, argv))
  {
    std::cerr << "ERROR: Failed to initialize Dist operations!\n";
    std::cerr.flush();
    Runtime::abort();
  }

  // get the MPI world comm rank of this process
  const int my_rank = Dist::Comm::world().rank();

  // initialize likwid marker api
  FEAT_MARKER_INIT;
  // initialize memory pool for main memory
  MemoryPool::initialize();

#ifdef FEAT_HAVE_CUDA
  // initialize memory pool for CUDA memory
  Util::cuda_initialize(my_rank, 1, 1, Util::cuda_get_device_count());
  Util::cuda_set_blocksize(256, 256, 256, 256, 256, 128);
#endif

  // check whether '---debug [<ranks...>]' option is given
  // if so, then trigger a breakpoint for the specified ranks
  for(int iarg = 1; iarg < argc; ++iarg)
  {
    if(strcmp(argv[iarg], "---build-info") == 0)
    {
      // dump some information about our build to the console
      if(my_rank == 0)
      {
        std::cout << "--- FEAT BUILD INFORMATION ---\n";
        //            123456789-123456789-123456789-
        std::cout << "__cplusplus...................: " << __cplusplus << "\n";
#ifdef __STDC_VERSION__
        std::cout << "__STDC_VERSION__..............: " << __STDC_VERSION__ << "\n";
#else
        std::cout << "__STDC_VERSION__..............: -N/A-\n";
#endif
#ifdef __STDC_HOSTED__
        std::cout << "__STDC_HOSTED__...............: " << __STDC_HOSTED__ << "\n";
#else
        std::cout << "__STDC_HOSTED__...............: -N/A-\n";
#endif
#ifdef __STDCPP_THREADS__
        std::cout << "__STDCPP_THREADS__............: " << __STDCPP_THREADS__ << "\n";
#else
        std::cout << "__STDCPP_THREADS__............: -N/A-\n";
#endif
#ifdef _OPENMP
        std::cout << "_OPENMP.......................: yes\n";
#else
        std::cout << "_OPENMP.......................: no\n";
#endif
        std::cout << "FEAT Version..................: " << FEAT::version_major << "." << FEAT::version_minor << "." << FEAT::version_patch << "\n";
#ifdef FEAT_GIT_SHA1
        std::cout << "FEAT_GIT_SHA1.................: " << FEAT_GIT_SHA1 << "\n";
#else
        std::cout << "FEAT_GIT_SHA1.................: -N/A-\n";
#endif
#ifdef FEAT_BUILD_ID
        std::cout << "FEAT_BUILD_ID.................: " << FEAT_BUILD_ID << "\n";
#else
        std::cout << "FEAT_BUILD_ID.................: -N/A-\n";
#endif
#ifdef FEAT_SOURCE_DIR
        std::cout << "FEAT_SOURCE_DIR...............: " << FEAT_SOURCE_DIR << "\n";
#else
        std::cout << "FEAT_SOURCE_DIR...............: -N/A-\n";
#endif
#ifdef FEAT_BUILD_DIR
        std::cout << "FEAT_BUILD_DIR................: " << FEAT_BUILD_DIR << "\n";
#else
        std::cout << "FEAT_BUILD_DIR................: -N/A-\n";
#endif
#ifdef FEAT_COMPILER
        std::cout << "FEAT_COMPILER.................: " << FEAT_COMPILER << "\n";
#else
        std::cout << "FEAT_COMPILER.................: -N/A-\n";
#endif
#ifdef FEAT_COMPILER_GNU
        std::cout << "FEAT_COMPILER_GNU.............: " << FEAT_COMPILER_GNU << "\n";
#else
        std::cout << "FEAT_COMPILER_GNU.............: -N/A-\n";
#endif
#ifdef FEAT_COMPILER_CLANG
        std::cout << "FEAT_COMPILER_CLANG...........: " << FEAT_COMPILER_CLANG << "\n";
#else
        std::cout << "FEAT_COMPILER_CLANG...........: -N/A-\n";
#endif
#ifdef FEAT_COMPILER_CRAY
        std::cout << "FEAT_COMPILER_CRAY............: " << FEAT_COMPILER_CRAY << "\n";
#else
        std::cout << "FEAT_COMPILER_CRAY............: -N/A-\n";
#endif
#ifdef FEAT_COMPILER_INTEL
        std::cout << "FEAT_COMPILER_INTEL...........: " << FEAT_COMPILER_INTEL << "\n";
#else
        std::cout << "FEAT_COMPILER_INTEL...........: -N/A-\n";
#endif
#ifdef FEAT_COMPILER_MICROSOFT
        std::cout << "FEAT_COMPILER_MICROSOFT.......: " << FEAT_COMPILER_MICROSOFT << "\n";
#else
        std::cout << "FEAT_COMPILER_MICROSOFT.......: -N/A-\n";
#endif
#ifdef FEAT_DEBUG_MODE
        std::cout << "FEAT_DEBUG_MODE...............: yes\n";
#else
        std::cout << "FEAT_DEBUG_MODE...............: no\n";
#endif
#ifdef FEAT_EICKT
        std::cout << "FEAT_EICKT....................: yes\n";
#else
        std::cout << "FEAT_EICKT....................: no\n";
#endif
#ifdef FEAT_INDEX_U32
        std::cout << "FEAT_INDEX_U32................: yes\n";
#else
        std::cout << "FEAT_INDEX_U32................: no\n";
#endif
#ifdef FEAT_MPI_THREAD_MULTIPLE
        std::cout << "FEAT_MPI_THREAD_MULTIPLE......: yes\n";
#else
        std::cout << "FEAT_MPI_THREAD_MULTIPLE......: no\n";
#endif
#ifdef FEAT_NO_CONFIG
        std::cout << "FEAT_NO_CONFIG................: yes\n";
#else
        std::cout << "FEAT_NO_CONFIG................: no\n";
#endif
#ifdef FEAT_OVERRIDE_MPI_OPS
        std::cout << "FEAT_OVERRIDE_MPI_OPS.........: yes\n";
#else
        std::cout << "FEAT_OVERRIDE_MPI_OPS.........: no\n";
#endif
#ifdef FEAT_USE_MKL_SPARSE_EXECUTOR
        std::cout << "FEAT_USE_MKL_SPARSE_EXECUTOR..: yes\n";
#else
        std::cout << "FEAT_USE_MKL_SPARSE_EXECUTOR..: no\n";
#endif
#ifdef FEAT_UNROLL_BANDED
        std::cout << "FEAT_UNROLL_BANDED............: yes\n";
#else
        std::cout << "FEAT_UNROLL_BANDED............: no\n";
#endif
#ifdef FEAT_HAVE_ALGLIB
        std::cout << "FEAT_HAVE_ALGLIB..............: yes\n";
#else
        std::cout << "FEAT_HAVE_ALGLIB..............: no\n";
#endif
#ifdef FEAT_HAVE_BOOST
        std::cout << "FEAT_HAVE_BOOST...............: yes\n";
#else
        std::cout << "FEAT_HAVE_BOOST...............: no\n";
#endif
#ifdef FEAT_HAVE_CGAL
        std::cout << "FEAT_HAVE_CGAL................: yes\n";
#else
        std::cout << "FEAT_HAVE_CGAL................: no\n";
#endif
#ifdef FEAT_HAVE_CUDA
        std::cout << "FEAT_HAVE_CUDA................: yes\n";
#else
        std::cout << "FEAT_HAVE_CUDA................: no\n";
#endif
#ifdef FEAT_HAVE_CUDSS
        std::cout << "FEAT_HAVE_CUDSS...............: yes\n";
#else
        std::cout << "FEAT_HAVE_CUDSS...............: no\n";
#endif
#ifdef FEAT_HAVE_DEATH_HANDLER
        std::cout << "FEAT_HAVE_DEATH_HANDLER.......: yes\n";
#else
        std::cout << "FEAT_HAVE_DEATH_HANDLER.......: no\n";
#endif
#ifdef FEAT_HAVE_FPARSER
        std::cout << "FEAT_HAVE_FPARSER.............: yes\n";
#else
        std::cout << "FEAT_HAVE_FPARSER.............: no\n";
#endif
#ifdef FEAT_HAVE_FLOATX
        std::cout << "FEAT_HAVE_FLOATX..............: yes\n";
#else
        std::cout << "FEAT_HAVE_FLOATX..............: no\n";
#endif
#ifdef FEAT_HAVE_HALFMATH
        std::cout << "FEAT_HAVE_HALFMATH............: yes\n";
#else
        std::cout << "FEAT_HAVE_HALFMATH............: no\n";
#endif
#ifdef FEAT_HAVE_HYPRE
        std::cout << "FEAT_HAVE_HYPRE...............: yes\n";
#else
        std::cout << "FEAT_HAVE_HYPRE...............: no\n";
#endif
#ifdef FEAT_HAVE_METIS
        std::cout << "FEAT_HAVE_METIS...............: yes\n";
#else
        std::cout << "FEAT_HAVE_METIS...............: no\n";
#endif
#ifdef FEAT_HAVE_MKL
        std::cout << "FEAT_HAVE_MKL.................: yes\n";
#else
        std::cout << "FEAT_HAVE_MKL.................: no\n";
#endif
#ifdef FEAT_HAVE_MPI
        std::cout << "FEAT_HAVE_MPI.................: yes\n";
#else
        std::cout << "FEAT_HAVE_MPI.................: no\n";
#endif
#ifdef FEAT_HAVE_OMP
        std::cout << "FEAT_HAVE_OMP.................: yes\n";
#else
        std::cout << "FEAT_HAVE_OMP.................: no\n";
#endif
#ifdef FEAT_HAVE_PARMETIS
        std::cout << "FEAT_HAVE_PARMETIS............: yes\n";
#else
        std::cout << "FEAT_HAVE_PARMETIS............: no\n";
#endif
#ifdef FEAT_HAVE_PMP
        std::cout << "FEAT_HAVE_PMP.................: yes\n";
#else
        std::cout << "FEAT_HAVE_PMP.................: no\n";
#endif
#ifdef FEAT_HAVE_QUADMATH
        std::cout << "FEAT_HAVE_QUADMATH............: yes\n";
#else
        std::cout << "FEAT_HAVE_QUADMATH............: no\n";
#endif
#ifdef FEAT_HAVE_SUITESPARSE
        std::cout << "FEAT_HAVE_SUITESPARSE.........: yes\n";
#else
        std::cout << "FEAT_HAVE_SUITESPARSE.........: no\n";
#endif
#ifdef FEAT_HAVE_SUPERLU_DIST
        std::cout << "FEAT_HAVE_SUPERLU_DIST........: yes\n";
#else
        std::cout << "FEAT_HAVE_SUPERLU_DIST........: no\n";
#endif
#ifdef FEAT_HAVE_TRIANGLE
        std::cout << "FEAT_HAVE_TRIANGLE............: yes\n";
#else
        std::cout << "FEAT_HAVE_TRIANGLE............: no\n";
#endif
#ifdef FEAT_HAVE_TRILINOS
        std::cout << "FEAT_HAVE_TRILINOS............: yes\n";
#else
        std::cout << "FEAT_HAVE_TRILINOS............: no\n";
#endif
#ifdef FEAT_HAVE_UMFPACK
        std::cout << "FEAT_HAVE_UMFPACK.............: yes\n";
#else
        std::cout << "FEAT_HAVE_UMFPACK.............: no\n";
#endif
#ifdef FEAT_HAVE_ZFP
        std::cout << "FEAT_HAVE_ZFP.................: yes\n";
#else
        std::cout << "FEAT_HAVE_ZFP.................: no\n";
#endif
#ifdef FEAT_HAVE_ZLIB
        std::cout << "FEAT_HAVE_ZLIB................: yes\n";
#else
        std::cout << "FEAT_HAVE_ZLIB................: no\n";
#endif
#ifdef FEAT_HAVE_ZOLTAN
        std::cout << "FEAT_HAVE_ZOLTAN..............: yes\n";
#else
        std::cout << "FEAT_HAVE_ZOLTAN..............: no\n";
#endif
        std::cout << "--- END OF FEAT BUILD INFORMATION ---\n";
        std::cout.flush();
      }
    }

    if(strcmp(argv[iarg], "---print-pid") == 0)
    {
#ifdef FEAT_HAVE_MPI
#if defined(_WIN32)
      std::cout << "Process ID " << std::setw(6) << Windows::get_current_process_id() << " runs rank " << my_rank << "\n";
#elif defined(__linux) || defined(__unix__)
      std::cout << "Process ID " << std::setw(6) << getpid() << " runs rank " << my_rank << "\n";
#endif // defined(_WIN32)
#else // no FEAT_HAVE_MPI
#if defined(_WIN32)
      std::cout << "Process ID " << std::setw(6) << Windows::get_current_process_id() << "\n";
#elif defined(__linux) || defined(__unix__)
      std::cout << "Process ID " << std::setw(6) << getpid() << "\n";
#endif // defined(_WIN32)
#endif // FEAT_HAVE_MPI
      std::cout.flush();
    }

#ifdef FEAT_COMPILER_MICROSOFT
    if(strcmp(argv[iarg], "---debug-break") == 0)
    {
#ifdef FEAT_HAVE_MPI
      // in an MPI case, the ranks to debug have to be specified
      for(++iarg; iarg < argc; ++iarg)
      {
        char* endptr = nullptr;
        int irank = int(std::strtol(argv[iarg], &endptr, 0));
        if((endptr != argv[iarg]) && (*endptr == '\0'))
        {
          if(irank == my_rank)
          {
            __debugbreak();
            break;
          }
        }
        else // argv[iarg] could not be parsed as a number
        {
          --iarg;
          break;
        }
      }
#else // no  FEAT_HAVE_MPI
      __debugbreak();
#endif // FEAT_HAVE_MPI
      continue;
    }
#endif // FEAT_COMPILER_MICROSOFT
  }

  _initialized = true;
  // finally initialize Syncguard
  SyncGuard::sync_on = true;
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
  }

#ifdef FEAT_HAVE_MPI
  ::MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::abort();
}

int Runtime::finalize()
{
  if (!_initialized)
  {
    std::cerr << "ERROR: Runtime::finalize called before Runtime::initialize!\n";
    std::cerr.flush();
    Runtime::abort();
  }
  if (_finalized)
  {
    std::cerr << "ERROR: Runtime::finalize called twice!\n";
    std::cerr.flush();
    Runtime::abort();
  }

  MemoryPool::finalize();
#ifdef FEAT_HAVE_CUDA
  Util::cuda_finalize();
#endif

  // finalize Dist operations
  Dist::finalize();
  // finalize Likwid markerAPI
  FEAT_MARKER_CLOSE;

  // reset device, which is/can be nesessary to have sanitizer and profiler work properly
  // this should not be called, if there are other cuda contexts live while FEAT is finalized,
  // for example if MPI is initilized by another library which shares the same process, or a feat app
  // is called by another app without reinitializing the cuda device
#if defined(FEAT_HAVE_CUDA) && defined(FEAT_FINALIZE_RESETS_DEVICE)
  Util::cuda_reset_device();
#endif
  _finalized = true;

  // return successful exit code
  return EXIT_SUCCESS;
}
