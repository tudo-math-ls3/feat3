#include <kernel/util/runtime.hpp>
#include <kernel/util/assertion.hpp>
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
PropertyMap Runtime::_global_property_map;
bool Runtime::_initialised = false;
bool Runtime::_finished = false;

PropertyMap * Runtime::global_property()
{
  if (!_initialised)
  {
    std::cerr << "ERROR: global_property_map not _initialised!" << std::endl;
    std::cerr << "       Call Runtime::initialise first!" << std::endl;
    std::cerr.flush();
    Runtime::abort();
  }
  return &_global_property_map;
}

void Runtime::initialise(int& argc, char**& argv)
{
  // Note:
  // On Windows, these two function calls MUST come before anything else,
  // otherwise one the following calls may cause the automated regression
  // test system to halt with an error prompt awaiting user interaction.
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

  int rank = 0;
#ifdef FEAT_HAVE_MPI
  // initialise MPI
  if(::MPI_Init(&argc, &argv) != MPI_SUCCESS)
    abort();
  // get rank for MemPool initialisation
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  (void)argc;
  (void)argv;
#endif

  // read in initial settings from provided ini file via system environment variable FEAT_INI_FILE
  if (const char* c_file_path = std::getenv("FEAT_INI_FILE"))
  {
    String property_file(c_file_path);
    if (std::ifstream(property_file).good())
    {
      _global_property_map.parse(property_file, true);
    }
    else
    {
      std::cerr << "WARNING: FEAT ini file " << property_file << " not found!" << std::endl;
    }
  }

  MemoryPool<Mem::Main>::initialise();

#ifdef FEAT_HAVE_CUDA
  MemoryPool<Mem::CUDA>::initialise(rank,
    atoi(_global_property_map.query("MPI.ranks_per_node", "1").c_str()),
    atoi(_global_property_map.query("MPI.ranks_per_uma", "1").c_str()),
    atoi(_global_property_map.query("MPI.gpus_per_node", "1").c_str())
    );
  // read in initial settings from Runtime and store them in MemoryPool<CUDA>
  Index misc = (Index)atoi(_global_property_map.query("CUDA.blocksize_misc", "256").c_str());
  Index reduction = (Index)atoi(_global_property_map.query("CUDA.blocksize_reduction", "256").c_str());
  Index spmv = (Index)atoi(_global_property_map.query("CUDA.blocksize_spmv", "256").c_str());
  Index axpy = (Index)atoi(_global_property_map.query("CUDA.blocksize_axpy", "256").c_str());
  MemoryPool<Mem::CUDA>::set_blocksize(misc, reduction, spmv, axpy);
#else
  (void)rank;
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

#ifdef FEAT_HAVE_MPI
  // finalise MPI
  ::MPI_Finalize();
#endif

  _finished = true;

  // return successful exit code
  return EXIT_SUCCESS;
}
