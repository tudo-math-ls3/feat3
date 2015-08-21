#include <kernel/util/runtime.hpp>
#include <kernel/util/assertion.hpp>

#include <cstdlib>

#ifdef FEAST_MPI
#include <mpi.h>
#endif

using namespace FEAST;

// static member initialisation
PropertyMap Runtime::_global_property_map;
bool Runtime::_initialised = false;
bool Runtime::_finished = false;
Index Runtime::_flops = Index(0);

PropertyMap & Runtime::global_property()
{
  ASSERT(_initialised == true, "global_property_map not _initialised! Call initialise first");
  return _global_property_map;
}

void Runtime::add_flops(Index flops)
{
  _flops += flops;
}

Index Runtime::get_flops()
{
  return _flops;
}

String Runtime::get_formated_flops(double seconds)
{
  double flops((double)_flops);
  flops /= seconds;
  flops /= 1000.; // kilo
  flops /= 1000.; // mega
  flops /= 1000.; // giga
  return stringify(flops) + " GFlop/s";
}

void Runtime::reset_flops()
{
  _flops = Index(0);
}

void Runtime::initialise(int& argc, char**& argv)
{
  int rank(0), nprocs(0);
  initialise(argc, argv, rank, nprocs);
}

void Runtime::initialise(int& argc, char**& argv, int& rank, int& nprocs)
{
  if (_initialised)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "Runtime::initialise called twice!");
  }
  if (_finished)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "Runtime::initialise called after finalise!");
  }

  // reset rank and nprocs
  rank = 0;
  nprocs = 0;

#ifdef FEAST_MPI
  // initialise MPI
  if(::MPI_Init(&argc, &argv) != MPI_SUCCESS)
    abort();
  // fetch rank and world size
  ::MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ::MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#else
  (void)argc;
  (void)argv;
  rank = 0;
  nprocs = 1;
#endif

  // read in initial settings from provided ini file via system environment variable FEAST_INI_FILE
  if (const char* c_file_path = std::getenv("FEAST_INI_FILE"))
  {
    String property_file(c_file_path);
    _global_property_map.parse(property_file, true);
  }

  _initialised = true;
}

void Runtime::abort()
{
#ifdef FEAST_MPI
  // abort MPI
  ::MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::abort();
}

int Runtime::finalise()
{
  if (!_initialised)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "Runtime::finalise called before initialise!");
  }
  if (_finished)
  {
    throw InternalError(__func__, __FILE__, __LINE__, "Runtime::finalise called twice!");
  }

#ifdef FEAST_MPI
  // finalise MPI
  ::MPI_Finalize();
#endif

  _finished = true;

  // return successful exit code
  return EXIT_SUCCESS;
}
