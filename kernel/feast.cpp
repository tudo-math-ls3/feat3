#include <kernel/feast.hpp>

#include <cstdlib>

#ifdef FEAST_MPI
#include <mpi.h>
#endif

namespace FEAST
{
  void initialise(int& argc, char**& argv)
  {
    int rank(0), nprocs(0);
    initialise(argc, argv, rank, nprocs);
  }

  void initialise(int& argc, char**& argv, int& rank, int& nprocs)
  {
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
#endif
  }

  void abort()
  {
#ifdef FEAST_MPI
    // abort MPI
    ::MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    std::abort();
  }

  int finalise()
  {
#ifdef FEAST_MPI
    // finalise MPI
    ::MPI_Finalize();
#endif
    // return successful exit code
    return EXIT_SUCCESS;
  }
} // namespace FEAST
