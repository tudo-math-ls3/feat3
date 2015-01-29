#include <kernel/base_header.hpp>

#ifdef FEAST_MPI

#include <mpi.h>


int main(int argc, char* argv[])
{
  int rank(0), nprocs(1);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Finalize();
  return 0;
}

#endif // FEAST_MPI
