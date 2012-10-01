#include <mpi.h>
#include <iostream>
#include <kernel/base_header.hpp>

int main(int argc, char* argv[])
{

#ifndef FEAST_SERIAL_MODE

  MPI_Init(&argc,&argv);

  int rank, numprocs;

  int tag(99), source(0), target(1), count(1);
  int buffer(5678);
  MPI_Status status;
  MPI_Request request(MPI_REQUEST_NULL);

  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(rank == source)
  {
    buffer = 5678;
    MPI_Isend(&buffer, count, MPI_INT, target, tag, MPI_COMM_WORLD, &request);
  }
  if(rank == target)
  {
    MPI_Irecv(&buffer, count, MPI_INT, source, tag, MPI_COMM_WORLD, &request);
  }

  MPI_Wait(&request, &status);

  if(rank == source)
  {
    std::cout << "Processor " << rank << " sent " << buffer << "!" << std::endl;
  }

  if(rank == target)
  {
    std::cout << "Processor " << rank << " got " << buffer << "!" << std::endl;
  }

  MPI_Finalize();

#endif

  return 0;
}
