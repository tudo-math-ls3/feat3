#include <iostream>
#include <kernel/scarc/solver_pattern.hpp>

#ifndef SERIAL
#include <mpi.h>
#endif
using namespace FEAST;
using namespace Foundation;

int main(int argc, char* argv[])
{

#ifndef SERIAL
  MPI_Init(&argc, &argv);
#endif


#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
