#ifndef SERIAL
#include <mpi.h>
#endif

#include <kernel/scarc/solver_pattern.hpp>
#include <iostream>

using namespace FEAST;
using namespace Foundation;

int main(int argc, char* argv[])
{

#ifndef SERIAL
  MPI_Init(&argc, &argv);
#endif
  (void)argc;
  (void)argv;
  std::cout<<"CTEST_FULL_OUTPUT"<<std::endl;


#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
