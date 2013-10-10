#include <kernel/base_header.hpp>
#include <kernel/foundation/control.hpp>
#include <kernel/foundation/data.hpp>
#include <kernel/foundation/topology.hpp>

#include <iostream>
#ifndef SERIAL
#  include <mpi.h>
#endif

using namespace FEAST;
using namespace Foundation;

int main(int argc, char* argv[])
{

#ifndef SERIAL
  MPI_Init(&argc, &argv);
#endif
  std::cout<<"CTEST_FULL_OUTPUT"<<std::endl;

  int me(0);
  //build communication structures
#ifndef SERIAL
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif


#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
