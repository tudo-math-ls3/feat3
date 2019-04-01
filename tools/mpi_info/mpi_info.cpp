// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <iostream>

#ifdef FEAT_HAVE_MPI
#include <mpi.h>

using namespace FEAT;

int main(int argc, char ** argv)
{
  int required = MPI_THREAD_MULTIPLE;
  int provided = MPI_THREAD_SINGLE;
  std::cout<<"Evaluating available thread support..."<<std::endl;
  if (::MPI_Init_thread(&argc, &argv, required, &provided) != MPI_SUCCESS)
  {
    std::cerr<<"MPI_Init_thread failed!"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  switch (provided)
  {
    case MPI_THREAD_MULTIPLE:
      std::cout<<"MPI_THREAD_MULTIPLE available"<<std::endl;
      std::cout<<"MPI_THREAD_SERIALIZED available"<<std::endl;
      std::cout<<"MPI_THREAD_FUNNELED available"<<std::endl;
      std::cout<<"MPI_THREAD_SINGLE available"<<std::endl;
      break;
    case MPI_THREAD_SERIALIZED:
      std::cout<<"MPI_THREAD_SERIALIZED available"<<std::endl;
      std::cout<<"MPI_THREAD_FUNNELED available"<<std::endl;
      std::cout<<"MPI_THREAD_SINGLE available"<<std::endl;
      break;
    case MPI_THREAD_FUNNELED:
      std::cout<<"MPI_THREAD_FUNNELED available"<<std::endl;
      std::cout<<"MPI_THREAD_SINGLE available"<<std::endl;
      break;
    case MPI_THREAD_SINGLE:
      std::cout<<"MPI_THREAD_SINGLE available"<<std::endl;
      break;
    default:
      std::cerr<<"invalid mpi provided value!"<<std::endl;
      break;
  }

  MPI_Finalize();
  return 0;
}
#else
int main(int , char **)
{
  std::cerr<<"no mpi support configured!"<<std::endl;
  return 0;
}
#endif
