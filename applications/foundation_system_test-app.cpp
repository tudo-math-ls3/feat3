#include <mpi.h>
#include <iostream>
#include <kernel/base_header.hpp>
#include <kernel/foundation/control.hpp>
#include <kernel/foundation/topology.hpp>

using namespace FEAST;
using namespace Foundation;

/*template<typename DT1_, typename DT2_, typename DT3_>
void test_check_equal_within_eps(DT1_ l, DT2_ r, DT3_ eps)
{
  if(abs(l - r) < eps)
    std::cout << "PASSED" << std::endl;
  else
    std::cout << "FAILED" << std::endl;
}*/

int main(int argc, char* argv[])
{

#ifndef FEAST_SERIAL_MODE

  MPI_Init(&argc,&argv);

  int rank, numprocs;

  ///TODO dedicated processes only
  //tell FEAST, how physical compute nodes and mesh patches are connected
  Topology<> network; //not needed for dummy LB
  Topology<> patches; //only needed for letting dummy LB know, how many patches actually exist
  patches.push_back();
  patches.push_back();
  patches.push_back();
  patches.push_back();
  LBConfig<Topology<> > lbconf(network, patches);

  Control<Parallel, SimpleLoadBalancingPolicy>::init(lbconf);

  //### output only ###
  int me;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  if(me == 0)
  {
    for(Index i(0) ; i < lbconf.patch_process_map.size() ; ++i)
      for(Index j(0) ; j < lbconf.patch_process_map.at(i).size() ; ++j)
        std::cout << "Process " << lbconf.patch_process_map.at(i).at(j) << " distributed to patch " << i << "." << std::endl;
  }
  //### end output only ###

  MPI_Finalize();

#else

  ///serial mode does exactly the same
  //tell FEAST, how physical compute nodes and mesh patches are connected
  Topology<> network; //not needed for dummy LB
  Topology<> patches; //only needed for letting dummy LB know, how many patches actually exist
  patches.push_back();
  patches.push_back();
  patches.push_back();
  patches.push_back();
  LBConfig<Topology<> > lbconf(network, patches);

  Control<Serial, SimpleLoadBalancingPolicy>::init(lbconf);

  //### output only ###
  for(Index i(0) ; i < lbconf.patch_process_map.size() ; ++i)
    for(Index j(0) ; j < lbconf.patch_process_map.at(i).size() ; ++j)
      std::cout << "Process " << lbconf.patch_process_map.at(i).at(j) << " distributed to patch " << i << "." << std::endl;
  //### end output only ###

#endif

  return 0;
}
