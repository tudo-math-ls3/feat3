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

  ///TODO dedicated processes only
  //tell FEAST, how physical compute nodes and mesh patches are connected
  Topology<> network; //not needed for dummy LB
  Topology<> patches; //only needed for letting dummy LB know, how many patches actually exist
  //patch topology: (will be done from file(s) )
  /*   0---1
   *   | X |
   *   2---3
   */
  patches.push_back();
  patches.at(0).push_back(1);
  patches.at(0).push_back(2);
  patches.at(0).push_back(3);
  patches.push_back();
  patches.at(1).push_back(0);
  patches.at(1).push_back(2);
  patches.at(1).push_back(3);
  patches.push_back();
  patches.at(2).push_back(0);
  patches.at(2).push_back(1);
  patches.at(2).push_back(3);
  patches.push_back();
  patches.at(3).push_back(0);
  patches.at(3).push_back(1);
  patches.at(3).push_back(2);

  CommStructures<Topology<> > lbconf(network, patches);

  //prepare local data for each process (i.t.m. 1 on 1 case only)
  PatchData<Mesh<rnt_2D>, Halo<0, pl_vertex, Mesh<rnt_2D, Topology<> > >, Topology<> > local_data;

  int me(0);
  //build communication structures
#ifndef SERIAL
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  Control<Parallel, SimpleLoadBalancingPolicy, SimpleDataFillPolicy>::init(lbconf, local_data, me);
#else
  Control<Serial, SimpleLoadBalancingPolicy, SimpleDataFillPolicy>::init(lbconf, local_data, me);
#endif

  //### begin test output only ###

  //test common data once
  if(me == 0)
  {
    for(Index i(0) ; i < lbconf.patch_process_map.size() ; ++i)
      for(Index j(0) ; j < lbconf.patch_process_map.at(i).size() ; ++j)
        std::cout << "Process " << lbconf.patch_process_map.at(i).at(j) << " distributed to patch " << i << "." << std::endl;

    for(Index i(0) ; i < lbconf.process_patch_map.size() ; ++i)
      for(Index j(0) ; j < lbconf.process_patch_map.at(i).size() ; ++j)
        std::cout << "Patch " << lbconf.process_patch_map.at(i).at(j) << " processed by process " << i << "." << std::endl;
  }

  //test local data
  std::cout << "Process " << me << " has: " << std::endl;
  std::cout << "    " << local_data.meshes_on_process_patch.size() << " meshes," << std::endl;
  for(Index i(0) ; i < local_data.meshes_on_process_patch.size() ; ++i)
  {
    std::cout << "        id: " << local_data.meshes_on_process_patch.at(i).get_id() << std::endl;
  }
  std::cout << "    " << local_data.halos_on_process_patch.size() << " halos," << std::endl;
  for(Index i(0) ; i < local_data.halos_on_process_patch.size() ; ++i)
  {
    std::cout << "        halo for comm with patch: " << local_data.halos_on_process_patch.at(i).get_other() << std::endl;
  }

  //### end test output only ###

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
