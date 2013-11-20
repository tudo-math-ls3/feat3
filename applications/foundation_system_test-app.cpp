#ifndef SERIAL
#  include <mpi.h>
#endif

#include <kernel/base_header.hpp>
#include <kernel/foundation/base.hpp>
#include <kernel/foundation/attribute.hpp>
#include <kernel/foundation/topology.hpp>
#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/foundation/partitioning.hpp>
#include <kernel/foundation/refinement.hpp>
#include <kernel/foundation/load_balancing.hpp>
#include <vector>
#include <iostream>

using namespace FEAST;
using namespace Foundation;
using namespace std;

int main(int argc, char* argv[])
{

#ifndef SERIAL
  MPI_Init(&argc, &argv);
#endif
  (void)argc;
  (void)argv;
  std::cout<<"CTEST_FULL_OUTPUT"<<std::endl;

  int me(0);
  int size(1);
#ifndef SERIAL
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  typedef double DT;
  typedef Attribute<DT, vector> AttrT;
  typedef vector<AttrT> AttrStoreT;
  typedef Topology<Index, vector, vector<Index> > TopoT;
  typedef Mesh<Dim2D, TopoT, vector> MeshT;
  typedef Halo<0, PLEdge, MeshT> BoundaryHaloT;
  typedef vector<BoundaryHaloT> BoundaryHaloStoreT;

  /*(0,1) (1,1)
   *  *----*
   *  |    |
   *  |    |
   *  *----*
   *(0,0) (1,0)
   */

  //create attributes for vertex coords
  AttrStoreT attrs;
  attrs.push_back(AttrT()); //vertex x-coords
  attrs.push_back(AttrT()); //vertex y-coords

  attrs.at(0).get_data().push_back(DT(0));
  attrs.at(1).get_data().push_back(DT(0));

  attrs.at(0).get_data().push_back(DT(1));
  attrs.at(1).get_data().push_back(DT(0));

  attrs.at(0).get_data().push_back(DT(0));
  attrs.at(1).get_data().push_back(DT(1));

  attrs.at(0).get_data().push_back(DT(1));
  attrs.at(1).get_data().push_back(DT(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  MeshT m(0);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);

  MeshT m_fine(m);

  BoundaryHaloStoreT boundaries;
  boundaries.push_back(BoundaryHaloT(m));
  boundaries.push_back(BoundaryHaloT(m));
  boundaries.push_back(BoundaryHaloT(m));
  boundaries.push_back(BoundaryHaloT(m));
  boundaries.at(0).push_back(0);
  boundaries.at(1).push_back(1);
  boundaries.at(2).push_back(2);
  boundaries.at(3).push_back(3);

  Index num_procs((Index)size);
  Index rank((Index)me);
  //Index level(4);

  std::cout << "p_" << me << " #procs: " << num_procs << std::endl;
  PData<Dim2D, TopoT, vector, Mesh, DT> p0(Partitioning<Mem::Main,
                                                        Algo::Generic,
                                                        Dim2D,
                                                        0,
                                                        pl_vertex>::execute(m,
                                                                            boundaries,
                                                                            num_procs,
                                                                            rank,
                                                                            attrs));

  std::cout << "p_" << me << " about to LB" << std::endl;
  LoadBalancing<LBPUniformCompScaledComm<DT> >::execute(p0);

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
