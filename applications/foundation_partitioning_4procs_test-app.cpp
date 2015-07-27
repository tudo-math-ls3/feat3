#include <kernel/base_header.hpp>

#include <kernel/foundation/comm_base.hpp>

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/archs.hpp>
#include <kernel/foundation/communication.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/foundation/attribute.hpp>
#include <kernel/foundation/topology.hpp>
#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/refinement.hpp>
#include <kernel/foundation/partitioning.hpp>
#include <kernel/foundation/mesh_control.hpp>
#include <kernel/foundation/halo_control.hpp>
#include <kernel/foundation/halo_interface.hpp>
#include <kernel/foundation/global_dot.hpp>
#include <kernel/foundation/global_synch_vec.hpp>
#include <kernel/foundation/global_product_mat_vec.hpp>
#include <kernel/foundation/global_defect.hpp>
#include <kernel/foundation/global_norm.hpp>
#include <kernel/foundation/gateway.hpp>
#include <kernel/foundation/aura.hpp>
#include <kernel/foundation/halo_frequencies.hpp>
#include <kernel/foundation/mesh_util.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/scarc/scarc_functor.hpp>
#include <kernel/scarc/matrix_conversion.hpp>

#include <iostream>
#include <limits>

using namespace FEAST;
using namespace Foundation;
using namespace Geometry;
using namespace ScaRC;

template<typename DT1_, typename DT2_, typename DT3_>
struct TestResult
{
  TestResult(DT1_ l, DT2_ r, DT3_ eps) :
    left(l),
    right(r),
    epsilon(eps)
  {
    //passed = std::abs(l - r) < eps ? true : false;
    passed = (l < r ? r - l : l - r) < eps ? true : false;
  }

  TestResult()
  {
  }

  DT1_ left;
  DT2_ right;
  DT3_ epsilon;
  bool passed;
};

template<typename DT1_, typename DT2_, typename DT3_>
TestResult<DT1_, DT2_, DT3_> test_check_equal_within_eps(DT1_ l, DT2_ r, DT3_ eps)
{
  return TestResult<DT1_, DT2_, DT3_>(l, r, eps);
}

void testmesh_hypercube_2D(Mesh<Dim2D>& target_mesh, std::vector<Attribute<double> >& attrs, std::vector<Halo<0, PLEdge, Mesh<Dim2D> > >& boundaries);

int main(int argc, char* argv[])
{
  int rank(0);
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  (void)argc;
  (void)argv;

  ///setup geometric problem data: mesh, its attributes and boundaries
  Mesh<Dim2D> mesh;
  std::vector<Attribute<double> > attrs;
  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries;
  testmesh_hypercube_2D(mesh, attrs, boundaries);

  if(!MeshUtil::iz_property_quad(mesh, attrs.at(0), attrs.at(1)))
    std::cout << "WARNING: rank " << rank << " " << " base mesh does not fulfill IZ-property!" << std::endl;

  ///provide memory for halos
  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  ///partitioning and initial loadbalancing
  auto p_i(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(mesh,
                                           boundaries,
                                           4, Index(rank),
                                           attrs
                                           ));

  //MeshUtil::establish_iz_property(p_i.basemesh, attrs.at(0), attrs.at(1));
  //MeshUtil::establish_iz_property(*((Mesh<Dim2D>*)(p_i.submesh.get())), *( (Attribute<double>*)(p_i.attrs.at(0).get()) ), *( (Attribute<double>*)(p_i.attrs.at(1).get()) ));

  if(!MeshUtil::iz_property_quad(*((Mesh<Dim2D>*)(p_i.submesh.get())), *( (Attribute<double>*)(p_i.attrs.at(0).get()) ), *( (Attribute<double>*)(p_i.attrs.at(1).get() ))))
    std::cout << "WARNING: rank " << rank << " " << " local mesh does not fulfill IZ-property!" << std::endl;

  if(!MeshUtil::iz_property_quad(p_i.basemesh, attrs.at(0), attrs.at(1)))
    std::cout << "WARNING: rank " << rank << " " << " global mesh (refined basemesh) does not fulfill IZ-property!" << std::endl;

  auto adjv_base(p_i.basemesh.get_adjacent_polytopes(pl_face, pl_vertex, p_i.submesh->get_map().at(0)));
  auto adjv_sub(p_i.submesh->get_adjacent_polytopes(pl_face, pl_vertex, 0));

  TestResult<double, double, double>* res = new TestResult<double, double, double>[adjv_base.size() * 2];
  for(Index i(0), iw(0) ; i < adjv_base.size(); ++i, iw+=2)
  {
    res[iw] = test_check_equal_within_eps(((Attribute<double>*)(p_i.attrs.at(0).get()))->at(adjv_sub.at(i)), attrs.at(0).at(adjv_base.at(i)), std::numeric_limits<double>::epsilon() * 1e8);
    res[iw + 1] = test_check_equal_within_eps(((Attribute<double>*)(p_i.attrs.at(1).get()))->at(adjv_sub.at(i)), attrs.at(1).at(adjv_base.at(i)), std::numeric_limits<double>::epsilon() * 1e8);
  }

  bool passed(true);
  for(Index i(0) ; i < adjv_base.size() * 2 ; ++i)
    if(!res[i].passed)
    {
      std::cout << "FAILED: rank " << rank << " " << res[i].left << " not within range (eps = " << res[i].epsilon << ") of " << res[i].right << "! (foundation partitioning test, FACE/VERTEX) " << std::endl;
    }

  if(passed)
    std::cout << "PASSED (rank " << rank <<"): (foundation partitioning test, FACE/VERTEX)" << std::endl;

  delete[] res;

  auto adjv_base_e(p_i.basemesh.get_adjacent_polytopes(pl_face, pl_edge, p_i.submesh->get_map().at(0)));
  auto adjv_sub_e(p_i.submesh->get_adjacent_polytopes(pl_face, pl_edge, 0));

  TestResult<double, double, double>* res_e = new TestResult<double, double, double>[adjv_sub_e.size() * 4];
  for(Index i(0), iw(0) ; i < adjv_base.size(); ++i, iw+=4)
  {
    auto adjv_base_ev(p_i.basemesh.get_adjacent_polytopes(pl_edge, pl_vertex, adjv_base_e.at(i)));
    auto adjv_sub_ev(p_i.submesh->get_adjacent_polytopes(pl_edge, pl_vertex, adjv_sub_e.at(i)));

    res_e[iw] =     test_check_equal_within_eps(((Attribute<double>*)(p_i.attrs.at(0).get()))->at(adjv_sub_ev.at(0)), attrs.at(0).at(adjv_base_ev.at(0)), std::numeric_limits<double>::epsilon() * 1e8);
    res_e[iw + 1] = test_check_equal_within_eps(((Attribute<double>*)(p_i.attrs.at(1).get()))->at(adjv_sub_ev.at(0)), attrs.at(1).at(adjv_base_ev.at(0)), std::numeric_limits<double>::epsilon() * 1e8);
    res_e[iw + 2] = test_check_equal_within_eps(((Attribute<double>*)(p_i.attrs.at(0).get()))->at(adjv_sub_ev.at(1)), attrs.at(0).at(adjv_base_ev.at(1)), std::numeric_limits<double>::epsilon() * 1e8);
    res_e[iw + 3] = test_check_equal_within_eps(((Attribute<double>*)(p_i.attrs.at(1).get()))->at(adjv_sub_ev.at(1)), attrs.at(1).at(adjv_base_ev.at(1)), std::numeric_limits<double>::epsilon() * 1e8);
  }

  bool passed_e(true);
  for(Index i(0) ; i < adjv_sub_e.size() * 4; ++i)
    if(!res_e[i].passed)
    {
      std::cout << "FAILED: rank " << rank << " " << res_e[i].left << " not within range (eps = " << res_e[i].epsilon << ") of " << res_e[i].right << "! (foundation partitioning test, EDGE/VERTEX) " << std::endl;
    }

  if(passed_e)
    std::cout << "PASSED (rank " << rank <<"): (foundation partitioning test, EDGE/VERTEX)" << std::endl;

  delete[] res_e;
#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}

void testmesh_hypercube_2D(Mesh<Dim2D>& target_mesh, std::vector<Attribute<double> >& attrs, std::vector<Halo<0, PLEdge, Mesh<Dim2D> > >& boundaries)
{
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));

  //setting up foundation mesh
  target_mesh.add_polytope(pl_vertex);
  target_mesh.add_polytope(pl_vertex);
  target_mesh.add_polytope(pl_vertex);
  target_mesh.add_polytope(pl_vertex);

  target_mesh.add_polytope(pl_edge);
  target_mesh.add_polytope(pl_edge);
  target_mesh.add_polytope(pl_edge);
  target_mesh.add_polytope(pl_edge);

  target_mesh.add_polytope(pl_face);

  target_mesh.add_adjacency(pl_face, pl_vertex, 0, 0);
  target_mesh.add_adjacency(pl_face, pl_vertex, 0, 1);
  target_mesh.add_adjacency(pl_face, pl_vertex, 0, 2);
  target_mesh.add_adjacency(pl_face, pl_vertex, 0, 3);

  target_mesh.add_adjacency(pl_edge, pl_vertex, 0, 0);
  target_mesh.add_adjacency(pl_edge, pl_vertex, 0, 1);

  target_mesh.add_adjacency(pl_edge, pl_vertex, 1, 2);
  target_mesh.add_adjacency(pl_edge, pl_vertex, 1, 3);

  target_mesh.add_adjacency(pl_edge, pl_vertex, 2, 0);
  target_mesh.add_adjacency(pl_edge, pl_vertex, 2, 2);

  target_mesh.add_adjacency(pl_edge, pl_vertex, 3, 1);
  target_mesh.add_adjacency(pl_edge, pl_vertex, 3, 3);

  /*target_mesh.add_adjacency(pl_vertex, pl_edge, 0, 0);
  target_mesh.add_adjacency(pl_vertex, pl_edge, 0, 2);
  target_mesh.add_adjacency(pl_vertex, pl_face, 0, 0);

  target_mesh.add_adjacency(pl_vertex, pl_edge, 1, 0);
  target_mesh.add_adjacency(pl_vertex, pl_edge, 1, 3);
  target_mesh.add_adjacency(pl_vertex, pl_face, 1, 0);

  target_mesh.add_adjacency(pl_vertex, pl_edge, 2, 1);
  target_mesh.add_adjacency(pl_vertex, pl_edge, 2, 2);
  target_mesh.add_adjacency(pl_vertex, pl_face, 2, 0);

  target_mesh.add_adjacency(pl_vertex, pl_edge, 3, 1);
  target_mesh.add_adjacency(pl_vertex, pl_edge, 3, 3);
  target_mesh.add_adjacency(pl_vertex, pl_face, 3, 0);*/

  boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D> >(target_mesh));
  boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D> >(target_mesh));
  boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D> >(target_mesh));
  boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D> >(target_mesh));
  boundaries.at(0).push_back(0);
  boundaries.at(1).push_back(1);
  boundaries.at(2).push_back(2);
  boundaries.at(3).push_back(3);
}
