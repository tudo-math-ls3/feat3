#include <kernel/base_header.hpp>

#include <kernel/foundation/comm_base.hpp>

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>
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
#include <kernel/foundation/gateway.hpp>
#include <kernel/foundation/aura.hpp>
#include <kernel/foundation/halo_frequencies.hpp>
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
#include <kernel/assembly/dirichlet_assembler.hpp>

#include <iostream>
#include <limits>

using namespace FEAST;
using namespace Foundation;
using namespace Geometry;

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

void check_global_dot_3D_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> mbuf2(mirrors.at(2).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));
  mirror_buffers.push_back(std::move(mbuf2));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  DenseVector<Mem::Main, double> fbuf2(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));
  freq_buffers.push_back(std::move(fbuf2));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  double result(0);
  GlobalDot<Mem::Main, Algo::Generic>::value(result, a, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 40., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_concave (dot)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_3D_concave (dot))" << std::endl;
}

void check_global_dot_3D_gateway_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> mbuf2(mirrors.at(2).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));
  mirror_buffers.push_back(std::move(mbuf2));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  DenseVector<Mem::Main, double> fbuf2(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));
  freq_buffers.push_back(std::move(fbuf2));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  GlobalDotGateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(a.dot(b, &gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 40., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_concave (dot by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_3D_concave (dot by gateway))" << std::endl;
}

void check_global_nrm2_3D_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> mbuf2(mirrors.at(2).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));
  mirror_buffers.push_back(std::move(mbuf2));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  DenseVector<Mem::Main, double> fbuf2(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));
  freq_buffers.push_back(std::move(fbuf2));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  double result(0);
  GlobalNorm2<Mem::Main, Algo::Generic>::value(result, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, sqrt(80.), std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_concave (norm2)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

void check_global_nrm2_3D_gateway_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> mbuf2(mirrors.at(2).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));
  mirror_buffers.push_back(std::move(mbuf2));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  DenseVector<Mem::Main, double> fbuf2(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));
  freq_buffers.push_back(std::move(fbuf2));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  GlobalNorm2Gateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, sqrt(80.), std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_cocave (norm2 by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

void check_global_synchvec0_3D_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  DenseVector<Mem::Main, double> sbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> sbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> sbuf2(mirrors.at(2).size());
  sendbufs.push_back(std::move(sbuf0));
  sendbufs.push_back(std::move(sbuf1));
  sendbufs.push_back(std::move(sbuf2));

  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> rbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> rbuf2(mirrors.at(2).size());
  recvbufs.push_back(std::move(rbuf0));
  recvbufs.push_back(std::move(rbuf1));
  recvbufs.push_back(std::move(rbuf2));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());
  other_ranks.push_back(p0.comm_halos.at(1)->get_other());
  other_ranks.push_back(p0.comm_halos.at(2)->get_other());

  GlobalSynchVec0<Mem::Main, Algo::Generic>::exec(a,
                                                  mirrors,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;
  TestResult<double, double, double> res4;
  TestResult<double, double, double> res5;
  TestResult<double, double, double> res6;
  TestResult<double, double, double> res7;
  res0 = test_check_equal_within_eps(a(0), rank == 0 ? 1. : ( rank == 1 ? 2. : ( rank == 2 ? 3. : 4.)), std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 0 ? 1. : ( rank == 1 ? 3. : ( rank == 2 ? 1. : 2.)), std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), rank == 0 ? 2. : ( rank == 1 ? 1. : ( rank == 2 ? 2. : 3.)), std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), rank == 0 ? 3. : ( rank == 1 ? 2. : ( rank == 2 ? 1. : 2.)), std::numeric_limits<double>::epsilon());
  res4 = test_check_equal_within_eps(a(4), rank == 0 ? 1. : ( rank == 1 ? 2. : ( rank == 2 ? 4. : 1.)), std::numeric_limits<double>::epsilon());
  res5 = test_check_equal_within_eps(a(5), rank == 0 ? 1. : ( rank == 1 ? 4. : ( rank == 2 ? 2. : 1.)), std::numeric_limits<double>::epsilon());
  res6 = test_check_equal_within_eps(a(6), rank == 0 ? 2. : ( rank == 1 ? 1. : ( rank == 2 ? 3. : 1.)), std::numeric_limits<double>::epsilon());
  res7 = test_check_equal_within_eps(a(7), rank == 0 ? 4. : ( rank == 1 ? 3. : ( rank == 2 ? 2. : 1.)), std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed && res4.passed && res5.passed && res6.passed && res7.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_concave (synch_vec0)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_3D_concave (synch_vec0) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_3D_concave (synch_vec0) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_3D_concave (synch_vec0) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_3D_concave (synch_vec0) ) " << std::endl;
  else if(!res4.passed)
    std::cout << "FAILED: " << res4.left << " not within range (eps = " << res4.epsilon << ") of " << res4.right << "! (foundation_global_ops_test_3D_concave (synch_vec0) ) " << std::endl;
  else if(!res5.passed)
    std::cout << "FAILED: " << res5.left << " not within range (eps = " << res5.epsilon << ") of " << res5.right << "! (foundation_global_ops_test_3D_concave (synch_vec0) ) " << std::endl;
  else if(!res6.passed)
    std::cout << "FAILED: " << res6.left << " not within range (eps = " << res6.epsilon << ") of " << res6.right << "! (foundation_global_ops_test_3D_concave (synch_vec0) ) " << std::endl;
  else if(!res7.passed)
    std::cout << "FAILED: " << res7.left << " not within range (eps = " << res7.epsilon << ") of " << res7.right << "! (foundation_global_ops_test_3D_concave (synch_vec0) ) " << std::endl;
}

void check_global_synchvec1_3D_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  DenseVector<Mem::Main, double> sbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> sbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> sbuf2(mirrors.at(2).size());
  sendbufs.push_back(std::move(sbuf0));
  sendbufs.push_back(std::move(sbuf1));
  sendbufs.push_back(std::move(sbuf2));

  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> rbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> rbuf2(mirrors.at(2).size());
  recvbufs.push_back(std::move(rbuf0));
  recvbufs.push_back(std::move(rbuf1));
  recvbufs.push_back(std::move(rbuf2));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());
  other_ranks.push_back(p0.comm_halos.at(1)->get_other());
  other_ranks.push_back(p0.comm_halos.at(2)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> mbuf2(mirrors.at(2).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));
  mirror_buffers.push_back(std::move(mbuf2));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  DenseVector<Mem::Main, double> fbuf2(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));
  freq_buffers.push_back(std::move(fbuf2));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  GlobalSynchVec1<Mem::Main, Algo::Generic>::exec(a,
                                                  mirrors,
                                                  frequencies,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;
  TestResult<double, double, double> res4;
  TestResult<double, double, double> res5;
  TestResult<double, double, double> res6;
  TestResult<double, double, double> res7;
  res0 = test_check_equal_within_eps(a(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), 1., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), 1., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), 1., std::numeric_limits<double>::epsilon());
  res4 = test_check_equal_within_eps(a(4), 1., std::numeric_limits<double>::epsilon());
  res5 = test_check_equal_within_eps(a(5), 1., std::numeric_limits<double>::epsilon());
  res6 = test_check_equal_within_eps(a(6), 1., std::numeric_limits<double>::epsilon());
  res7 = test_check_equal_within_eps(a(7), 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed && res4.passed && res5.passed && res6.passed && res7.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_concave (synch_vec1)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_3D_concave (synch_vec1) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_3D_concave (synch_vec1) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_3D_concave (synch_vec1) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_3D_concave (synch_vec1) ) " << std::endl;
  else if(!res4.passed)
    std::cout << "FAILED: " << res4.left << " not within range (eps = " << res4.epsilon << ") of " << res4.right << "! (foundation_global_ops_test_3D_concave (synch_vec1) ) " << std::endl;
  else if(!res5.passed)
    std::cout << "FAILED: " << res5.left << " not within range (eps = " << res5.epsilon << ") of " << res5.right << "! (foundation_global_ops_test_3D_concave (synch_vec1) ) " << std::endl;
  else if(!res6.passed)
    std::cout << "FAILED: " << res6.left << " not within range (eps = " << res6.epsilon << ") of " << res6.right << "! (foundation_global_ops_test_3D_concave (synch_vec1) ) " << std::endl;
  else if(!res7.passed)
    std::cout << "FAILED: " << res7.left << " not within range (eps = " << res7.epsilon << ") of " << res7.right << "! (foundation_global_ops_test_3D_concave (synch_vec1) ) " << std::endl;
}

void check_global_nrm2sqr_3D_gateway_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> mbuf2(mirrors.at(2).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));
  mirror_buffers.push_back(std::move(mbuf2));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  DenseVector<Mem::Main, double> fbuf2(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));
  freq_buffers.push_back(std::move(fbuf2));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  GlobalNorm2SquaredGateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2sqr(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 80., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_concave (norm2sqr by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

void check_global_synchvec0_3D_gateway_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  DenseVector<Mem::Main, double> sbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> sbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> sbuf2(mirrors.at(2).size());
  sendbufs.push_back(std::move(sbuf0));
  sendbufs.push_back(std::move(sbuf1));
  sendbufs.push_back(std::move(sbuf2));

  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> rbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> rbuf2(mirrors.at(2).size());
  recvbufs.push_back(std::move(rbuf0));
  recvbufs.push_back(std::move(rbuf1));
  recvbufs.push_back(std::move(rbuf2));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());
  other_ranks.push_back(p0.comm_halos.at(1)->get_other());
  other_ranks.push_back(p0.comm_halos.at(2)->get_other());

  GlobalSynchVec0Gateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double>, VectorMirror<Mem::Main, double> > gate(
                                                                                                                          mirrors,
                                                                                                                          other_ranks,
                                                                                                                          sendbufs,
                                                                                                                          recvbufs);
  a.synch0(&gate);


  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;
  TestResult<double, double, double> res4;
  TestResult<double, double, double> res5;
  TestResult<double, double, double> res6;
  TestResult<double, double, double> res7;
  res0 = test_check_equal_within_eps(a(0), rank == 0 ? 1. : ( rank == 1 ? 2. : ( rank == 2 ? 3. : 4.)), std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 0 ? 1. : ( rank == 1 ? 3. : ( rank == 2 ? 1. : 2.)), std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), rank == 0 ? 2. : ( rank == 1 ? 1. : ( rank == 2 ? 2. : 3.)), std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), rank == 0 ? 3. : ( rank == 1 ? 2. : ( rank == 2 ? 1. : 2.)), std::numeric_limits<double>::epsilon());
  res4 = test_check_equal_within_eps(a(4), rank == 0 ? 1. : ( rank == 1 ? 2. : ( rank == 2 ? 4. : 1.)), std::numeric_limits<double>::epsilon());
  res5 = test_check_equal_within_eps(a(5), rank == 0 ? 1. : ( rank == 1 ? 4. : ( rank == 2 ? 2. : 1.)), std::numeric_limits<double>::epsilon());
  res6 = test_check_equal_within_eps(a(6), rank == 0 ? 2. : ( rank == 1 ? 1. : ( rank == 2 ? 3. : 1.)), std::numeric_limits<double>::epsilon());
  res7 = test_check_equal_within_eps(a(7), rank == 0 ? 4. : ( rank == 1 ? 3. : ( rank == 2 ? 2. : 1.)), std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed && res4.passed && res5.passed && res6.passed && res7.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_concave (synch_vec0 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_3D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_3D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_3D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_3D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res4.passed)
    std::cout << "FAILED: " << res4.left << " not within range (eps = " << res4.epsilon << ") of " << res4.right << "! (foundation_global_ops_test_3D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res5.passed)
    std::cout << "FAILED: " << res5.left << " not within range (eps = " << res5.epsilon << ") of " << res5.right << "! (foundation_global_ops_test_3D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res6.passed)
    std::cout << "FAILED: " << res6.left << " not within range (eps = " << res6.epsilon << ") of " << res6.right << "! (foundation_global_ops_test_3D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res7.passed)
    std::cout << "FAILED: " << res7.left << " not within range (eps = " << res7.epsilon << ") of " << res7.right << "! (foundation_global_ops_test_3D_concave (synch_vec0 by gateway) ) " << std::endl;
}

void check_global_synchvec1_3D_gateway_concave(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(2));
  attrs.at(2).get_data().push_back(double(2));

  //creating foundation mesh
  Mesh<Dim3D> m(0);
  for (Index i(0); i < 20; i++)
    m.add_polytope(pl_vertex);

  for (Index i(0); i < 36; i++)
    m.add_polytope(pl_edge);

  for (Index i(0); i < 21; i++)
    m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);
  m.add_polytope(pl_polyhedron);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 5);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 6);
  m.add_adjacency(pl_vertex, pl_edge, 1, 11);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 3);
  m.add_adjacency(pl_vertex, pl_face, 1, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 5);
  m.add_adjacency(pl_vertex, pl_edge, 2, 7);
  m.add_adjacency(pl_vertex, pl_edge, 2, 12);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 1);
  m.add_adjacency(pl_vertex, pl_face, 2, 4);
  m.add_adjacency(pl_vertex, pl_face, 2, 8);
  m.add_adjacency(pl_vertex, pl_face, 2, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 1);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 2);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_edge, 3, 8);
  m.add_adjacency(pl_vertex, pl_edge, 3, 13);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);
  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 4);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);
  m.add_adjacency(pl_vertex, pl_face, 3, 9);
  m.add_adjacency(pl_vertex, pl_face, 3, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 2);

  m.add_adjacency(pl_vertex, pl_edge, 4, 2);
  m.add_adjacency(pl_vertex, pl_edge, 4, 9);
  m.add_adjacency(pl_vertex, pl_edge, 4, 14);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 5);
  m.add_adjacency(pl_vertex, pl_face, 4, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 2);

  m.add_adjacency(pl_vertex, pl_edge, 5, 3);
  m.add_adjacency(pl_vertex, pl_edge, 5, 7);
  m.add_adjacency(pl_vertex, pl_edge, 5, 15);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 6);
  m.add_adjacency(pl_vertex, pl_face, 5, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 3);
  m.add_adjacency(pl_vertex, pl_edge, 6, 4);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_edge, 6, 16);
  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);
  m.add_adjacency(pl_vertex, pl_face, 6, 6);
  m.add_adjacency(pl_vertex, pl_face, 6, 7);
  m.add_adjacency(pl_vertex, pl_face, 6, 11);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 4);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 17);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);
  m.add_adjacency(pl_vertex, pl_face, 7, 7);
  m.add_adjacency(pl_vertex, pl_face, 7, 12);
  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 2);

  m.add_adjacency(pl_vertex, pl_edge, 8, 10);
  m.add_adjacency(pl_vertex, pl_edge, 8, 18);
  m.add_adjacency(pl_vertex, pl_edge, 8, 23);
  m.add_adjacency(pl_vertex, pl_face, 8, 13);
  m.add_adjacency(pl_vertex, pl_face, 8, 3);
  m.add_adjacency(pl_vertex, pl_face, 8, 8);
  m.add_adjacency(pl_vertex, pl_polyhedron, 8, 0);

  m.add_adjacency(pl_vertex, pl_edge, 9, 11);
  m.add_adjacency(pl_vertex, pl_edge, 9, 18);
  m.add_adjacency(pl_vertex, pl_edge, 9, 24);
  m.add_adjacency(pl_vertex, pl_face, 9, 13);
  m.add_adjacency(pl_vertex, pl_face, 9, 3);
  m.add_adjacency(pl_vertex, pl_face, 9, 9);
  m.add_adjacency(pl_vertex, pl_polyhedron, 9, 0);

  m.add_adjacency(pl_vertex, pl_edge, 10, 12);
  m.add_adjacency(pl_vertex, pl_edge, 10, 19);
  m.add_adjacency(pl_vertex, pl_edge, 10, 23);
  m.add_adjacency(pl_vertex, pl_edge, 10, 25);
  m.add_adjacency(pl_vertex, pl_face, 10, 13);
  m.add_adjacency(pl_vertex, pl_face, 10, 14);
  m.add_adjacency(pl_vertex, pl_face, 10, 4);
  m.add_adjacency(pl_vertex, pl_face, 10, 8);
  m.add_adjacency(pl_vertex, pl_face, 10, 10);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 10, 1);

  m.add_adjacency(pl_vertex, pl_edge, 11, 13);
  m.add_adjacency(pl_vertex, pl_edge, 11, 19);
  m.add_adjacency(pl_vertex, pl_edge, 11, 20);
  m.add_adjacency(pl_vertex, pl_edge, 11, 24);
  m.add_adjacency(pl_vertex, pl_edge, 11, 26);
  m.add_adjacency(pl_vertex, pl_edge, 11, 28);
  m.add_adjacency(pl_vertex, pl_face, 11, 4);
  m.add_adjacency(pl_vertex, pl_face, 11, 5);
  m.add_adjacency(pl_vertex, pl_face, 11, 9);
  m.add_adjacency(pl_vertex, pl_face, 11, 11);
  m.add_adjacency(pl_vertex, pl_face, 11, 13);
  m.add_adjacency(pl_vertex, pl_face, 11, 14);
  m.add_adjacency(pl_vertex, pl_face, 11, 15);
  m.add_adjacency(pl_vertex, pl_face, 11, 16);
  m.add_adjacency(pl_vertex, pl_face, 11, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 0);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 11, 3);

  m.add_adjacency(pl_vertex, pl_edge, 12, 14);
  m.add_adjacency(pl_vertex, pl_edge, 12, 20);
  m.add_adjacency(pl_vertex, pl_edge, 12, 27);
  m.add_adjacency(pl_vertex, pl_edge, 12, 29);
  m.add_adjacency(pl_vertex, pl_face, 12, 5);
  m.add_adjacency(pl_vertex, pl_face, 12, 12);
  m.add_adjacency(pl_vertex, pl_face, 12, 15);
  m.add_adjacency(pl_vertex, pl_face, 12, 16);
  m.add_adjacency(pl_vertex, pl_face, 12, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 12, 3);

  m.add_adjacency(pl_vertex, pl_edge, 13, 15);
  m.add_adjacency(pl_vertex, pl_edge, 13, 21);
  m.add_adjacency(pl_vertex, pl_edge, 13, 25);
  m.add_adjacency(pl_vertex, pl_face, 13, 6);
  m.add_adjacency(pl_vertex, pl_face, 13, 10);
  m.add_adjacency(pl_vertex, pl_face, 13, 14);
  m.add_adjacency(pl_vertex, pl_polyhedron, 13, 1);

  m.add_adjacency(pl_vertex, pl_edge, 14, 16);
  m.add_adjacency(pl_vertex, pl_edge, 14, 21);
  m.add_adjacency(pl_vertex, pl_edge, 14, 22);
  m.add_adjacency(pl_vertex, pl_edge, 14, 26);
  m.add_adjacency(pl_vertex, pl_edge, 14, 30);
  m.add_adjacency(pl_vertex, pl_face, 14, 6);
  m.add_adjacency(pl_vertex, pl_face, 14, 7);
  m.add_adjacency(pl_vertex, pl_face, 14, 11);
  m.add_adjacency(pl_vertex, pl_face, 14, 14);
  m.add_adjacency(pl_vertex, pl_face, 14, 15);
  m.add_adjacency(pl_vertex, pl_face, 14, 17);
  m.add_adjacency(pl_vertex, pl_face, 14, 18);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 1);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 14, 3);

  m.add_adjacency(pl_vertex, pl_edge, 15, 17);
  m.add_adjacency(pl_vertex, pl_edge, 15, 22);
  m.add_adjacency(pl_vertex, pl_edge, 15, 27);
  m.add_adjacency(pl_vertex, pl_edge, 15, 31);
  m.add_adjacency(pl_vertex, pl_face, 15, 7);
  m.add_adjacency(pl_vertex, pl_face, 15, 12);
  m.add_adjacency(pl_vertex, pl_face, 15, 15);
  m.add_adjacency(pl_vertex, pl_face, 15, 17);
  m.add_adjacency(pl_vertex, pl_face, 15, 19);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 2);
  m.add_adjacency(pl_vertex, pl_polyhedron, 15, 3);

  m.add_adjacency(pl_vertex, pl_edge, 16, 28);
  m.add_adjacency(pl_vertex, pl_edge, 16, 32);
  m.add_adjacency(pl_vertex, pl_edge, 16, 34);
  m.add_adjacency(pl_vertex, pl_face, 16, 16);
  m.add_adjacency(pl_vertex, pl_face, 16, 18);
  m.add_adjacency(pl_vertex, pl_face, 16, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 16, 3);

  m.add_adjacency(pl_vertex, pl_edge, 17, 29);
  m.add_adjacency(pl_vertex, pl_edge, 17, 32);
  m.add_adjacency(pl_vertex, pl_edge, 17, 35);
  m.add_adjacency(pl_vertex, pl_face, 17, 16);
  m.add_adjacency(pl_vertex, pl_face, 17, 19);
  m.add_adjacency(pl_vertex, pl_face, 17, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 17, 3);

  m.add_adjacency(pl_vertex, pl_edge, 18, 30);
  m.add_adjacency(pl_vertex, pl_edge, 18, 33);
  m.add_adjacency(pl_vertex, pl_edge, 18, 34);
  m.add_adjacency(pl_vertex, pl_face, 18, 17);
  m.add_adjacency(pl_vertex, pl_face, 18, 18);
  m.add_adjacency(pl_vertex, pl_face, 18, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 18, 3);

  m.add_adjacency(pl_vertex, pl_edge, 19, 31);
  m.add_adjacency(pl_vertex, pl_edge, 19, 33);
  m.add_adjacency(pl_vertex, pl_edge, 19, 35);
  m.add_adjacency(pl_vertex, pl_face, 19, 17);
  m.add_adjacency(pl_vertex, pl_face, 19, 19);
  m.add_adjacency(pl_vertex, pl_face, 19, 20);
  m.add_adjacency(pl_vertex, pl_polyhedron, 19, 3);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(8);
  DenseVector<Mem::Main, double> b(8);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1),  attrs.at(2));

  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(1).get()));
  auto cell_sub_set2(HaloInterface<0, Dim3D>::convert(p0.comm_halos.at(2).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);
  VectorMirror<Mem::Main, double> target_mirror2;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror2, space, cell_sub_set2);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));
  mirrors.push_back(std::move(target_mirror2));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  DenseVector<Mem::Main, double> sbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> sbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> sbuf2(mirrors.at(2).size());
  sendbufs.push_back(std::move(sbuf0));
  sendbufs.push_back(std::move(sbuf1));
  sendbufs.push_back(std::move(sbuf2));

  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> rbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> rbuf2(mirrors.at(2).size());
  recvbufs.push_back(std::move(rbuf0));
  recvbufs.push_back(std::move(rbuf1));
  recvbufs.push_back(std::move(rbuf2));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());
  other_ranks.push_back(p0.comm_halos.at(1)->get_other());
  other_ranks.push_back(p0.comm_halos.at(2)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> mbuf2(mirrors.at(2).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));
  mirror_buffers.push_back(std::move(mbuf2));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  DenseVector<Mem::Main, double> fbuf2(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));
  freq_buffers.push_back(std::move(fbuf2));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  GlobalSynchVec1Gateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double>, VectorMirror<Mem::Main, double> > gate(
                                                                                                                          mirrors,
                                                                                                                          frequencies,
                                                                                                                          other_ranks,
                                                                                                                          sendbufs,
                                                                                                                          recvbufs);
  a.synch1(&gate);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;
  TestResult<double, double, double> res4;
  TestResult<double, double, double> res5;
  TestResult<double, double, double> res6;
  TestResult<double, double, double> res7;
  res0 = test_check_equal_within_eps(a(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), 1., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), 1., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), 1., std::numeric_limits<double>::epsilon());
  res4 = test_check_equal_within_eps(a(4), 1., std::numeric_limits<double>::epsilon());
  res5 = test_check_equal_within_eps(a(5), 1., std::numeric_limits<double>::epsilon());
  res6 = test_check_equal_within_eps(a(6), 1., std::numeric_limits<double>::epsilon());
  res7 = test_check_equal_within_eps(a(7), 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed && res4.passed && res5.passed && res6.passed && res7.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D_concave (synch_vec1 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_3D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_3D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_3D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_3D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res4.passed)
    std::cout << "FAILED: " << res4.left << " not within range (eps = " << res4.epsilon << ") of " << res4.right << "! (foundation_global_ops_test_3D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res5.passed)
    std::cout << "FAILED: " << res5.left << " not within range (eps = " << res5.epsilon << ") of " << res5.right << "! (foundation_global_ops_test_3D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res6.passed)
    std::cout << "FAILED: " << res6.left << " not within range (eps = " << res6.epsilon << ") of " << res6.right << "! (foundation_global_ops_test_3D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res7.passed)
    std::cout << "FAILED: " << res7.left << " not within range (eps = " << res7.epsilon << ") of " << res7.right << "! (foundation_global_ops_test_3D_concave (synch_vec1 by gateway) ) " << std::endl;
}

int main(int argc, char* argv[])
{
  int me(0);
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
  (void)argc;
  (void)argv;
  std::cout<<"CTEST_FULL_OUTPUT"<<std::endl;

#ifndef SERIAL
  check_global_dot_3D_concave((Index)me);
  check_global_dot_3D_gateway_concave((Index)me);
  check_global_nrm2_3D_concave((Index)me);
  check_global_nrm2_3D_gateway_concave((Index)me);
  check_global_synchvec0_3D_concave((Index)me);
  check_global_synchvec1_3D_concave((Index)me);
  check_global_nrm2sqr_3D_gateway_concave((Index)me);
  check_global_synchvec0_3D_gateway_concave((Index)me);
  check_global_synchvec1_3D_gateway_concave((Index)me);
#else
  std::cout << "Parallel tests unavailable on sole process " << me << std::endl;
#endif

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
