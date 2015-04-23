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
#include <kernel/foundation/global_dot.hpp>
#include <kernel/foundation/global_synch_vec.hpp>
#include <kernel/foundation/global_product_mat_vec.hpp>
#include <kernel/foundation/global_defect.hpp>
#include <kernel/foundation/global_norm.hpp>
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
using namespace LAFEM;

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

void check_global_dot_2D(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  double result(0);
  GlobalDot<Mem::Main>::value(result, a, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 18., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (dot)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_2D (dot))" << std::endl;
}

void check_global_dot_2D_gateway(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalDotGateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(a.dot(b, &gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 18., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (dot by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_2D (dot))" << std::endl;
}

void check_global_nrm2_2D(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  double result(0);
  GlobalNorm2<Mem::Main>::value(result, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 6., std::numeric_limits<double>::epsilon()*10);
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (norm2)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

void check_global_nrm2_2D_gateway(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalNorm2Gateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 6., std::numeric_limits<double>::epsilon()*10.);
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (norm2 by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

void check_global_synchvec0_2D(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec0<Mem::Main>::exec(a,
                                                  mirrors,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs,
                                                  tags);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;
  res0 = test_check_equal_within_eps(a(0), rank == 3 ? 2. : 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 3 ? 1. : 2., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), rank == 3 ? 2. : 2., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), rank == 3 ? 4. : 4., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (synch_vec0)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;

}

void check_global_synchvec1_2D(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec1<Mem::Main>::exec(a,
                                                  mirrors,
                                                  frequencies,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs,
                                                  tags);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;
  res0 = test_check_equal_within_eps(a(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), 1., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), 1., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (synch_vec1)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;

}

void check_global_synchvec0_2D_gateway(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec0Gateway<Mem::Main, DenseVector<Mem::Main, double>, VectorMirror<Mem::Main, double> > gate(
                                                                                                                          mirrors,
                                                                                                                          other_ranks,
                                                                                                                          sendbufs,
                                                                                                                          recvbufs,
                                                                                                                          tags);
  a.synch0(&gate);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;

  res0 = test_check_equal_within_eps(a(0), rank == 3 ? 2. : 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 3 ? 1. : 2., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), rank == 3 ? 2. : 2., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), rank == 3 ? 4. : 4., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (synch_vec0 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D (synch_vec0by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway) ) " << std::endl;

}

void check_global_synchvec1_2D_gateway(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec1Gateway<Mem::Main, DenseVector<Mem::Main, double>, VectorMirror<Mem::Main, double> > gate(
                                                                                                                          mirrors,
                                                                                                                          frequencies,
                                                                                                                          other_ranks,
                                                                                                                          sendbufs,
                                                                                                                          recvbufs,
                                                                                                                          tags);
  a.synch1(&gate);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;

  res0 = test_check_equal_within_eps(a(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), 1., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), 1., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (synch_vec1 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;

}

void check_global_nrm2sqr_2D_gateway(Index rank)
{

  //create attributes for vertex coords
  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
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

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           4, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  DenseVector<Mem::Main, double> b(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  Index* polytopes_in_subset0 = new Index[3];
  Index* polytopes_in_subset1 = new Index[3];
  Index* polytopes_in_subset2 = new Index[3];

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }
  else
  {
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset0);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), polytopes_in_subset1);
    HaloControl<dim_2D>::fill_sizes(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), polytopes_in_subset2);
  }

  CellSubSet<Shape::Hypercube<2> > cell_sub_set0(polytopes_in_subset0);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set1(polytopes_in_subset1);
  CellSubSet<Shape::Hypercube<2> > cell_sub_set2(polytopes_in_subset2);

  if(rank == 0)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else if(rank == 1 || rank == 2)
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }
  else
  {
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set0);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(1).get())), cell_sub_set1);
    HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLEdge, Mesh<Dim2D> >*)(p0.comm_halos.at(2).get())), cell_sub_set2);
  }

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
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

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalNorm2SquaredGateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2sqr(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 36., std::numeric_limits<double>::epsilon()*10.);
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (norm2sqr by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
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
  check_global_dot_2D((Index)me);
  check_global_dot_2D_gateway((Index)me);
  check_global_nrm2_2D((Index)me);
  check_global_nrm2_2D_gateway((Index)me);
  check_global_synchvec0_2D((Index)me);
  check_global_synchvec1_2D((Index)me);
  check_global_synchvec0_2D_gateway((Index)me);
  check_global_synchvec1_2D_gateway((Index)me);
  check_global_nrm2sqr_2D_gateway((Index)me);
#else
  std::cout << "Parallel tests unavailable on sole process " << me << std::endl;
#endif

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
