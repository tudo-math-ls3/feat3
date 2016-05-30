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
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>

#include <iostream>
#include <limits>

using namespace FEAT;
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

void check_global_dot_1D(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  //Foundation::HaloControl<Foundation::dim_1D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  //MeshPart<ConformalMesh<Shape::Hypercube<1> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_1D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  double result(0);
  GlobalDot<Mem::Main>::value(result, a, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 6., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (dot)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_1D (dot) ) " << std::endl;

  delete[] size_set;
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
                                           2, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(6);
  DenseVector<Mem::Main, double> b(6);
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

  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set0);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

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

  delete[] size_set;
}

void check_global_dot_3D(Index rank)
{
  /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;

  Refinement<Mem::Main,
             mrt_standard>::execute(m_fine, &halos, attrs);

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(18);
  DenseVector<Mem::Main, double> b(18);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1), attrs.at(2));

  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set0(HaloInterface<0, Dim3D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set0);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  double result(0);
  GlobalDot<Mem::Main>::value(result, a, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 54., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (dot)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_3D (dot))" << std::endl;

  delete[] size_set;
}

void check_global_synchvec0_1D(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec0<Mem::Main>::exec(a,
                                                  mirrors,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs,
                                                  tags);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  res0 = test_check_equal_within_eps(a(0), rank == 0 ? 1. : 2., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 0 ? 2. : 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (synch_vec0)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_1D (synch_vec0) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_1D (synch_vec0) ) " << std::endl;

  delete[] size_set;
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

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(6);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1));

  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim2D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

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
  TestResult<double, double, double> res4;
  TestResult<double, double, double> res5;
  res0 = test_check_equal_within_eps(a(0), rank == 0 ? 1. : 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 0 ? 1. : 1., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), rank == 0 ? 2. : 1., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), rank == 0 ? 1. : 2., std::numeric_limits<double>::epsilon());
  res4 = test_check_equal_within_eps(a(4), rank == 0 ? 2. : 2., std::numeric_limits<double>::epsilon());
  res5 = test_check_equal_within_eps(a(5), rank == 0 ? 2. : 2., std::numeric_limits<double>::epsilon());

  if(res0.passed && res1.passed && res2.passed && res3.passed && res4.passed && res5.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (synch_vec0)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (epsll = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;
  else if(!res4.passed)
    std::cout << "FAILED: " << res4.left << " not within range (eps = " << res4.epsilon << ") of " << res4.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;
  else if(!res5.passed)
    std::cout << "FAILED: " << res5.left << " not within range (eps = " << res5.epsilon << ") of " << res5.right << "! (foundation_global_ops_test_2D (synch_vec0) ) " << std::endl;

  delete[] size_set;
}

void check_global_synchvec0_3D(Index rank)
{
   /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;


  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(18);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1), attrs.at(2));

  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim3D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec0<Mem::Main>::exec(a,
                                                  mirrors,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs,
                                                  tags);

  TestResult<double, double, double> res[18];
  const double ref[18][2] =
      {
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 2., 2.},
        { 1., 2.},
        { 2., 1.},
        { 2., 2.},
        { 2., 2.},
        { 1., 1.},
        { 2., 2.},
        { 2., 2.},
        { 2., 2.},
        { 2., 2.},
        { 2., 2.}
      };

  for( Index i(0); i < 18; i++)
  {
  res[i] = test_check_equal_within_eps(a(i), rank == 0 ? ref[i][0] : ref[i][1], std::numeric_limits<double>::epsilon());
  }
  if( res[0].passed && res[1].passed && res[2].passed && res[3].passed && res[4].passed && res[5].passed && res[6].passed && res[7].passed && res[8].passed && res[9].passed && res[10].passed && res[11].passed && res[12].passed && res[13].passed && res[14].passed && res[15].passed && res[16].passed && res[17].passed )
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (synch_vec0)" << std::endl;
  else
  {
    for( Index j(0); j < 18; j++)
      if(!res[j].passed)
        std::cout << "FAILED: " << res[j].left << " not within range (eps = " << res[j].epsilon << ") of " << res[j].right << "! (foundation_global_ops_test_3D (synch_vec0) ) " << std::endl;
  }

  delete[] size_set;
}

void check_global_synchvec1_1D(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

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
  res0 = test_check_equal_within_eps(a(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (synch_vec1)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_1D (synch_vec1) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_1D (synch_vec1) ) " << std::endl;

  delete[] size_set;
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

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(6);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1));

  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim2D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

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

  TestResult<double, double, double> res0, res1, res2, res3, res4, res5;
  res0 = test_check_equal_within_eps(a(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), 1., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), 1., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), 1., std::numeric_limits<double>::epsilon());
  res4 = test_check_equal_within_eps(a(4), 1., std::numeric_limits<double>::epsilon());
  res5 = test_check_equal_within_eps(a(5), 1., std::numeric_limits<double>::epsilon());

  if(res0.passed && res1.passed && res2.passed && res3.passed && res4.passed && res5.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (synch_vec1)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;
  else if(!res4.passed)
    std::cout << "FAILED: " << res4.left << " not within range (eps = " << res4.epsilon << ") of " << res4.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;
  else if(!res5.passed)
    std::cout << "FAILED: " << res5.left << " not within range (eps = " << res5.epsilon << ") of " << res5.right << "! (foundation_global_ops_test_2D (synch_vec1) ) " << std::endl;

  delete[] size_set;
}

void check_global_synchvec1_3D(Index rank)
{
   /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;


  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(18);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1), attrs.at(2));

  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim3D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

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

  TestResult<double, double, double> res[18];

  for( Index i(0); i < 18; i++)
  {
  res[i] = test_check_equal_within_eps(a(i), 1., std::numeric_limits<double>::epsilon());
  }
  if( res[0].passed && res[1].passed && res[2].passed && res[3].passed && res[4].passed && res[5].passed && res[6].passed && res[7].passed && res[8].passed && res[9].passed && res[10].passed && res[11].passed && res[12].passed && res[13].passed && res[14].passed && res[15].passed && res[16].passed && res[17].passed )
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (synch_vec1)" << std::endl;
  else
  {
    for( Index j(0); j < 18; j++)
      if(!res[j].passed)
        std::cout << "FAILED: " << res[j].left << " not within range (eps = " << res[j].epsilon << ") of " << res[j].right << "! (foundation_global_ops_test_3D (synch_vec1) ) " << std::endl;
  }

  delete[] size_set;
}

void check_global_product_mat_vec_1D(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  SparseMatrixCSR<Mem::Main, double> mat_sys;
  Assembly::SymbolicAssembler::assemble_matrix_std1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  DenseVector<Mem::Main, double> result(a.size(), double(0));
  auto tags(HaloTags::value(p0.comm_halos));
  GlobalProductMat0Vec1<Mem::Main>::exec(
                                                        result,
                                                        mat_sys,
                                                        a,
                                                        mirrors,
                                                        other_ranks,
                                                        sendbufs,
                                                        recvbufs,
                                                        tags);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  res0 = test_check_equal_within_eps(result(0), 0., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(result(1), 0., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (product_mat0vec1)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_1D (product_mat0vec1) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_1D (product_mat0vec1) ) " << std::endl;

  delete[] size_set;
}

void check_global_defect_1D(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  SparseMatrixCSR<Mem::Main, double> mat_sys;
  Assembly::SymbolicAssembler::assemble_matrix_std1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  DenseVector<Mem::Main, double> result(a.size(), double(0));
  auto tags(HaloTags::value(p0.comm_halos));
  GlobalDefect<Mem::Main>::exec(
                                               result,
                                               b,
                                               mat_sys,
                                               a,
                                               mirrors,
                                               other_ranks,
                                               sendbufs,
                                               recvbufs,
                                               tags);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  res0 = test_check_equal_within_eps(result(0), 2., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(result(1), 2., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (defect)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_1D (defect) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_1D (defect) ) " << std::endl;

  delete[] size_set;
}

void check_global_dot_1D_gateway(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  //Foundation::HaloControl<Foundation::dim_1D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<1> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_1D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalDotGateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(a.dot(b, &gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 6., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (dot by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_1D (dot by gateway) ) " << std::endl;

  delete[] size_set;
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

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(6);
  DenseVector<Mem::Main, double> b(6);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1));

  //Foundation::HaloControl<Foundation::dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim2D>::convert(h_new));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalDotGateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(a.dot(b, &gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 18., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (dot by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_2D (dot by gateway) ) " << std::endl;

  delete[] size_set;
}

void check_global_dot_3D_gateway(Index rank)
{
   /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;


  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(18);
  DenseVector<Mem::Main, double> b(18);
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

  //Foundation::HaloControl<Foundation::dim_3D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim3D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim3D>::convert(h_new));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<3> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_3D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim3D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalDotGateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(a.dot(b, &gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 54., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (dot by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_3D (dot by gateway) ) " << std::endl;

  delete[] size_set;
}

void check_global_nrm2_1D(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  //Foundation::HaloControl<Foundation::dim_1D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<1> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_1D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  double result(0);
  GlobalNorm2<Mem::Main>::value(result, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, std::sqrt(12.), std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (nrm2)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_1D (nrm2) ) " << std::endl;

  delete[] size_set;
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

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(6);
  DenseVector<Mem::Main, double> b(6);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1));

  //Foundation::HaloControl<Foundation::dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim2D>::convert(h_new));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  double result(0);
  GlobalNorm2<Mem::Main>::value(result, b, frequencies);


  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 6., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (nrm2)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_2D (nrm2) ) " << std::endl;

  delete[] size_set;
}

void check_global_nrm2_3D(Index rank)
{
   /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;


  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(18);
  DenseVector<Mem::Main, double> b(18);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1), attrs.at(2));

  //Foundation::HaloControl<Foundation::dim_3D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim3D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim3D>::convert(h_new));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_3D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim3D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  double result(0);
  GlobalNorm2<Mem::Main>::value(result, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, std::sqrt(108.), std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (nrm2)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_3D (nrm2) ) " << std::endl;

  delete[] size_set;
}

void check_global_nrm2_1D_gateway(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < b.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  //Foundation::HaloControl<Foundation::dim_1D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<1> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_1D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(b.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalNorm2Gateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, std::sqrt(12.), std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (nrm2 by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_1D (nrm2 by gateway) ) " << std::endl;

  delete[] size_set;
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

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> b(6);
  for(Index i(0) ; i < b.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1));

  //Foundation::HaloControl<Foundation::dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim2D>::convert(h_new));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(b.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalNorm2Gateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 6., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (nrm2 by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_2D (nrm2 by gateway) ) " << std::endl;

  delete[] size_set;
}

void check_global_nrm2_3D_gateway(Index rank)
{
   /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;


  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> b(18);
  for(Index i(0) ; i < b.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1), attrs.at(2));

  //Foundation::HaloControl<Foundation::dim_3D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim3D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim3D>::convert(h_new));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<3> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_3D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim3D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(b.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalNorm2Gateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, std::sqrt(108.), std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (nrm2 by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_3D (nrm2 by gateway) ) " << std::endl;

  delete[] size_set;
}

void check_global_nrm2sqr_1D_gateway(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < b.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  //Foundation::HaloControl<Foundation::dim_1D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<1> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_1D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim1D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(b.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalNorm2SquaredGateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2sqr(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 12., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (nrm2sqr by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_1D (nrm2sqr by gateway) ) " << std::endl;

  delete[] size_set;
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

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> b(6);
  for(Index i(0) ; i < b.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1));

  //Foundation::HaloControl<Foundation::dim_2D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim2D>::convert(h_new));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_2D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim2D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(b.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalNorm2SquaredGateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2sqr(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 36., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (nrm2sqr by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_2D (nrm2sqr by gateway) ) " << std::endl;

  delete[] size_set;
}

void check_global_nrm2sqr_3D_gateway(Index rank)
{
   /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;


  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> b(18);
  for(Index i(0) ; i < b.size() ; ++i)
  {
    b(i, 2.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1), attrs.at(2));

  //Foundation::HaloControl<Foundation::dim_3D>::fill_sizes(*((Halo<0, PLVertex, Mesh<Dim3D> >*)(p0.comm_halos.at(0).get())), polytopes_in_subset);
  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim3D>::convert(h_new));

  //MeshPart<Geometry::ConformalMesh<Shape::Hypercube<3> > > cell_sub_set(polytopes_in_subset);
  //HaloControl<dim_3D>::fill_target_set(*((Halo<0, PLVertex, Mesh<Dim3D> >*)(p0.comm_halos.at(0).get())), cell_sub_set);

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(b.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  GlobalNorm2SquaredGateway<Mem::Main, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2sqr(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 108., std::numeric_limits<double>::epsilon());
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (nrm2sqr by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_3D (nrm2sqr by gateway) ) " << std::endl;

  delete[] size_set;
}

void check_global_synchvec0_1D_gateway(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

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
  res0 = test_check_equal_within_eps(a(0), rank == 0 ? 1. : 2., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 0 ? 2. : 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (synch_vec0 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_1D (synch_vec0 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_1D (synch_vec0 by gateway) ) " << std::endl;

  delete[] size_set;
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

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(6);
  DenseVector<Mem::Main, double> b(6);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1));

  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim2D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

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
  TestResult<double, double, double> res4;
  TestResult<double, double, double> res5;
  res0 = test_check_equal_within_eps(a(0), rank == 0 ? 1. : 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 0 ? 1. : 1., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), rank == 0 ? 2. : 1., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), rank == 0 ? 1. : 2., std::numeric_limits<double>::epsilon());
  res4 = test_check_equal_within_eps(a(4), rank == 0 ? 2. : 2., std::numeric_limits<double>::epsilon());
  res5 = test_check_equal_within_eps(a(5), rank == 0 ? 2. : 2., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed && res2.passed && res3.passed && res4.passed && res5.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (synch_vec0 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (epsll = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway)) ) " << std::endl;
  else if(!res4.passed)
    std::cout << "FAILED: " << res4.left << " not within range (eps = " << res4.epsilon << ") of " << res4.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway) ) " << std::endl;
  else if(!res5.passed)
    std::cout << "FAILED: " << res5.left << " not within range (eps = " << res5.epsilon << ") of " << res5.right << "! (foundation_global_ops_test_2D (synch_vec0 by gateway) ) " << std::endl;

  delete[] size_set;
}

void check_global_synchvec0_3D_gateway(Index rank)
{
   /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;


  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(18);
  DenseVector<Mem::Main, double> b(18);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1), attrs.at(2));

  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim3D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec0Gateway<Mem::Main, DenseVector<Mem::Main, double>, VectorMirror<Mem::Main, double> > gate(
                                                                                                                          mirrors,
                                                                                                                          other_ranks,
                                                                                                                          sendbufs,
                                                                                                                          recvbufs,
                                                                                                                          tags);
  a.synch0(&gate);

  TestResult<double, double, double> res[18];
  const double ref[18][2] =
      {
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 1., 1.},
        { 2., 2.},
        { 1., 2.},
        { 2., 1.},
        { 2., 2.},
        { 2., 2.},
        { 1., 1.},
        { 2., 2.},
        { 2., 2.},
        { 2., 2.},
        { 2., 2.},
        { 2., 2.}
      };

  for( Index i(0); i < 18; i++)
  {
  res[i] = test_check_equal_within_eps(a(i), rank == 0 ? ref[i][0] : ref[i][1], std::numeric_limits<double>::epsilon());
  }
  if( res[0].passed && res[1].passed && res[2].passed && res[3].passed && res[4].passed && res[5].passed && res[6].passed && res[7].passed && res[8].passed && res[9].passed && res[10].passed && res[11].passed && res[12].passed && res[13].passed && res[14].passed && res[15].passed && res[16].passed && res[17].passed )
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (synch_vec0 by gateway)" << std::endl;
  else
  {
    for( Index j(0); j < 18; j++)
      if(!res[j].passed)
        std::cout << "FAILED: " << res[j].left << " not within range (eps = " << res[j].epsilon << ") of " << res[j].right << "! (foundation_global_ops_test_3D (synch_vec0 by gateway) ) " << std::endl;
  }

  delete[] size_set;
}

void check_global_synchvec1_1D_gateway(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

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
  res0 = test_check_equal_within_eps(a(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (synch_vec1 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_1D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_1D (synch_vec1 by gateway) ) " << std::endl;

  delete[] size_set;
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

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim2D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(6);
  DenseVector<Mem::Main, double> b(6);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1));

  auto h_new(Aura<Mem::Main, Halo<0, PLVertex, Mesh<Dim2D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim2D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

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

  TestResult<double, double, double> res0, res1, res2, res3, res4, res5;
  res0 = test_check_equal_within_eps(a(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), 1., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), 1., std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), 1., std::numeric_limits<double>::epsilon());
  res4 = test_check_equal_within_eps(a(4), 1., std::numeric_limits<double>::epsilon());
  res5 = test_check_equal_within_eps(a(5), 1., std::numeric_limits<double>::epsilon());

  if(res0.passed && res1.passed && res2.passed && res3.passed && res4.passed && res5.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D (synch_vec1 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res4.passed)
    std::cout << "FAILED: " << res4.left << " not within range (eps = " << res4.epsilon << ") of " << res4.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;
  else if(!res5.passed)
    std::cout << "FAILED: " << res5.left << " not within range (eps = " << res5.epsilon << ") of " << res5.right << "! (foundation_global_ops_test_2D (synch_vec1 by gateway) ) " << std::endl;

  delete[] size_set;
}

void check_global_synchvec1_3D_gateway(Index rank)
{
   /*
   * (0,1,1)  (1,1,1)
   *      *----*
   *     /    /|
   *(0,1,0)(1,1,0)
   *   *----*  *(1,0,1)
   *   | /  | /
   *   |/   |/
   *   *----*
   *(0,0,0) (1,0,0)
   */

  //create attributes for vertex coords
  std::vector<Attribute<double> > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords
  attrs.push_back(Attribute<double>()); //vertex y-coords
  attrs.push_back(Attribute<double>()); //vertex z-coords

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

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(0));
  attrs.at(2).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(1));
  attrs.at(2).get_data().push_back(double(1));

  /*  2    3
   *  *-1--*
   *  2    |
   *  |    3
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim3D> m;
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);
  m.add_polytope(pl_edge);

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_polytope(pl_polyhedron);

  //vertex 0
  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 10);

  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 3);
  m.add_adjacency(pl_vertex, pl_face, 0, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 0, 0);

  //vertex 1
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 4);

  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);
  m.add_adjacency(pl_vertex, pl_face, 1, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 1, 0);

  //vertex 2
  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_edge, 2, 11);

  m.add_adjacency(pl_vertex, pl_face, 2, 0);
  m.add_adjacency(pl_vertex, pl_face, 2, 3);
  m.add_adjacency(pl_vertex, pl_face, 2, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 2, 0);

  //vertex 3
  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 5);

  m.add_adjacency(pl_vertex, pl_face, 3, 2);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 3, 0);

  //vertex 4
  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_edge, 4, 7);

  m.add_adjacency(pl_vertex, pl_face, 4, 1);
  m.add_adjacency(pl_vertex, pl_face, 4, 2);
  m.add_adjacency(pl_vertex, pl_face, 4, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 4, 0);

  //vertex 5
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_edge, 5, 6);
  m.add_adjacency(pl_vertex, pl_edge, 5, 8);

  m.add_adjacency(pl_vertex, pl_face, 5, 1);
  m.add_adjacency(pl_vertex, pl_face, 5, 2);
  m.add_adjacency(pl_vertex, pl_face, 5, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 5, 0);

  //vertex 6
  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 9);
  m.add_adjacency(pl_vertex, pl_edge, 6, 10);

  m.add_adjacency(pl_vertex, pl_face, 6, 1);
  m.add_adjacency(pl_vertex, pl_face, 6, 3);
  m.add_adjacency(pl_vertex, pl_face, 6, 4);

  m.add_adjacency(pl_vertex, pl_polyhedron, 6, 0);

  //vertex 7
  m.add_adjacency(pl_vertex, pl_edge, 7, 8);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_edge, 7, 11);

  m.add_adjacency(pl_vertex, pl_face, 7, 1);
  m.add_adjacency(pl_vertex, pl_face, 7, 3);
  m.add_adjacency(pl_vertex, pl_face, 7, 5);

  m.add_adjacency(pl_vertex, pl_polyhedron, 7, 0);

  Mesh<Dim3D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim3D>, double> > > > halos;


  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);


  std::vector<Halo<0, PLFace, Mesh<Dim3D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim3D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(18);
  DenseVector<Mem::Main, double> b(18);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<3> > confmeshtype_;

  Index* size_set(new Index[4]);
  MeshControl<dim_3D>::fill_sizes(*((Mesh<Dim3D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_3D>::fill_adjacencies(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_3D>::fill_vertex_sets(*((Mesh<Dim3D>*)(p0.submesh.get())), confmesh, attrs.at(0), attrs.at(1), attrs.at(2));


  auto h_new(Aura<Mem::Main, Halo<0, PLFace, Mesh<Dim3D> > >::value(p0.comm_halos)); //we can do this since we know there is only one other process
  auto cell_sub_set(HaloInterface<0, Dim3D>::convert(h_new));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<3> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

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

  TestResult<double, double, double> res[18];

  for( Index i(0); i < 18; i++)
  {
  res[i] = test_check_equal_within_eps(a(i), 1., std::numeric_limits<double>::epsilon());
  }
  if( res[0].passed && res[1].passed && res[2].passed && res[3].passed && res[4].passed && res[5].passed && res[6].passed && res[7].passed && res[8].passed && res[9].passed && res[10].passed && res[11].passed && res[12].passed && res[13].passed && res[14].passed && res[15].passed && res[16].passed && res[17].passed )
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_3D (synch_vec1 by gateway)" << std::endl;
  else
  {
    for( Index j(0); j < 18; j++)
      if(!res[j].passed)
        std::cout << "FAILED: " << res[j].left << " not within range (eps = " << res[j].epsilon << ") of " << res[j].right << "! (foundation_global_ops_test_3D (synch_vec1 by gateway) ) " << std::endl;
  }

  delete[] size_set;
}

void check_global_product_mat_vec_1D_gateway(Index rank)
{
  /* (0)  (1)
   *  *----*
   */

  std::vector<Attribute<double>, std::allocator<Attribute<double> > > attrs;
  attrs.push_back(Attribute<double>()); //vertex x-coords

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(0).get_data().push_back(double(1));
  /*
   *  *--0-*
   *  0    1
   */

  //creating foundation mesh
  Mesh<Dim1D> m(0);

  m.add_polytope(pl_vertex);
  m.add_polytope(pl_vertex);

  m.add_polytope(pl_edge);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 0);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim1D>, double> > > > halos;

  Mesh<Dim1D> m_fine(m);

  //set up halos

  Refinement<Mem::Main,
             mrt_standard,
             hrt_refine>::execute(m_fine, &halos, attrs);

  /*  *----*----*
   *      (2)
   */

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           2, rank,
                                           attrs
                                           ));

  DenseVector<Mem::Main, double> a(2);
  DenseVector<Mem::Main, double> b(2);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));

  auto cell_sub_set(HaloInterface<0, Dim1D>::convert(p0.comm_halos.at(0).get()));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, cell_sub_set);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
  sendbufs.push_back(std::move(sbuf));
  recvbufs.push_back(std::move(rbuf));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf(mirrors.at(0).size());
  mirror_buffers.push_back(std::move(mbuf));

  DenseVector<Mem::Main, double> fbuf(a.size());

  auto frequencies(HaloFrequencies<Mem::Main>::value(mirrors, mirror_buffers, fbuf));

  SparseMatrixCSR<Mem::Main, double> mat_sys;
  Assembly::SymbolicAssembler::assemble_matrix_std1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  DenseVector<Mem::Main, double> result(a.size(), double(0));
  auto tags(HaloTags::value(p0.comm_halos));
  GlobalProductMat0Vec1Gateway<Mem::Main,
                               DenseVector<Mem::Main, double>,
                               SparseMatrixCSR<Mem::Main, double>,
                               VectorMirror<Mem::Main, double>
                               > gate(
                                      mirrors,
                                      other_ranks,
                                      sendbufs,
                                      recvbufs,
                                      tags);

  mat_sys.apply(result, a, &gate);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  res0 = test_check_equal_within_eps(result(0), 0., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(result(1), 0., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_1D (product_mat0vec1 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_1D (product_mat0vec1 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_1D (product_mat0vec1 by gateway) ) " << std::endl;

  delete[] size_set;
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
  check_global_dot_1D((Index)me);
  check_global_dot_2D((Index)me);
  check_global_dot_3D((Index)me);
  check_global_synchvec0_1D((Index)me);
  check_global_synchvec0_2D((Index)me);
  check_global_synchvec0_3D((Index)me);
  check_global_synchvec1_1D((Index)me);
  check_global_synchvec1_2D((Index)me);
  check_global_synchvec1_3D((Index)me);
  check_global_product_mat_vec_1D((Index)me);
  check_global_defect_1D((Index)me);
  check_global_dot_1D_gateway((Index)me);
  check_global_dot_2D_gateway((Index)me);
  check_global_dot_3D_gateway((Index)me);
  check_global_nrm2_1D((Index)me);
  check_global_nrm2_2D((Index)me);
  check_global_nrm2_3D((Index)me);
  check_global_nrm2_1D_gateway((Index)me);
  check_global_nrm2_2D_gateway((Index)me);
  check_global_nrm2_3D_gateway((Index)me);
  check_global_nrm2sqr_1D_gateway((Index)me);
  check_global_nrm2sqr_2D_gateway((Index)me);
  check_global_nrm2sqr_3D_gateway((Index)me);
  check_global_synchvec0_1D_gateway((Index)me);
  check_global_synchvec0_2D_gateway((Index)me);
  check_global_synchvec0_3D_gateway((Index)me);
  check_global_synchvec1_1D_gateway((Index)me);
  check_global_synchvec1_2D_gateway((Index)me);
  check_global_synchvec1_3D_gateway((Index)me);
  check_global_product_mat_vec_1D_gateway((Index)me);
#else
  std::cout << "Parallel tests unavailable on sole process " << me << std::endl;
#endif

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
