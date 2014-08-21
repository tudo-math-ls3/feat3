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

void check_global_dot_2D_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
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

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  double result(0);
  GlobalDot<Mem::Main, Algo::Generic>::value(result, a, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 16., std::numeric_limits<double>::epsilon()*10);
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D_concave (dot)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

void check_global_synchvec0_2D_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  DenseVector<Mem::Main, double> sbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> sbuf1(mirrors.at(1).size());
  sendbufs.push_back(std::move(sbuf0));
  sendbufs.push_back(std::move(sbuf1));
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> rbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf1(mirrors.at(1).size());
  recvbufs.push_back(std::move(rbuf0));
  recvbufs.push_back(std::move(rbuf1));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());
  other_ranks.push_back(p0.comm_halos.at(1)->get_other());

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec0<Mem::Main, Algo::Generic>::exec(a,
                                                  mirrors,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs,
                                                  tags);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  TestResult<double, double, double> res3;

  res0 = test_check_equal_within_eps(a(0), rank == 0 ? 2. : ( rank == 1 ? 3. : 2.), std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 0 ? 3. : ( rank == 1 ? 2. : 3.), std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), rank == 0 ? 1. : ( rank == 1 ? 1. : 1.), std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), rank == 0 ? 2. : ( rank == 1 ? 1. : 1.), std::numeric_limits<double>::epsilon());

  if(res0.passed && res1.passed && res2.passed && res3.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D_concave (synch_vec0)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (epsll = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D_concave (synch_vec0) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D_concave (synch_vec0) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D_concave (synch_vec0) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D_concave (synch_vec0) ) " << std::endl;
}

void check_global_synchvec1_2D_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> sbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> rbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf1(mirrors.at(1).size());
  sendbufs.push_back(std::move(sbuf0));
  sendbufs.push_back(std::move(sbuf1));
  recvbufs.push_back(std::move(rbuf0));
  recvbufs.push_back(std::move(rbuf1));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());
  other_ranks.push_back(p0.comm_halos.at(1)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec1<Mem::Main, Algo::Generic>::exec(a,
                                                  mirrors,
                                                  frequencies,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs,
                                                  tags);

  TestResult<double, double, double> res0, res1, res2, res3;
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

void check_global_synchvec0_2D_gateway_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  DenseVector<Mem::Main, double> sbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> sbuf1(mirrors.at(1).size());
  sendbufs.push_back(std::move(sbuf0));
  sendbufs.push_back(std::move(sbuf1));
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> rbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf1(mirrors.at(1).size());
  recvbufs.push_back(std::move(rbuf0));
  recvbufs.push_back(std::move(rbuf1));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());
  other_ranks.push_back(p0.comm_halos.at(1)->get_other());

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec0Gateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double>, VectorMirror<Mem::Main, double> > gate(
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

  res0 = test_check_equal_within_eps(a(0), rank == 0 ? 2. : ( rank == 1 ? 3. : 2.), std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(a(1), rank == 0 ? 3. : ( rank == 1 ? 2. : 3.), std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(a(2), rank == 0 ? 1. : ( rank == 1 ? 1. : 1.), std::numeric_limits<double>::epsilon());
  res3 = test_check_equal_within_eps(a(3), rank == 0 ? 2. : ( rank == 1 ? 1. : 1.), std::numeric_limits<double>::epsilon());

  if(res0.passed && res1.passed && res2.passed && res3.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D_concave (synch_vec0 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (epsll = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D_concave (synch_vec0by gateway) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D_concave (synch_vec0 by gateway) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D_concave (synch_vec0 by gateway) ) " << std::endl;
}

void check_global_synchvec1_2D_gateway_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
                                           attrs
                                          ));

  DenseVector<Mem::Main, double> a(4);
  for(Index i(0) ; i < a.size() ; ++i)
  {
    a(i, 1.);
  }

  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;

  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)p0.submesh.get()), size_set);

  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)p0.submesh.get()), confmesh, attrs.at(0), attrs.at(1));

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  DenseVector<Mem::Main, double> sbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> sbuf1(mirrors.at(1).size());
  DenseVector<Mem::Main, double> rbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> rbuf1(mirrors.at(1).size());
  sendbufs.push_back(std::move(sbuf0));
  sendbufs.push_back(std::move(sbuf1));
  recvbufs.push_back(std::move(rbuf0));
  recvbufs.push_back(std::move(rbuf1));

  std::vector<Index> other_ranks;
  other_ranks.push_back(p0.comm_halos.at(0)->get_other());
  other_ranks.push_back(p0.comm_halos.at(1)->get_other());

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  auto tags(HaloTags::value(p0.comm_halos));
  GlobalSynchVec1Gateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double>, VectorMirror<Mem::Main, double> > gate(
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
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D_concave (synch_vec1 by gateway)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (epsll = " << res0.epsilon << ") of " << res0.right << "! (foundation_global_ops_test_2D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (foundation_global_ops_test_2D_concave (synch_vec1by gateway) ) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (foundation_global_ops_test_2D_concave (synch_vec1 by gateway) ) " << std::endl;
  else if(!res3.passed)
    std::cout << "FAILED: " << res3.left << " not within range (eps = " << res3.epsilon << ") of " << res3.right << "! (foundation_global_ops_test_2D_concave (synch_vec1 by gateway) ) " << std::endl;
}

void check_global_dot_2D_gateway_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
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

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  GlobalDotGateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(a.dot(b, &gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 16., std::numeric_limits<double>::epsilon()*10.);
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D_concave (dot by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_2D_concave (dot by gateway) ) " << std::endl;
}

void check_global_nrm2_2D_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
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

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  double result(0);
  GlobalNorm2<Mem::Main, Algo::Generic>::value(result, b, frequencies);

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, std::sqrt(32.), std::numeric_limits<double>::epsilon()*10);
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D_concave (norm2)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

void check_global_nrm2_2D_gateway_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
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

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  GlobalNorm2Gateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, std::sqrt(32.), std::numeric_limits<double>::epsilon()*10.);
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D_concave (norm2 by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "!" << std::endl;
}

void check_global_nrm2sqr_2D_gateway_concave(Index rank)
{
  /*(0,1) (1,1) (2,1)
   *  *----*----*
   *  |    |    |
   *  |    |    |
   *  *----*----*
   *(0,0) (1,0) (2,0)
   *  |    |
   *  *----*
   *(0,-1) (1,-1)
   */

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

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(1));

  attrs.at(0).get_data().push_back(double(2));
  attrs.at(1).get_data().push_back(double(0));

  attrs.at(0).get_data().push_back(double(0));
  attrs.at(1).get_data().push_back(double(-1));

  attrs.at(0).get_data().push_back(double(1));
  attrs.at(1).get_data().push_back(double(-1));

  /*  2    3    4
   *  *-1--*-6--*
   *  2  0 | 1  |
   *  |    3    4
   *  *--0-*-5--*
   *  0 2  1    5
   *  8    9
   *  *--7-*
   *  6    7
   */

  //creating foundation mesh
  Mesh<Dim2D> m(0);
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

  m.add_polytope(pl_face);
  m.add_polytope(pl_face);
  m.add_polytope(pl_face);

  m.add_adjacency(pl_vertex, pl_edge, 0, 0);
  m.add_adjacency(pl_vertex, pl_edge, 0, 2);
  m.add_adjacency(pl_vertex, pl_edge, 0, 8);
  m.add_adjacency(pl_vertex, pl_face, 0, 0);
  m.add_adjacency(pl_vertex, pl_face, 0, 2);

  m.add_adjacency(pl_vertex, pl_edge, 1, 0);
  m.add_adjacency(pl_vertex, pl_edge, 1, 3);
  m.add_adjacency(pl_vertex, pl_edge, 1, 9);
  m.add_adjacency(pl_vertex, pl_face, 1, 0);
  m.add_adjacency(pl_vertex, pl_face, 1, 1);
  m.add_adjacency(pl_vertex, pl_face, 1, 2);

  m.add_adjacency(pl_vertex, pl_edge, 2, 1);
  m.add_adjacency(pl_vertex, pl_edge, 2, 2);
  m.add_adjacency(pl_vertex, pl_face, 2, 0);

  m.add_adjacency(pl_vertex, pl_edge, 3, 1);
  m.add_adjacency(pl_vertex, pl_edge, 3, 3);
  m.add_adjacency(pl_vertex, pl_edge, 3, 6);
  m.add_adjacency(pl_vertex, pl_face, 3, 0);
  m.add_adjacency(pl_vertex, pl_face, 3, 1);

  m.add_adjacency(pl_vertex, pl_edge, 4, 4);
  m.add_adjacency(pl_vertex, pl_edge, 4, 6);
  m.add_adjacency(pl_vertex, pl_face, 4, 1);

  m.add_adjacency(pl_vertex, pl_edge, 5, 4);
  m.add_adjacency(pl_vertex, pl_edge, 5, 5);
  m.add_adjacency(pl_vertex, pl_face, 5, 1);

  m.add_adjacency(pl_vertex, pl_edge, 6, 7);
  m.add_adjacency(pl_vertex, pl_edge, 6, 8);
  m.add_adjacency(pl_vertex, pl_face, 6, 2);

  m.add_adjacency(pl_vertex, pl_edge, 7, 7);
  m.add_adjacency(pl_vertex, pl_edge, 7, 9);
  m.add_adjacency(pl_vertex, pl_face, 7, 2);

  Mesh<Dim2D> m_fine(m);

  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries; //not needed

  auto p0(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(
                                           m_fine,
                                           boundaries,
                                           3, rank,
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

  auto cell_sub_set0(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(0).get()));
  auto cell_sub_set1(HaloInterface<0, Dim2D>::convert(p0.comm_halos.at(1).get()));
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  VectorMirror<Mem::Main, double> target_mirror0;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror0, space, cell_sub_set0);
  VectorMirror<Mem::Main, double> target_mirror1;
  Assembly::MirrorAssembler::assemble_mirror(target_mirror1, space, cell_sub_set1);

  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  mirrors.push_back(std::move(target_mirror0));
  mirrors.push_back(std::move(target_mirror1));

  std::vector<DenseVector<Mem::Main, double> > mirror_buffers;
  DenseVector<Mem::Main, double> mbuf0(mirrors.at(0).size());
  DenseVector<Mem::Main, double> mbuf1(mirrors.at(1).size());
  mirror_buffers.push_back(std::move(mbuf0));
  mirror_buffers.push_back(std::move(mbuf1));

  std::vector<DenseVector<Mem::Main, double> > freq_buffers;
  DenseVector<Mem::Main, double> fbuf0(a.size());
  DenseVector<Mem::Main, double> fbuf1(a.size());
  freq_buffers.push_back(std::move(fbuf0));
  freq_buffers.push_back(std::move(fbuf1));

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mirror_buffers, freq_buffers));

  GlobalNorm2SquaredGateway<Mem::Main, Algo::Generic, DenseVector<Mem::Main, double> > gate(frequencies);
  double result(b.norm2sqr(&gate));

  TestResult<double, double, double> res;
  res = test_check_equal_within_eps(result, 32., std::numeric_limits<double>::epsilon()*20.);
  if(res.passed)
    std::cout << "PASSED (rank " << rank <<"): foundation_global_ops_test_2D_concave (nrm2sqr by gateway)" << std::endl;
  else
    std::cout << "FAILED: " << res.left << " not within range (eps = " << res.epsilon << ") of " << res.right << "! (foundation_global_ops_test_2D_concave (nrm2sqr by gateway) ) " << std::endl;
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
  check_global_dot_2D_concave((Index)me);
  check_global_synchvec0_2D_concave((Index)me);
  check_global_synchvec1_2D_concave((Index)me);
  check_global_synchvec0_2D_gateway_concave((Index)me);
  check_global_synchvec1_2D_gateway_concave((Index)me);
  check_global_dot_2D_gateway_concave((Index)me);
  check_global_nrm2_2D_concave((Index)me);
  check_global_nrm2_2D_gateway_concave((Index)me);
  check_global_nrm2sqr_2D_gateway_concave((Index)me);
#else
  std::cout << "Parallel tests unavailable on sole process " << me << std::endl;
#endif

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
