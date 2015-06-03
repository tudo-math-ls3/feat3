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
#include <kernel/scarc/gateway_creation.hpp>
#include <kernel/scarc/scarc_data.hpp>
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
using namespace ScaRC;
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

void check_dot_gateway_1D(Index rank)
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
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  SynchronisedPreconditionedScaRCData<double,
                                       Mem::Main,
                                       DenseVector<Mem::Main, double>,
                                       VectorMirror<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double> > pdata(std::move(mat_sys), std::move(mat_sys), std::move(a), std::move(b));

  pdata.halo_frequencies() = std::move(frequencies);

  auto dot_gate(GatewayCreation<got_dot>::value(pdata));

  double res(pdata.sol().dot(pdata.rhs(), &dot_gate));

  TestResult<double, double, double> res0;
  res0 = test_check_equal_within_eps(res, 6., std::numeric_limits<double>::epsilon());
  if(res0.passed)
    std::cout << "PASSED (rank " << rank <<"): scarc_gateway_test_1D (dot)" << std::endl;
  else
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (scarc_gateway_test_1D (dot) ) " << std::endl;
}

void check_nrm2_gateway_1D(Index rank)
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
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  SynchronisedPreconditionedScaRCData<double,
                                       Mem::Main,
                                       DenseVector<Mem::Main, double>,
                                       VectorMirror<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double> > pdata(std::move(mat_sys), std::move(mat_sys), std::move(a), std::move(b));

  pdata.halo_frequencies() = std::move(frequencies);

  auto nrm2_gate(GatewayCreation<got_nrm2>::value(pdata));

  double res(pdata.rhs().norm2(&nrm2_gate));

  TestResult<double, double, double> res0;
  res0 = test_check_equal_within_eps(res, std::sqrt(12.), std::numeric_limits<double>::epsilon());
  if(res0.passed)
    std::cout << "PASSED (rank " << rank <<"): scarc_gateway_test_1D (nrm2)" << std::endl;
  else
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (scarc_gateway_test_1D (nrm2) ) " << std::endl;
}

void check_nrm2sqr_gateway_1D(Index rank)
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
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  SynchronisedPreconditionedScaRCData<double,
                                       Mem::Main,
                                       DenseVector<Mem::Main, double>,
                                       VectorMirror<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double> > pdata(std::move(mat_sys), std::move(mat_sys), std::move(a), std::move(b));

  pdata.halo_frequencies() = std::move(frequencies);
  auto tags(HaloTags::value(p0.comm_halos));
  pdata.tags() = std::move(tags);

  auto nrm2sqr_gate(GatewayCreation<got_nrm2sqr>::value(pdata));

  double res(pdata.rhs().norm2sqr(&nrm2sqr_gate));

  TestResult<double, double, double> res0;
  res0 = test_check_equal_within_eps(res, 12., std::numeric_limits<double>::epsilon());
  if(res0.passed)
    std::cout << "PASSED (rank " << rank <<"): scarc_gateway_test_1D (nrm2sqr)" << std::endl;
  else
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (scarc_gateway_test_1D (nrm2sqr) ) " << std::endl;
}

void check_synch0_gateway_1D(Index rank)
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
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  SynchronisedPreconditionedScaRCData<double,
                                       Mem::Main,
                                       DenseVector<Mem::Main, double>,
                                       VectorMirror<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double> > pdata(std::move(mat_sys), std::move(mat_sys), std::move(a), std::move(b));

  pdata.vector_mirrors() = std::move(mirrors);
  pdata.vector_mirror_sendbufs() = std::move(sendbufs);
  pdata.vector_mirror_recvbufs() = std::move(recvbufs);
  pdata.dest_ranks().push_back(rank == 0 ? 1 : 0);

  pdata.halo_frequencies() = std::move(frequencies);
  auto tags(HaloTags::value(p0.comm_halos));
  pdata.tags() = std::move(tags);

  auto synch0_gate(GatewayCreation<got_synch_vec0>::value(pdata));

  pdata.sol().synch0(&synch0_gate);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  res0 = test_check_equal_within_eps(pdata.sol()(0), rank == 0 ? 1. : 2., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(pdata.sol()(1), rank == 0 ? 2. : 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): scarc_gateway_test_1D (synch0)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (scarc_gateway_test_1D (synch0) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (scarc_gateway_test_1D (synch0) ) " << std::endl;
}

void check_synch1_gateway_1D(Index rank)
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
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  SynchronisedPreconditionedScaRCData<double,
                                       Mem::Main,
                                       DenseVector<Mem::Main, double>,
                                       VectorMirror<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double> > pdata(std::move(mat_sys), std::move(mat_sys), std::move(a), std::move(b));

  pdata.vector_mirrors() = std::move(mirrors);
  pdata.vector_mirror_sendbufs() = std::move(sendbufs);
  pdata.vector_mirror_recvbufs() = std::move(recvbufs);
  pdata.dest_ranks().push_back(rank == 0 ? 1 : 0);

  pdata.halo_frequencies() = std::move(frequencies);
  auto tags(HaloTags::value(p0.comm_halos));
  pdata.tags() = std::move(tags);

  auto synch1_gate(GatewayCreation<got_synch_vec1>::value(pdata));

  pdata.sol().synch1(&synch1_gate);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  res0 = test_check_equal_within_eps(pdata.sol()(0), 1., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(pdata.sol()(1), 1., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): scarc_gateway_test_1D (synch1)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (scarc_gateway_test_1D (synch1) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (scarc_gateway_test_1D (synch1) ) " << std::endl;
}

void check_prodmat0vec1_gateway_1D(Index rank)
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
  a(rank == 0 ? 1 : 0, 5.);

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
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  SynchronisedPreconditionedScaRCData<double,
                                       Mem::Main,
                                       DenseVector<Mem::Main, double>,
                                       VectorMirror<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double>,
                                       SparseMatrixCSR<Mem::Main, double> > pdata(std::move(mat_sys), std::move(mat_sys), std::move(a), std::move(b));

  pdata.vector_mirrors() = std::move(mirrors);
  pdata.vector_mirror_sendbufs() = std::move(sendbufs);
  pdata.vector_mirror_recvbufs() = std::move(recvbufs);
  pdata.dest_ranks().push_back(rank == 0 ? 1 : 0);

  pdata.halo_frequencies() = std::move(frequencies);
  auto tags(HaloTags::value(p0.comm_halos));
  pdata.tags() = std::move(tags);

  auto prod_gate(GatewayCreation<got_product_mat0_vec1>::value(pdata));

  DenseVector<Mem::Main, double> res(pdata.sol().size());
  pdata.sys().apply(res, pdata.sol(), &prod_gate);

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;

  res0 = test_check_equal_within_eps(res(0), rank == 0 ? -4. : 8., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(res(1), rank == 0 ? 8. : -4., std::numeric_limits<double>::epsilon());
  if(res0.passed && res1.passed)
    std::cout << "PASSED (rank " << rank <<"): scarc_gateway_test_1D (prodmat0vec1)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (scarc_gateway_test_1D (prodmat0vec1) ) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (scarc_gateway_test_1D (prodmat0vec1) ) " << std::endl;
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
  check_dot_gateway_1D((Index)me);
  check_nrm2_gateway_1D((Index)me);
  check_nrm2sqr_gateway_1D((Index)me);
  check_synch0_gateway_1D((Index)me);
  check_synch1_gateway_1D((Index)me);
  check_prodmat0vec1_gateway_1D((Index)me);
#else
  std::cout << "Parallel tests unavailable on sole process " << me << std::endl;
#endif

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
