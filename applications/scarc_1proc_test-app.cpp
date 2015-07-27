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
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/scarc/scarc_functor.hpp>
#include <kernel/scarc/matrix_conversion.hpp>

#include <iostream>
#include <limits>

using namespace FEAST;
using namespace Foundation;
using namespace Geometry;
using namespace ScaRC;
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

void check_scarc_rich_1D(Index rank)
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

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries;
  boundaries.push_back(std::move(Halo<0, PLVertex, Mesh<Dim1D> >(m_fine)));
  boundaries.push_back(std::move(Halo<0, PLVertex, Mesh<Dim1D> >(m_fine)));
  boundaries.at(0).push_back(0);
  boundaries.at(1).push_back(1);

  /*auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           1, rank,
                                           attrs
                                           ));*/

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  //MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);
  MeshControl<dim_1D>::fill_sizes(m_fine, size_set);

  confmeshtype_ confmesh(size_set);
  //MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  //MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));
  MeshControl<dim_1D>::fill_adjacencies(m_fine, confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(m_fine, confmesh, attrs.at(0));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);

  SparseMatrixCSR<Mem::Main, double> mat_sys;
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  DenseVector<Mem::Main, double> vec_rhs(space.get_num_dofs(), double(0));
  Assembly::Common::ConstantFunction rhs_func(1.0);
  Assembly::Common::ForceFunctional<Assembly::Common::ConstantFunction> rhs_functional(rhs_func);
  Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, rhs_functional, space, cubature_factory);

  Assembly::UnitFilterAssembler<Geometry::ConformalMesh<Shape::Hypercube<1> > > dirichlet;
  auto bound_sub_set_0(HaloInterface<0, Dim1D>::convert(boundaries.at(0)));
  auto bound_sub_set_1(HaloInterface<0, Dim1D>::convert(boundaries.at(1)));
  dirichlet.add_mesh_part(bound_sub_set_0);
  dirichlet.add_mesh_part(bound_sub_set_1);

  DenseVector<Mem::Main, double> vec_sol(space.get_num_dofs(), double(0));

  UnitFilter<Mem::Main, double> filter(space.get_num_dofs());
  dirichlet.assemble(filter, space);

  SparseMatrixCSR<Mem::Main, double> mat_localsys;
  mat_localsys.clone(mat_sys);

  SparseMatrixCOO<Mem::Main, double> mat_precon_temp(mat_localsys.rows(), mat_localsys.columns());
  for(Index i(0) ; i < mat_localsys.rows() ; ++i)
    mat_precon_temp(i, i, double(0.75) * (double(1)/mat_localsys(i, i)));

  SparseMatrixCSR<Mem::Main, double> mat_precon(mat_precon_temp);

  ///filter system
  filter.filter_mat(mat_sys);
  filter.filter_mat(mat_localsys);
  filter.filter_rhs(vec_rhs);
  filter.filter_sol(vec_sol);
  filter.filter_mat(mat_precon); //TODO: check if -> NO! we do this in the solver program when applying the correction filter after preconditioning

  SynchronisedPreconditionedFilteredScaRCData<double,
                                              Mem::Main,
                                              DenseVector<Mem::Main, double>,
                                              VectorMirror<Mem::Main, double>,
                                              SparseMatrixCSR<Mem::Main, double>,
                                              SparseMatrixCSR<Mem::Main, double>,
                                              UnitFilter<Mem::Main, double> > data(std::move(mat_sys), std::move(mat_precon), std::move(vec_sol), std::move(vec_rhs), std::move(filter));

#ifndef SERIAL
  Communicator c(MPI_COMM_WORLD);
#else
  Communicator c(0);
#endif
  data.communicators().push_back(std::move(c));

  //data.source_ranks() = std::move(sourceranks);
  data.localsys() = std::move(mat_localsys);

  ///layer 0 (local layer)
  std::shared_ptr<ScaRCFunctorBase<double,
                                   Mem::Main,
                                   DenseVector<Mem::Main, double>,
                                   VectorMirror<Mem::Main, double>,
                                   SparseMatrixCSR<Mem::Main, double>,
                                   SparseMatrixCSR<Mem::Main, double>,
                                   UnitFilter<Mem::Main, double>,
                                   std::vector,
                                   Index> > local_solver(new ScaRCFunctorRichardson1<double,
                                                                                       Mem::Main,
                                                                                       DenseVector<Mem::Main, double>,
                                                                                       VectorMirror<Mem::Main, double>,
                                                                                       SparseMatrixCSR<Mem::Main, double>,
                                                                                       SparseMatrixCSR<Mem::Main, double>,
                                                                                       UnitFilter<Mem::Main, double>,
                                                                                       std::vector,
                                                                                       Index>(data) );

  ///layer 0 (local layer), preconditioner
  std::shared_ptr<ScaRCFunctorBase<double,
                                   Mem::Main,
                                   DenseVector<Mem::Main, double>,
                                   VectorMirror<Mem::Main, double>,
                                   SparseMatrixCSR<Mem::Main, double>,
                                   SparseMatrixCSR<Mem::Main, double>,
                                   UnitFilter<Mem::Main, double>,
                                   std::vector,
                                   Index> > local_precon(new ScaRCFunctorPreconSpM1V1<double,
                                                                                              Mem::Main,
                                                                                              DenseVector<Mem::Main, double>,
                                                                                              VectorMirror<Mem::Main, double>,
                                                                                              SparseMatrixCSR<Mem::Main, double>,
                                                                                              SparseMatrixCSR<Mem::Main, double>,
                                                                                              UnitFilter<Mem::Main, double>,
                                                                                              std::vector,
                                                                                              Index>(data) );

  if(rank == 0)
  {
    std::cout << "A0 " << data.sys() << std::endl;
    std::cout << "A1 " << data.localsys() << std::endl;
    std::cout << "P " << data.precon() << std::endl;
    std::cout << "b " << data.rhs() << std::endl;
    std::cout << "x_0 " << data.sol() << std::endl;
  }

  local_solver->reset_preconditioner(local_precon);
  local_solver->execute();
  std::cout << rank << ", #iters local Rich: " << local_solver->iterations() << std::endl;

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  res0 = test_check_equal_within_eps(data.sol()(0), 0., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(data.sol()(1), 0., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(data.sol()(2), 0.125, std::numeric_limits<double>::epsilon()*1e8);

  std::cout << "SOL " << data.sol() << " SOL" << std::endl;
  std::cout << "s0 " << data.sol()(0) << std::endl;
  std::cout << "s1 " << data.sol()(1) << std::endl;
  std::cout << "s2 " << data.sol()(2) << std::endl;

  std::cout << "p0 " << res0.passed << std::endl;
  std::cout << "p1 " << res1.passed << std::endl;
  std::cout << "p2 " << res2.passed << std::endl;

  if(res0.passed && res1.passed && res2.passed)
    std::cout << "PASSED (rank " << rank <<"): scarc_test_1D (rich)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " (0) not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (scarc_test_1D (rich)) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " (1) not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (scarc_test_1D (rich)) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " (2) not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (scarc_test_1D (rich)) " << std::endl;
}

void check_scarc_pcg_1D(Index rank)
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

  std::vector<Halo<0, PLVertex, Mesh<Dim1D> > > boundaries;
  boundaries.push_back(std::move(Halo<0, PLVertex, Mesh<Dim1D> >(m_fine)));
  boundaries.push_back(std::move(Halo<0, PLVertex, Mesh<Dim1D> >(m_fine)));
  boundaries.at(0).push_back(0);
  boundaries.at(1).push_back(1);

  /*auto p0(Partitioning<Mem::Main,
                       Dim1D,
                       0,
                       pl_vertex>::execute(m_fine,
                                           boundaries,
                                           1, rank,
                                           attrs
                                           ));*/

  typedef ConformalMesh<Shape::Hypercube<1> > confmeshtype_;

  Index* size_set(new Index[2]);
  //MeshControl<dim_1D>::fill_sizes(*((Mesh<Dim1D>*)(p0.submesh.get())), size_set);
  MeshControl<dim_1D>::fill_sizes(m_fine, size_set);

  confmeshtype_ confmesh(size_set);
  //MeshControl<dim_1D>::fill_adjacencies(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh);
  //MeshControl<dim_1D>::fill_vertex_sets(*((Mesh<Dim1D>*)(p0.submesh.get())), confmesh, attrs.at(0));
  MeshControl<dim_1D>::fill_adjacencies(m_fine, confmesh);
  MeshControl<dim_1D>::fill_vertex_sets(m_fine, confmesh, attrs.at(0));

  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<1> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<1> > > > space(trafo);

  SparseMatrixCSR<Mem::Main, double> mat_sys;
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  DenseVector<Mem::Main, double> vec_rhs(space.get_num_dofs(), double(0));
  Assembly::Common::ConstantFunction rhs_func(1.0);
  Assembly::Common::ForceFunctional<Assembly::Common::ConstantFunction> rhs_functional(rhs_func);
  Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, rhs_functional, space, cubature_factory);

  Assembly::UnitFilterAssembler<Geometry::ConformalMesh<Shape::Hypercube<1> > > dirichlet;
  auto bound_sub_set_0(HaloInterface<0, Dim1D>::convert(boundaries.at(0)));
  auto bound_sub_set_1(HaloInterface<0, Dim1D>::convert(boundaries.at(1)));
  dirichlet.add_mesh_part(bound_sub_set_0);
  dirichlet.add_mesh_part(bound_sub_set_1);

  DenseVector<Mem::Main, double> vec_sol(space.get_num_dofs(), double(0));

  UnitFilter<Mem::Main, double> filter(space.get_num_dofs());
  dirichlet.assemble(filter, space);

  SparseMatrixCSR<Mem::Main, double> mat_localsys;
  mat_localsys.clone(mat_sys);

  SparseMatrixCOO<Mem::Main, double> mat_precon_temp(mat_localsys.rows(), mat_localsys.columns());
  for(Index i(0) ; i < mat_localsys.rows() ; ++i)
    mat_precon_temp(i, i, double(0.75) * (double(1)/mat_localsys(i, i)));

  SparseMatrixCSR<Mem::Main, double> mat_precon(mat_precon_temp);

  ///filter system
  filter.filter_mat(mat_sys);
  filter.filter_mat(mat_localsys);
  filter.filter_rhs(vec_rhs);
  filter.filter_sol(vec_sol);
  filter.filter_mat(mat_precon); //TODO: check if -> NO! we do this in the solver program when applying the correction filter after preconditioning

  SynchronisedPreconditionedFilteredScaRCData<double,
                                              Mem::Main,
                                              DenseVector<Mem::Main, double>,
                                              VectorMirror<Mem::Main, double>,
                                              SparseMatrixCSR<Mem::Main, double>,
                                              SparseMatrixCSR<Mem::Main, double>,
                                              UnitFilter<Mem::Main, double> > data(std::move(mat_sys), std::move(mat_precon), std::move(vec_sol), std::move(vec_rhs), std::move(filter));

#ifndef SERIAL
  Communicator c(MPI_COMM_WORLD);
#else
  Communicator c(0);
#endif
  data.communicators().push_back(std::move(c));

  //data.source_ranks() = std::move(sourceranks);
  data.localsys() = std::move(mat_localsys);

  ///layer 0 (local layer)
  std::shared_ptr<ScaRCFunctorBase<double,
                                   Mem::Main,
                                   DenseVector<Mem::Main, double>,
                                   VectorMirror<Mem::Main, double>,
                                   SparseMatrixCSR<Mem::Main, double>,
                                   SparseMatrixCSR<Mem::Main, double>,
                                   UnitFilter<Mem::Main, double>,
                                   std::vector,
                                   Index> > local_solver(new ScaRCFunctorPCG1<double,
                                                                                       Mem::Main,
                                                                                       DenseVector<Mem::Main, double>,
                                                                                       VectorMirror<Mem::Main, double>,
                                                                                       SparseMatrixCSR<Mem::Main, double>,
                                                                                       SparseMatrixCSR<Mem::Main, double>,
                                                                                       UnitFilter<Mem::Main, double>,
                                                                                       std::vector,
                                                                                       Index>(data) );

  ///layer 0 (local layer), preconditioner
  std::shared_ptr<ScaRCFunctorBase<double,
                                   Mem::Main,
                                   DenseVector<Mem::Main, double>,
                                   VectorMirror<Mem::Main, double>,
                                   SparseMatrixCSR<Mem::Main, double>,
                                   SparseMatrixCSR<Mem::Main, double>,
                                   UnitFilter<Mem::Main, double>,
                                   std::vector,
                                   Index> > local_precon(new ScaRCFunctorPreconSpM1V1<double,
                                                                                              Mem::Main,
                                                                                              DenseVector<Mem::Main, double>,
                                                                                              VectorMirror<Mem::Main, double>,
                                                                                              SparseMatrixCSR<Mem::Main, double>,
                                                                                              SparseMatrixCSR<Mem::Main, double>,
                                                                                              UnitFilter<Mem::Main, double>,
                                                                                              std::vector,
                                                                                              Index>(data) );

  if(rank == 0)
  {
    std::cout << "A0 " << data.sys() << std::endl;
    std::cout << "A1 " << data.localsys() << std::endl;
    std::cout << "P " << data.precon() << std::endl;
    std::cout << "b " << data.rhs() << std::endl;
    std::cout << "x_0 " << data.sol() << std::endl;
  }

  local_solver->reset_preconditioner(local_precon);
  local_solver->execute();
  std::cout << rank << ", #iters local PCG: " << local_solver->iterations() << std::endl;

  TestResult<double, double, double> res0;
  TestResult<double, double, double> res1;
  TestResult<double, double, double> res2;
  res0 = test_check_equal_within_eps(data.sol()(0), 0., std::numeric_limits<double>::epsilon());
  res1 = test_check_equal_within_eps(data.sol()(1), 0., std::numeric_limits<double>::epsilon());
  res2 = test_check_equal_within_eps(data.sol()(2), 0.125, std::numeric_limits<double>::epsilon()*1e8);

  std::cout << "SOL " << data.sol() << " SOL" << std::endl;
  std::cout << "s0 " << data.sol()(0) << std::endl;
  std::cout << "s1 " << data.sol()(1) << std::endl;
  std::cout << "s2 " << data.sol()(2) << std::endl;

  std::cout << "p0 " << res0.passed << std::endl;
  std::cout << "p1 " << res1.passed << std::endl;
  std::cout << "p2 " << res2.passed << std::endl;

  if(res0.passed && res1.passed && res2.passed)
    std::cout << "PASSED (rank " << rank <<"): scarc_test_1D (rich)" << std::endl;
  else if(!res0.passed)
    std::cout << "FAILED: " << res0.left << " (0) not within range (eps = " << res0.epsilon << ") of " << res0.right << "! (scarc_test_1D (rich)) " << std::endl;
  else if(!res1.passed)
    std::cout << "FAILED: " << res1.left << " (1) not within range (eps = " << res1.epsilon << ") of " << res1.right << "! (scarc_test_1D (rich)) " << std::endl;
  else if(!res2.passed)
    std::cout << "FAILED: " << res2.left << " (2) not within range (eps = " << res2.epsilon << ") of " << res2.right << "! (scarc_test_1D (rich)) " << std::endl;
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
  check_scarc_rich_1D((Index)me);
  check_scarc_pcg_1D((Index)me);
#else
  std::cout << "Parallel tests unavailable on sole process " << me << std::endl;
#endif

#ifndef SERIAL
  MPI_Finalize();
#endif

  return 0;
}
