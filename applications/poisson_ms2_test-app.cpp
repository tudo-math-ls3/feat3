#ifndef SERIAL
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
#include <kernel/scarc/scarc_functor.hpp>
#include <kernel/scarc/matrix_conversion.hpp>
#include <kernel/scarc/scarc_log.hpp>

#include <iostream>
#include <limits>

using namespace FEAST;
using namespace Foundation;
using namespace Geometry;
using namespace ScaRC;

void testmesh_hypercube_2D(Mesh<Dim2D>& target_mesh, std::vector<Attribute<double> >& attrs, std::vector<Halo<0, PLEdge, Mesh<Dim2D> > >& boundaries);

#ifndef SERIAL
int main(int argc, char* argv[])
#else
int main()
#endif
{
#ifndef SERIAL
  int rank(0), nprocs(1);
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
  (void)argc;
  (void)argv;

  ScaRCLog<> log;

  ///setup geometric problem data: mesh, its attributes and boundaries
  Mesh<Dim2D> mesh;
  std::vector<Attribute<double> > attrs;
  std::vector<Halo<0, PLEdge, Mesh<Dim2D> > > boundaries;
  testmesh_hypercube_2D(mesh, attrs, boundaries);

  ///provide memory for halos
  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > halos;

  ///minimal level 1
  std::vector<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >, std::allocator<std::shared_ptr<HaloBase<Mesh<Dim2D>, double> > > > boundaries_copy;
  for(Index i(0) ; i < boundaries.size() ; ++i)
    boundaries_copy.push_back(std::shared_ptr<HaloBase<Mesh<Dim2D>, double> >(new Halo<0, PLEdge, Mesh<Dim2D>, double>(boundaries.at(i))));

  Refinement<Mem::Main,
    Algo::Generic,
    mrt_standard>::execute(mesh, &boundaries_copy, attrs);
  boundaries.clear();

  for(auto& bc_i : boundaries_copy)
    boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D>, double>(*((Halo<0, PLEdge, Mesh<Dim2D>, double>*)(bc_i.get()))));

  ///partitioning and initial loadbalancing
  auto p_i(Partitioning<Mem::Main,
                       Algo::Generic,
                       Dim2D,
                       0,
                       pl_vertex>::execute(mesh,
                                           boundaries,
                                           Index(nprocs), Index(rank),
                                           attrs
                                           ));
  if(!MeshUtil::iz_property(p_i.basemesh, attrs.at(0), attrs.at(1)))
    log.checkin_line(" ...iz prop hurt at basemesh!");

  if(!MeshUtil::iz_property(*((Mesh<Dim2D>*)(p_i.submesh.get())), *( (Attribute<double>*)(p_i.attrs.at(0).get()) ), *( (Attribute<double>*)(p_i.attrs.at(1).get()) )))
    log.checkin_line(" ...iz prop hurt at submesh!");



  ///convert mesh into geometry mesh
  ///TODO mesh interface analogously to HaloInterface
  typedef ConformalMesh<Shape::Hypercube<2> > confmeshtype_;
  Index* size_set(new Index[3]);
  MeshControl<dim_2D>::fill_sizes(*((Mesh<Dim2D>*)(p_i.submesh.get())), size_set);
  confmeshtype_ confmesh(size_set);
  MeshControl<dim_2D>::fill_adjacencies(*((Mesh<Dim2D>*)(p_i.submesh.get())), confmesh);
  MeshControl<dim_2D>::fill_vertex_sets(*((Mesh<Dim2D>*)(p_i.submesh.get())), confmesh, *( (Attribute<double>*)(p_i.attrs.at(0).get()) ), *( (Attribute<double>*)(p_i.attrs.at(1).get()) ));
  delete[] size_set;


  ///convert comm_halos
  ///TODO: Peter-> why is CSS's copy assignment private, while copy-ctor is not?
  std::vector<std::shared_ptr<CellSubSet<Shape::Hypercube<2> > > > halo_subsets;
  for(auto& ch_i : p_i.comm_halos)
  {
    CellSubSet<Shape::Hypercube<2> > cell_sub_set(HaloInterface<0, Dim2D>::convert(ch_i.get()));
    halo_subsets.push_back(std::shared_ptr<CellSubSet<Shape::Hypercube<2> > >(new CellSubSet<Shape::Hypercube<2> >(cell_sub_set)));
  }

  ///Trafo, Space, and boundary assembler
  Trafo::Standard::Mapping<ConformalMesh<Shape::Hypercube<2> > > trafo(confmesh);
  Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > space(trafo);
  Assembly::DirichletAssembler<Space::Lagrange1::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2> > > > > dirichlet(space);

  ///convert boundaries and add to dirichlet assembler
  std::vector<std::shared_ptr<CellSubSet<Shape::Hypercube<2> > > > boundary_subsets;
  for(auto& b_i : p_i.boundaries)
  {
    CellSubSet<Shape::Hypercube<2> > cell_sub_set(HaloInterface<0, Dim2D>::convert(&b_i));
    boundary_subsets.push_back(std::shared_ptr<CellSubSet<Shape::Hypercube<2> > >(new CellSubSet<Shape::Hypercube<2> >(cell_sub_set)));
    dirichlet.add_cell_set(cell_sub_set);
  }
  ///vector mirrors
  std::vector<VectorMirror<Mem::Main, double> > mirrors;
  for(auto& hs_i : halo_subsets)
  {
    VectorMirror<Mem::Main, double> target_mirror;
    Assembly::MirrorAssembler::assemble_mirror(target_mirror, space, *hs_i);
    mirrors.push_back(std::move(target_mirror));
  }

  ///buffers
  std::vector<DenseVector<Mem::Main, double> > sendbufs;
  std::vector<DenseVector<Mem::Main, double> > recvbufs;
  std::vector<DenseVector<Mem::Main, double> > mbufs;
  for(auto& m_i : mirrors)
  {
    DenseVector<Mem::Main, double> smbuf(m_i.size());
    DenseVector<Mem::Main, double> rmbuf(m_i.size());
    DenseVector<Mem::Main, double> mmbuf(m_i.size());
    sendbufs.push_back(std::move(smbuf));
    recvbufs.push_back(std::move(rmbuf));
    mbufs.push_back(std::move(mmbuf));
  }
  if(mirrors.size() > 0)
  {
    DenseVector<Mem::Main, double> sbuf(mirrors.at(0).size());
    DenseVector<Mem::Main, double> rbuf(mirrors.at(0).size());
    sendbufs.push_back(std::move(sbuf));
    recvbufs.push_back(std::move(rbuf));
  }

  ///destination ranks
  std::vector<Index> other_ranks;
  for(auto& ch_i : p_i.comm_halos)
    other_ranks.push_back(ch_i->get_other());

  ///assembly
  SparseMatrixCSR<Mem::Main, double> mat_sys;
  Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);
  mat_sys.format();
  Cubature::DynamicFactory cubature_factory("gauss-legendre:2");
  Assembly::Common::LaplaceOperator laplace;
  Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, laplace, space, cubature_factory);

  DenseVector<Mem::Main, double> fbuf(mat_sys.rows());

  auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mbufs, fbuf));

  DenseVector<Mem::Main, double> vec_rhs(space.get_num_dofs(), double(0));
  Assembly::Common::ConstantFunction rhs_func(1.0);
  Assembly::Common::ForceFunctional<Assembly::Common::ConstantFunction> rhs_functional(rhs_func);
  Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, rhs_functional, space, cubature_factory);

  DenseVector<Mem::Main, double> vec_sol(space.get_num_dofs(), double(0));

  UnitFilter<Mem::Main, double> filter(space.get_num_dofs());
  dirichlet.assemble(filter);

  //synching for type-0 -> type-1
  auto tags(HaloTags::value(p_i.comm_halos));
  auto mat_localsys(MatrixConversion<Mem::Main, double, Index, SparseMatrixCSR>::value(mat_sys, mirrors, other_ranks, tags));
  GlobalSynchVec0<Mem::Main, Algo::Generic>::exec(vec_rhs,
                                                  mirrors,
                                                  other_ranks,
                                                  sendbufs,
                                                  recvbufs,
                                                  tags);

  SparseMatrixCOO<Mem::Main, double> mat_precon_temp(mat_localsys.rows(), mat_localsys.columns());
  for(Index i(0) ; i < mat_localsys.rows() ; ++i)
    mat_precon_temp(i, i, double(0.75) * (double(1)/mat_localsys(i, i)));

  SparseMatrixCSR<Mem::Main, double> mat_precon(mat_precon_temp);

  ///filter system
  filter.filter_mat<Algo::Generic>(mat_sys);
  filter.filter_mat<Algo::Generic>(mat_localsys);
  filter.filter_rhs<Algo::Generic>(vec_rhs);
  filter.filter_sol<Algo::Generic>(vec_sol);
  filter.filter_mat<Algo::Generic>(mat_precon); //TODO: check if -> NO! we do this in the solver program when applying the correction filter after preconditioning

  SynchronisedPreconditionedFilteredScaRCData<double,
                                              Mem::Main,
                                              DenseVector<Mem::Main, double>,
                                              VectorMirror<Mem::Main, double>,
                                              SparseMatrixCSR<Mem::Main, double>,
                                              SparseMatrixCSR<Mem::Main, double>,
                                              UnitFilter<Mem::Main, double> > data(std::move(mat_sys), std::move(mat_precon), std::move(vec_sol), std::move(vec_rhs), std::move(filter));

  data.vector_mirrors() = std::move(mirrors);
  data.vector_mirror_sendbufs() = std::move(sendbufs);
  data.vector_mirror_recvbufs() = std::move(recvbufs);
  data.dest_ranks() = std::move(other_ranks);
  data.tags() = std::move(tags);

#ifndef SERIAL
  Communicator c(MPI_COMM_WORLD);
#else
  Communicator c(0);
#endif
  data.communicators().push_back(std::move(c));

  //data.source_ranks() = std::move(sourceranks);
  data.localsys() = std::move(mat_localsys);

  data.halo_frequencies() = std::move(frequencies);
  /*log.checkin_line(" basemesh: ");
  for(Index i(0) ; i < p_i.basemesh.num_polytopes(pl_face) ; ++i)
  {
    log.checkin_line("   QUAD: ", i);
    auto verts(p_i.basemesh.get_adjacent_polytopes(pl_face, pl_vertex, i));
    log.checkin_line("     v=(" , attrs.at(0).at(verts.at(0)) , "," , attrs.at(1).at(verts.at(0)) , ")");
    log.checkin_line("     v=(" , attrs.at(0).at(verts.at(1)) , "," , attrs.at(1).at(verts.at(1)) , ")");
    log.checkin_line("     v=(" , attrs.at(0).at(verts.at(2)) , "," , attrs.at(1).at(verts.at(2)) , ")");
    log.checkin_line("     v=(" , attrs.at(0).at(verts.at(3)) , "," , attrs.at(1).at(verts.at(3)) , ")");

    auto edges(p_i.basemesh.get_adjacent_polytopes(pl_face, pl_edge, i));
    for(auto e_i : edges)
    {
      log.checkin_line("   EDGE: ", e_i);
      auto verts_e(p_i.basemesh.get_adjacent_polytopes(pl_edge, pl_vertex, e_i));
      log.checkin_line("     v=(" , attrs.at(0).at(verts_e.at(0)) , "," , attrs.at(1).at(verts_e.at(0)) , ")");
      log.checkin_line("     v=(" , attrs.at(0).at(verts_e.at(1)) , "," , attrs.at(1).at(verts_e.at(1)) , ")");
    }
  }*/
  /*log.checkin_line(" submesh: ");
  for(Index i(0) ; i < ((Mesh<Dim2D>*)(p_i.submesh.get()))->num_polytopes(pl_face) ; ++i)
  {
    log.checkin_line("   QUAD: ", i);
    auto verts(((Mesh<Dim2D>*)(p_i.submesh.get()))->get_adjacent_polytopes(pl_face, pl_vertex, i));
    log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(verts.at(0)) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(verts.at(0)) , ")");
    log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(verts.at(1)) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(verts.at(1)) , ")");
    log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(verts.at(2)) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(verts.at(2)) , ")");
    log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(verts.at(3)) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(verts.at(3)) , ")");

    auto edges(((Mesh<Dim2D>*)(p_i.submesh.get()))->get_adjacent_polytopes(pl_face, pl_edge, i));
    for(auto e_i : edges)
    {
      log.checkin_line("   EDGE: ", e_i);
      auto verts_e(((Mesh<Dim2D>*)(p_i.submesh.get()))->get_adjacent_polytopes(pl_edge, pl_vertex, e_i));
      log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(verts_e.at(0)) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(verts_e.at(0)) , ")");
      log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(verts_e.at(1)) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(verts_e.at(1)) , ")");
    }
  }*/
  /*log.checkin_line(" comm-halos: ");
  for(Index i(0) ; i < p_i.comm_halos.size(); ++i)
  {
    Index other(p_i.comm_halos.at(i)->get_other());
    log.checkin_line("   HALO: ", i, " links to ", other);
    if(p_i.comm_halos.at(i)->get_level() == pl_vertex)
    {
      for(Index j(0) ; j < p_i.comm_halos.at(i)->size() ; ++j)
      {
        log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(Index(p_i.comm_halos.at(i)->get_element(j))) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(Index(p_i.comm_halos.at(i)->get_element(j))) , ")");
      }
    }
    else
    {
      for(Index j(0) ; j < p_i.comm_halos.at(i)->size() ; ++j)
      {
        auto verts(((Mesh<Dim2D>*)(p_i.submesh.get()))->get_adjacent_polytopes(pl_edge, pl_vertex, p_i.comm_halos.at(i)->get_element(j)));
        log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(Index(verts.at(0))) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(Index(verts.at(0))) , ")");
        log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(Index(verts.at(1))) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(Index(verts.at(1))) , ")");
      }
    }
  }
  log.checkin_line(" boundaries: ");
  for(Index i(0) ; i < p_i.boundaries.size(); ++i)
  {
    log.checkin_line("   BC: ", i);
    for(Index j(0) ; j < p_i.boundaries.at(i).size() ; ++j)
    {
      auto verts(((Mesh<Dim2D>*)(p_i.submesh.get()))->get_adjacent_polytopes(pl_edge, pl_vertex, p_i.boundaries.at(i).get_element(j)));
      log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(Index(verts.at(0))) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(Index(verts.at(0))) , ")");
      log.checkin_line("     v=(" , ( (Attribute<double>*)(p_i.attrs.at(0).get()) )->at(Index(verts.at(1))) , "," , ( (Attribute<double>*)(p_i.attrs.at(1).get()) )->at(Index(verts.at(1))) , ")");
    }
  }*/
  /*log.checkin_line("converted submesh:");
  auto submesh_vertex_set(confmesh.get_vertex_set());
  for(Index i(0) ; i < 4 ; ++i)
  {
    log.checkin_line("    v=(", submesh_vertex_set[i][0], "," , submesh_vertex_set[i][1], ")");
  }*/

  log.checkin_line(" A0: ", data.sys());
  log.checkin_line(" A1: ", data.localsys());
  log.checkin_line(" P1: ", data.precon());
  log.checkin_line(" b1: ", data.rhs());
  log.checkin_line(" x1: ", data.sol());

  if(Comm::size() > 1)
  {
    ///layer 1 (global layer)
    std::shared_ptr<ScaRCFunctorBase<double,
      Mem::Main,
      DenseVector<Mem::Main, double>,
      VectorMirror<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      UnitFilter<Mem::Main, double>,
      std::vector,
      Index,
      Algo::Generic> > solver(new ScaRCFunctorPCG0<double,
          Mem::Main,
          DenseVector<Mem::Main, double>,
          VectorMirror<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          UnitFilter<Mem::Main, double>,
          std::vector,
          Index,
          Algo::Generic>(data) );

    ///layer 1 (global layer), preconditioner
    std::shared_ptr<ScaRCFunctorBase<double,
      Mem::Main,
      DenseVector<Mem::Main, double>,
      VectorMirror<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      UnitFilter<Mem::Main, double>,
      std::vector,
      Index,
      Algo::Generic> > block_smoother(new ScaRCFunctorPreconBlock<double,
          Mem::Main,
          DenseVector<Mem::Main, double>,
          VectorMirror<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          UnitFilter<Mem::Main, double>,
          std::vector,
          Index,
          Algo::Generic>(data) );

    ///layer 0 (local layer)
    std::shared_ptr<ScaRCFunctorBase<double,
      Mem::Main,
      DenseVector<Mem::Main, double>,
      VectorMirror<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      UnitFilter<Mem::Main, double>,
      std::vector,
      Index,
      Algo::Generic> > local_solver(new ScaRCFunctorPCG1<double,
          Mem::Main,
          DenseVector<Mem::Main, double>,
          VectorMirror<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          UnitFilter<Mem::Main, double>,
          std::vector,
          Index,
          Algo::Generic>(data) );

    ///layer 0 (local layer), preconditioner
    std::shared_ptr<ScaRCFunctorBase<double,
      Mem::Main,
      DenseVector<Mem::Main, double>,
      VectorMirror<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      UnitFilter<Mem::Main, double>,
      std::vector,
      Index,
      Algo::Generic> > local_precon(new ScaRCFunctorPreconSpM1V1<double,
          Mem::Main,
          DenseVector<Mem::Main, double>,
          VectorMirror<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          UnitFilter<Mem::Main, double>,
          std::vector,
          Index,
          Algo::Generic>(data) );


    solver->reset_preconditioner(block_smoother);
    block_smoother->reset_preconditioner(local_solver);
    local_solver->reset_preconditioner(local_precon);

    local_solver->conv_check() = true;
    local_solver->max_iters() = 1000;

    solver->execute();
    log.checkin_line("#iters global CG: ", solver->iterations());
    log.checkin_line("#iters local CG: ", local_solver->iterations());
    log.checkin_line("sol: ", data.sol());
  }
  else
  {
    ///layer 0 (local layer)
    std::shared_ptr<ScaRCFunctorBase<double,
      Mem::Main,
      DenseVector<Mem::Main, double>,
      VectorMirror<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      UnitFilter<Mem::Main, double>,
      std::vector,
      Index,
      Algo::Generic> > local_solver(new ScaRCFunctorPCG1<double,
          Mem::Main,
          DenseVector<Mem::Main, double>,
          VectorMirror<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          UnitFilter<Mem::Main, double>,
          std::vector,
          Index,
          Algo::Generic>(data) );

    ///layer 0 (local layer), preconditioner
    std::shared_ptr<ScaRCFunctorBase<double,
      Mem::Main,
      DenseVector<Mem::Main, double>,
      VectorMirror<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      SparseMatrixCSR<Mem::Main, double>,
      UnitFilter<Mem::Main, double>,
      std::vector,
      Index,
      Algo::Generic> > local_precon(new ScaRCFunctorPreconSpM1V1<double,
          Mem::Main,
          DenseVector<Mem::Main, double>,
          VectorMirror<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          SparseMatrixCSR<Mem::Main, double>,
          UnitFilter<Mem::Main, double>,
          std::vector,
          Index,
          Algo::Generic>(data) );


    local_solver->reset_preconditioner(local_precon);

    local_solver->conv_check() = true;
    local_solver->max_iters() = 1000;

    local_solver->execute();
    log.checkin_line("#iters local CG: ", local_solver->iterations());
    log.checkin_line("sol: ", data.sol());
  }
  /*std::cout << rank << ", #iters global CG: " << solver->iterations() << std::endl;
  std::cout << rank << ", #iters local CG: " << local_solver->iterations() << std::endl;
  std::cout << rank << " sol: " << data.sol() << std::endl;*/

  log.synch();

  Comm::barrier();
  if(rank == 0)
    std::cout << log.msg << std::endl;

#ifndef SERIAL
  MPI_Finalize();
#endif

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

  target_mesh.add_adjacency(pl_vertex, pl_edge, 0, 0);
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
  target_mesh.add_adjacency(pl_vertex, pl_face, 3, 0);

  boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D> >(target_mesh));
  boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D> >(target_mesh));
  boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D> >(target_mesh));
  boundaries.push_back(Halo<0, PLEdge, Mesh<Dim2D> >(target_mesh));
  boundaries.at(0).push_back(0);
  boundaries.at(1).push_back(1);
  boundaries.at(2).push_back(2);
  boundaries.at(3).push_back(3);
}
#endif
