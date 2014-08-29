#include <test_system/test_system.hpp>

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
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;

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

#ifdef SERIAL
template<typename Tag_, typename Algo_, typename DataType_>
class ScaRCFunctorTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    ScaRCFunctorTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("ScaRCFunctorTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      int rank(0), nprocs(1);
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

      std::vector<DenseVector<Mem::Main, double> > freq_buffers;
      for(Index i(0) ; i < mbufs.size() ; ++i)
      {
        DenseVector<Mem::Main, double> fbuf(mat_sys.rows());
        freq_buffers.push_back(std::move(fbuf));
      }

      auto frequencies(HaloFrequencies<Mem::Main, Algo::Generic>::value(mirrors, mbufs, freq_buffers));

      DenseVector<Mem::Main, double> vec_rhs(space.get_num_dofs(), double(0));
      Assembly::Common::ConstantFunction rhs_func(1.0);
      Assembly::Common::ForceFunctional<Assembly::Common::ConstantFunction> rhs_functional(rhs_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, rhs_functional, space, cubature_factory);

      DenseVector<Mem::Main, double> vec_sol(space.get_num_dofs(), double(0));

      UnitFilter<Mem::Main, double> filter(space.get_num_dofs());
      dirichlet.assemble(filter);

      auto tags(HaloTags::value(p_i.comm_halos));
      //synching for type-0 -> type-1
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

      //data.source_ranks() = std::move(sourceranks);
      data.localsys() = std::move(mat_localsys);

      data.halo_frequencies() = std::move(frequencies);
      data.tags() = std::move(tags);
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

      //std::cout << data.sol() << std::endl;
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(0), DataType_(0), std::numeric_limits<DataType_>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(1), DataType_(0), std::numeric_limits<DataType_>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(2), DataType_(0), std::numeric_limits<DataType_>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(3), DataType_(0), std::numeric_limits<DataType_>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(4), DataType_(0), std::numeric_limits<DataType_>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(5), DataType_(0), std::numeric_limits<DataType_>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(6), DataType_(0), std::numeric_limits<DataType_>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(7), DataType_(0), std::numeric_limits<DataType_>::epsilon());
      TEST_CHECK_EQUAL_WITHIN_EPS(data.sol()(8), DataType_(0.09375), std::numeric_limits<DataType_>::epsilon());
    }
};
ScaRCFunctorTest<Mem::Main, Algo::Generic,  double> sf_cpu_double("ELL double");
#endif
