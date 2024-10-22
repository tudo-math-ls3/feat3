// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/voxel_amavanka.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/assembly/function_integral_jobs.hpp>
#include <kernel/voxel_assembly/defo_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/geometry/unit_cube_patch_generator.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/statistics.hpp>

#ifdef FEAT_HAVE_OMP
#include "omp.h"
#else
// define omp function if omp is not found
inline void omp_set_num_threads(int){}
#endif

namespace SaddlePointAssemblyBench
{
  using namespace FEAT;

  template<typename Trafo_>
  struct RTQ0
  {
    typedef Space::CroRavRanTur::Element<Trafo_> V;
    typedef Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<0>> P;
    static const char* name() {return "RT/Q0";}
  };

  template<typename Trafo_>
  struct Q2P1
  {
    typedef Space::Lagrange2::Element<Trafo_> V;
    typedef Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<1>> P;
    static const char* name() {return "Q2/P1";}
  };

  template<typename DT_, typename IT_, int dim_>
  struct SpaceHelper
  {
    static constexpr int dim = dim_;
    typedef DT_ DataType;
    typedef IT_ IndexType;
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  };

  struct BenchResults
  {
    int num_worker_threads;
    /// coloring time
    double coloring_time;
    /// total symbolic matrix assembly time in seconds
    double sym_mat_asm_time;
    /// total symbolic vanka assembly time in seconds
    double sym_vanka_asm_time;
    /// total numeric matrix assembly time in seconds
    double num_asm_mat_time;
    /// total numeric vanka assembly time in seconds
    double num_asm_vanka_time;
    /// total matrix-vector product time in seconds
    double solver_apply_time;
    double all_time;
    /// number of numeric assemblies performed
    int num_elements;
    /// number of numeric assemblies performed
    int matrix_nnz;
    /// number of matrix-vector products performed
    int vanka_nnz;

    BenchResults() :
      num_worker_threads(1),
      coloring_time(0.0),
      sym_mat_asm_time(0.0), sym_vanka_asm_time(0.0), num_asm_mat_time(0.0), num_asm_vanka_time(0.0),
      solver_apply_time(0.0), all_time(0.0), num_elements(0), matrix_nnz(0), vanka_nnz(0)
    {
    }
  }; // struct BenchResults

  constexpr double omega = 0.5;
  constexpr int num_outer_iterations = 5;
  constexpr int num_solver_iterations = 100;


  template<template<typename> class Space_, typename DT_, typename IT_, int dim_>
  void bench_amavanka_inner(Geometry::RootMeshNode<typename SpaceHelper<DT_, IT_, dim_>::MeshType>& mesh_node, BenchResults& bres, PreferredBackend backend_asm, PreferredBackend backend_calc)
  {
    constexpr int dim = dim_;
    typedef SpaceHelper<DT_, IT_, dim_> SPHelp;
    typedef typename SPHelp::DataType DataType;
    typedef typename SPHelp::IndexType IndexType;
    typedef typename SPHelp::MeshType MeshType;
    typedef typename SPHelp::TrafoType TrafoType;
    typedef typename Space_<TrafoType>::V SpaceVeloType;
    typedef typename Space_<TrafoType>::P SpacePresType;

    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, dim_, dim_> MatrixTypeA;
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, dim_, 1> MatrixTypeB;
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, dim_> MatrixTypeD;
    typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> MatrixType;
    typedef typename MatrixType::VectorTypeL VectorType;
    typedef LAFEM::UnitFilterBlocked<DT_, IT_, dim> VeloFilterType;
    typedef LAFEM::NoneFilter<DT_, IT_> PresFilterType;
    typedef LAFEM::TupleFilter<VeloFilterType, PresFilterType> FilterType;

    StopWatch asm_sym_matrix_watch, asm_sym_vanka_watch, asm_num_matrix_watch, asm_num_vanka_watch;
    StopWatch coloring_watch, solver_apply_watch, all_watch;

    all_watch.start();
    // get our mesh
    MeshType& mesh = *mesh_node.get_mesh();
    TrafoType trafo(mesh);
    SpaceVeloType space_v(trafo);
    SpacePresType space_p(trafo);

    //create coloring from our mesh
    coloring_watch.start();
    const auto& verts_at_elem = mesh.template get_index_set<dim, 0>();
    Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);
    Adjacency::Graph elems_at_elem(Adjacency::RenderType::injectify, verts_at_elem, elems_at_vert);

    // create coloring
    Adjacency::Coloring col(elems_at_elem);
    coloring_watch.stop();

    MatrixType matrix;
    MatrixTypeA& mat_a = matrix.block_a();
    MatrixTypeB& mat_b = matrix.block_b();
    MatrixTypeD& mat_d = matrix.block_d();
    FilterType filter;

    Cubature::DynamicFactory cubature("auto-degree:5");


    matrix.format();
    Backend::set_preferred_backend(backend_asm);
    // assemble A, B and D
    asm_sym_matrix_watch.start();
    Assembly::SymbolicAssembler::assemble_matrix_std1(mat_a, space_v);
    VoxelAssembly::VoxelDefoAssembler<SpaceVeloType, DataType, IndexType> voxel_defo(space_v, col);
    voxel_defo.nu = DataType(2.3);
    asm_sym_matrix_watch.stop();

    asm_num_matrix_watch.start();
    // Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a, dudv_op, space_v, cubature);
    voxel_defo.assemble_matrix1(mat_a, space_v, Cubature::DynamicFactory(cubature));
    Assembly::GradPresDivVeloAssembler::assemble(mat_b, mat_d, space_v, space_p, cubature);
    asm_num_matrix_watch.stop();


    // assemble filter
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:0"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:1"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:2"));

    // finally, assemble the filters
    if constexpr(dim == 2)
    {
      Analytic::Common::ParProfileVector<DT_> inflow_func(0.0, 0.0, 0.0, 1.0, 1.0);
      unit_asm.assemble(filter.template at<0>(), space_v, inflow_func);
    }
    else
    {
      Analytic::Common::PoiseuillePipeFlow<DT_, dim> inflow_func(Tiny::Vector<DT_, dim>{0.0, 0.5, 0.5}, Tiny::Vector<DT_, dim>{1.0, 0., 0.}, DT_(0.5));
      unit_asm.assemble(filter.template at<0>(), space_v, inflow_func);
    }

    // filter matrices
    filter.template at<0>().filter_mat(mat_a);
    filter.template at<0>().filter_offdiag_row_mat(mat_b);

    // create RHS and SOL
    VectorType vec_rhs = matrix.create_vector_l();
    VectorType vec_sol = matrix.create_vector_l();
    VectorType vec_def = matrix.create_vector_l();
    vec_rhs.format();
    vec_sol.format();
    vec_def.format();

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    asm_sym_vanka_watch.start();
    auto vanka = Solver::new_amavanka(matrix, filter, omega);
    auto fgmres = Solver::new_gmres(matrix, filter, Index(std::min(num_solver_iterations,16)), DataType(0), vanka);
    fgmres->set_min_iter(num_solver_iterations);
    fgmres->set_max_iter(num_solver_iterations);
    fgmres->skip_defect_calc(true);
    fgmres->set_plot_mode(Solver::PlotMode::none);
    fgmres->init_symbolic();
    asm_sym_vanka_watch.stop();
    bres.matrix_nnz = int(matrix.template used_elements<LAFEM::Perspective::pod>());
    bres.vanka_nnz = int(vanka->data_size());
    bres.num_elements = int(mesh.get_num_elements());

    for(int k = 0; k < num_outer_iterations; ++k)
    {
      Backend::set_preferred_backend(backend_asm);
      asm_num_matrix_watch.start();
      mat_a.format();
      voxel_defo.assemble_matrix1(mat_a, space_v, Cubature::DynamicFactory(cubature));
      asm_num_matrix_watch.stop();
      asm_num_vanka_watch.start();
      fgmres->init_numeric();
      asm_num_vanka_watch.stop();
      Backend::set_preferred_backend(backend_calc);
      solver_apply_watch.start();
      vec_sol.format();
      fgmres->apply(vec_sol, vec_rhs);
      vec_def.axpy(vec_sol, DataType(1));
      solver_apply_watch.stop();
      Backend::set_preferred_backend(backend_asm);
      asm_num_vanka_watch.start();
      fgmres->done_numeric();
      asm_num_vanka_watch.stop();

    }

    fgmres->done_symbolic();
    all_watch.stop();

    bres.coloring_time = coloring_watch.elapsed();
    bres.num_asm_mat_time = asm_num_matrix_watch.elapsed();
    bres.num_asm_vanka_time = asm_num_vanka_watch.elapsed();
    bres.sym_vanka_asm_time = asm_sym_vanka_watch.elapsed();
    bres.sym_mat_asm_time = asm_sym_vanka_watch.elapsed();
    bres.solver_apply_time = solver_apply_watch.elapsed();
    bres.all_time = all_watch.elapsed();
  }

  template<template<typename> class Space_, typename DT_, typename IT_, int dim_>
  void bench_voxel_amavanka_inner(Geometry::RootMeshNode<typename SpaceHelper<DT_, IT_, dim_>::MeshType>& mesh_node, BenchResults& bres, PreferredBackend backend_asm, PreferredBackend backend_calc)
  {
    constexpr int dim = dim_;
    typedef SpaceHelper<DT_, IT_, dim_> SPHelp;
    typedef typename SPHelp::DataType DataType;
    typedef typename SPHelp::IndexType IndexType;
    typedef typename SPHelp::MeshType MeshType;
    typedef typename SPHelp::TrafoType TrafoType;
    typedef typename Space_<TrafoType>::V SpaceVeloType;
    typedef typename Space_<TrafoType>::P SpacePresType;

    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, dim_, dim_> MatrixTypeA;
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, dim_, 1> MatrixTypeB;
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, dim_> MatrixTypeD;
    typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> MatrixType;
    typedef typename MatrixType::VectorTypeL VectorType;
    typedef LAFEM::UnitFilterBlocked<DT_, IT_, dim> VeloFilterType;
    typedef LAFEM::NoneFilter<DT_, IT_> PresFilterType;
    typedef LAFEM::TupleFilter<VeloFilterType, PresFilterType> FilterType;

    StopWatch asm_sym_matrix_watch, asm_sym_vanka_watch, asm_num_matrix_watch, asm_num_vanka_watch;
    StopWatch coloring_watch, solver_apply_watch, all_watch;

    all_watch.start();
    // get our mesh
    MeshType& mesh = *mesh_node.get_mesh();
    TrafoType trafo(mesh);
    SpaceVeloType space_v(trafo);
    SpacePresType space_p(trafo);

    //create coloring from our mesh
    coloring_watch.start();
    const auto& verts_at_elem = mesh.template get_index_set<dim, 0>();
    Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);
    Adjacency::Graph elems_at_elem(Adjacency::RenderType::injectify, verts_at_elem, elems_at_vert);

    // create coloring
    Adjacency::Coloring col(elems_at_elem);
    coloring_watch.stop();

    MatrixType matrix;
    MatrixTypeA& mat_a = matrix.block_a();
    MatrixTypeB& mat_b = matrix.block_b();
    MatrixTypeD& mat_d = matrix.block_d();
    FilterType filter;

    Cubature::DynamicFactory cubature("auto-degree:5");


    matrix.format();
    Backend::set_preferred_backend(backend_asm);
    // assemble A, B and D
    asm_sym_matrix_watch.start();
    Assembly::SymbolicAssembler::assemble_matrix_std1(mat_a, space_v);
    VoxelAssembly::VoxelDefoAssembler<SpaceVeloType, DataType, IndexType> voxel_defo(space_v, col);
    voxel_defo.nu = DataType(2.3);
    asm_sym_matrix_watch.stop();

    asm_num_matrix_watch.start();
    // Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_a, dudv_op, space_v, cubature);
    voxel_defo.assemble_matrix1(mat_a, space_v, Cubature::DynamicFactory(cubature));
    Assembly::GradPresDivVeloAssembler::assemble(mat_b, mat_d, space_v, space_p, cubature);
    asm_num_matrix_watch.stop();


    // assemble filter
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:0"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:1"));
    unit_asm.add_mesh_part(*mesh_node.find_mesh_part("bnd:2"));

    // finally, assemble the filters
    if constexpr(dim == 2)
    {
      Analytic::Common::ParProfileVector<DT_> inflow_func(0.0, 0.0, 0.0, 1.0, 1.0);
      unit_asm.assemble(filter.template at<0>(), space_v, inflow_func);
    }
    else
    {
      Analytic::Common::PoiseuillePipeFlow<DT_, dim> inflow_func(Tiny::Vector<DT_, dim>{0.0, 0.5, 0.5}, Tiny::Vector<DT_, dim>{1.0, 0., 0.}, DT_(0.5));
      unit_asm.assemble(filter.template at<0>(), space_v, inflow_func);
    }

    // filter matrices
    filter.template at<0>().filter_mat(mat_a);
    filter.template at<0>().filter_offdiag_row_mat(mat_b);

    // create RHS and SOL
    VectorType vec_rhs = matrix.create_vector_l();
    VectorType vec_sol = matrix.create_vector_l();
    VectorType vec_def = matrix.create_vector_l();
    vec_rhs.format();
    vec_sol.format();
    vec_def.format();

    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    asm_sym_vanka_watch.start();
    auto vanka = Solver::new_voxel_amavanka<MatrixType, FilterType, Adjacency::Coloring, FEAT::Intern::VankaAssemblyPolicy::batchedAssembly, FEAT::Intern::VankaMacroPolicy::uniformMacros>(matrix, filter, col, omega);
    auto fgmres = Solver::new_gmres(matrix, filter, Index(std::min(num_solver_iterations, 16)), DataType(0), vanka);
    fgmres->set_min_iter(num_solver_iterations);
    fgmres->set_max_iter(num_solver_iterations);
    fgmres->skip_defect_calc(true);
    fgmres->set_plot_mode(Solver::PlotMode::none);
    fgmres->init_symbolic();
    asm_sym_vanka_watch.stop();
    bres.matrix_nnz = int(matrix.template used_elements<LAFEM::Perspective::pod>());
    bres.vanka_nnz = int(vanka->data_size());
    bres.num_elements = int(mesh.get_num_elements());

    for(int k = 0; k < num_outer_iterations; ++k)
    {
      Backend::set_preferred_backend(backend_asm);
      asm_num_matrix_watch.start();
      mat_a.format();
      voxel_defo.assemble_matrix1(mat_a, space_v, Cubature::DynamicFactory(cubature));
      asm_num_matrix_watch.stop();
      asm_num_vanka_watch.start();
      fgmres->init_numeric();
      asm_num_vanka_watch.stop();
      Backend::set_preferred_backend(backend_calc);
      solver_apply_watch.start();
      vec_sol.format();
      fgmres->apply(vec_sol, vec_rhs);
      vec_def.axpy(vec_sol, DataType(1));
      solver_apply_watch.stop();
      Backend::set_preferred_backend(backend_asm);
      asm_num_vanka_watch.start();
      fgmres->done_numeric();
      asm_num_vanka_watch.stop();

    }

    fgmres->done_symbolic();
    all_watch.stop();

    bres.coloring_time = coloring_watch.elapsed();
    bres.num_asm_mat_time = asm_num_matrix_watch.elapsed();
    bres.num_asm_vanka_time = asm_num_vanka_watch.elapsed();
    bres.sym_vanka_asm_time = asm_sym_vanka_watch.elapsed();
    bres.sym_mat_asm_time = asm_sym_vanka_watch.elapsed();
    bres.solver_apply_time = solver_apply_watch.elapsed();
    bres.all_time = all_watch.elapsed();
  }

  template<template<typename> class Space_, typename DT_, typename IT_, int dim_>
  void run_bench_amavanka(const int level, BenchResults& bres,
            PreferredBackend backend_asm, PreferredBackend backend_calc)
  {
    typedef typename SpaceHelper<DT_, IT_, dim_>::MeshType MeshType;
    typedef Geometry::RootMeshNode<MeshType> MeshNode;
    // create mesh node
    std::vector<int> ranks;
    std::unique_ptr<MeshNode> mesh_node;
    Geometry::UnitCubePatchGenerator<MeshType>::create_unique(0, 1, mesh_node, ranks);

    // refine a few times
    for(int i = 0; i < level; ++i)
    {
      mesh_node = mesh_node->refine_unique();
    }

    bench_amavanka_inner<Space_, DT_, IT_, dim_>(*mesh_node, bres, backend_asm, backend_calc);

  }

  template<template<typename> class Space_, typename DT_, typename IT_, int dim_>
  void run_bench_voxel_amavanka(const int level, BenchResults& bres,
            PreferredBackend backend_asm, PreferredBackend backend_calc)
  {
    typedef typename SpaceHelper<DT_, IT_, dim_>::MeshType MeshType;
    typedef Geometry::RootMeshNode<MeshType> MeshNode;
    // create mesh node
    std::vector<int> ranks;
    std::unique_ptr<MeshNode> mesh_node;
    Geometry::UnitCubePatchGenerator<MeshType>::create_unique(0, 1, mesh_node, ranks);

    // refine a few times
    for(int i = 0; i < level; ++i)
    {
      mesh_node = mesh_node->refine_unique();
    }

    bench_voxel_amavanka_inner<Space_, DT_, IT_, dim_>(*mesh_node, bres, backend_asm, backend_calc);

  }

  template<template<typename> class Space_, typename DT_, typename IT_, int dim_>
  void run(SimpleArgParser& args)
  {
    // parse levels
    Index lvl_min(0);
    Index lvl_max(1);
    args.parse("level", lvl_max, lvl_min);
    const bool generic_run = args.check("no-generic-apply") < 0;


    // parse threads
    std::vector<std::size_t> num_threads;
    {
      auto* p = args.query("threads");
      if(p != nullptr)
      {
        for(std::size_t i(0); i < p->second.size(); ++i)
        {
          std::size_t t(0u);
          if(!p->second.at(i).parse(t))
          {
            std::cout << "ERROR: Failed to parse '" << p->second.at(i) << "' as thread count\n";
            return;
          }
          num_threads.push_back(t);
        }
      }
    }
    if(num_threads.empty())
    {
      num_threads.push_back(1u);
    }

    #ifdef FEAT_HAVE_CUDA
    const std::size_t run_size = 3*(lvl_max-lvl_min+1);
    #else
    const std::size_t run_size = (lvl_max-lvl_min+1);
    #endif
    std::deque<std::deque<std::array<BenchResults, 2>>> results(run_size);

    std::cout << "Starting benchmark...\n";
    for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
    {
      std::cout << "Starting level " << lvl << "\n";
      const Index cnt = lvl-lvl_min;
      if(generic_run)
      {
        std::cout << "Starting generic generic test...\n";
        auto& result = results.at(cnt);
        result.resize(num_threads.size());
        for(Index i = 0; i < num_threads.size(); ++i)
        {
          PreferredBackend backend_asm = PreferredBackend::generic;
          PreferredBackend backend_calc = PreferredBackend::mkl;
          omp_set_num_threads(int(num_threads[i]));
          run_bench_amavanka<Space_, DT_ , IT_, dim_>(int(lvl), result.at(i).at(0), backend_asm, backend_calc);
          run_bench_voxel_amavanka<Space_, DT_ , IT_, dim_>(int(lvl), result.at(i).at(1), backend_asm, backend_calc);
        }
      }
      #ifdef FEAT_HAVE_CUDA
      {
        std::cout << "Starting generic cuda test...\n";
        auto& result = results.at(run_size/3 + cnt);
        result.resize(num_threads.size());
        for(Index i = 0; i < num_threads.size(); ++i)
        {
          PreferredBackend backend_asm = PreferredBackend::generic;
          PreferredBackend backend_calc = PreferredBackend::cuda;
          omp_set_num_threads(int(num_threads[i]));
          run_bench_amavanka<Space_, DT_ , IT_, dim_>(int(lvl), result.at(i).at(0), backend_asm, backend_calc);
          run_bench_voxel_amavanka<Space_, DT_ , IT_, dim_>(int(lvl), result.at(i).at(1), backend_asm, backend_calc);
        }
      }
      {
        std::cout << "Starting cuda cuda test...\n";
        auto& result = results.at(2*run_size/3 + cnt);
        result.resize(num_threads.size());
        for(Index i = 0; i < num_threads.size(); ++i)
        {
          PreferredBackend backend_asm = PreferredBackend::cuda;
          PreferredBackend backend_calc = PreferredBackend::cuda;
          omp_set_num_threads(int(num_threads[i]));
          run_bench_amavanka<Space_, DT_ , IT_, dim_>(int(lvl), result.at(i).at(0), backend_asm, backend_calc);
          run_bench_voxel_amavanka<Space_, DT_ , IT_, dim_>(int(lvl), result.at(i).at(1), backend_asm, backend_calc);
        }
      }
      #endif

    }
    std::cout << "Finished benchmark...\n";
    {
      std::cout << "\n------------------------------------\n";
      std::cout << "\n Sizes \n";
      {
        std::cout << " LVL       Num Elements         NNZ Matrix         NNZ Vanka  \n";
        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ":    ";

          if(generic_run)
          {
            const auto& b = results.at(lvl-lvl_min).at(0).at(0);
            std::cout << stringify(b.num_elements).pad_front(10).pad_back(20);
            std::cout << stringify(b.matrix_nnz).pad_front(10).pad_back(20);
            std::cout << stringify(b.vanka_nnz).pad_front(10).pad_back(20) << "\n";
          }
          else
          {
            #ifdef FEAT_HAVE_CUDA
            const auto& b = results.at((run_size/3)+lvl-lvl_min).at(0).at(0);
            std::cout << stringify(b.num_elements).pad_front(10).pad_back(20);
            std::cout << stringify(b.matrix_nnz).pad_front(10).pad_back(20);
            std::cout << stringify(b.vanka_nnz).pad_front(10).pad_back(20) << "\n";
            #endif
          }
        }

      }
      std::cout << "\n------------------------------------\n";
      if(generic_run)
      {
        std::cout << "Generic Generic All Timings:\n";
        String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                    + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
        std::cout << line_h << "\n";
        String line_v = String(" LVL   ");
        for(auto h : num_threads)
        {
          line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
        }
        for(auto h : num_threads)
        {
          line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
        }
        std::cout << line_v << "\n";


        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ":    ";

          for(std::size_t i = 0; i < num_threads.size(); ++i)
          {
            const auto& b = results.at(lvl-lvl_min).at(i).at(0);
            double t = b.all_time;
            std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
          }
          for(std::size_t i = 0; i < num_threads.size(); ++i)
          {
            const auto& b = results.at(lvl-lvl_min).at(i).at(1);
            double t = b.all_time;
            std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
          }
          std::cout << "\n";
        }
      }
    }

#ifdef FEAT_HAVE_CUDA
    {
      // std::cout << "\n------------------------------------\n";
      std::cout << "\nGeneric CUDA All Timings:\n";
      String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                   + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
      std::cout << line_h << "\n";
      String line_v = String(" LVL   ");
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      std::cout << line_v << "\n";


      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ":    ";

        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at((run_size/3)+lvl-lvl_min).at(i).at(0);
          double t = b.all_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at((run_size/3) + lvl-lvl_min).at(i).at(1);
          double t = b.all_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        std::cout << "\n";
      }
    }

    {
      // std::cout << "\n------------------------------------\n";
      std::cout << "\nCUDA CUDA All Timings:\n";
      String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                   + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
      std::cout << line_h << "\n";
      String line_v = String(" LVL   ");
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      std::cout << line_v << "\n";


      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ":    ";

        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at(2*(run_size/3)+lvl-lvl_min).at(i).at(0);
          double t = b.all_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at(2*(run_size/3) + lvl-lvl_min).at(i).at(1);
          double t = b.all_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        std::cout << "\n";
      }
    }
#endif // FEAT_HAVE_CUDA

    {
      std::cout << "\n------------------------------------\n";
      if(generic_run)
      {
        std::cout << "Generic Generic Asm Vanka Num Timings:\n";
        String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                    + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
        std::cout << line_h << "\n";
        String line_v = String(" LVL   ");
        for(auto h : num_threads)
        {
          line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
        }
        for(auto h : num_threads)
        {
          line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
        }
        std::cout << line_v << "\n";


        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ":    ";

          for(std::size_t i = 0; i < num_threads.size(); ++i)
          {
            const auto& b = results.at(lvl-lvl_min).at(i).at(0);
            double t = b.num_asm_vanka_time;
            std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
          }
          for(std::size_t i = 0; i < num_threads.size(); ++i)
          {
            const auto& b = results.at(lvl-lvl_min).at(i).at(1);
            double t = b.num_asm_vanka_time;
            std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
          }
          std::cout << "\n";
        }
      }
    }
#ifdef FEAT_HAVE_CUDA
    {
      // std::cout << "\n------------------------------------\n";
      std::cout << "\nGeneric CUDA Asm Vanka Num Timings:\n";
      String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                   + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
      std::cout << line_h << "\n";
      String line_v = String(" LVL   ");
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      std::cout << line_v << "\n";


      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ":    ";

        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at((run_size/3)+lvl-lvl_min).at(i).at(0);
          double t = b.num_asm_vanka_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at((run_size/3) + lvl-lvl_min).at(i).at(1);
          double t = b.num_asm_vanka_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        std::cout << "\n";
      }
    }

    {
      // std::cout << "\n------------------------------------\n";
      std::cout << "\nCUDA CUDA Asm Vanka Num Timings:\n";
      String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                   + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
      std::cout << line_h << "\n";
      String line_v = String(" LVL   ");
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      std::cout << line_v << "\n";


      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ":    ";

        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at(2*(run_size/3)+lvl-lvl_min).at(i).at(0);
          double t = b.num_asm_vanka_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at(2*(run_size/3) + lvl-lvl_min).at(i).at(1);
          double t = b.num_asm_vanka_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        std::cout << "\n";
      }
    }
#endif // FEAT_HAVE_CUDA

    {
      std::cout << "\n------------------------------------\n";
      if(generic_run)
      {
        std::cout << "Generic Generic Asm Matrix Num Timings:\n";
        String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                    + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
        std::cout << line_h << "\n";
        String line_v = String(" LVL   ");
        for(auto h : num_threads)
        {
          line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
        }
        for(auto h : num_threads)
        {
          line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
        }
        std::cout << line_v << "\n";


        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ":    ";

          for(std::size_t i = 0; i < num_threads.size(); ++i)
          {
            const auto& b = results.at(lvl-lvl_min).at(i).at(0);
            double t = b.num_asm_mat_time;
            std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
          }
          for(std::size_t i = 0; i < num_threads.size(); ++i)
          {
            const auto& b = results.at(lvl-lvl_min).at(i).at(1);
            double t = b.num_asm_mat_time;
            std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
          }
          std::cout << "\n";
        }
      }
    }
#ifdef FEAT_HAVE_CUDA
    {
      // std::cout << "\n------------------------------------\n";
      std::cout << "\nGeneric CUDA Asm Matrix Num Timings:\n";
      String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                   + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
      std::cout << line_h << "\n";
      String line_v = String(" LVL   ");
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      std::cout << line_v << "\n";


      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ":    ";

        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at((run_size/3)+lvl-lvl_min).at(i).at(0);
          double t = b.num_asm_mat_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at((run_size/3) + lvl-lvl_min).at(i).at(1);
          double t = b.num_asm_mat_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        std::cout << "\n";
      }
    }

    {
      // std::cout << "\n------------------------------------\n";
      std::cout << "\nCUDA CUDA Asm Matrix Num Timings:\n";
      String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                   + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
      std::cout << line_h << "\n";
      String line_v = String(" LVL   ");
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      std::cout << line_v << "\n";


      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ":    ";

        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at(2*(run_size/3)+lvl-lvl_min).at(i).at(0);
          double t = b.num_asm_mat_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at(2*(run_size/3) + lvl-lvl_min).at(i).at(1);
          double t = b.num_asm_mat_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        std::cout << "\n";
      }
    }
#endif // FEAT_HAVE_CUDA

    {
      std::cout << "\n------------------------------------\n";
      if(generic_run)
      {
        std::cout << "Generic Generic Solver Apply Timings:\n";
        String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                    + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
        std::cout << line_h << "\n";
        String line_v = String(" LVL   ");
        for(auto h : num_threads)
        {
          line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
        }
        for(auto h : num_threads)
        {
          line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
        }
        std::cout << line_v << "\n";


        for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
        {
          std::cout << stringify(lvl).pad_front(2) << ":    ";

          for(std::size_t i = 0; i < num_threads.size(); ++i)
          {
            const auto& b = results.at(lvl-lvl_min).at(i).at(0);
            double t = b.solver_apply_time;
            std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
          }
          for(std::size_t i = 0; i < num_threads.size(); ++i)
          {
            const auto& b = results.at(lvl-lvl_min).at(i).at(1);
            double t = b.solver_apply_time;
            std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
          }
          std::cout << "\n";
        }
      }
    }
#ifdef FEAT_HAVE_CUDA
    {
      // std::cout << "\n------------------------------------\n";
      std::cout << "\nGeneric CUDA Solver Apply Timings:\n";
      String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                   + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
      std::cout << line_h << "\n";
      String line_v = String(" LVL   ");
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      std::cout << line_v << "\n";


      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ":    ";

        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at((run_size/3)+lvl-lvl_min).at(i).at(0);
          double t = b.solver_apply_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at((run_size/3) + lvl-lvl_min).at(i).at(1);
          double t = b.solver_apply_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        std::cout << "\n";
      }
    }

    {
      // std::cout << "\n------------------------------------\n";
      std::cout << "\nCUDA CUDA Solver Apply Timings:\n";
      String line_h = String("       ") + String("Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*20)
                   + String("Voxel-Vanka").pad_front(num_threads.size()*10).pad_back(num_threads.size()*10);
      std::cout << line_h << "\n";
      String line_v = String(" LVL   ");
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      for(auto h : num_threads)
      {
        line_v += String(stringify(h) + "threads").pad_front(10).pad_back(20);
      }
      std::cout << line_v << "\n";


      for(Index lvl(lvl_min); lvl <= lvl_max; ++lvl)
      {
        std::cout << stringify(lvl).pad_front(2) << ":    ";

        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at(2*(run_size/3)+lvl-lvl_min).at(i).at(0);
          double t = b.solver_apply_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        for(std::size_t i = 0; i < num_threads.size(); ++i)
        {
          const auto& b = results.at(2*(run_size/3) + lvl-lvl_min).at(i).at(1);
          double t = b.solver_apply_time;
          std::cout << stringify_fp_fix(t, 6, 10).pad_front(10).pad_back(20);
        }
        std::cout << "\n";
      }
    }
#endif // FEAT_HAVE_CUDA
  }

  void main(int argc, char** argv)
  {
    SimpleArgParser args(argc, argv);

    args.support("level");
    args.support("threads");
    args.support("dim");
    args.support("no-generic-apply");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if( !unsupported.empty() )
    {
      // print all unsupported options to cerr
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option '--" << (*it).second << "'\n";
      return;
    }
    //parse dimension
    int dim(2);
    args.parse("dim", dim);


    if(dim == 2)
    {
      run<Q2P1, double, Index, 2>(args);
    }
    else if(dim == 3)
    {
      run<Q2P1, double, Index, 3>(args);
    }
    else
    {
      std::cerr << "ERROR: No valid dimension given!\n";
    }
  }

}

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  SaddlePointAssemblyBench::main(argc, argv);
  return 0;
}
