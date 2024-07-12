// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/umfpack.hpp>
#include <kernel/solver/superlu.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/geometry/atlas/bezier.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/direct_stokes_solver.hpp>

#include <control/domain/voxel_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/voxel_transfer_assembler.hpp>
#include <control/statistics.hpp>

#include <vector>

namespace StokesVoxelYPipe
{
  using namespace FEAT;

  struct Counts
  {
    static constexpr std::size_t num_ranks = 0u;
    static constexpr std::size_t num_elems = 1u;
    static constexpr std::size_t num_dofs_l = 2u; // local dofs
    static constexpr std::size_t num_dofs_g = 3u; // global dofs [ != sum(num_dofs_l) ]
    static constexpr std::size_t num_nze = 4u;
    static constexpr std::size_t bytes_domain = 5u;
    static constexpr std::size_t bytes_system = 6u;
    static constexpr std::size_t elems_mirror = 7u;
    static constexpr std::size_t count = 8u;
  };

  struct Times
  {
    static constexpr std::size_t asm_total = 0u;
    static constexpr std::size_t asm_gate = 1u;
    static constexpr std::size_t asm_muxer = 2u;
    static constexpr std::size_t asm_transfer = 3u;
    static constexpr std::size_t asm_matrix = 4u;
    static constexpr std::size_t count = 5u;
  };

  template<typename T_>
  inline T_ sum(const std::vector<T_>& dv)
  {
    T_ t = T_(0);
    for(auto x : dv) t += x;
    return t;
  }

  template<typename T_, std::size_t n_>
  inline T_ sum(const std::vector<std::array<T_, n_>>& dv, std::size_t i)
  {
    XASSERT(i < n_);
    T_ t = T_(0);
    for(const auto& x : dv) t += x[i];
    return t;
  }

  struct BenchStats
  {
    // per [level][index]
    std::vector<std::array<unsigned long long, Counts::count>> counts, counts_sum, counts_max;

    // per [level][index]
    std::vector<std::array<double, Times::count>> times;

    // (physical, virtual)
    std::array<unsigned long long, 2> mem_use, mem_use_sum, mem_use_max;

    double toe_asm_rhs;

    explicit BenchStats(std::size_t virt_size) :
      counts(virt_size),
      counts_sum(virt_size),
      counts_max(virt_size),
      times(virt_size),
      toe_asm_rhs(0.0)
    {
      for(std::size_t i(0u); i < virt_size; ++i)
      {
        for(std::size_t j(0u); j < Counts::count; ++j)
          counts[i][j] = 0ull;
        for(std::size_t j(0u); j < Times::count; ++j)
          times[i][j] = 0.0;
      }
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(counts.data(), counts_sum.data(), counts.size()*Counts::count, Dist::dt_unsigned_long_long, Dist::op_sum);
      comm.allreduce(counts.data(), counts_max.data(), counts.size()*Counts::count, Dist::dt_unsigned_long_long, Dist::op_max);

      comm.allreduce(mem_use.data(), mem_use_sum.data(), 2u, Dist::op_sum);
      comm.allreduce(mem_use.data(), mem_use_max.data(), 2u, Dist::op_max);
    }

    String format() const
    {
      String s;

      s += "\nAssembly Timings:\n";
      s += "                Gate        Muxer     Transfer       Matrix        Total        RHS\n";
      s += "Overall : " +
        stringify_fp_fix(sum(times, Times::asm_gate), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_muxer), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_transfer), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_matrix), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_total), 6, 10) + " / " +
        stringify_fp_fix(toe_asm_rhs, 6) + "\n";
      for(std::size_t i(0); i < times.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(times[i][Times::asm_gate], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_muxer], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_transfer], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_matrix], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_total], 6, 10) + "\n";
      }

      s += "\nBasic Statistics:\n";
      s += "           Ranks       Elements [    per Patch ]           Dofs [    per Patch ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify(counts_max[i][Counts::num_ranks]).pad_front(6) + " / " +
          stringify(counts_sum[i][Counts::num_elems]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_elems]).pad_front(12) + " ] / " +
          stringify(counts    [i][Counts::num_dofs_g]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_dofs_l]).pad_front(12) + " ]\n";
      }

      s += "\nMirror #elems Statistics:\n";
      s += "Overall : " + stringify(sum(counts_sum, Counts::elems_mirror)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::elems_mirror)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::elems_mirror]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::elems_mirror]).pad_front(15) + " ]\n";
      }

      s += "\nMemory Usage Statistics:\n";
      s += String("Peak Physical") +
        ": " + stringify_fp_fix(double(mem_use_sum[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB " +
        "[ " + stringify_fp_fix(double(mem_use_max[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB ]\n";
      //": " + stringify(mem_use_sum[0]).pad_front(15) +
      //" [" + stringify(mem_use_max[0]).pad_front(15) + " ]\n";
      s += String("Peak Virtual.") +
        ": " + stringify_fp_fix(double(mem_use_sum[1])/(1024.0*1024.0*1024.0), 6, 15) + " GiB " +
        "[ " + stringify_fp_fix(double(mem_use_max[1])/(1024.0*1024.0*1024.0), 6, 15) + " GiB ]\n";
      //": " + stringify(mem_use_sum[1]).pad_front(15) +
      //" [" + stringify(mem_use_max[1]).pad_front(15) + " ]\n";

      s += "\nDomain Level Bytes Statistics:\n";
      s += "Overall : " + stringify(sum(counts_sum, Counts::bytes_domain)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::bytes_domain)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::bytes_domain]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::bytes_domain]).pad_front(15) + " ]\n";
      }

      s += "\nSystem Level Bytes Statistics:\n";
      s += "Overall : " + stringify(sum(counts_sum, Counts::bytes_system)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::bytes_system)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::bytes_system]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::bytes_system]).pad_front(15) + " ]\n";
      }

      return s;
    }
  };

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());
    comm.print(String(100u, '*'));

    // dump system call
    {
      String s("Arguments: ");
      s.append(argv[0]);
      for(int i(1); i < argc; ++i)
        s.append(" ").append(argv[i]);
      comm.print(s);
    }

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    args.support("vtk");
    args.support("level");
    args.support("dump");
    args.support("iters");
    args.support("steps");
    args.support("cycle");
    args.support("space");
    args.support("debug");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());

      // abort
      FEAT::Runtime::abort();
    }

    // define our types
    static constexpr int shape_dim = 2;
    typedef Shape::Hypercube<shape_dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceTypeVelo;
    typedef Space::Discontinuous::ElementP1<TrafoType> SpaceTypePres;

    typedef Tiny::Vector<Real, shape_dim> WorldPoint;

    // create Y-pipe chart
    Geometry::Atlas::Bezier<MeshType> ypipe_chart;
    ypipe_chart.push_vertex(WorldPoint{-1.0, 0.0});
    ypipe_chart.push_vertex(WorldPoint{-1.0, -0.4});
    ypipe_chart.push_vertex(WorldPoint{0.35, -0.45});
    ypipe_chart.push_control(WorldPoint{1.25, -0.45});
    ypipe_chart.push_control(WorldPoint{1.25, -0.9});
    ypipe_chart.push_vertex(WorldPoint{2.15, -0.9});
    ypipe_chart.push_vertex(WorldPoint{4.0, -0.9});
    ypipe_chart.push_vertex(WorldPoint{4.0, -0.45});
    ypipe_chart.push_vertex(WorldPoint{2.15, -0.45});
    ypipe_chart.push_control(WorldPoint{1.7, -0.45});
    ypipe_chart.push_control(WorldPoint{1.25, -0.18});
    ypipe_chart.push_vertex(WorldPoint{1.25, 0.0});
    ypipe_chart.push_control(WorldPoint{1.25, 0.18});
    ypipe_chart.push_control(WorldPoint{1.7, 0.45});
    ypipe_chart.push_vertex(WorldPoint{2.15, 0.45});
    ypipe_chart.push_vertex(WorldPoint{4.0, 0.45});
    ypipe_chart.push_vertex(WorldPoint{4.0, 0.9});
    ypipe_chart.push_vertex(WorldPoint{2.15, 0.9});
    ypipe_chart.push_control(WorldPoint{1.25, 0.9});
    ypipe_chart.push_control(WorldPoint{1.25, 0.45});
    ypipe_chart.push_vertex(WorldPoint{0.35, 0.45});
    ypipe_chart.push_vertex(WorldPoint{-1.0, 0.45});
    ypipe_chart.push_close();


    TimeStamp time_stamp;

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceTypeVelo, SpaceTypePres> DomainLevelType;
    typedef Control::Domain::VoxelDomainLevelWrapper<DomainLevelType> VoxelDomainLevelType;
    Control::Domain::VoxelDomainControl<VoxelDomainLevelType> domain(comm, true);

    domain.create_base_mesh_2d(3, 2, 0.0, 3.0, -1.0, 1.0);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);

    domain.create_voxel_map_from_chart(ypipe_chart, true, 0.0);
    domain.create_hierarchy();

    auto hit_l = [](auto p) {return p[0] < 0.0001;};
    auto hit_r = [](auto p) {return p[0] > 2.9999;};

    // create boundary meshparts on each level
    for(Index i(0); i < domain.size_physical(); ++i)
    {
      Geometry::RootMeshNode<MeshType>& mesh_node = *domain.at(i)->get_mesh_node();
      Geometry::HitTestFactory<decltype(hit_l), MeshType> bnd_factory_l(hit_l, *mesh_node.get_mesh());
      Geometry::HitTestFactory<decltype(hit_r), MeshType> bnd_factory_r(hit_r, *mesh_node.get_mesh());

      auto bnd_l = bnd_factory_l.make_unique();
      auto bnd_r = bnd_factory_r.make_unique();

      Geometry::GlobalMaskedBoundaryFactory<MeshType> boundary_factory(*mesh_node.get_mesh());
      boundary_factory.add_mask_meshpart(*bnd_l);
      boundary_factory.add_mask_meshpart(*bnd_r);
      for(const auto& v : mesh_node.get_halo_map())
        boundary_factory.add_halo(v.first, *v.second);
      boundary_factory.compile(domain.at(i).layer().comm());
      mesh_node.add_mesh_part("bnd:w", boundary_factory.make_unique());
      mesh_node.add_mesh_part("bnd:l", std::move(bnd_l));
      mesh_node.add_mesh_part("bnd:r", std::move(bnd_r));
    }

    Statistics::toe_partition = time_stamp.elapsed_now();

    comm.print("\nDesired Levels: " + domain.format_desired_levels());
    comm.print(  "Chosen  Levels: " + domain.format_chosen_levels());

    comm.print("\nChosen Partitioning Info:\n" + domain.get_chosen_parti_info());

    comm.print("Base-Mesh Creation Time: " + domain.get_watch_base_mesh().elapsed_string().pad_front(7));
    comm.print("Voxel-Map Creation Time: " + domain.get_watch_voxel_map().elapsed_string().pad_front(7));
    comm.print("Hierarchy Creation Time: " + domain.get_watch_hierarchy().elapsed_string().pad_front(7));

    // dump domain info if desired
    if(args.check("dump") >= 0)
    {
      comm.print("\nDomain Ancestry:");
      comm.allprint(domain.dump_ancestry());
      comm.print("\nDomain Layers:");
      comm.allprint(domain.dump_layers());
      comm.print("\nDomain Layer Levels:");
      comm.allprint(domain.dump_layer_levels());
      comm.print("\nDomain Virtual Levels:");
      comm.allprint(domain.dump_virt_levels());
    }

    // define our arch types
    typedef double DataType;
    typedef Index IndexType;

    // define our system level
    typedef Control::StokesBlockedCombinedSystemLevel<shape_dim, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = domain.size_physical();

    // enable solver expressions
    Statistics::enable_solver_expressions = true;

    BenchStats stats(domain.size_virtual());

    {
      std::size_t n = domain.size_virtual() > domain.size_physical() ? domain.size_physical()+1 : domain.size_physical();
      for(std::size_t i(0); i < n; ++i)
      {
        const auto& vdl = domain.at(i);
        if(vdl.is_parent())
          stats.counts[i][Counts::bytes_domain] = vdl.level_c().bytes() + vdl.level_p().bytes();
        else if(vdl.is_child())
          stats.counts[i][Counts::bytes_domain] = vdl.level_c().bytes();
        else
          stats.counts[i][Counts::bytes_domain] = vdl.level().bytes();
      }
    }

    Index iters = 3; // mg iterations
    Index steps = 8; // smoothing steps
    args.parse("iters", iters);
    args.parse("steps", steps);

    Solver::MultiGridCycle cycle(Solver::MultiGridCycle::V);
    args.parse("cycle", cycle);

    comm.print(String("Smooth Steps.: ") + stringify(steps));
    //comm.print(String("MG Iterations: ") + stringify(iters));
    comm.print(String("MG Cycle.....: ") + stringify(cycle));

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    const String cubature("gauss-legendre:5");// + stringify(SpaceTypeVelo::local_degree+2));

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      domain.at(i)->domain_asm.compile_all_elements();
      system_levels.at(i)->assemble_gates(domain.at(i));
      stats.times[i][Times::asm_gate] += ts.elapsed_now();
      if((i+1) < domain.size_virtual())
      {
        TimeStamp ts2;
        system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
        stats.times[i][Times::asm_muxer] += ts2.elapsed_now();
      }
      stats.times[i][Times::asm_total] += ts.elapsed_now();
    }

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts3;
      if((i+1) < num_levels)
        Control::VoxelTransferAssembler::assemble_stokes_blocked_transfers(*system_levels.at(i), *system_levels.at(i+1), domain.at(i), domain.at(i+1), cubature, true, true, true);
      else if((i+1) < domain.size_virtual())
        Control::VoxelTransferAssembler::assemble_stokes_blocked_transfers(*system_levels.at(i), domain.at(i), domain.at(i+1), cubature, true, true, true);
      stats.times[i][Times::asm_transfer] += ts3.elapsed_now();
    }

    /* ***************************************************************************************** */

    Assembly::Common::LaplaceOperatorBlocked<shape_dim> laplace_op;

    std::deque<typename SystemLevelType::LocalSystemMatrix> local_matrices(num_levels);

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
      system_levels.at(i)->assemble_pres_struct(domain.at(i)->space_pres);
      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->domain_asm, domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
      Assembly::assemble_bilinear_operator_matrix_1(domain.at(i)->domain_asm, system_levels.at(i)->matrix_a.local(), laplace_op, domain.at(i)->space_velo, cubature);
      system_levels.at(i)->compile_system_matrix();

      double tt = ts.elapsed_now();
      stats.times[i][Times::asm_total] += tt;
      stats.times[i][Times::asm_matrix] += tt;
    }

    /* ***************************************************************************************** */

    Analytic::Common::ParProfileVector<DataType> inflow_profile(0.0, -0.45, 0.0, 0.45);

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;

      Assembly::UnitFilterAssembler<MeshType> unit_asm_l, unit_asm_w;

      auto& filter_v_l = system_levels.at(i)->filter_velo.local().template at<1>().find_or_add("bnd:l");
      auto& filter_v_w = system_levels.at(i)->filter_velo.local().template at<1>().find_or_add("bnd:w");

      // get inflow meshpart
      auto* part_l = domain.at(i)->get_mesh_node()->find_mesh_part("bnd:l");
      if(part_l)
        unit_asm_l.add_mesh_part(*part_l);
      auto* part_w = domain.at(i)->get_mesh_node()->find_mesh_part("bnd:w");
      if(part_w)
        unit_asm_w.add_mesh_part(*part_w);

      // assemble filter
      unit_asm_l.assemble(filter_v_l, domain.at(i)->space_velo, inflow_profile);
      unit_asm_w.assemble(filter_v_w, domain.at(i)->space_velo);

      system_levels.at(i)->compile_system_filter();

      stats.times[i][Times::asm_total] += ts.elapsed_now();

      stats.counts[i][Counts::num_ranks] = Index(domain.at(i).layer().comm().size());
      stats.counts[i][Counts::num_elems] = domain.at(i)->get_mesh().get_num_elements();
      stats.counts[i][Counts::num_dofs_g] = system_levels.at(i)->matrix_sys.rows();
      stats.counts[i][Counts::num_dofs_l] = system_levels.at(i)->matrix_sys.local().rows();
      stats.counts[i][Counts::num_nze] = system_levels.at(i)->matrix_sys.local().used_elements();
      stats.counts[i][Counts::bytes_system] = system_levels.at(i)->bytes();
      stats.counts[i][Counts::elems_mirror] = 0;
      auto& sys_mirrors =  system_levels.at(i)->gate_sys._mirrors;
      for (auto& mirror : sys_mirrors)
      {
        stats.counts[i][Counts::elems_mirror] += mirror.template at<0>().num_indices()
          + mirror.template at<1>().num_indices();
      }
    }


    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;

      local_matrices.at(i).block_a() = system_levels.at(i)->matrix_a.convert_to_1();
      local_matrices.at(i).block_b() = system_levels.at(i)->matrix_b.local().clone(LAFEM::CloneMode::Weak);
      local_matrices.at(i).block_d() = system_levels.at(i)->matrix_d.local().clone(LAFEM::CloneMode::Weak);

      for(const auto& fv : system_levels.at(i)->filter_velo.local().at<1>())
      {
        fv.second.filter_mat(local_matrices.at(i).block_a());
        fv.second.filter_offdiag_row_mat(local_matrices.at(i).block_b());
      }

      double tt = ts.elapsed_now();
      stats.times[i][Times::asm_total] += tt;
      stats.times[i][Times::asm_matrix] += tt;
    }

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

    // fetch our finest levels
    //DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    TimeStamp ts_rhs;

    // create new vector
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

    vec_sol.format();
    vec_rhs.format();

    // and filter it
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_rhs);

    stats.toe_asm_rhs = ts_rhs.elapsed_now();

    Statistics::toe_assembly = stamp_ass.elapsed_now();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    auto multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalSystemMatrix,
      typename SystemLevelType::GlobalSystemFilter,
      typename SystemLevelType::GlobalSystemTransfer
      > >(domain.size_virtual());

    // push all levels except the coarse most one
    for (Index i(0); i < num_levels; ++i)
    {
      const SystemLevelType& lvl = *system_levels.at(i);

      if((i+1) < domain.size_virtual())
      {
        auto vanka = Solver::new_amavanka(local_matrices.at(i), lvl.filter_sys.local());
        auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
        //auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, 0.5, schwarz);
        auto smoother = Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, steps, 0.0, schwarz);
        smoother->set_min_iter(steps);
        smoother->set_max_iter(steps);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
#ifdef FEAT_HAVE_UMFPACK
      else if(domain.back_layer().comm().size() == 1)
      {
        auto cgsolver = Solver::new_direct_stokes_solver(lvl.matrix_sys, lvl.filter_sys);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
#endif //  FEAT_HAVE_UMFPACK
      else
      {
        auto cgsolver = Solver::new_jacobi_precond(lvl.matrix_sys, lvl.filter_sys, 1.0);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
    }

    auto mgv = Solver::new_multigrid(multigrid_hierarchy, cycle);
    //mgv->set_adapt_cgc(Solver::MultiGridAdaptCGC::MinEnergy);
    auto solver = Solver::new_richardson(the_system_level.matrix_sys, the_system_level.filter_sys, 1.0, mgv);
    //auto solver = Solver::new_pcg(the_system_level.matrix_sys, the_system_level.filter_sys, mgv);

    // enable plotting
    solver->set_plot_mode(Solver::PlotMode::iter);

    // set tolerance
    solver->set_plot_name("Multigrid");
    //solver->set_min_iter(iters);
    //solver->set_max_iter(50);
    solver->set_tol_rel(1E-8);
    solver->set_max_iter(1000);

    // initialize
    multigrid_hierarchy->init();
    solver->init();

    Statistics::reset();

    comm.print("");

    TimeStamp at;
    auto result = Solver::solve(*solver, vec_sol, vec_rhs, the_system_level.matrix_sys, the_system_level.filter_sys);

    const double solver_toe(at.elapsed_now());

    if (!Solver::status_success(result))
    {
      comm.print("Solver execution FAILED, with status: " + stringify(result));
    }

    // release solver
    solver->done();
    multigrid_hierarchy->done();

    if (args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./stokes-voxel-ypipe");
      vtk_name += "-lvl" + stringify(domain.front()->get_level_index());
      vtk_name += "-n" + stringify(comm.size());

      auto node = domain.front()->get_mesh_node()->refine_unique();

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(*node->get_mesh());

      // write velocity
      exporter.add_vertex_vector("velo", vec_sol.local().template at<0>());

      // loop over all mesh parts
      auto part_names = node->get_mesh_part_names(true);
      std::vector<double> vtx_data(node->get_mesh()->get_num_entities(0), 0.0);
      for(auto it = part_names.begin(); it != part_names.end(); ++it)
      {
        // get the mesh part
        auto* part = node->find_mesh_part(*it);
        if(part == nullptr)
          continue;

        // get the vertex target set
        Geometry::TargetSet& trg = part->template get_target_set<0>();

        // mark all indexes vertices
        for(Index i(0); i < trg.get_num_entities(); ++i)
          vtx_data[trg[i]] = 1.0;

        // add variable
        exporter.add_vertex_scalar(*it, vtx_data.data());

        // unmark vertices
        for(Index i(0); i < trg.get_num_entities(); ++i)
          vtx_data[trg[i]] = 0.0;
      }


      // finally, write the VTK file
      exporter.write(vtk_name, comm);
    }

    // get memory info
    {
      MemoryUsage meminfo;
      stats.mem_use[0] = meminfo.get_peak_physical();
      stats.mem_use[1] = meminfo.get_peak_virtual();
    }

    stats.sync(comm);

    comm.print("\nMultigrid Timings:");
    comm.print("                Defect /     Smoother /     Transfer /       Coarse /        Total");
    comm.print("Overall : " +
      stringify_fp_fix(multigrid_hierarchy->get_time_defect(), 6, 12) + " / " +
      stringify_fp_fix(multigrid_hierarchy->get_time_smooth(), 6, 12) + " / " +
      stringify_fp_fix(multigrid_hierarchy->get_time_transfer(), 6, 12) + " / " +
      stringify_fp_fix(multigrid_hierarchy->get_time_coarse(), 6, 12) + " / " +
      stringify_fp_fix(multigrid_hierarchy->get_time_defect()+multigrid_hierarchy->get_time_smooth()
        +multigrid_hierarchy->get_time_transfer()+multigrid_hierarchy->get_time_coarse(), 6, 12));
    for(int i(0); i < int(multigrid_hierarchy->size_physical()); ++i)
    {
      comm.print("Level " + stringify(i).pad_front(2) + ": " +
        stringify_fp_fix(multigrid_hierarchy->get_time_defect(i), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_smooth(i), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_transfer(i), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_coarse(i), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_defect(i)+multigrid_hierarchy->get_time_smooth(i)
          +multigrid_hierarchy->get_time_transfer(i)+multigrid_hierarchy->get_time_coarse(i), 6, 12));
    }

    comm.print(stats.format());

    FEAT::Control::Statistics::report(solver_toe, 0, MeshType::ShapeType::dimension, system_levels, domain);

    comm.print(FEAT::Statistics::get_formatted_solver_internals());
    comm.print(FEAT::Statistics::get_formatted_solver_tree().trim());
  }
} // namespace StokesVoxelYPipe

int main(int argc, char* argv [])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  try
  {
    StokesVoxelYPipe::main(argc, argv);
  }
  catch (const std::exception& exc)
  {
    std::cerr << "ERROR: unhandled exception: " << exc.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return 0;
}
