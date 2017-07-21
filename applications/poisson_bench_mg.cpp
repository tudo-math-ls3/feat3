//
// Geometric Multigrid on Recursive Partitions Benchmarking Application
// --------------------------------------------------------------------
// This is an application is used as a testing ground and for benchmarking
// our proof-of-concept prototype implementation of geometric multigrid
// solvers on recursive/hierarchic partitions in the context of next-gen
// high-performance <insert-fancy-buzzword-here> EXA-scale computing.
//
// WARNING:
// Do NOT use this application as a base for your own applications unless you have
// some seriously masochistic tendencies and/or you really know what you're doing.
//
// \author Peter Zajac
//
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/umfpack.hpp>
#include <kernel/util/dist.hpp>

#include <control/domain/unit_cube_domain_control.hpp>
#include <control/scalar_basic.hpp>
#include <control/statistics.hpp>

#include <vector>

namespace PoissonDirichlet
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
    std::array<unsigned long, 2> mem_use, mem_use_sum, mem_use_max;

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

  template<typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our arch types
    typedef Mem::Main MemType;
    typedef double DataType;
    typedef Index IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef typename DomainLevelType::SpaceType SpaceType;

    // fetch our mesh and shape types
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainControlType::ShapeType ShapeType;

    // choose our desired analytical solution
    Analytic::Common::ExpBubbleFunction<ShapeType::dimension> sol_func;

    // define our system level
    typedef Control::ScalarUnitFilterSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = domain.size_physical();

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

#ifdef FEAT_HAVE_UMFPACK
    const bool umf_cgs = (domain.back_layer().comm().size() == 1);
#else
    const bool umf_cgs = false;
#endif

    Index iters = 3; // mg iterations
    Index steps = 4; // smoothing steps
    args.parse("iters", iters);
    args.parse("steps", steps);

    Solver::MultiGridCycle cycle(Solver::MultiGridCycle::V);
    args.parse("cycle", cycle);

    comm.print(String("Space........: ") + SpaceType::name());
    comm.print(String("Coarse Solver: ") + (umf_cgs ? "UMFPACK" : "Jacobi"));
    comm.print(String("Smooth Steps.: ") + stringify(steps));
    comm.print(String("MG Iterations: ") + stringify(iters));
    comm.print(String("MG Cycle.....: ") + stringify(cycle));

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("gauss-legendre:" + stringify(SpaceType::local_degree+1));

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system_levels.at(i)->assemble_gate(domain.at(i));
      stats.times[i][Times::asm_gate] += ts.elapsed_now();
      if((i+1) < domain.size_virtual())
      {
        TimeStamp ts2;
        system_levels.at(i)->assemble_coarse_muxer(domain.at(i+1));
        stats.times[i][Times::asm_muxer] += ts2.elapsed_now();
        TimeStamp ts3;
        system_levels.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);
        stats.times[i][Times::asm_transfer] += ts3.elapsed_now();
      }
      stats.times[i][Times::asm_total] += ts.elapsed_now();
    }

    /* ***************************************************************************************** */

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system_levels.at(i)->assemble_laplace_matrix(domain.at(i)->space, cubature);
      double tt = ts.elapsed_now();
      stats.times[i][Times::asm_total] += tt;
      stats.times[i][Times::asm_matrix] += tt;
    }

    /* ***************************************************************************************** */

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system_levels.at(i)->assemble_homogeneous_unit_filter(*domain.at(i), domain.at(i)->space);

      // apply filter
      system_levels.at(i)->filter_sys.local().filter_mat(system_levels.at(i)->matrix_sys.local());
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
        stats.counts[i][Counts::elems_mirror] += mirror.num_indices();
      }
    }

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    TimeStamp ts_rhs;

    // create new vector
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

    vec_sol.format();
    vec_rhs.format();

    {
      // get the local vector
      typename SystemLevelType::LocalSystemVector& vec_f = vec_rhs.local();

      // assemble the force
      Assembly::Common::LaplaceFunctional<decltype(sol_func)> force_func(sol_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_f, force_func, the_domain_level.space, cubature);

      // sync the vector
      vec_rhs.sync_0();
    }

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
        auto jacobi = Solver::new_jacobi_precond(lvl.matrix_sys, lvl.filter_sys, 0.7);
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, 1.0, jacobi);
        smoother->set_min_iter(steps);
        smoother->set_max_iter(steps);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
#ifdef FEAT_HAVE_UMFPACK
      else if(umf_cgs)
      {
        // create UMFPACK coarse grid solver
        auto umfpack = std::make_shared<Solver::Umfpack>(lvl.matrix_sys.local());
        auto cgsolver = Solver::new_schwarz_precond(umfpack, lvl.filter_sys);
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
    auto solver = Solver::new_richardson(the_system_level.matrix_sys, the_system_level.filter_sys, 1.0, mgv);

    // enable plotting
    if(comm.rank() == 0)
    {
      solver->set_plot_mode(Solver::PlotMode::iter);
    }

    // set tolerance
    solver->set_plot_name("Multigrid");
    solver->set_min_iter(iters);
    solver->set_max_iter(iters);
    //solver->set_tol_rel(1E-8);
    //solver->set_max_iter(1000);

    // initialise
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

    // get memory info
    {
      auto meminfo = Util::get_memory_usage();
      stats.mem_use[0] = meminfo.peak_physical;
      stats.mem_use[1] = meminfo.peak_virtual;
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

  template<template<typename> class TSpace_>
  void run_space(Dist::Comm& comm, SimpleArgParser& args)
  {
    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef TSpace_<TrafoType> SpaceType;

    // create a time-stamp
    TimeStamp time_stamp;

    // let's create our domain

    // create our domain control
    typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
    Control::Domain::HierarchUnitCubeDomainControl2<DomainLevelType> domain(comm, args.query("level")->second);

    Statistics::toe_partition = time_stamp.elapsed_now();

    {
      std::deque<std::pair<int,int>> ilvs = domain.get_level_indices();
      auto it = ilvs.begin();
      String slvls("Levels: ");
      ((slvls += stringify(it->first)) += ":") += stringify(it->second);
      for(++it; it != ilvs.end(); ++it)
        (((slvls += ", ") += stringify(it->first)) += ":") += stringify(it->second);
      comm.print(slvls);
    }

    // dump domain info if desired
    if(args.check("dump") >= 0)
    {
      comm.allprint(domain.dump_layers());
      comm.allprint(domain.dump_layer_levels());
      comm.allprint(domain.dump_virt_levels());
    }

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }

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
    args.support("level");
    args.support("dump");
    args.support("iters");
    args.support("steps");
    args.support("cycle");
    args.support("space");

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

    int space_degree = 1;
    args.parse("space", space_degree);

    switch(space_degree)
    {
    case 1:
      run_space<Space::Lagrange1::Element>(comm, args);
      break;

    case 2:
      run_space<Space::Lagrange2::Element>(comm, args);
      break;

    case 3:
      run_space<Space::Lagrange3::Element>(comm, args);
      break;

    default:
      comm.print("ERROR: invalid space argument");
      FEAT::Runtime::abort();
      return;
    }
    comm.print(String(100u, '#'));
  }
} // namespace PoissonDirichlet

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialise(argc, argv);
  try
  {
    PoissonDirichlet::main(argc, argv);
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
  return FEAT::Runtime::finalise();
}
