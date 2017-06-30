#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/solver/legacy_preconditioners.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/util/dist.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>
#include <control/statistics.hpp>

namespace PoissonScaRC
{
  using namespace FEAT;

  template<typename MemType_, typename DataType_, typename IndexType_>
  class PoissonScarcSystemLevel :
    public Control::ScalarUnitFilterSystemLevel<MemType_, DataType_, IndexType_>
  {
  public:
    typedef Control::ScalarUnitFilterSystemLevel<MemType_, DataType_, IndexType_> BaseClass;

    /// the local system matrix
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    /// the local system filter
    typename BaseClass::LocalSystemFilter local_filter_sys;

    /// the local transfer operator
    typename BaseClass::LocalSystemTransfer local_transfer_sys;

    void compile_local_transfer(const PoissonScarcSystemLevel& coarse_level)
    {
      // clone transfer matrix
      auto locmat = this->transfer_sys.local().get_mat_prol().clone(LAFEM::CloneMode::Weak);

      Global::synch_matrix(locmat, *this->gate_sys._comm, this->gate_sys._ranks, this->gate_sys._mirrors, coarse_level.gate_sys._mirrors);

      this->local_transfer_sys.get_mat_prol() = locmat.clone(LAFEM::CloneMode::Shallow);
      this->local_transfer_sys.get_mat_rest() = locmat.transpose();
    }

    void compile_local_system()
    {
      // convert system matrix
      this->local_matrix_sys = this->matrix_sys.convert_to_1();

      // clone local filter
      this->local_filter_sys = this->filter_sys.local().clone();

      // filter local system matrix
      this->local_filter_sys.filter_mat(this->local_matrix_sys);
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

    // fetch our mesh and shape types
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainControlType::ShapeType ShapeType;

    // choose our desired analytical solution
    Analytic::Common::ExpBubbleFunction<ShapeType::dimension> sol_func;

    // define our system level
    typedef PoissonScarcSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = domain.size_physical();

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("auto-degree:5");

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    comm.print("Assembling gates, muxers and transfers...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gate(domain.at(i));
      if((i+1) < domain.size_virtual())
      {
        system_levels.at(i)->assemble_coarse_muxer(domain.at(i+1));
        system_levels.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);
      }
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system matrices...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_laplace_matrix(domain.at(i)->space, cubature);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_homogeneous_unit_filter(*domain.at(i), domain.at(i)->space);
    }

    /* ***************************************************************************************** */

    comm.print("Compiling local systems...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->compile_local_system();
      if((i+1) < domain.size_virtual())
      {
        system_levels.at(i)->compile_local_transfer(*system_levels.at(i+1));
      }
    }

    Statistics::toe_assembly = stamp_ass.elapsed_now();

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

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

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("Setting up ScaRC solver...");

    // the local multigrid hierarchy
    auto multigrid_hierarchy_local = std::make_shared<
      Solver::MultiGridHierarchy<
      typename SystemLevelType::LocalSystemMatrix,
      typename SystemLevelType::LocalSystemFilter,
      typename SystemLevelType::LocalSystemTransfer
        > >(domain.size_virtual());

    // the global multigrid hierarchy
    auto multigrid_hierarchy_global = std::make_shared<
      Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalSystemMatrix,
      typename SystemLevelType::GlobalSystemFilter,
      typename SystemLevelType::GlobalSystemTransfer
        > >(domain.size_virtual());

    // setup local multigrid hierarchy
    for (Index i(0); i < num_levels; ++i)
    {
      const SystemLevelType& lvl = *system_levels.at(i);

      if((i+1) < domain.size_virtual())
      {
        // create a smoother
        auto jacobi = Solver::new_jacobi_precond(lvl.local_matrix_sys, lvl.local_filter_sys, 0.7);
        auto smoother = Solver::new_richardson(lvl.local_matrix_sys, lvl.local_filter_sys, 1.0, jacobi);
        smoother->set_min_iter(4);
        smoother->set_max_iter(4);
        multigrid_hierarchy_local->push_level(lvl.local_matrix_sys, lvl.local_filter_sys, lvl.local_transfer_sys, smoother, smoother, smoother);
      }
      else
      {
        // create a local coarse grid solver
        auto cgsolver = Solver::new_pcg(lvl.local_matrix_sys, lvl.local_filter_sys);
        multigrid_hierarchy_local->push_level(lvl.local_matrix_sys, lvl.local_filter_sys, cgsolver);
      }
    }

    // setup global multigrid hierarchy
    for (Index i(0); i < num_levels; ++i)
    {
      const SystemLevelType& lvl = *system_levels.at(i);

      if((i+1) < domain.size_virtual())
      {
        // create a local multigrid solver for this level
        auto local_mg = Solver::new_multigrid(multigrid_hierarchy_local, Solver::MultiGridCycle::V, int(i));

        // set adaptive coarse grid correction
        local_mg->set_adapt_cgc(Solver::MultiGridAdaptCGC::MinEnergy);

        // put it into a schwarz preconditioner
        auto schwarz = Solver::new_schwarz_precond(local_mg, lvl.filter_sys);

        // and create a Richardson for the Schwarz
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, 1.0, schwarz);
        smoother->set_min_iter(4);
        smoother->set_max_iter(4);
        multigrid_hierarchy_global->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
      else
      {
        // create a global coarse grid solver
        auto cgsolver = Solver::new_pcg(lvl.matrix_sys, lvl.filter_sys);
        multigrid_hierarchy_global->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
    }

    // create a global multigrid
    auto global_mg = Solver::new_multigrid(multigrid_hierarchy_global, Solver::MultiGridCycle::V);
    auto solver = Solver::new_richardson(the_system_level.matrix_sys, the_system_level.filter_sys, 1.0, global_mg);

    // enable plotting
    if(comm.rank() == 0)
    {
      solver->set_plot_name("ScaRC");
      solver->set_plot_mode(Solver::PlotMode::iter);
    }

    // set tolerance
    solver->set_tol_rel(1E-8);
    solver->set_max_iter(100);

    // initialise
    multigrid_hierarchy_local->init();
    multigrid_hierarchy_global->init();
    solver->init();

    Statistics::reset();

    TimeStamp at;

    // solve
    auto result = Solver::solve(*solver, vec_sol, vec_rhs, the_system_level.matrix_sys, the_system_level.filter_sys);

    if (!Solver::status_success(result))
    {
      comm.print("Solver execution FAILED, with status: " + stringify(result));
    }

    const double solver_toe(at.elapsed_now());


    // release solver
    solver->done();
    multigrid_hierarchy_global->done();
    multigrid_hierarchy_local->done();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("no-err") < 0)
    {
      // Compute and print the H0-/H1-errors
      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute
        (vec_sol.local(), sol_func, the_domain_level.space, cubature);

      // synchronise all local errors
      errors.synchronise(comm);

      // print errors
      comm.print("");
      comm.print(errors.format_string());
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./poisson-scarc");
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(comm.size());

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

      // project velocity and pressure
      typename SystemLevelType::LocalSystemVector vtx_sol, vtx_rhs;
      Assembly::DiscreteVertexProjector::project(vtx_sol, vec_sol.local(), the_domain_level.space);
      Assembly::DiscreteVertexProjector::project(vtx_rhs, vec_rhs.local(), the_domain_level.space);

      // write velocity
      exporter.add_vertex_scalar("sol", vtx_sol.elements());
      exporter.add_vertex_scalar("rhs", vtx_rhs.elements());

      // finally, write the VTK file
      exporter.write(vtk_name, comm);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if(args.check("statistics") >= 0)
    {
      comm.print("\nGlobal Multigrid Timings:");
      comm.print("              Defect /   Smoother /   Transfer /     Coarse");
      comm.print("Overall : " +
          stringify_fp_fix(multigrid_hierarchy_global->get_time_defect(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_global->get_time_smooth(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_global->get_time_transfer(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_global->get_time_coarse(), 3, 10));
      for(int i(0); i < int(multigrid_hierarchy_global->size_physical()); ++i)
      {
        comm.print("Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(multigrid_hierarchy_global->get_time_defect(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_global->get_time_smooth(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_global->get_time_transfer(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_global->get_time_coarse(i), 3, 10));
      }
      comm.print("\nLocal Multigrid Timings:");
      comm.print("              Defect /   Smoother /   Transfer /     Coarse");
      comm.print("Overall : " +
          stringify_fp_fix(multigrid_hierarchy_local->get_time_defect(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_local->get_time_smooth(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_local->get_time_transfer(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_local->get_time_coarse(), 3, 10));
      for(int i(0); i < int(multigrid_hierarchy_local->size_physical()); ++i)
      {
        comm.print("Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(multigrid_hierarchy_local->get_time_defect(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_local->get_time_smooth(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_local->get_time_transfer(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy_local->get_time_coarse(i), 3, 10));
      }

      FEAT::Control::Statistics::report(solver_toe, 0, MeshType::ShapeType::dimension, system_levels, domain);

      comm.print("\n");
      comm.print(FEAT::Statistics::get_formatted_flops(solver_toe, comm.size()));
      comm.print(FEAT::Statistics::get_formatted_times(solver_toe));
      comm.print(FEAT::Statistics::get_formatted_solver_internals("default"));
      comm.print("\n");
      comm.print(FEAT::Statistics::get_formatted_solver_tree("default").trim());
    }
  }

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

#ifdef FEAT_HAVE_MPI
    comm.print("NUM-PROCS: " + stringify(comm.size()));
#endif

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("mesh");
    args.support("level");
    args.support("no-err");
    args.support("vtk");
    args.support("statistics");
    args.support("dump");

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

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
      FEAT::Runtime::abort();
    }
    if(args.check("level") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
      FEAT::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // create a time-stamp
    TimeStamp time_stamp;

    // let's create our domain
    comm.print("Preparing domain...");

    // create our domain control
    typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, false);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(args.query("mesh")->second);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("LVL-MAX: " + stringify(domain.max_level_index()) + " [" + stringify(domain.get_desired_level_max()) + "]");
    //comm.print("LVL-MED: " + stringify(domain.med_level_index()) + " [" + stringify(domain.get_desired_level_med()) + "]");
    comm.print("LVL-MIN: " + stringify(domain.min_level_index()) + " [" + stringify(domain.get_desired_level_min()) + "]");

    // dump domain info if desired
    if(args.check("dump") >= 0)
    {
      domain.dump_layers();
      domain.dump_layer_levels();
      domain.dump_virt_levels();
    }

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }
} // namespace PoissonScaRC

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialise(argc, argv);
  try
  {
    PoissonScaRC::main(argc, argv);
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
