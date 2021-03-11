// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/solver/legacy_preconditioners.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>

namespace PoissonNeumann
{
  using namespace FEAT;

  template<typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our arch types
    typedef Mem::Main MemType;
    typedef Real DataType;
    typedef Index IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh and shape types
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainControlType::ShapeType ShapeType;

    // choose our desired analytical solution
    Analytic::Common::CosineWaveFunction<ShapeType::dimension> sol_func;

    // define our system level
    typedef Control::ScalarMeanFilterSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = domain.size_physical();

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("auto-degree:5");

    /* ***************************************************************************************** */

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
      system_levels.at(i)->assemble_mean_filter(domain.at(i)->space, cubature);
    }

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // fetch matrix and filter
    GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
    GlobalSystemFilter& filter = the_system_level.filter_sys;

    // create new vectors
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

    comm.print("Setting up solver...");

    // create a multigrid cycle
    auto multigrid_hierarchy = std::make_shared<Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalSystemMatrix,
      typename SystemLevelType::GlobalSystemFilter,
      typename SystemLevelType::GlobalSystemTransfer
      > >(domain.size_virtual());

    // scaling factor
    DataType omega = DataType(0.2);

    // push levels into MGV
    for(std::size_t i(0); i < system_levels.size(); ++i)
    {
      // get a reference to the system level
      const SystemLevelType& lvl = *system_levels.at(i);

      // Is this the virtual coarse level?
      if((i+1) < domain.size_virtual())
      {
        // No, so add a smoothing level
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, omega);
        smoother->set_max_iter(4);
        smoother->set_min_iter(4);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
      else
      {
        // Yes, that's the virtual coarse level
        auto coarse_solver = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, omega);
        coarse_solver->set_max_iter(4);
        coarse_solver->set_min_iter(4);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, coarse_solver);
      }
    }

    // create our solver
    auto mgv = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

    auto solver = Solver::new_pcg(matrix, filter, mgv);

    // enable plotting
    solver->set_plot_mode(Solver::PlotMode::iter);

    // set tolerance
    solver->set_tol_rel(1E-8);
    solver->set_max_iter(100);

    comm.print("Solving...");

    // initialize
    multigrid_hierarchy->init();
    solver->init();

    // solve
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // release solver
    solver->done();
    multigrid_hierarchy->done();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("no-err") < 0)
    {
      // Compute and print the H0-/H1-errors
      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute
        (vec_sol.local(), sol_func, the_domain_level.space, cubature);

      // synhronise all local errors
      errors.synchronize(comm);

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
      String vtk_name = String("./poisson-neumann");
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
  }

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("mesh");
    args.support("level");
    args.support("no-err");
    args.support("vtk");

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
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(args.query("mesh")->second);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }
} // namespace PoissonNeumann

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialize(argc, argv);
  try
  {
    PoissonNeumann::main(argc, argv);
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
  return FEAT::Runtime::finalize();
}
