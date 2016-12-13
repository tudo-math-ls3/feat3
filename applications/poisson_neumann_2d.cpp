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

#include <control/domain/unit_cube_domain_control.hpp>
#include <control/scalar_basic.hpp>

namespace PoissonNeumann2D
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

    // choose our desired analytical solution
    Analytic::Common::CosineWaveFunction<2> sol_func;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;

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

    comm.print("Assembling gates...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gate(domain.at(i));
    }

    /* ***************************************************************************************** */

    comm.print("Assembling transfers...");

    for (Index i(0); (i+1) < domain.size_virtual(); ++i)
    {
      system_levels.at(i)->assemble_coarse_muxer(domain.at(i+1));
      system_levels.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);
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
      > >();

    // scaling factor
    DataType omega = DataType(0.2);

    // push levels into MGV
    auto it_end = --system_levels.end();
    for (auto it = system_levels.begin(); it != it_end; ++it)
    {
      auto smoother = Solver::new_richardson((*it)->matrix_sys, (*it)->filter_sys, omega);
      smoother->set_max_iter(4);
      smoother->set_min_iter(4);
      multigrid_hierarchy->push_level((*it)->matrix_sys, (*it)->filter_sys, (*it)->transfer_sys, smoother, smoother, smoother);
    }

    // create coarse grid solver
    {
      auto coarse_solver = Solver::new_richardson(system_levels.back()->matrix_sys, system_levels.back()->filter_sys, omega);
      coarse_solver->set_max_iter(4);
      coarse_solver->set_min_iter(4);
      multigrid_hierarchy->push_level(system_levels.back()->matrix_sys, system_levels.back()->filter_sys, coarse_solver);
    }

    // create our solver
    auto mgv = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);
    auto solver = Solver::new_pcg(matrix, filter, mgv);

    // enable plotting
    solver->set_plot(comm.rank() == 0);

    // set tolerance
    solver->set_tol_rel(1E-8);
    solver->set_max_iter(100);

    comm.print("Solving...");

    // initialise
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
      String vtk_name = String("./poisson-neumann-2d");
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

    // initialise
#ifdef FEAT_HAVE_MPI
    comm.print("NUM-PROCS: " + stringify(comm.size()));
#endif

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
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
      // abort
      FEAT::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // parse levels
    std::deque<int> lvls;
    {
      auto p = args.query("level");
      if(p == nullptr)
      {
        lvls.push_back(5);
      }
      else
      {
        XASSERTM(!p->second.empty(), "no levels given to --level option");

        lvls.resize(p->second.size());
        for(std::size_t i(0); i < p->second.size(); ++i)
        {
          if(!p->second.at(i).parse(lvls.at(i)))
          {
            comm.print(std::cerr, "ERROR: failed to parse '" + p->second.at(i) + "' as level");
            FEAT::Runtime::abort();
          }
        }
      }
      if(lvls.size() < std::size_t(2))
        lvls.push_back(0);
    }

#ifndef DEBUG
    try
#endif
    {
      TimeStamp stamp1;

      // let's create our domain
      typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
      Control::Domain::HierarchUnitCubeDomainControl<DomainLevelType> domain(comm, lvls);

      // plot our levels
      comm.print("LVL-MAX: " + stringify(domain.max_level_index()) + " [" + stringify(lvls.front()) + "]");
      comm.print("LVL-MIN: " + stringify(domain.min_level_index()) + " [" + stringify(lvls.back()) + "]");

      domain.dump_layers();
      domain.dump_layer_levels();
      domain.dump_virt_levels();

      // run our application
      run(args, domain);

      TimeStamp stamp2;

      // get times
      long long time1 = stamp2.elapsed_micros(stamp1);

      // accumulate times over all processes
      long long time2 = time1 * (long long) comm.size();

      // print time
      comm.print("Run-Time: "
        + TimeStamp::format_micros(time1, TimeFormat::m_s_m) + " ["
        + TimeStamp::format_micros(time2, TimeFormat::m_s_m) + "]");
    }
#ifndef DEBUG
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
#endif // DEBUG
  }
} // namespace PoissonNeumann2D

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialise(argc, argv);
  PoissonNeumann2D::main(argc, argv);
  return FEAT::Runtime::finalise();
}
