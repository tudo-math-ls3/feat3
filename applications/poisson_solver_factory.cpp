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
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/precon_wrapper.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/solver/matrix_stock.hpp>
#include <kernel/solver/solver_factory.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>
#include <control/statistics.hpp>


namespace PoissonSolverFactory
{
  using namespace FEAT;

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

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename MeshType::ShapeType ShapeType;

    // define our system level
    typedef Control::ScalarUnitFilterSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = Index(domain.size_physical());

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
      Analytic::Common::ConstantFunction<ShapeType::dimension> one_func(1.0);
      Assembly::Common::ForceFunctional<decltype(one_func)> force_func(one_func);
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

    comm.print("Creating solver tree");

    ////////// MATRIX STOCK
    Solver::MatrixStock<
      typename SystemLevelType::GlobalSystemMatrix,
      typename SystemLevelType::GlobalSystemFilter,
      typename SystemLevelType::GlobalSystemTransfer> matrix_stock(domain.size_virtual());
    for (auto& system_level : system_levels)
    {
      matrix_stock.systems.push_back(system_level->matrix_sys.clone(LAFEM::CloneMode::Shallow));
      matrix_stock.gates_row.push_back(&system_level->gate_sys);
      matrix_stock.gates_col.push_back(&system_level->gate_sys);
      matrix_stock.filters.push_back(system_level->filter_sys.clone(LAFEM::CloneMode::Shallow));
      matrix_stock.muxers.push_back(&system_level->coarse_muxer_sys);
      matrix_stock.transfers.push_back(system_level->transfer_sys.clone(LAFEM::CloneMode::Shallow));
    }

    comm.print("Solving linear system...");

    String solver_ini_name;
    args.parse("solver-ini", solver_ini_name);
    PropertyMap property_map;
    property_map.parse(solver_ini_name, true);
    auto solver = Solver::SolverFactory::create_scalar_solver(matrix_stock, &property_map, "linsolver");

    matrix_stock.hierarchy_init();

    // initialise
    solver->init();

    Statistics::reset();

    TimeStamp at;

    // solve
    auto result = Solver::solve(*solver, vec_sol, vec_rhs, the_system_level.matrix_sys, the_system_level.filter_sys);

    if (!Solver::status_success(result))
    {
      comm.print("Solver execution FAILED, with status: " + stringify(result));
    }

    double solver_toe(at.elapsed_now());

    FEAT::Control::Statistics::report(solver_toe, args.check("statistics"), MeshType::ShapeType::dimension,
      system_levels, domain);

    if (args.check("statistics") >= 0)
        comm.print("\n" + FEAT::Statistics::get_formatted_solver_internals("default"));

    // release solver
    solver->done();
    matrix_stock.hierarchy_done();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./poisson-solver-factory");
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

    if (args.check("test-iter") >= 0)
    {
      int num_iter = (int)((Solver::IterativeSolver<GlobalSystemVector>*)solver.get())->get_num_iter();
      int iter_target(0);
      args.parse("test-iter", iter_target);
      if (num_iter < iter_target - 1 || num_iter > iter_target + 1)
      {
        comm.print("FAILED");
        throw InternalError(__func__, __FILE__, __LINE__, "iter count deviation! " + stringify(num_iter) + " vs " + stringify(iter_target));
      }
    }
  }

  template<typename Shape_>
  void run_shape(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
    comm.print(String("Shape-Type: ") + Shape_::name());

    // define our mesh type
    typedef Shape_ ShapeType;
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
    domain.create(mesh_reader);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("LVL-MAX: " + stringify(domain.max_level_index()) + " [" + stringify(domain.get_desired_level_max()) + "]");
    comm.print("LVL-MIN: " + stringify(domain.min_level_index()) + " [" + stringify(domain.get_desired_level_min()) + "]");

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
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
    args.support("statistics");
    args.support("test-iter");
    args.support("solver-ini");
    args.support("vtk");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unknown option '--" << (*it).second << "'" << std::endl;

      comm.print(std::cerr, "Supported options are:");
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
    if(args.check("solver-ini") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--solver-ini <ini-file>' is missing!");
      FEAT::Runtime::abort();
    }

    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    comm.print(String("Mesh File Type: ") + mesh_type);

    // run 2D or 3D ?
    if(mesh_type == "conformal:hypercube:1:1")
      run_shape<Shape::Hypercube<1>>(args, comm, mesh_reader);
    else if(mesh_type == "conformal:hypercube:2:2")
      run_shape<Shape::Hypercube<2>>(args, comm, mesh_reader);
    else if(mesh_type == "conformal:hypercube:3:3")
      run_shape<Shape::Hypercube<3>>(args, comm, mesh_reader);
    else if(mesh_type == "conformal:simplex:2:2")
      run_shape<Shape::Simplex<2>>(args, comm, mesh_reader);
    else if(mesh_type == "conformal:simplex:3:3")
      run_shape<Shape::Simplex<3>>(args, comm, mesh_reader);
    else
    {
      comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
      FEAT::Runtime::abort();
    }
  }
} // namespace PoissonSolverFactory

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialise(argc, argv);
  try
  {
    PoissonSolverFactory::main(argc, argv);
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
