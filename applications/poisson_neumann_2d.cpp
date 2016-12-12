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
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/precon_wrapper.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>

#include <control/domain/unit_cube_domain_control.hpp>
#include <control/scalar_basic.hpp>

namespace PoissonNeumann2D
{
  using namespace FEAT;

  template<typename MeshType_, typename TargetMatrixSolve_>
  void run(const Dist::Comm& comm, SimpleArgParser& args, Control::Domain::DomainControl<MeshType_>& domain)
  {
    // define our mesh type
    typedef MeshType_ MeshType;

    // define our arch types
    typedef typename Mem::Main MemType;
    typedef typename TargetMatrixSolve_::DataType DataType;
    typedef typename TargetMatrixSolve_::IndexType IndexType;

    // choose our desired analytical solution
    Analytic::Common::CosineWaveFunction<2> sol_func;

    // define our domain type
    typedef Control::Domain::DomainControl<MeshType_> DomainControlType;

    // define our system level
    typedef Control::ScalarMeanFilterSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    // define our trafo and FE spaces
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // define our assembler level
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef Control::ScalarBasicAssemblerLevel<SpaceType> AssemblerLevelType;

    // get our domain level and layer
    typedef typename DomainControlType::LayerType DomainLayerType;
    const DomainLayerType& layer = *domain.get_layers().back();
    const std::deque<DomainLevelType*>& domain_levels = domain.get_levels();

    std::deque<SystemLevelType*> system_levels;
    std::deque<AssemblerLevelType*> asm_levels;

    const Index num_levels = domain_levels.size();

    //Lin-Solve phase related typedefs
    //Main-CSR or CUDA-ELL
    typedef typename TargetMatrixSolve_::MemType MemTypeSolve;
    typedef Control::ScalarMeanFilterSystemLevel<MemTypeSolve, DataType, IndexType, TargetMatrixSolve_> SystemLevelTypeSolve;

    std::deque<SystemLevelTypeSolve*> system_levels_solve;

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      asm_levels.push_back(new AssemblerLevelType(*domain_levels.at(i)));
      system_levels.push_back(new SystemLevelType());
    }

    /* ***************************************************************************************** */

    comm.print("Creating gates...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gate(layer, *domain_levels.at(i), asm_levels.at(i)->space);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system matrices...");

    Cubature::DynamicFactory cubature("auto-degree:5");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_laplace_matrix(asm_levels.at(i)->space, cubature);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_mean_filter(asm_levels.at(i)->space, cubature);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling transfer operators...");

    for (Index i(1); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_transfer(asm_levels.at(i)->space, asm_levels.at(i-1)->space, cubature);
    }

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    //typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    //typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain_levels.back();
    SystemLevelType& the_system_level = *system_levels.back();
    AssemblerLevelType& the_asm_level = *asm_levels.back();

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
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_f, force_func, the_asm_level.space, cubature);

      // sync the vector
      vec_rhs.sync_0();
    }

    // and filter it
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_rhs);

    ////////////////// solver type conversion ////////////////////////

    // get our global solver system types
    typedef typename SystemLevelTypeSolve::GlobalSystemVector GlobalSystemVectorSolve;
    typedef typename SystemLevelTypeSolve::GlobalSystemMatrix GlobalSystemMatrixSolve;
    typedef typename SystemLevelTypeSolve::GlobalSystemFilter GlobalSystemFilterSolve;

    comm.print("Converting assembled linear system from " + SystemLevelType::LocalScalarMatrix::name() +
      ", Mem:" + MemType::name() + " to " + SystemLevelTypeSolve::LocalScalarMatrix::name() + ", Mem:" +
      MemTypeSolve::name() + "...");

    //convert system and transfer levels
    for (Index i(0); i < num_levels; ++i)
    {
      //system levels must be converted first, because transfer levels use their converted gates
      system_levels_solve.push_back(new SystemLevelTypeSolve());
      system_levels_solve.back()->convert(*system_levels.at(i));
    }

    // get our global solve matrix and filter
    GlobalSystemMatrixSolve& matrix_solve = (*system_levels_solve.back()).matrix_sys;
    GlobalSystemFilterSolve& filter_solve = (*system_levels_solve.back()).filter_sys;

    //convert rhs and sol vectors
    GlobalSystemVectorSolve vec_rhs_solve;
    vec_rhs_solve.convert(&system_levels_solve.back()->gate_sys, vec_rhs);
    GlobalSystemVectorSolve vec_sol_solve;
    vec_sol_solve.convert(&system_levels_solve.back()->gate_sys, vec_sol);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a multigrid cycle
    auto multigrid_hierarchy = std::make_shared<Solver::MultiGridHierarchy<
      typename SystemLevelTypeSolve::GlobalSystemMatrix,
      typename SystemLevelTypeSolve::GlobalSystemFilter,
      typename SystemLevelTypeSolve::GlobalSystemTransfer
      > >();

    // scaling factor
    DataType omega = DataType(0.2);

    // push levels into MGV
    auto it_end = --system_levels_solve.rend();
    for (auto it = system_levels_solve.rbegin(); it != it_end; ++it)
    {
      auto smoother = Solver::new_scale_precond((*it)->filter_sys, omega);
      multigrid_hierarchy->push_level((*it)->matrix_sys, (*it)->filter_sys, (*it)->transfer_sys, smoother, smoother, smoother);
    }

    // create coarse grid solver
    auto coarse_solver = Solver::new_scale_precond(system_levels_solve.front()->filter_sys, omega);
    multigrid_hierarchy->push_level(system_levels_solve.front()->matrix_sys, system_levels_solve.front()->filter_sys, coarse_solver);

    // create our solver
    auto mgv = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);
    auto solver = Solver::new_pcg(matrix_solve, filter_solve, mgv);
    //auto solver = Solver::new_bicgstab(matrix_solve, filter_solve, mgv);
    //auto solver = Solver::new_fgmres(matrix_solve, filter_solve, 8, 0.0, mgv);
    //auto solver = Solver::new_richardson(matrix_solve, filter_solve, 1.0, mgv);

    // enable plotting
    solver->set_plot(comm.rank() == 0);

    // set tolerance
    solver->set_tol_rel(1E-8);
    solver->set_max_iter(1000);

    // initialise
    multigrid_hierarchy->init();
    solver->init();

    // solve
    Solver::solve(*solver, vec_sol_solve, vec_rhs_solve, matrix_solve, filter_solve);

    // release solver
    solver->done();
    multigrid_hierarchy->done();

    // download solution
    vec_sol.convert(&system_levels.back()->gate_sys, vec_sol_solve);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("no-err") < 0)
    {
      // Compute and print the H0-/H1-errors
      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute
        (vec_sol.local(), sol_func, the_asm_level.space, cubature);

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
      Assembly::DiscreteVertexProjector::project(vtx_sol, vec_sol.local(), the_asm_level.space);
      Assembly::DiscreteVertexProjector::project(vtx_rhs, vec_rhs.local(), the_asm_level.space);

      // write velocity
      exporter.add_vertex_scalar("sol", vtx_sol.elements());
      exporter.add_vertex_scalar("rhs", vtx_rhs.elements());

      // finally, write the VTK file
      exporter.write(vtk_name, comm);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // clean up
    while (!system_levels.empty())
    {
      delete system_levels.back();
      system_levels.pop_back();
    }
    while (!asm_levels.empty())
    {
      delete asm_levels.back();
      asm_levels.pop_back();
    }

    while (!system_levels_solve.empty())
    {
      delete system_levels_solve.back();
      system_levels_solve.pop_back();
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
    args.support("mem");

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

    int lvl_max = 3;
    int lvl_min = 0;
    args.parse("level", lvl_max, lvl_min);

    FEAT::String mem_string = "main";
    args.parse("mem", mem_string);

#ifndef DEBUG
    try
#endif
    {
      TimeStamp stamp1;

      // let's create our domain
      Control::Domain::UnitCubeDomainControl<MeshType> domain(&comm, lvl_max, lvl_min);

      // plot our levels
      comm.print("LVL-MIN: " + stringify(domain.get_levels().front()->get_level_index()) + " [" + stringify(lvl_min) + "]");
      comm.print("LVL-MAX: " + stringify(domain.get_levels().back()->get_level_index()) + " [" + stringify(lvl_max) + "]");

      // run our application
      if (mem_string == "main")
      {
        run<MeshType, LAFEM::SparseMatrixCSR<Mem::Main, double, Index> >(comm, args, domain);
      }
#ifdef FEAT_HAVE_CUDA
      else if(mem_string == "cuda")
      {
        run<MeshType, LAFEM::SparseMatrixELL<Mem::CUDA, double, Index> >(comm, args, domain);
      }
#endif
      else
      {
        throw InternalError("Memory type " + mem_string + " not known!");
      }

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
