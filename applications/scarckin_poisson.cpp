#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/dist.hpp>
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
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/spai_precond.hpp>
#include <kernel/global/synch_mat.hpp>
#include <kernel/solver/multigrid.hpp>
#include <control/statistics.hpp>
#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>

#include <kernel/util/whistle_punk.hpp>

namespace PoissonDirichlet2D
{
  using namespace FEAT;

  template<typename TargetMatrixSolve_, typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our arch types
    typedef typename Mem::Main MemType;
    typedef typename TargetMatrixSolve_::DataType DataType;
    typedef typename TargetMatrixSolve_::IndexType IndexType;

    // choose our desired analytical solution
    Analytic::Common::ExpBubbleFunction<2> sol_func;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;

    // define our system level
    typedef Control::ScalarUnitFilterSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = domain.size_physical();

    //Lin-Solve phase related typedefs
    //Main-CSR or CUDA-ELL
    typedef typename TargetMatrixSolve_::MemType MemTypeSolve;
    typedef Control::ScalarUnitFilterSystemLevel<MemTypeSolve, DataType, IndexType, TargetMatrixSolve_> SystemLevelTypeSolve;

    std::deque<std::shared_ptr<SystemLevelTypeSolve>> system_levels_solve;

    // create stokes and system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }


    Cubature::DynamicFactory cubature("auto-degree:5");

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

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
      Assembly::Common::LaplaceFunctional<decltype(sol_func)> force_func(sol_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_f, force_func, the_domain_level.space, cubature);

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
    typedef typename SystemLevelTypeSolve::GlobalSystemTransfer GlobalSystemTransferSolve;
    typedef typename SystemLevelTypeSolve::LocalSystemMatrix LocalSystemMatrixSolve;
    typedef typename SystemLevelTypeSolve::LocalSystemFilter LocalSystemFilterSolve;
    typedef typename SystemLevelTypeSolve::LocalSystemTransfer LocalSystemTransferSolve;

    comm.print("Converting assembled linear system from " + SystemLevelType::LocalScalarMatrix::name() +
      ", Mem:" + MemType::name() + " to " + SystemLevelTypeSolve::LocalScalarMatrix::name() + ", Mem:" +
      MemTypeSolve::name() + "...");

    //convert system and transfer levels
    for (Index i(0); i < num_levels; ++i)
    {
      //system levels must be converted first, because transfer levels use their converted gates
      system_levels_solve.push_back(std::make_shared<SystemLevelTypeSolve>());
      system_levels_solve.back()->convert(*system_levels.at(i));
    }

    // get our global solve matrix and filter
    GlobalSystemMatrixSolve& matrix_solve = (*system_levels_solve.front()).matrix_sys;
    GlobalSystemFilterSolve& filter_solve = (*system_levels_solve.front()).filter_sys;

    //convert rhs and sol vectors
    GlobalSystemVectorSolve vec_rhs_solve;
    vec_rhs_solve.convert(&system_levels_solve.front()->gate_sys, vec_rhs);
    GlobalSystemVectorSolve vec_sol_solve;
    vec_sol_solve.convert(&system_levels_solve.front()->gate_sys, vec_sol);

    Util::WhistlePunk<> wp;

    ///Create local matrices
    auto kt(system_levels.begin());
    std::deque<LocalSystemMatrixSolve> local_matrices_solve;
    for (auto it(system_levels_solve.begin()); it != system_levels_solve.end(); ++it, ++kt)
    {
      auto& kt_ranks((*kt)->gate_sys._ranks);
      auto& kt_mirrors((*kt)->gate_sys._mirrors);

      wp.log(" A0: ", *((*kt)->matrix_sys));
      auto it_matrix = (*((*kt)->matrix_sys)).clone();
      Global::synch_matrix(it_matrix, *(*kt)->gate_sys.get_comm(), kt_ranks, kt_mirrors, kt_mirrors);

      wp.log(" A1: ", it_matrix);
      LocalSystemMatrixSolve it_matrix_solve;
      it_matrix_solve.convert(it_matrix);
      local_matrices_solve.push_back(std::move(it_matrix_solve));
    }

    ///Create a local multigrid hierarchy that can be used in the smoothers
    //#### INNER MG ####
    auto inner_multigrid_hierarchy(std::make_shared<Solver::MultiGridHierarchy<LocalSystemMatrixSolve, LocalSystemFilterSolve, LocalSystemTransferSolve> >());

    //all other levels
    Index level(0);
    const auto it_end = --system_levels_solve.end();
    for (auto it = system_levels_solve.begin(); it != it_end; ++it, ++level)
    {
      auto local_precond = Solver::new_jacobi_precond(local_matrices_solve.at(level), *((*it)->filter_sys), 0.7);
      auto local_solver = Solver::new_richardson(local_matrices_solve.at(level), *((*it)->filter_sys), 1.0, local_precond);
      //local_solver->set_max_iter(2);
      local_solver->set_max_iter(4);
      local_solver->set_min_iter(4);
      inner_multigrid_hierarchy->push_level(
        local_matrices_solve.at(level),     // the system matrix for this level
        (*it)->filter_sys.local(),          // the system filter for this level
        (*it)->transfer_sys.local(),        // the transfer operator for this level
        local_solver,       // pre-smoother
        local_solver,       // post
        local_solver,
        nullptr,
        double(1./double(Index(local_matrices_solve.size()) - level))
      );
    }

    //coarse
    {
      auto inner_coarse_solver = Solver::new_pcg(local_matrices_solve.back(), *(system_levels_solve.back()->filter_sys));

      inner_multigrid_hierarchy->push_level(
        local_matrices_solve.back(),       // the coarse-level system matrix
        *(system_levels_solve.back()->filter_sys),       // the coarse-level system filter
        inner_coarse_solver     // the coarse-level solver
      );
    }

    inner_multigrid_hierarchy->init();
    //#### END INNER MG ####

    //#### OUTER MG ####
    auto outer_multigrid_hierarchy(std::make_shared<Solver::MultiGridHierarchy<GlobalSystemMatrixSolve, GlobalSystemFilterSolve, GlobalSystemTransferSolve> >());

    //all other levels
    level = 0;
    for (auto it = system_levels_solve.begin(); it != it_end; ++it, ++level)
    {
      auto inner_multigrid = Solver::new_scarcmultigrid(
          inner_multigrid_hierarchy,
          Solver::MultiGridCycle::W,
          -int(level+1),
          0);

      auto inner_solver = Solver::new_richardson(local_matrices_solve.at(level), *((*it)->filter_sys), DataType(1), inner_multigrid);
      inner_solver->init();
      //inner_solver->set_max_iter(2);
      inner_solver->set_tol_rel(1e-1);
      inner_solver->set_plot(comm.rank() == 0);

      //auto local_precond = Solver::new_jacobi_precond(local_matrices_solve.at(level), *((*it)->filter_sys), 0.7);
      //auto local_solver = Solver::new_richardson(local_matrices_solve.at(level), *((*it)->filter_sys), 1.0, local_precond);
      //local_solver->set_max_iter(20);
      auto smoother = Solver::new_schwarz_precond(inner_solver, (*it)->filter_sys);
      //smoother->set_max_iter(4);
      //smoother->set_min_iter(4);
      outer_multigrid_hierarchy->push_level(
        (*it)->matrix_sys,     // the system matrix for this level
        (*it)->filter_sys,     // the system filter for this level
        (*it)->transfer_sys,   // the transfer operator for this level
        smoother,       // pre-smoother
        smoother,       // post
        smoother        // the peak-smoother
      );
    }

    //coarse
    {
      auto outer_coarse_solver = Solver::new_pcg(system_levels_solve.back()->matrix_sys, system_levels_solve.back()->filter_sys);

      outer_multigrid_hierarchy->push_level(
        system_levels_solve.back()->matrix_sys,       // the coarse-level system matrix
        system_levels_solve.back()->filter_sys,       // the coarse-level system filter
        outer_coarse_solver     // the coarse-level solver
      );
    }


    auto outer_multigrid = Solver::new_multigrid(
      outer_multigrid_hierarchy,
      Solver::MultiGridCycle::W);

    auto outer_solver = Solver::new_richardson(system_levels_solve.front()->matrix_sys, system_levels_solve.front()->filter_sys, DataType(1), outer_multigrid);
    outer_multigrid_hierarchy->init();
    //#### OUTER MG ####

    // enable plotting
    outer_solver->set_plot(comm.rank() == 0);

    //auto outer_krylov(Solver::new_pcg(system_levels_solve.back()->matrix_sys, system_levels_solve.back()->filter_sys, outer_solver));

    // set tolerance
    outer_solver->set_tol_rel(1E-8);
    //outer_krylov->set_tol_rel(1E-8);
    //outer_solver->set_max_iter(1000);
    outer_solver->set_max_iter(100);
    //outer_krylov->set_max_iter(1000);

    // initialise
    //outer_krylov->set_plot(rank == 0);
    //outer_krylov->init();
    outer_solver->init();

    Statistics::reset_flops();
    Statistics::reset_times();
    Statistics::reset_solver_statistics();

    TimeStamp at;

    // solve
    //Solver::solve(*outer_solver, vec_sol_solve, vec_rhs_solve, matrix_solve, filter_solve);
    //Solver::solve(*outer_krylov, vec_sol_solve, vec_rhs_solve, matrix_solve, filter_solve);
    Solver::solve(*outer_solver, vec_sol_solve, vec_rhs_solve, matrix_solve, filter_solve);

    double solver_toe(at.elapsed_now());
    FEAT::Control::Statistics::report(solver_toe, args.check("statistics"), MeshType::ShapeType::dimension,
      system_levels, domain);

    // release solver
    //outer_krylov->done();
    outer_solver->done();
    //outer_multigrid_hierarchy->done();
    //inner_multigrid_hierarchy->done();

    // download solution
    vec_sol.convert(&system_levels.front()->gate_sys, vec_sol_solve);

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
      String vtk_name = String("./scarckin-poisson-dirichlet-2d");
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

    //logger ouput
    wp.synch();

    //if(Util::Comm::rank() == 0)
    //  std::cout << wp.msg;
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
    args.support("level");
    args.support("no-err");
    args.support("vtk");
    args.support("statistics");
    args.support("mem");
    args.support("mesh");
    args.support("parti-type");
    args.support("parti-name");
    args.support("parti-rank-elems");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unknown option '--" << (*it).second << "'" << std::endl;
      // abort
      FEAT::Runtime::abort();
    }

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
      FEAT::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

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
      comm.print("Preparing domain...");

      // query mesh filename list
      const std::deque<String>& mesh_filenames = args.query("mesh")->second;

      // create our domain control
      typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
      Control::Domain::PartiDomainControl<DomainLevelType> domain(comm);

      // let the controller parse its arguments
      if(!domain.parse_args(args))
      {
        FEAT::Runtime::abort();
      }

      // read the base-mesh
      domain.read_mesh(mesh_filenames);

      // try to create the partition
      domain.create_partition();

      comm.print("Creating mesh hierarchy...");

      // create the level hierarchy
      domain.create_hierarchy(lvl_max, lvl_min);

      // plot our levels
      comm.print("LVL-MAX: " + stringify(domain.max_level_index()) + " [" + stringify(lvl_max) + "]");
      comm.print("LVL-MIN: " + stringify(domain.min_level_index()) + " [" + stringify(lvl_min) + "]");

      // run our application
      if (mem_string == "main")
      {
        run<LAFEM::SparseMatrixCSR<Mem::Main, double, Index> >(args, domain);
      }
#ifdef FEAT_HAVE_CUDA
      else if(mem_string == "cuda")
      {
        run<LAFEM::SparseMatrixELL<Mem::CUDA, double, Index> >(args, domain);
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
      comm.print("Run-Time: " + stringify(TimeStamp::format_micros(time1, TimeFormat::m_s_m)) + " [" +
        stringify(TimeStamp::format_micros(time2, TimeFormat::m_s_m)) + "]");
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
} // namespace PoissonDirichlet2D

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialise(argc, argv);
  PoissonDirichlet2D::main(argc, argv);
  return FEAT::Runtime::finalise();
}
