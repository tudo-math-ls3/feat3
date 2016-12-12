#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/solver/legacy_preconditioners.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/schur_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>

#include <control/domain/unit_cube_domain_control.hpp>
#include <control/stokes_basic.hpp>

namespace StokesDriCav2D
{
  using namespace FEAT;

  template<typename T_>
  struct VeloFuncX
  {
    static T_ eval (T_ x, T_ y)
    {
      if((y > T_(0.99)) && (x > T_(0)) && (x < T_(1)))
        return T_(1);
      else
        return T_(0);
    }
  };

  template<typename MeshType_>
  void run(const Dist::Comm& comm, SimpleArgParser& args, Control::Domain::DomainControl<MeshType_>& domain)
  {
    // define our mesh type
    typedef MeshType_ MeshType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // define our arch types
    typedef Mem::Main MemType;
    typedef double DataType;
    typedef Index IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<MeshType_> DomainControlType;

    // define our system level
    typedef Control::StokesUnitVeloMeanPresSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    // define our trafo and FE spaces
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // define our assembler level
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef Control::StokesBasicAssemblerLevel<SpaceVeloType, SpacePresType> AssemblerLevelType;

    // get our domain level and layer
    typedef typename DomainControlType::LayerType DomainLayerType;
    const DomainLayerType& layer = *domain.get_layers().back();
    const std::deque<DomainLevelType*>& domain_levels = domain.get_levels();

    std::deque<SystemLevelType*> system_levels;
    std::deque<AssemblerLevelType*> asm_levels;

    const Index num_levels = Index(domain_levels.size());

    // create stokes and system levels
    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.push_back(new AssemblerLevelType(*domain_levels.at(i)));
      system_levels.push_back(new SystemLevelType());
    }

    /* ***************************************************************************************** */

    comm.print("Creating gates...");

    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gates(layer, *domain_levels.at(i), asm_levels.at(i)->space_velo, asm_levels.at(i)->space_pres);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system matrices...");

    Cubature::DynamicFactory cubature("auto-degree:5");

    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_velocity_laplace_matrix(asm_levels.at(i)->space_velo, cubature);
      system_levels.at(i)->assemble_grad_div_matrices(asm_levels.at(i)->space_velo, asm_levels.at(i)->space_pres, cubature);
      system_levels.at(i)->compile_system_matrix();
    }

    // assemble Schur-matrix on finest level
    {
      // get the local matrix S
      auto& mat_loc_s = system_levels.back()->matrix_s.local();

      // assemble matrix structure?
      Assembly::SymbolicAssembler::assemble_matrix_std1(mat_loc_s, asm_levels.back()->space_pres);

      // assemble schur matrix
      mat_loc_s.format();
      Assembly::Common::IdentityOperator id_op;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_loc_s, id_op, asm_levels.back()->space_pres, cubature, -DataType(1));
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local velocity filter
      auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // create unit-filter assemblers
      Assembly::UnitFilterAssembler<MeshType> unit_asm;

      // loop over all boundary components
      for(int k(0); k < 4; ++k)
      {
        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain_levels.at(i)->get_mesh_node()->find_mesh_part_node(String("bnd:") + stringify(k));
        XASSERT(mesh_part_node != nullptr);

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        auto* mesh_part = mesh_part_node->get_mesh();
        if(mesh_part == nullptr)
          continue;

        unit_asm.add_mesh_part(*mesh_part);
      }

      Analytic::StaticWrapperFunction<2, VeloFuncX> velo_x_func;

      // assemble the filters
      unit_asm.assemble(fil_loc_v.get(0), asm_levels.at(i)->space_velo, velo_x_func);
      unit_asm.assemble(fil_loc_v.get(1), asm_levels.at(i)->space_velo);

      // assemble pressure mean filter
      system_levels.at(i)->assemble_pressure_mean_filter(asm_levels.at(i)->space_pres, cubature);

      // compile system filter
      system_levels.at(i)->compile_system_filter();
    }

    /* ***************************************************************************************** */

    comm.print("Assembling transfer matrices...");

    for(Index i(1); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_transfers(
        asm_levels.at(i)->space_velo, asm_levels.at(i)->space_pres,
        asm_levels.at(i-1)->space_velo, asm_levels.at(i-1)->space_pres, cubature);
    }

    /* ***************************************************************************************** */

    // get our global system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain_levels.back();
    SystemLevelType& the_system_level = *system_levels.back();
    AssemblerLevelType& the_asm_level = *asm_levels.back();

    // get our global matrix and filter
    GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
    GlobalSystemFilter& filter = the_system_level.filter_sys;

    // create our RHS and SOL vectors
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

    // format the vectors
    vec_sol.format();
    vec_rhs.format();

    // and filter it
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_rhs);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // our A/S block solvers
    std::shared_ptr<Solver::SolverBase<GlobalVeloVector>> solver_a(nullptr);
    std::shared_ptr<Solver::SolverBase<GlobalPresVector>> solver_s(nullptr);

    // create a multigrid cycle A-solver
    auto multigrid_hierarchy_a = std::make_shared<
      Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalMatrixBlockA,
      typename SystemLevelType::GlobalVeloFilter,
      typename SystemLevelType::GlobalVeloTransfer
        > >();

    {
      // push levels into MGV
      auto it_end = --system_levels.rend();
      for (auto it = system_levels.rbegin(); it != it_end; ++it)
      {
        auto smoother = Solver::new_jacobi_precond((*it)->matrix_a, (*it)->filter_velo);
        multigrid_hierarchy_a->push_level((*it)->matrix_a, (*it)->filter_velo, (*it)->transfer_velo, smoother, smoother, smoother);
      }

      // create coarse grid solver
      auto coarse_solver = Solver::new_jacobi_precond(system_levels.front()->matrix_a, system_levels.front()->filter_velo);
      multigrid_hierarchy_a->push_level(system_levels.front()->matrix_a, system_levels.front()->filter_velo, coarse_solver);

      // set our A-solver
      solver_a = Solver::new_multigrid(multigrid_hierarchy_a, Solver::MultiGridCycle::V);
    }

    // create S-solver
    {
      // create a local ILU(0) for S
      /// \todo do not use global pressure filter here...
      auto loc_ilu = Solver::new_ilu_precond(*the_system_level.matrix_s, *the_system_level.filter_pres, Index(0));

      // make it Schwarz...
      auto glob_ilu = Solver::new_schwarz_precond(loc_ilu, the_system_level.filter_pres);

      // set our S-solver
      solver_s = glob_ilu;
    }

    // create a global Schur-Complement preconditioner
    auto schur = Solver::new_schur_precond(
        the_system_level.matrix_a,
        the_system_level.matrix_b,
        the_system_level.matrix_d,
        the_system_level.filter_velo,
        the_system_level.filter_pres,
        solver_a,
        solver_s,
        Solver::SchurType::full
        );

    // create our solver
    auto solver = Solver::new_pcg(matrix, filter, schur);

    // enable plotting
    solver->set_plot(comm.rank() == 0);

    solver->set_max_iter(1000);

    // initialise
    multigrid_hierarchy_a->init();
    solver->init();

    // solve
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // release solver
    solver->done();
    multigrid_hierarchy_a->done();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if(args.check("no-err") < 0)
    {
      // define reference solution functions
      Analytic::Common::ConstantFunction<2> zero_func;

      // compute local errors
      auto vi = Assembly::VelocityAnalyser::compute(vec_sol.local().template at<0>(), the_asm_level.space_velo, cubature);
      auto pi = Assembly::ScalarErrorComputer<0>::compute(
        vec_sol.local().template at<1>(), zero_func, the_asm_level.space_pres, cubature);

      // synhronise all local errors
      vi.synchronise(comm);
      pi.synchronise(comm);

      // print errors
      comm.print("");
      comm.print(vi.format_string());
      comm.print("Pressure..: " + stringify_fp_sci(pi.norm_h0));
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    //*
    if(args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./stokes-dricav-2d");
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(comm.size());

      // write VTK file
      the_asm_level.write_vtk(vtk_name, *vec_sol, comm);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // clean up
    while(!system_levels.empty())
    {
      delete system_levels.back();
      system_levels.pop_back();
    }
    while(!asm_levels.empty())
    {
      delete asm_levels.back();
      asm_levels.pop_back();
    }
  }

  void main(int argc, char* argv[])
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
      run(comm, args, domain);

      TimeStamp stamp2;

      // get times
      long long time1 = stamp2.elapsed_micros(stamp1);

      // accumulate times over all processes
      long long time2 = time1 * (long long)comm.size();

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

} // namespace StokesDriCav2D

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialise(argc, argv);
  StokesDriCav2D::main(argc, argv);
  return FEAT::Runtime::finalise();
}
