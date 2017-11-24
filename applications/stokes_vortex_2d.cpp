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
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/uzawa_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>

#include <control/domain/unit_cube_domain_control.hpp>
#include <control/stokes_basic.hpp>

namespace StokesVortex2D
{
  using namespace FEAT;

  template<typename T_>
  struct VeloFuncX
  {
    static T_ eval (T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return Math::sin(pi*x) * Math::cos(pi*y);
    }

    static T_ der_x(T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return pi * Math::cos(pi*x) * Math::cos(pi*y);
    }

    static T_ der_y(T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return -pi * Math::sin(pi*x) * Math::sin(pi*y);
    }
  };

  template<typename T_>
  struct VeloFuncY
  {
    static T_ eval (T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return -Math::cos(pi*x) * Math::sin(pi*y);
    }
    static T_ der_x(T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return pi * Math::sin(pi*x) * Math::sin(pi*y);
    }
    static T_ der_y(T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return -pi * Math::cos(pi*x) * Math::cos(pi*y);
    }
  };

  template<typename T_>
  struct PresFunc
  {
    static T_ eval (T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return Math::cos(pi*x) * Math::cos(pi*y);
    }
  };

  template<typename T_>
  struct ForceFuncX
  {
    static T_ eval (T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return -pi * (T_(1) - T_(2)*pi) * Math::sin(pi*x) * Math::cos(pi*y);
    }
  };

  template<typename T_>
  struct ForceFuncY
  {
    static T_ eval (T_ x, T_ y)
    {
      const T_ pi(Math::pi<T_>());
      return -pi * (T_(1) + T_(2)*pi) * Math::cos(pi*x) * Math::sin(pi*y);
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

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // define our system level
    typedef Control::StokesUnitVeloMeanPresSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

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
      system_levels.at(i)->assemble_gates(domain.at(i));
      if((i+1) < domain.size_virtual())
      {
        system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
        system_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature);
      }
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system matrices...");

    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_velocity_laplace_matrix(domain.at(i)->space_velo, cubature);
      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
      system_levels.at(i)->compile_system_matrix();
    }

    // assemble Schur-matrix on finest level
    {
      // get the local matrix S
      auto& mat_loc_s = system_levels.front()->matrix_s.local();

      // assemble matrix structure?
      Assembly::SymbolicAssembler::assemble_matrix_std1(mat_loc_s, domain.front()->space_pres);

      // assemble schur matrix
      mat_loc_s.format();
      Assembly::Common::IdentityOperator id_op;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_loc_s, id_op, domain.front()->space_pres, cubature, -DataType(1));
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local velocity filter
      auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // create unit-filter assemblers
      Assembly::UnitFilterAssembler<MeshType> unit_asm_x, unit_asm_y;

      // loop over all boundary components
      for(int k(0); k < 4; ++k)
      {
        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(String("bnd:") + stringify(k));
        XASSERT(mesh_part_node != nullptr);

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        auto* mesh_part = mesh_part_node->get_mesh();
        if(mesh_part == nullptr)
          continue;

        // boundary regions 0 and 1 are the bottom and top edges; 2 and 3 are the left and right edges
        if(k < 2)
          unit_asm_y.add_mesh_part(*mesh_part);
        else
          unit_asm_x.add_mesh_part(*mesh_part);
      }

      // assemble the filters
      unit_asm_x.assemble(fil_loc_v.get(0), domain.at(i)->space_velo);
      unit_asm_y.assemble(fil_loc_v.get(1), domain.at(i)->space_velo);

      // assemble pressure mean filter
      system_levels.at(i)->assemble_pressure_mean_filter(domain.at(i)->space_pres, cubature);

      // compile system filter
      system_levels.at(i)->compile_system_filter();
    }

    /* ***************************************************************************************** */

    // get our global system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // get our global matrix and filter
    GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
    GlobalSystemFilter& filter = the_system_level.filter_sys;

    // create our RHS and SOL vectors
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

    // format the vectors
    vec_sol.format();
    vec_rhs.format();

    // assemble RHS vector
    {
      // get the Fx,Fy subvectors
      typename SystemLevelType::LocalScalarVector& vec_fx = vec_rhs.local().template at<0>().template at<0>();
      typename SystemLevelType::LocalScalarVector& vec_fy = vec_rhs.local().template at<0>().template at<1>();

      // assemble the two forces
      Analytic::StaticWrapperFunction<2, ForceFuncX> force_x;
      Analytic::StaticWrapperFunction<2, ForceFuncY> force_y;
      Assembly::Common::ForceFunctional<decltype(force_x)> force_func_x(force_x);
      Assembly::Common::ForceFunctional<decltype(force_y)> force_func_y(force_y);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_fx, force_func_x, the_domain_level.space_velo, cubature);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_fy, force_func_y, the_domain_level.space_velo, cubature);

      // sync the vector
      vec_rhs.sync_0();
    }

    // and filter it
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_rhs);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // our A/S block solvers
    std::shared_ptr<Solver::SolverBase<GlobalVeloVector>> solver_a(nullptr);
    std::shared_ptr<Solver::SolverBase<GlobalPresVector>> solver_s(nullptr);

    // compute the initial defect norm
    const DataType def_init = vec_rhs.norm2();

    // choose a relative tolerance
    const DataType tol_rel = 1E-8;

    // create a multigrid cycle A-solver
    auto multigrid_hierarchy_a = std::make_shared<
      Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalMatrixBlockA,
      typename SystemLevelType::GlobalVeloFilter,
      typename SystemLevelType::GlobalVeloTransfer
        > >(domain.size_virtual());
    {
      // push levels into MGV
      for (std::size_t i(0); i < system_levels.size(); ++i)
      {
        const SystemLevelType& lvl = *system_levels.at(i);

        // Is this the virtual coarse level?
        if((i+1) < domain.size_virtual())
        {
          auto smoother = Solver::new_jacobi_precond(lvl.matrix_a, lvl.filter_velo);
          multigrid_hierarchy_a->push_level(lvl.matrix_a, lvl.filter_velo, lvl.transfer_velo, smoother, smoother, smoother);
        }
        else
        {
          auto coarse_solver = Solver::new_jacobi_precond(lvl.matrix_a, lvl.filter_velo);
          multigrid_hierarchy_a->push_level(lvl.matrix_a, lvl.filter_velo, coarse_solver);
        }
      }

      auto mgv = Solver::new_multigrid(multigrid_hierarchy_a, Solver::MultiGridCycle::V);

      // create a PCG solver
      auto pcg = Solver::new_pcg(the_system_level.matrix_a, the_system_level.filter_velo, mgv);

      // set its tolerances
      pcg->set_tol_rel(1E-2);
      pcg->set_tol_abs(def_init * tol_rel * 1E-1);
      pcg->set_max_iter(1000);
      //if(comm.rank() == 0)
      //{
      //  pcg->set_plot_mode(rank == 0);
      //}

      // set our A-solver
      solver_a = pcg;
    }

    // create S-solver
    {
      // create a local ILU(0) for S
      auto loc_ilu = Solver::new_ilu_precond(the_system_level.matrix_s.local(), the_system_level.filter_pres.local(), Index(0));

      // make it Schwarz...
      auto glob_ilu = Solver::new_schwarz_precond(loc_ilu, the_system_level.filter_pres);

      // set our S-solver
      solver_s = glob_ilu;
    }

    // create a global Uzawa preconditioner
    auto uzawa = Solver::new_uzawa_precond(
      the_system_level.matrix_a,
      the_system_level.matrix_b,
      the_system_level.matrix_d,
      the_system_level.filter_velo,
      the_system_level.filter_pres,
      solver_a,
      solver_s,
      Solver::UzawaType::upper
    );

    // create our solver
    auto solver = Solver::new_richardson(matrix, filter, DataType(1), uzawa);

    // enable plotting
    if(comm.rank() == 0)
    {
      solver->set_plot_mode(Solver::PlotMode::iter);
    }

    // set tolerance
    solver->set_tol_rel(tol_rel);

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
      Analytic::StaticWrapperFunction<2, VeloFuncX, true, true> velo_x_func;
      Analytic::StaticWrapperFunction<2, VeloFuncY, true, true> velo_y_func;
      Analytic::StaticWrapperFunction<2, PresFunc> pres_func;

      // fetch our vector components
      const auto& vx = vec_sol.local().template at<0>().get(0);
      const auto& vy = vec_sol.local().template at<0>().get(1);
      const auto& vp = vec_sol.local().template at<1>();

      // compute local errors
      Assembly::ScalarErrorInfo<DataType> vxerr = Assembly::ScalarErrorComputer<1>::compute(
        vx, velo_x_func, the_domain_level.space_velo, cubature);
      Assembly::ScalarErrorInfo<DataType> vyerr = Assembly::ScalarErrorComputer<1>::compute(
        vy, velo_y_func, the_domain_level.space_velo, cubature);
      Assembly::ScalarErrorInfo<DataType> vperr = Assembly::ScalarErrorComputer<0>::compute(
        vp, pres_func, the_domain_level.space_pres, cubature);

      // synhronise all local errors
      vxerr.synchronise(comm);
      vyerr.synchronise(comm);

      // compute field errors
      DataType vv_h0 = Math::sqrt(Math::sqr(vxerr.norm_h0) + Math::sqr(vyerr.norm_h0));
      DataType vv_h1 = Math::sqrt(Math::sqr(vxerr.norm_h1) + Math::sqr(vyerr.norm_h1));

      // print errors
      if (comm.rank() == 0)
      {
        std::cout << "Velocity H0-Error: " << stringify_fp_sci(vv_h0, 12) << " [ ";
        std::cout << stringify_fp_sci(vxerr.norm_h0, 12) << " , " << stringify_fp_sci(vyerr.norm_h0, 12) << " ]" << std::endl;
        std::cout << "Velocity H1-Error: " << stringify_fp_sci(vv_h1, 12) << " [ ";
        std::cout << stringify_fp_sci(vxerr.norm_h1, 12) << " , " << stringify_fp_sci(vyerr.norm_h1, 12) << " ]" << std::endl;
        std::cout << "Pressure H0-Error: " << stringify_fp_sci(vperr.norm_h0, 12) << std::endl;
      }
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    //*
    if(args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./stokes-vortex-2d");
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(comm.size());

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

      // project velocity and pressure
      LAFEM::DenseVector<Mem::Main, double, Index> vtx_vx, vtx_vy;
      Assembly::DiscreteVertexProjector::project(vtx_vx, (*vec_sol).template at<0>().get(0), the_domain_level.space_velo);
      Assembly::DiscreteVertexProjector::project(vtx_vy, (*vec_sol).template at<0>().get(1), the_domain_level.space_velo);
      exporter.add_vertex_vector("velocity", vtx_vx.elements(), vtx_vy.elements());

      // project pressure
      Cubature::DynamicFactory cub("auto-degree:2");
      LAFEM::DenseVector<Mem::Main, double, Index> vtx_p;
      Assembly::DiscreteCellProjector::project(vtx_p, (*vec_sol).template at<1>(), the_domain_level.space_pres, cub);

      // write pressure
      exporter.add_cell_scalar("pressure", vtx_p.elements());

      // finally, write the VTK file
      exporter.write(vtk_name, comm);
    }
  }

  void main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

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
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // parse levels
    int lvl_max = 3;
    int lvl_min = 0;
    args.parse("level", lvl_max, lvl_min);

    // create a time-stamp
    TimeStamp time_stamp;

    // let's create our domain
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    Control::Domain::UnitCubeDomainControl<DomainLevelType> domain(comm, lvl_max, lvl_min);

    // plot our levels
    comm.print("LVL-MAX: " + stringify(domain.max_level_index()) + " [" + stringify(lvl_max) + "]");
    comm.print("LVL-MIN: " + stringify(domain.min_level_index()) + " [" + stringify(lvl_min) + "]");

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }
} // namespace StokesVortex2D

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialise(argc, argv);
  try
  {
    StokesVortex2D::main(argc, argv);
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
