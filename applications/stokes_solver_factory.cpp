#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/solver/legacy_preconditioners.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcr.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/schur_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/solver/matrix_stock.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_basic.hpp>
#include <control/statistics.hpp>
#include <control/solver_factory.hpp>

namespace StokesPoiseuille2D
{
  using namespace FEAT;

  template<typename T_>
  struct VeloFuncX
  {
    static T_ eval (T_, T_ y) {return y * (T_(1) - y);}
    static T_ der_x(T_, T_  ) {return T_(0);}
    static T_ der_y(T_, T_ y) {return T_(1) - T_(2)*y;}
  };

  template<typename T_>
  struct VeloFuncY
  {
    static T_ eval (T_, T_) {return T_(0);}
    static T_ der_x(T_, T_) {return T_(0);}
    static T_ der_y(T_, T_) {return T_(0);}
  };

  template<typename T_>
  struct PresFunc
  {
    static T_ eval (T_ x, T_) {return T_(2)*(T_(1) - x);}
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
    typedef Control::StokesUnitVeloNonePresSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

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

    // our inflow BC function
    Analytic::StaticWrapperFunction<2, VeloFuncX> inflow_func;

    // the names of the mesh parts on which to assemble
    std::deque<String> part_names;
    part_names.push_back("bnd:l");
    part_names.push_back("bnd:t");
    part_names.push_back("bnd:b");

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local velocity filter
      auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> unit_asm;

      // loop over all boundary parts except for the right one, which is outflow
      for(const auto& name : part_names)
      {
        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
        XASSERT(mesh_part_node != nullptr);

        auto* mesh_part = mesh_part_node->get_mesh();
        if (mesh_part != nullptr)
        {
          // add to boundary assembler
          unit_asm.add_mesh_part(*mesh_part);
        }
      }

      // assemble the filters
      unit_asm.assemble(fil_loc_v.template at<0>(), domain.at(i)->space_velo, inflow_func);
      unit_asm.assemble(fil_loc_v.template at<1>(), domain.at(i)->space_velo);

      // finally, compile the system filter
      system_levels.at(i)->compile_system_filter();
    }

    Statistics::toe_assembly = stamp_ass.elapsed_now();

    /* ***************************************************************************************** */

    // get our global solver system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    //typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    //typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    //typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    //typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // create new vectors
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

    // format the vectors
    vec_sol.format();
    vec_rhs.format();

    // and filter it
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_rhs);

    /* ***************************************************************************************** */

    comm.print("Creating solver tree");

    ////////// MATRIX STOCK
    Solver::MatrixStock<
      typename SystemLevelType::GlobalMatrixBlockA,
      typename SystemLevelType::GlobalVeloFilter,
      typename SystemLevelType::GlobalVeloTransfer> matrix_stock_a(domain.size_virtual());
    for (auto& system_level : system_levels)
    {
      matrix_stock_a.systems.push_back(system_level->matrix_a.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_a.gates_row.push_back(&system_level->gate_velo);
      matrix_stock_a.gates_col.push_back(&system_level->gate_velo);
      matrix_stock_a.filters.push_back(system_level->filter_velo.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_a.muxers.push_back(&system_level->coarse_muxer_velo);
      matrix_stock_a.transfers.push_back(system_level->transfer_velo.clone(LAFEM::CloneMode::Shallow));
    }

    Solver::MatrixStock<
      typename SystemLevelType::GlobalSchurMatrix,
      typename SystemLevelType::GlobalPresFilter,
      typename SystemLevelType::GlobalPresTransfer> matrix_stock_s(domain.size_virtual());
    for (auto& system_level : system_levels)
    {
      matrix_stock_s.systems.push_back(system_level->matrix_s.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_s.gates_row.push_back(&system_level->gate_pres);
      matrix_stock_s.gates_col.push_back(&system_level->gate_pres);
      matrix_stock_s.filters.push_back(system_level->filter_pres.clone(LAFEM::CloneMode::Shallow));
      matrix_stock_s.muxers.push_back(&system_level->coarse_muxer_pres);
      matrix_stock_s.transfers.push_back(system_level->transfer_pres.clone(LAFEM::CloneMode::Shallow));
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // our A/S block solvers
    String solver_ini_name;
    args.parse("solver-ini", solver_ini_name);
    PropertyMap property_map;
    property_map.parse(solver_ini_name, true);
    auto solver_a = Control::SolverFactory::create_scalar_solver(matrix_stock_a, &property_map, "linsolver_a");
    auto solver_s = Control::SolverFactory::create_scalar_solver(matrix_stock_s, &property_map, "linsolver_s");

    matrix_stock_a.hierarchy_init();
    matrix_stock_s.hierarchy_init();

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
    auto solver = Solver::new_pcr(the_system_level.matrix_sys, the_system_level.filter_sys, schur);

    // enable plotting
    if(comm.rank() == 0)
    {
      solver->set_plot_mode(Solver::PlotMode::iter);
    }

    solver->set_max_iter(1000);

    // initialise
    solver->init();

    Statistics::reset();

    TimeStamp at;

    // solve
    Solver::solve(*solver, vec_sol, vec_rhs, the_system_level.matrix_sys, the_system_level.filter_sys);

    double solver_toe(at.elapsed_now());

    FEAT::Control::Statistics::report(solver_toe, args.check("statistics"), MeshType::ShapeType::dimension,
      system_levels, domain);

    // release solver
    solver->done();
    matrix_stock_a.hierarchy_done();
    matrix_stock_s.hierarchy_done();

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
      const auto& vx = vec_sol.local().template at<0>().template at<0>();
      const auto& vy = vec_sol.local().template at<0>().template at<1>();
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
      String vtk_name = String("./stokes-solver-factory");
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

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("test-iter") >= 0)
    {
      int num_iter = (int)solver->get_num_iter();
      int iter_target(0);
      args.parse("test-iter", iter_target);
      if (num_iter < iter_target - 1 || num_iter > iter_target + 1)
      {
        comm.print("FAILED");
        throw InternalError(__func__, __FILE__, __LINE__, "iter count deviation! " + stringify(num_iter) + " vs " + stringify(iter_target));
      }
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
    Control::Domain::add_supported_pdc_args(args);
    args.support("mesh");
    args.support("level");
    args.support("no-err");
    args.support("vtk");
    args.support("statistics");
    args.support("solver-ini");
    args.support("test-iter");

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

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // create a time-stamp
    TimeStamp time_stamp;

    // let's create our domain
    comm.print("Preparing domain...");

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
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

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }
} // namespace StokesPoiseuille2D

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialise(argc, argv);
  try
  {
    StokesPoiseuille2D::main(argc, argv);
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
