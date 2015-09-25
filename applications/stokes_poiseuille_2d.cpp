#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/solver/basic_vcycle.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/schur_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>

#include <control/domain/unit_cube_domain_control.hpp>
#include <control/stokes_basic.hpp>

namespace StokesPoiseuille2D
{
  using namespace FEAST;

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

  template<
    typename SpaceVelo_,
    typename SpacePres_>
  class StokesUnitSquarePoiseuilleAssemblerLevel :
    public Control::StokesBasicAssemblerLevel<SpaceVelo_, SpacePres_>
  {
  public:
    typedef Control::StokesBasicAssemblerLevel<SpaceVelo_, SpacePres_> BaseClass;
    typedef typename SpaceVelo_::MeshType MeshType;

  public:
    explicit StokesUnitSquarePoiseuilleAssemblerLevel(typename BaseClass::DomainLevelType& dom_lvl) :
      BaseClass(dom_lvl)
    {
    }

    template<typename SystemLevel_>
    void assemble_velocity_filter(SystemLevel_& sys_level)
    {
      // get our global velocity filter
      typename SystemLevel_::GlobalVeloFilter& fil_glob_v = sys_level.filter_velo;

      // get our local velocity filter
      typename SystemLevel_::LocalVeloFilter& fil_loc_v = *fil_glob_v;

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> unit_asm;

      // loop over the first three boundary components,
      // the fourth boundary component is the outflow region
      for(int i(0); i < 3; ++i)
      {
        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = this->domain_level.get_mesh_node()->find_mesh_part_node(String("bnd:") + stringify(i));

        // found it?
        if(mesh_part_node == nullptr)
          throw InternalError("Mesh Part Node 'bnd:" + stringify(i) + "' not found!");

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        auto* mesh_part = mesh_part_node->get_mesh();
        if(mesh_part == nullptr)
          continue;

        // add our boundary mesh part
        unit_asm.add_mesh_part(*mesh_part);
      }

      // our inflow BC function
      Analytic::StaticWrapperFunction<2, VeloFuncX> inflow_func;

      // finally, assemble the filters
      unit_asm.assemble(fil_loc_v.template at<0>(), this->space_velo, inflow_func);
      unit_asm.assemble(fil_loc_v.template at<1>(), this->space_velo);
    }

    template<typename SystemLevel_>
    void assemble_pressure_filter(SystemLevel_& /*sys_level*/)
    {
      // nothing to do
    }

    template<typename SystemLevel_>
    void assemble_system_filter(SystemLevel_& sys_level)
    {
      assemble_velocity_filter(sys_level);
      assemble_pressure_filter(sys_level);

      // clone into system filter
      (*sys_level.filter_sys).template at<0>() = (*sys_level.filter_velo).clone();
      (*sys_level.filter_sys).template at<1>() = (*sys_level.filter_pres).clone();
    }

    template<typename SystemLevel_>
    typename SystemLevel_::GlobalSystemVector assemble_rhs_vector(SystemLevel_& sys_level)
    {
      // create new vector, format and filter it
      typename SystemLevel_::GlobalSystemVector vec_rhs = sys_level.matrix_sys.create_vector_r();
      vec_rhs.format();
      sys_level.filter_sys.filter_rhs(vec_rhs);
      return vec_rhs;
    }

    template<typename SystemLevel_>
    typename SystemLevel_::GlobalSystemVector assemble_sol_vector(SystemLevel_& sys_level)
    {
      typename SystemLevel_::GlobalSystemVector vec_sol = sys_level.matrix_sys.create_vector_r();
      vec_sol.format();
      sys_level.filter_sys.filter_sol(vec_sol);
      return vec_sol;
    }

    template<typename SystemLevel_>
    void analyse_sol_vector(bool plot, SystemLevel_& sys_level, const typename SystemLevel_::GlobalSystemVector& vec_sol)
    {
      typedef typename SystemLevel_::DataType DataType;

      // define reference solution functions
      Analytic::StaticWrapperFunction<2, VeloFuncX, true, true> velo_x_func;
      Analytic::StaticWrapperFunction<2, VeloFuncY, true, true> velo_y_func;
      Analytic::StaticWrapperFunction<2, PresFunc> pres_func;

      // fetch our vector components
      const auto& vx = (*vec_sol).template at<0>().template at<0>();
      const auto& vy = (*vec_sol).template at<0>().template at<1>();
      const auto& vp = (*vec_sol).template at<1>();

      // compute local errors
      Assembly::ScalarErrorInfo<DataType> vxerr = Assembly::ScalarErrorComputer<1>::compute(
        vx, velo_x_func, this->space_velo, this->cubature);
      Assembly::ScalarErrorInfo<DataType> vyerr = Assembly::ScalarErrorComputer<1>::compute(
        vy, velo_y_func, this->space_velo, this->cubature);
      Assembly::ScalarErrorInfo<DataType> vperr = Assembly::ScalarErrorComputer<0>::compute(
        vp, pres_func, this->space_pres, this->cubature);

      // synhronise all local errors
      vxerr.norm_h0 = sys_level.gate_sys.norm2(vxerr.norm_h0);
      vyerr.norm_h0 = sys_level.gate_sys.norm2(vyerr.norm_h0);
      vxerr.norm_h1 = sys_level.gate_sys.norm2(vxerr.norm_h1);
      vyerr.norm_h1 = sys_level.gate_sys.norm2(vyerr.norm_h1);
      vperr.norm_h0 = sys_level.gate_sys.norm2(vperr.norm_h0);

      // compute field errors
      DataType vv_h0 = Math::sqrt(Math::sqr(vxerr.norm_h0) + Math::sqr(vyerr.norm_h0));
      DataType vv_h1 = Math::sqrt(Math::sqr(vxerr.norm_h1) + Math::sqr(vyerr.norm_h1));

      // print errors
      if (plot)
      {
        std::cout << "Velocity H0-Error: " << stringify_fp_sci(vv_h0, 12) << " [ ";
        std::cout << stringify_fp_sci(vxerr.norm_h0, 12) << " , " << stringify_fp_sci(vyerr.norm_h0, 12) << " ]" << std::endl;
        std::cout << "Velocity H1-Error: " << stringify_fp_sci(vv_h1, 12) << " [ ";
        std::cout << stringify_fp_sci(vxerr.norm_h1, 12) << " , " << stringify_fp_sci(vyerr.norm_h1, 12) << " ]" << std::endl;
        std::cout << "Pressure H0-Error: " << stringify_fp_sci(vperr.norm_h0, 12) << std::endl;
      }
    }
  }; // class StokesUnitSquarePoiseuilleAssemblerLevel<...>


  template<typename MeshType_, typename TargetMatrixSolve_>
  void run(const int rank, const int nprocs, SimpleArgParser& args, Control::Domain::DomainControl<MeshType_>& domain)
  {
    // define our mesh type
    typedef MeshType_ MeshType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // define our arch types
    typedef Mem::Main MemType;
    typedef typename TargetMatrixSolve_::DataType DataType;
    typedef typename TargetMatrixSolve_::IndexType IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<MeshType_> DomainControlType;

    // define our system level
    typedef Control::StokesUnitVeloNonePresSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    // define our transfer level
    typedef Control::StokesBasicTransferLevel<SystemLevelType> TransferLevelType;

    // define our trafo and FE spaces
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // define our assembler level
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef StokesUnitSquarePoiseuilleAssemblerLevel<SpaceVeloType, SpacePresType> AssemblerLevelType;

    // get our domain level and layer
    typedef typename DomainControlType::LayerType DomainLayerType;
    const DomainLayerType& layer = *domain.get_layers().back();
    const std::deque<DomainLevelType*>& domain_levels = domain.get_levels();

    std::deque<SystemLevelType*> system_levels;
    std::deque<AssemblerLevelType*> asm_levels;
    std::deque<TransferLevelType*> transfer_levels;

    const Index num_levels = Index(domain_levels.size());

    //Lin-Solve phase related typedefs
    //Main-CSR or CUDA-ELL
    typedef typename TargetMatrixSolve_::MemType MemTypeSolve;
    typedef Control::StokesUnitVeloNonePresSystemLevel<dim, MemTypeSolve, DataType, IndexType, TargetMatrixSolve_> SystemLevelTypeSolve;
    typedef Control::StokesBasicTransferLevel<SystemLevelTypeSolve> TransferLevelTypeSolve;

    std::deque<SystemLevelTypeSolve*> system_levels_solve;
    std::deque<TransferLevelTypeSolve*> transfer_levels_solve;

    // create stokes and system levels
    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.push_back(new AssemblerLevelType(*domain_levels.at(i)));
      system_levels.push_back(new SystemLevelType());
      if(i > 0)
      {
        transfer_levels.push_back(new TransferLevelType(*system_levels.at(i-1), *system_levels.at(i)));
      }
    }

    /* ***************************************************************************************** */

    if(rank == 0)
    {
      std::cout << "Creating gates..." << std::endl;
    }

    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_gates(layer, *system_levels.at(i));
    }

    /* ***************************************************************************************** */

    if(rank == 0)
    {
      std::cout << "Assembling system matrices..." << std::endl;
    }

    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_system_matrix(*system_levels.at(i));
    }

    // assemble Schur-matrix on finest level
    asm_levels.back()->assemble_schur_matrix(*system_levels.back());

    /* ***************************************************************************************** */

    if(rank == 0)
    {
      std::cout << "Assembling system filters..." << std::endl;
    }

    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_system_filter(*system_levels.at(i));
    }

    /* ***************************************************************************************** */

    if(rank == 0)
    {
      std::cout << "Assembling transfer matrices..." << std::endl;
    }

    for(Index i(0); (i+1) < num_levels; ++i)
    {
      asm_levels.at(i+1)->assemble_system_transfer(*transfer_levels.at(i), *asm_levels.at(i));
    }

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain_levels.back();
    SystemLevelType& the_system_level = *system_levels.back();
    AssemblerLevelType& the_asm_level = *asm_levels.back();

    // create our RHS and SOL vectors
    GlobalSystemVector vec_rhs = the_asm_level.assemble_rhs_vector(the_system_level);
    GlobalSystemVector vec_sol = the_asm_level.assemble_sol_vector(the_system_level);

    ////////////////// solver type conversion ////////////////////////

    // get our global solver system types
    typedef typename SystemLevelTypeSolve::GlobalSystemVector GlobalSystemVectorSolve;
    typedef typename SystemLevelTypeSolve::GlobalSystemMatrix GlobalSystemMatrixSolve;
    typedef typename SystemLevelTypeSolve::GlobalSystemFilter GlobalSystemFilterSolve;
    typedef typename SystemLevelTypeSolve::GlobalVeloVector GlobalVeloVectorSolve;
    typedef typename SystemLevelTypeSolve::GlobalPresVector GlobalPresVectorSolve;

    if (rank == 0)
    {
      std::cout << "Converting assembled linear system from " + SystemLevelType::LocalScalarMatrix::name() <<
        ", Mem:" << MemType::name() << " to " << SystemLevelTypeSolve::LocalScalarMatrix::name() << ", Mem:" <<
        MemTypeSolve::name() << "..." << std::endl;
    }

    //convert system and transfer levels
    for (Index i(0); i < num_levels; ++i)
    {
      //system levels must be converted first, because transfer levels use their converted gates
      system_levels_solve.push_back(new SystemLevelTypeSolve());
      system_levels_solve.back()->convert(*system_levels.at(i));
      if (i > 0)
      {
        transfer_levels_solve.push_back(new TransferLevelTypeSolve());
        transfer_levels_solve.back()->convert(*system_levels_solve.at(i-1), *system_levels_solve.at(i), *transfer_levels.at(i-1));
      }
    }

    SystemLevelTypeSolve& the_system_level_solve = *system_levels_solve.back();

    // get our global solve matrix and filter
    GlobalSystemMatrixSolve& matrix_solve = the_system_level_solve.matrix_sys;
    GlobalSystemFilterSolve& filter_solve = the_system_level_solve.filter_sys;

    //convert rhs and sol vectors
    GlobalSystemVectorSolve vec_rhs_solve;
    vec_rhs_solve.convert(&system_levels_solve.back()->gate_sys, vec_rhs);
    GlobalSystemVectorSolve vec_sol_solve;
    vec_sol_solve.convert(&system_levels_solve.back()->gate_sys, vec_sol);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // our A/S block solvers
    std::shared_ptr<Solver::SolverBase<GlobalVeloVectorSolve>> solver_a(nullptr);
    std::shared_ptr<Solver::SolverBase<GlobalPresVectorSolve>> solver_s(nullptr);

    // create a multigrid cycle A-solver
    {
      auto mgv = std::make_shared<
        Solver::BasicVCycle<
          typename SystemLevelTypeSolve::GlobalMatrixBlockA,
          typename SystemLevelTypeSolve::GlobalVeloFilter,
          typename TransferLevelTypeSolve::GlobalVeloTransferMatrix
        > > ();

      // create coarse grid solver
      auto coarse_solver = Solver::new_jacobi_precond(system_levels_solve.front()->matrix_a, system_levels_solve.front()->filter_velo);
      mgv->set_coarse_level(system_levels_solve.front()->matrix_a, system_levels_solve.front()->filter_velo, coarse_solver);

      // push levels into MGV
      auto jt = transfer_levels_solve.begin();
      for(auto it = ++system_levels_solve.begin(); it != system_levels_solve.end(); ++it, ++jt)
      {
        auto smoother = Solver::new_jacobi_precond((*it)->matrix_a, (*it)->filter_velo);
        mgv->push_level((*it)->matrix_a, (*it)->filter_velo, (*jt)->prol_velo, (*jt)->rest_velo, smoother, smoother);
      }

      // set our A-solver
      solver_a = mgv;
    }

    // create S-solver
    {
      // create a local ILU(0) for S
      //auto loc_ilu = Solver::new_ilu_precond(*the_system_level_solve.matrix_s, *the_system_level_solve.filter_pres, Index(0));
      auto loc_jac = Solver::new_jacobi_precond(*the_system_level_solve.matrix_s, *the_system_level_solve.filter_pres, 0.7);
      auto loc_pcg = Solver::new_pcg(*the_system_level_solve.matrix_s, *the_system_level_solve.filter_pres, loc_jac);

      // make it Schwarz...
      //auto glob_ilu = Solver::new_schwarz_precond(loc_ilu, the_system_level_solve.filter_pres);
      auto glob_ilu = Solver::new_schwarz_precond(loc_pcg, the_system_level_solve.filter_pres);

      // set our S-solver
      solver_s = glob_ilu;
    }

    // create a global Schur-Complement preconditioner
    auto schur = Solver::new_schur_precond(
        the_system_level_solve.matrix_a,
        the_system_level_solve.matrix_b,
        the_system_level_solve.matrix_d,
        the_system_level_solve.filter_velo,
        the_system_level_solve.filter_pres,
        solver_a,
        solver_s,
        Solver::SchurType::full
      );

    // create our solver
    auto solver = Solver::new_pcg(matrix_solve, filter_solve, schur);

    // enable plotting
    solver->set_plot(rank == 0);

    solver->set_max_iter(1000);

    // initialise
    solver->init();

    // solve
    Solver::solve(*solver, vec_sol_solve, vec_rhs_solve, matrix_solve, filter_solve);

    // release solver
    solver->done();

    // download solution
    vec_sol.convert(&system_levels.back()->gate_sys, vec_sol_solve);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if(args.check("no-err") < 0)
    {
      the_asm_level.analyse_sol_vector(rank == 0, the_system_level, vec_sol);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    //*
    if(args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./stokes-poiseuille-2d");
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(nprocs);

      // write VTK file
      the_asm_level.write_vtk(vtk_name, *vec_sol, rank, nprocs);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // clean up
    while(!transfer_levels.empty())
    {
      delete transfer_levels.back();
      transfer_levels.pop_back();
    }
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

    while(!transfer_levels_solve.empty())
    {
      delete transfer_levels_solve.back();
      transfer_levels_solve.pop_back();
    }
    while(!system_levels_solve.empty())
    {
      delete system_levels_solve.back();
      system_levels_solve.pop_back();
    }
  }

  int main(int argc, char* argv[])
  {
    int rank = 0;
    int nprocs = 0;

    // initialise
    FEAST::Runtime::initialise(argc, argv, rank, nprocs);
#ifndef SERIAL
    if (rank == 0)
    {
      std::cout << "NUM-PROCS: " << nprocs << std::endl;
    }
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
      if (rank == 0)
      {
        for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
          std::cerr << "ERROR: unknown option '--" << (*it).second << "'" << std::endl;
      }
      // abort
      FEAST::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;

    int lvl_max = 3;
    int lvl_min = -1;
    args.parse("level", lvl_max, lvl_min);

    FEAST::String mem_string = "main";
    args.parse("mem", mem_string);

#ifndef DEBUG
    try
#endif
    {
      TimeStamp stamp1;

      // let's create our domain
      Control::Domain::UnitCubeDomainControl<MeshType> domain(rank, nprocs, lvl_max, lvl_min);


      // plot our levels
      if(rank == 0)
      {
        std::cout << "LVL-MIN: " << domain.min_level_index() << " [" << lvl_min << "]" << std::endl;
        std::cout << "LVL-MAX: " << domain.max_level_index() << " [" << lvl_max << "]" << std::endl;
      }

      // run our application
      if (mem_string == "main")
      {
        run<MeshType, LAFEM::SparseMatrixCSR<Mem::Main, double, Index> >(rank, nprocs, args, domain);
      }
      else if(mem_string == "cuda")
      {
        run<MeshType, LAFEM::SparseMatrixELL<Mem::CUDA, double, Index> >(rank, nprocs, args, domain);
      }

      TimeStamp stamp2;

      // get times
      long long time1 = stamp2.elapsed_micros(stamp1);

      // accumulate times over all processes
      long long time2 = time1 * (long long)nprocs;

      // print time
      if(rank == 0)
      {
        std::cout << "Run-Time: "
          << TimeStamp::format_micros(time1, TimeFormat::m_s_m) << " ["
          << TimeStamp::format_micros(time2, TimeFormat::m_s_m) << "]" << std::endl;
      }
    }
#ifndef DEBUG
    catch (const std::exception& exc)
    {
      std::cerr << "ERROR: unhandled exception: " << exc.what() << std::endl;
      FEAST::Runtime::abort();
    }
    catch (...)
    {
      std::cerr << "ERROR: unknown exception" << std::endl;
      FEAST::Runtime::abort();
    }
#endif // DEBUG

    // okay
    return FEAST::Runtime::finalise();
  }

} // namespace StokesPoiseuille2D

int main(int argc, char* argv[])
{
  return StokesPoiseuille2D::main(argc, argv);
}
