#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
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

namespace StokesVortex2D
{
  using namespace FEAST;

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

  template<
    typename SpaceVelo_,
    typename SpacePres_>
  class StokesUnitSquareVortexAssemblerLevel :
    public Control::StokesBasicAssemblerLevel<SpaceVelo_, SpacePres_>
  {
  public:
    typedef Control::StokesBasicAssemblerLevel<SpaceVelo_, SpacePres_> BaseClass;
    typedef typename SpaceVelo_::MeshType MeshType;

  public:
    explicit StokesUnitSquareVortexAssemblerLevel(typename BaseClass::DomainLevelType& dom_lvl) :
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
      Assembly::UnitFilterAssembler<MeshType> unit_asm_x;
      Assembly::UnitFilterAssembler<MeshType> unit_asm_y;

      // loop over all boundary components
      for(int i(0); i < 4; ++i)
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

        // boundary regions 0 and 1 are the bottom and top edges; 2 and 3 are the left and right edges
        if(i < 2)
          unit_asm_y.add_mesh_part(*mesh_part);
        else
          unit_asm_x.add_mesh_part(*mesh_part);
      }

      // finally, assemble the filters
      unit_asm_x.assemble(fil_loc_v.template at<0>(), this->space_velo);
      unit_asm_y.assemble(fil_loc_v.template at<1>(), this->space_velo);
    }

    //template<typename SystemLevel_> void assemble_pressure_filter(SystemLevel_&) {}

    template<typename SystemLevel_>
    void assemble_pressure_filter(SystemLevel_& sys_level)
    {
      // get our global pressure filter
      typename SystemLevel_::GlobalPresFilter& fil_glob_p = sys_level.filter_pres;

      // get our local pressure filter
      typedef typename SystemLevel_::LocalPresFilter MeanFilterType;
      MeanFilterType& fil_loc_p = *fil_glob_p;

      // create two global vectors
      typename SystemLevel_::GlobalPresVector vec_glob_v(&sys_level.gate_pres), vec_glob_w(&sys_level.gate_pres);

      // fetch the local vectors
      typename SystemLevel_::LocalPresVector& vec_loc_v = *vec_glob_v;
      typename SystemLevel_::LocalPresVector& vec_loc_w = *vec_glob_w;

      // fetch the frequency vector of the pressure gate
      typename SystemLevel_::LocalPresVector& vec_loc_f = sys_level.gate_pres._freqs;

      // assemble the mean filter
      Assembly::MeanFilterAssembler::assemble(vec_loc_v, vec_loc_w, this->space_pres, this->cubature);

      // synchronise the vectors
      vec_glob_v.sync_1();
      vec_glob_w.sync_0();

      // build the mean filter
      fil_loc_p = MeanFilterType(vec_loc_v.clone(), vec_loc_w.clone(), vec_loc_f.clone());
    }

    template<typename SystemLevel_>
    void assemble_system_filter(SystemLevel_& sys_level)
    {
      assemble_velocity_filter(sys_level);
      assemble_pressure_filter(sys_level);

      // clone into system filter
      (*sys_level.filter_sys).template at<0>() = (*sys_level.filter_velo).clone(LAFEM::CloneMode::Shallow);
      (*sys_level.filter_sys).template at<1>() = (*sys_level.filter_pres).clone(LAFEM::CloneMode::Shallow);
    }

    template<typename SystemLevel_>
    typename SystemLevel_::GlobalSystemVector assemble_rhs_vector(SystemLevel_& sys_level)
    {
      // create new vector
      typename SystemLevel_::GlobalSystemVector vec_rhs = sys_level.matrix_sys.create_vector_r();
      vec_rhs.format();

      // get the Fx,Fy subvectors
      typename SystemLevel_::LocalScalarVector& vec_fx = (*vec_rhs).template at<0>().template at<0>();
      typename SystemLevel_::LocalScalarVector& vec_fy = (*vec_rhs).template at<0>().template at<1>();

      // assemble the two forces
      Analytic::StaticWrapperFunction<2, ForceFuncX> force_x;
      Analytic::StaticWrapperFunction<2, ForceFuncY> force_y;
      Assembly::Common::ForceFunctional<decltype(force_x)> force_func_x(force_x);
      Assembly::Common::ForceFunctional<decltype(force_y)> force_func_y(force_y);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_fx, force_func_x, this->space_velo, this->cubature);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_fy, force_func_y, this->space_velo, this->cubature);

      // sync the vector
      vec_rhs.sync_0();

      // and filter it
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
  };

  template<typename MeshType_>
  void run(const int rank, const int nprocs, SimpleArgParser& args, Control::Domain::DomainControl<MeshType_>& domain)
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

    // define our transfer level
    typedef Control::StokesBasicTransferLevel<SystemLevelType> TransferLevelType;

    // define our trafo and FE spaces
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // define our assembler level
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef StokesUnitSquareVortexAssemblerLevel<SpaceVeloType, SpacePresType> AssemblerLevelType;

    // get our domain level and layer
    typedef typename DomainControlType::LayerType DomainLayerType;
    const DomainLayerType& layer = *domain.get_layers().back();
    const std::deque<DomainLevelType*>& domain_levels = domain.get_levels();

    std::deque<SystemLevelType*> system_levels;
    std::deque<AssemblerLevelType*> asm_levels;
    std::deque<TransferLevelType*> transfer_levels;

    const Index num_levels = Index(domain_levels.size());

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

    // create our RHS and SOL vectors
    GlobalSystemVector vec_rhs = the_asm_level.assemble_rhs_vector(the_system_level);
    GlobalSystemVector vec_sol = the_asm_level.assemble_sol_vector(the_system_level);

    // get our global matrix and filter
    GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
    GlobalSystemFilter& filter = the_system_level.filter_sys;

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
    {
      auto mgv = std::make_shared<
        Solver::BasicVCycle<
          typename SystemLevelType::GlobalMatrixBlockA,
          typename SystemLevelType::GlobalVeloFilter,
          typename TransferLevelType::GlobalVeloTransferMatrix
        > > ();

      // create coarse grid solver
      auto coarse_solver = Solver::new_jacobi_precond(system_levels.front()->matrix_a, system_levels.front()->filter_velo);
      mgv->set_coarse_level(system_levels.front()->matrix_a, system_levels.front()->filter_velo, coarse_solver);

      // push levels into MGV
      auto jt = transfer_levels.begin();
      for(auto it = ++system_levels.begin(); it != system_levels.end(); ++it, ++jt)
      {
        auto smoother = Solver::new_jacobi_precond((*it)->matrix_a, (*it)->filter_velo);
        mgv->push_level((*it)->matrix_a, (*it)->filter_velo, (*jt)->prol_velo, (*jt)->rest_velo, smoother, smoother);
      }

      // create a PCG solver
      auto pcg = Solver::new_pcg
      /* std::make_shared<
        Solver::PCG<
          typename SystemLevelType::GlobalMatrixBlockA,
          typename SystemLevelType::GlobalVeloFilter
        > > */(system_levels.back()->matrix_a, system_levels.back()->filter_velo, mgv);

      // set its tolerances
      pcg->set_tol_rel(1E-2);
      pcg->set_tol_abs(def_init * tol_rel * 1E-1);
      pcg->set_max_iter(1000);
      //pcg->set_plot(rank == 0);

      // set our A-solver
      solver_a = pcg;
    }

    // create S-solver
    {
      // create a local ILU(0) for S
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
        Solver::SchurType::upper
      );

    // create our solver
    auto solver = Solver::new_richardson(matrix, filter, DataType(1), schur);

    // enable plotting
    solver->set_plot(rank == 0);

    // set tolerance
    solver->set_tol_rel(tol_rel);

    solver->set_max_iter(1000);

    // initialise
    solver->init();

    // solve
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // release solver
    solver->done();

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
      String vtk_name = String("./stokes-vortex-2d");
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
    int lvl_min = 0;
    args.parse("level", lvl_max, lvl_min);

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
        std::cout << "LVL-MIN: " << domain.get_levels().front()->get_level_index() << " [" << lvl_min << "]" << std::endl;
        std::cout << "LVL-MAX: " << domain.get_levels().back ()->get_level_index() << " [" << lvl_max << "]" << std::endl;
      }

      // run our application
      run(rank, nprocs, args, domain);

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

} // namespace StokesVortex2D

int main(int argc, char* argv[])
{
  return StokesVortex2D::main(argc, argv);
}
