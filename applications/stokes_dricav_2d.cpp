#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/solver/basic_vcycle.hpp>
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

  template<
    typename SpaceVelo_,
    typename SpacePres_>
  class StokesUnitSquareDriCavAssemblerLevel :
    public Control::StokesBasicAssemblerLevel<SpaceVelo_, SpacePres_>
  {
  public:
    typedef Control::StokesBasicAssemblerLevel<SpaceVelo_, SpacePres_> BaseClass;
    typedef typename SpaceVelo_::MeshType MeshType;

  public:
    explicit StokesUnitSquareDriCavAssemblerLevel(typename BaseClass::DomainLevelType& dom_lvl) :
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
        unit_asm.add_mesh_part(*mesh_part);
      }

      Analytic::StaticWrapperFunction<2, VeloFuncX> velo_x_func;

      // finally, assemble the filters
      unit_asm.assemble(fil_loc_v.template at<0>(), this->space_velo, velo_x_func);
      unit_asm.assemble(fil_loc_v.template at<1>(), this->space_velo);
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
      fil_loc_p = MeanFilterType(vec_loc_v.clone(), vec_loc_w.clone(), vec_loc_f.clone(), sys_level.gate_pres.get_comm());
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
      Analytic::Common::ConstantFunction<2> zero_func;

      // compute local errors
      auto vi = Assembly::VelocityAnalyser::compute((*vec_sol).template at<0>(), this->space_velo, this->cubature);
      DataType vp_h0 = Assembly::ScalarErrorComputer<0>::compute(
        (*vec_sol).template at<1>(), zero_func, this->space_pres, this->cubature).norm_h0;

      // synhronise all local errors
      vi.norm_h0 = sys_level.gate_sys.norm2(vi.norm_h0);
      vi.norm_h1 = sys_level.gate_sys.norm2(vi.norm_h1);
      vi.divergence = sys_level.gate_sys.norm2(vi.divergence);
      vi.vorticity = sys_level.gate_sys.norm2(vi.vorticity);
      vi.norm_h0_comp[0] = sys_level.gate_sys.norm2(vi.norm_h0_comp[0]);
      vi.norm_h0_comp[1] = sys_level.gate_sys.norm2(vi.norm_h0_comp[1]);
      vi.norm_h1_comp[0] = sys_level.gate_sys.norm2(vi.norm_h1_comp[0]);
      vi.norm_h1_comp[1] = sys_level.gate_sys.norm2(vi.norm_h1_comp[1]);
      vp_h0 = sys_level.gate_sys.norm2(vp_h0);

      // print errors
      if (plot)
      {
        std::cout << vi;
        std::cout << "Pressure..: " << vp_h0 << std::endl;
      }
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

    // define our transfer level
    typedef Control::StokesBasicTransferLevel<SystemLevelType> TransferLevelType;

    // define our trafo and FE spaces
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // define our assembler level
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef StokesUnitSquareDriCavAssemblerLevel<SpaceVeloType, SpacePresType> AssemblerLevelType;

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

    comm.print("Creating gates...");

    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_gates(layer, *system_levels.at(i));
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system matrices...");

    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_system_matrix(*system_levels.at(i));
    }

    // assemble Schur-matrix on finest level
    asm_levels.back()->assemble_schur_matrix(*system_levels.back());

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_system_filter(*system_levels.at(i));
    }

    /* ***************************************************************************************** */

    comm.print("Assembling transfer matrices...");

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

    // create a multigrid cycle A-solver
    {
      auto mgv = std::make_shared<
        Solver::BasicVCycle<
        typename SystemLevelType::GlobalMatrixBlockA,
        typename SystemLevelType::GlobalVeloFilter,
        typename TransferLevelType::GlobalVeloTransferMatrix
        > >();

      // create coarse grid solver
      auto coarse_solver = Solver::new_jacobi_precond(system_levels.front()->matrix_a, system_levels.front()->filter_velo);
      mgv->set_coarse_level(system_levels.front()->matrix_a, system_levels.front()->filter_velo, coarse_solver);

      // push levels into MGV
      auto jt = transfer_levels.begin();
      for (auto it = ++system_levels.begin(); it != system_levels.end(); ++it, ++jt)
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
      the_asm_level.analyse_sol_vector(comm.rank() == 0, the_system_level, vec_sol);
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
