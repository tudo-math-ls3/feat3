//
// The 2D Nonsteady Navier-Stokes CP-Q2/Q1 Toy-Code Solver (TM)
// ------------------------------------------------------------
// This application implements a "simple" parallel non-steady Navier-Stokes solver using
// the "CP" approach with Q2/Q1 space and Crank-Nicolson time discretisation.
//
// ---------------
// !!! WARNING !!!
// ---------------
// This application is a "toy code" solver, i.e. it is meant as a playground for the
// HPC guys to tweak their parallel poisson solvers for more interesting scenarious than
// poisson on the unit-square. You can furthermore generate fancy videos of vortex
// streets to impress your girlfriend or to show your parents what you are being paid for.
// But: Do not expect this application to be accurate in time and/or space,
// so do *NOT* use it to perform any serious PDE/FEM analysis work !!!
// </Warning>
//
// This application is a "flow-through-a-domain" solver, i.e. it handles Navier-Stokes
// equations with an inflow and and outflow region without any external forces, moving
// boundaries or any other fancy stuff.
//
// This application has three pre-configured benchmark problems, which can be launched
// by specifying the "--setup <config>" command line arguments, where <config> specifies
// one of the following:
//
// --setup square
// Loads the "Poiseuille-Flow-On-Unit-Square" problem.
// This is the most simple of the three pre-confgured problems, where the time-dependent
// solution converges to a steady-state Poiseuille-Flow.
//
// --setup nozzle
// Loads the "Jet-Flow-Through-A-Nozzle" problem.
// Similar to the previous problem, but on a slightly more complex domain.
// The solution also converges to a steady-state flow.
//
// --setup bench1
// Loads the famous "Flow-Around-A-Cylinder" problem (non-steady version).
// This is the problem which generates the fancy "Von-Karman vortex shedding".
// In contrast to the previous problems, this solution is periodic.
//
// Moreover, this application can be configured by specifying further options, which can
// be used to define new (non-preconfigured) problems or override the pre-defined settings.
//
// IMPORTANT #1:
// In any case, you will need to specify the path to the mesh directory, as the application
// will fail to find the required mesh-file otherwise. You can either use the '--mesh-path'
// option (see below) for this or you can specify the mesh-path by defining the
// "FEAT3_PATH_MESHES" environment variable to point to the "data/meshes" directory of
// your FEAT3 checkout.
//
// IMPORTANT #2:
// If you adjust the minimum and/or maximum mesh levels for one of the pre-configured
// problems, then it may be necessary to increase the number of time-steps in some cases,
// as the non-linear solver may run amok otherwise...
//
//
// In the following, the options are categorised:
//
//
// Domain / Mesh Specification Options
// -----------------------------------
// The input mesh file is specified by 2 options, namely '--mesh-path' and '--mesh-file'.
//
// The '--mesh-path <dir>' option specifies the path to the directory which contains
// the mesh files. This usually points to the 'data/meshes' sub-directory of the FEAT3
// root directory.
//
// The '--mesh-file <file>' option specifies the filename of the mesh-file. Note that
// the filename is relative to the mesh-path specified by the previous option.
//
// The '--level <max> [<min>]' option specifies the desired minimum and maximum refinement
// levels to be used. The <min> parameter is optional and is set to 0 if not given.
//
// The '--rank-elems <N>' option specifies the minimum number of elements per rank for
// the partitioner. Default: 4
//
//
// Operator Specification Options
// ------------------------------
// The '--nu <nu>' option specifies the viscosity of the fluid.
//
// The '--deformation' option (without parameters) switches to the "deformation tensor"
// (aka Du:Dv) formulation. Without this option, the "gradient tensor" formulation is used
// for the assembly of the diffusive term.
//
//
// Boundary Condition Specification Options
// ----------------------------------------
// This application supports only very limited customisation of boundary conditions, which
// is limited to specifiying
//  1) the mesh-part for the parabolic inflow region,
//  2) the mesh-part for the "do-nothing" outflow region
// All other mesh-parts are treated as "no-flow" boundaries.
//
// The '--part-in <name>' and '--part-out <name>' options specify the names of the mesh-parts
// that serve as the inflow and outflow regions, respectively.
//
// The option '--profile <x0> <y0> <x1> <y1>' specifies the four coordinates that define
// the line segment of the parabolic inflow profile.
//
// The option '--vmax <V>' specifies the maximum inflow velocity, which is set to 1 by default.
//
//
// Time Interval/Stepping Options
// ------------------------------
// The Crank-Nicolson time-stepping scheme can be configured by three options.
//
// The '--time-max <T>' option specifies the end of the desired time interval [0,T].
//
// The '--time-steps <N>' option specifies the total number of equidistant time-steps
// for the whole time interval [0,T]. The "mesh width" of the time discretisation is
// then given by T/N.
//
// The '--max-time-steps <N>' sets the maximum number of time-steps to perform.
// This can be used if one is only interested in performing a fixed number of time steps
// for performance benchmarks or preconditioner testing.
// Note: This option does NOT affect the time-stepping "mesh width".
//
//
// Non-Linear Solver Options
// -------------------------
// This application implements a simple variant of the "CP" (Coupled-solver, Projection-preconditioner)
// approach. The non-linear solver as well as its nested "DPM" (Discrete Projection Method) always
// perform a fixed number of iterations without any convergence control.
//
// The '--nl-steps <N>' option specifies the number of non-linear iterations to perform per
// time-step. By default, only 1 step is performed, which yields a solver with semi-implicit
// treatment of the non-linear convection.
//
// The '--dpm-steps <N>' option specifies the number of DPM iterations to perform per
// non-linear itation. By default, only 1 step is perfomed.
//
//
// Linear Solver/Preconditioner Options
// ------------------------------------
// This application uses 2 multigrid solvers as linear preconditioners for the DPM:
//  1) a Richardson-Multigrid for the linearised Burgers system in velocity space (A-solver)
//  2) a PCG-Multigrid for the Poisson problem in pressure space (S-solver)
//
// Both multigrids use a damped Jacobi smoother as well as a Jacobi "smoother" as the coarse-grid solver.
//
// These two multigrid solvers are configured by the same set of options with slightly
// different names: options for the A-solver are postfixed by "-a", whereas options for
// the S-solver are postfixed by '-s'.
//
// The '--max-iter-[a|s] <N>' option sets the maximum allowed multigrid iterations.
// Default: 25 for A-solver, 50 for S-solver
//
// The '--tol-rel-[a|s] <eps>' options sets the relative tolerance for the multigrid.
// Default: 1E-5 for both A- and S-solver
//
// The '--smooth-[a|s] <N>' option sets the number of pre-/post-smoothing steps.
// Default: 4 for both A- and S-solver
//
// The '--damp-[a|s] <omega>' option sets the damping parameter for the smoother.
// Default: 0.7 for both A- and S-solver
//
// Furthermore, it is possible to use a simple (one-grid) damped Jacobi-Iteration instead
// of multigrid as the A-solver. This can be achieved by supplying the '--no-multigrid-a' option.
//
//
// VTK-Export Options
// ------------------
// The option '--vtk <name> [<step>]' can be used to export the solutions for the
// time-steps. The <name> parameter specifies the base filename of the exported (P)VTU files,
// which is postfixed with the number of ranks (if MPI is used) and the time-step index.
// The optional <step> parameter specifies the stepping of the VTU export, i.e. only
// those time steps are exported if the index of the time-step is a multiple of the <step>
// parameter. Example: the option '--vtk <name> 50' will write every 50-th time-step to
// a corresponding VTU file.
//
//
// \author Peter Zajac
//
#include <kernel/util/runtime.hpp>
#include <kernel/util/mpi_cout.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/jacobi_precond.hpp>

#include <control/domain/partitioner_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/statistics.hpp>

#include <deque>

namespace NaverStokesCP2D
{
  using namespace FEAT;

  // helper functions for padded console output
  static inline void dump_line(String s, String t)
  {
    Util::mpi_cout(s.pad_back(30, '.') + ": " + t + "\n");
  }

  template<typename T_>
  static inline void dump_line(String s, T_ t)
  {
    Util::mpi_cout(s.pad_back(30, '.') + ": " + stringify(t) + "\n");
  }

  static inline void dump_time(String s, double t, double total)
  {
    Util::mpi_cout(s.pad_back(30, '.') + ": " + stringify_fp_fix(t, 3, 10)
      + " (" + stringify_fp_fix(100.0*t/total,3,7) + "%)\n");
  }

  /**
   * \brief Configuration auxiliary class
   *
   * This class is responsible for storing the various application parameters
   * which are set by the user from the command line.
   */
  class Config
  {
  public:
    /// filename of the mesh file
    String mesh_file;
    /// path to the mesh directory
    String mesh_path;

    /// minimum and maximum levels (as configured)
    Index level_min_in, level_max_in;

    /// minimum and maximum levels (after partitioning)
    Index level_min, level_max;

    /// base-name of VTK files
    String vtk_name;
    /// stepping of VTK output
    Index vtk_step;

    // name of inflow mesh-part
    String part_name_in;

    // name of outflow mesh-part
    String part_name_out;

    // -------------------------------

    /// use deformation tensor?
    bool deformation;

    /// viscosity
    Real nu;

    /// inflow profile line segment coordinates
    Real ix0, iy0, ix1, iy1;

    /// maximum inflow velocity
    Real vmax;

    // -------------------------------

    /// maximum simulation time
    Real time_max;

    /// number of time-steps for the total simulation time
    Index time_steps;

    /// maximum number of time-steps to perform
    // this may be < time_steps to enforce premature stop
    Index max_time_steps;

    // -------------------------------

    // number of non-linear steps per time-step
    Index nonlin_steps;

    // number of linear DPM steps per non-linear step
    Index dpm_steps;

    // -------------------------------

    // use multigrid for A-solver ?
    bool multigrid_a;

    // maximum number of iterations for velocity mg
    Index max_iter_a;

    // relative tolerance for velocity mg
    Real tol_rel_a;

    // smoothing steps for velocity mg
    Index smooth_steps_a;

    // damping parameter for velocity smoother
    Real smooth_damp_a;

    // -------------------------------

    // maximum number of iterations for pressure mg
    Index max_iter_s;

    // relative tolerance for pressure mg
    Real tol_rel_s;

    // smoothing steps for pressure mg
    Index smooth_steps_s;

    // damping parameter for pressure smoother
    Real smooth_damp_s;

  public:
    Config() :
      level_min_in(0),
      level_max_in(0),
      level_min(0),
      level_max(0),
      vtk_step(0),
      deformation(false),
      nu(1.0),
      ix0(0.0),
      iy0(0.0),
      ix1(0.0),
      iy1(0.0),
      vmax(1.0),
      time_max(0.0),
      time_steps(0),
      max_time_steps(0),
      nonlin_steps(1),
      dpm_steps(1),
      multigrid_a(false),
      max_iter_a(25),
      tol_rel_a(1E-5),
      smooth_steps_a(4),
      smooth_damp_a(0.7),
      max_iter_s(50),
      tol_rel_s(1E-5),
      smooth_steps_s(4),
      smooth_damp_s(0.7)
    {
      const char* mpath = getenv("FEAT3_PATH_MESHES");
      if(mpath != nullptr)
        mesh_path = mpath;
    }

    bool parse_args(SimpleArgParser& args)
    {
      String s;
      if(args.parse("setup", s) > 0)
      {
        if(s.compare_no_case("square") == 0)
          setup_square();
        else if(s.compare_no_case("nozzle") == 0)
          setup_nozzle();
        else if(s.compare_no_case("bench1") == 0)
          setup_bench1();
        else
        {
          Util::mpi_cout("ERROR: unknown setup '" + s + "'");
          return false;
        }
      }
      deformation = (args.check("deformation") >= 0);
      args.parse("mesh-path", mesh_path);
      args.parse("mesh-file", mesh_file);
      if(args.parse("vtk", vtk_name, vtk_step) == 1)
        vtk_step = 1; // vtk-name given, but not vtk-step, so set to 1
      args.parse("level", level_max_in, level_min_in);
      level_max = level_max_in;
      level_min = level_min_in;
      args.parse("nu", nu);
      args.parse("part-in", part_name_in);
      args.parse("part-out", part_name_out);
      args.parse("profile", ix0, iy0, ix1, iy1, vmax);
      args.parse("time-max", time_max);
      args.parse("time-steps", time_steps);
      if(args.parse("max-time-steps", max_time_steps) < 1)
        max_time_steps = time_steps;
      args.parse("nl-steps", nonlin_steps);
      args.parse("dpm-steps", dpm_steps);
      multigrid_a = (args.check("no-multigrid-a") < 0);
      args.parse("max-iter-a", max_iter_a);
      args.parse("tol-rel-a", tol_rel_a);
      args.parse("smooth-a", smooth_steps_a);
      args.parse("damp-a", smooth_damp_a);
      args.parse("max-iter-s", max_iter_s);
      args.parse("tol-rel-s", tol_rel_s);
      args.parse("smooth-s", smooth_steps_s);
      args.parse("damp-s", smooth_damp_s);

      return true;
    }

    void dump()
    {
      Util::mpi_cout("Configuration Summary:\n");
      dump_line("Mesh File", mesh_file);
      dump_line("Mesh Path", mesh_path);
      dump_line("Level-Min", stringify(level_min) + " [" + stringify(level_min_in) + "]");
      dump_line("Level-Max", stringify(level_max) + " [" + stringify(level_max_in) + "]");
      dump_line("VTK-Name", vtk_name);
      dump_line("VTK-Step", vtk_step);
      dump_line("Inflow-Part", part_name_in);
      dump_line("Outflow-Part", part_name_out);
      dump_line("Inflow-Profile", "( " + stringify(ix0) + " , " + stringify(iy0) + " ) - ( "
        + stringify(ix1) + " , " + stringify(iy1) + " )");
      dump_line("V-Max", vmax);
      dump_line("Tensor", (deformation ? "Deformation" : "Gradient"));
      dump_line("Nu", nu);
      dump_line("Time-Max", time_max);
      dump_line("Time-Steps", time_steps);
      dump_line("Max Time-Steps", max_time_steps);
      dump_line("Non-Linear Steps", nonlin_steps);
      dump_line("Linear DPM Steps", dpm_steps);
      dump_line("A: Solver", (multigrid_a ? "Multigrid" : "Jacobi"));
      dump_line("A: Max-Iter", max_iter_a);
      dump_line("A: Tol-Rel", tol_rel_a);
      dump_line("A: Smooth Steps", smooth_steps_a);
      dump_line("A: Smooth Damp", smooth_damp_a);
      dump_line("S: Max-Iter", max_iter_s);
      dump_line("S: Tol-Rel", tol_rel_s);
      dump_line("S: Smooth Steps", smooth_steps_s);
      dump_line("S: Smooth Damp", smooth_damp_s);
    }

    // Setup: Poiseuille-Flow on unit-square
    void setup_square()
    {
      mesh_file = "unit-square-quad.xml";
      part_name_in = "left";
      part_name_out = "right";
      level_min = level_min_in = Index(0);
      level_max = level_max_in = Index(7);
      nu = 1E-3;
      ix0 = 0.0;
      iy0 = 0.0;
      ix1 = 0.0;
      iy1 = 1.0;
      vmax = 1.0;
      time_max = 3.0;
      time_steps = max_time_steps = 1200;
    }

    // Setup: nozzle-jet simulation
    void setup_nozzle()
    {
      mesh_file = "nozzle-2-quad.xml";
      part_name_in = "bnd:l";
      part_name_out = "bnd:r";
      level_min = level_min_in = Index(0);
      level_max = level_max_in = Index(6);
      nu = 1E-3;
      ix0 = 0.0;
      iy0 = -0.5;
      ix1 = 0.0;
      iy1 = 0.5;
      vmax = 1.0;
      time_max = 7.0;
      time_steps = max_time_steps = 3500;
    }

    // Setup: flow around a cylinder
    void setup_bench1()
    {
      mesh_file = "bench1-quad.xml";
      part_name_in = "left";
      part_name_out = "right";
      level_min = level_min_in = Index(0);
      level_max = level_max_in = Index(4);
      nu = 1E-3;
      ix0 = 0.0;
      iy0 = 0.0;
      ix1 = 0.0;
      iy1 = 0.41;
      vmax = 1.5;
      time_max = 3.0;
      time_steps = max_time_steps = 4500;
    }
  }; // class Config

  /**
   * \brief Navier-Stokes System Level class
   *
   * This extends the StokesBlockedSystemLevel by the corresponding filters for
   * the velocity and pressure sub-systems.
   */
  template<
    int dim_,
    typename MemType_ = Mem::Main,
    typename DataType_ = Real,
    typename IndexType_ = Index,
    typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
    typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
    typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
    typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>>
  class NavierStokesBlockedSystemLevel :
    public Control::StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
  {
  public:
    typedef Control::StokesBlockedSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;

    // define local filter types
    typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalVeloFilter;
    typedef LAFEM::NoneFilter<MemType_, DataType_, IndexType_> LocalPresNoneFilter;
    typedef LAFEM::UnitFilter<MemType_, DataType_, IndexType_> LocalPresUnitFilter;
    typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresNoneFilter> LocalSystemFilter;

    // define global filter types
    typedef Global::Filter<LocalVeloFilter, typename BaseClass::VeloMirror> GlobalVeloFilter;
    typedef Global::Filter<LocalPresNoneFilter, typename BaseClass::PresMirror> GlobalPresNoneFilter;
    typedef Global::Filter<LocalPresUnitFilter, typename BaseClass::PresMirror> GlobalPresUnitFilter;
    typedef Global::Filter<LocalSystemFilter, typename BaseClass::SystemMirror> GlobalSystemFilter;

    // (global) filters
    GlobalSystemFilter filter_sys;
    GlobalVeloFilter filter_velo;
    GlobalPresUnitFilter filter_pres_unit;
  }; // class NavierStokesBlockedSystemLevel

  /**
   * \brief Navier-Stokes Assembler Level class
   *
   * This extends the StokesBlockedAssemblerLevel by the assembly of the filters
   * as well as the burgers matrix and the RHS vector.
   */
  template<
    typename SpaceVelo_,
    typename SpacePres_>
  class NavierStokesBlockedAssemblerLevel :
    public Control::StokesBlockedAssemblerLevel<SpaceVelo_, SpacePres_>
  {
  public:
    typedef Control::StokesBlockedAssemblerLevel<SpaceVelo_, SpacePres_> BaseClass;
    typedef typename BaseClass::MeshType MeshType;

  public:
    explicit NavierStokesBlockedAssemblerLevel(typename BaseClass::DomainLevelType& dom_lvl) :
      BaseClass(dom_lvl)
    {
    }

    template<typename VeloSystemLevel_>
    void assemble_velo_filter(const Config& cfg, VeloSystemLevel_& velo_sys_level)
    {
      // get our global system filter
      typename VeloSystemLevel_::GlobalVeloFilter& fil_glob = velo_sys_level.filter_velo;

      // get our local system filter
      typename VeloSystemLevel_::LocalVeloFilter& fil_loc = *fil_glob;

      // create unit-filter assemblers
      Assembly::UnitFilterAssembler<MeshType> unit_asm, unit_asm_inflow;
      bool have_inflow(false);

      // loop over all boundary parts
      std::deque<String> part_names = this->domain_level.get_mesh_node()->get_mesh_part_names();
      for(const auto& name : part_names)
      {
        // skip internal meshparts
        if(name.starts_with('_'))
          continue;

        // skip outflow part
        if(name == cfg.part_name_out)
          continue;

        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = this->domain_level.get_mesh_node()->find_mesh_part_node(name);

        // found it?
        XASSERT(mesh_part_node != nullptr);

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        auto* mesh_part = mesh_part_node->get_mesh();
        if (mesh_part != nullptr)
        {
          // add to boundary assembler
          if(name == cfg.part_name_in)
          {
            unit_asm_inflow.add_mesh_part(*mesh_part);
            have_inflow = true;
          }
          else
          {
            unit_asm.add_mesh_part(*mesh_part);
          }
        }
      }

      // assemble the filter
      unit_asm.assemble(fil_loc, this->space_velo);

      if(!have_inflow)
        return;

      // create parabolic inflow profile
      Analytic::Common::ParProfileVector inflow(cfg.ix0, cfg.iy0, cfg.ix1, cfg.iy1, cfg.vmax);

      // assemble inflow BC
      unit_asm_inflow.assemble(fil_loc, this->space_velo, inflow);
    }

    template<typename PresSystemLevel_>
    void assemble_pres_filter(const Config& cfg, PresSystemLevel_& pres_sys_level)
    {
      // get our global system filter
      typename PresSystemLevel_::GlobalPresUnitFilter& fil_glob = pres_sys_level.filter_pres_unit;

      // get our local system filter
      typename PresSystemLevel_::LocalPresUnitFilter& fil_loc = *fil_glob;

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> unit_asm;

      // try to fetch the corresponding mesh part node
      auto* mesh_part_node = this->domain_level.get_mesh_node()->find_mesh_part_node(cfg.part_name_out);

      // found it?
      XASSERT(mesh_part_node != nullptr);

      // let's see if we have that mesh part
      // if it is nullptr, then our patch is not adjacent to that boundary part
      auto* mesh_part = mesh_part_node->get_mesh();
      if (mesh_part != nullptr)
      {
        unit_asm.add_mesh_part(*mesh_part);
      }

      // assemble the filter
      unit_asm.assemble(fil_loc, this->space_pres);
    }

    template<typename GlobalVeloVector_>
    void assemble_rhs_vector(const Config& cfg, const Real delta_t, GlobalVeloVector_& vec_rhs_v, const GlobalVeloVector_& vec_sol_v)
    {
      typedef typename GlobalVeloVector_::DataType DataType;
      typedef typename GlobalVeloVector_::IndexType IndexType;
      Assembly::BurgersAssembler<DataType, IndexType, 2> burgers_rhs;
      burgers_rhs.deformation = cfg.deformation;
      burgers_rhs.nu = -cfg.nu;
      burgers_rhs.beta = -DataType(1);
      burgers_rhs.theta = DataType(1) / delta_t;

      // assemble RHS vector
      (*vec_rhs_v).format();
      burgers_rhs.assemble(this->space_velo, this->cubature, (*vec_sol_v), nullptr, &(*vec_rhs_v));

      // synchronise RHS vector
      vec_rhs_v.sync_0();
    }

    template<typename GlobalMatrixBlockA_, typename GlobalVeloVector_>
    void assemble_burgers_matrix(const Config& cfg, const Real delta_t, GlobalMatrixBlockA_& matrix_a, const GlobalVeloVector_& vec_conv)
    {
      typedef typename GlobalVeloVector_::DataType DataType;
      typedef typename GlobalVeloVector_::IndexType IndexType;
      Assembly::BurgersAssembler<DataType, IndexType, 2> burgers_mat;
      burgers_mat.deformation = cfg.deformation;
      burgers_mat.nu = cfg.nu;
      burgers_mat.beta = DataType(1);
      burgers_mat.theta = DataType(1) / delta_t;

      // "restrict" our convection vector onto that level;
      // this exploits the 2-level numbering of the Q2 convection vector
      typename GlobalVeloVector_::LocalVectorType vec_cv(*vec_conv, (*matrix_a).rows(), IndexType(0));

      // format and assemble the matrix
      (*matrix_a).format();
      burgers_mat.assemble(this->space_velo, this->cubature, vec_cv, &(*matrix_a), nullptr);
    }
  }; // class NavierStokesBlockedAssemblerLevel


  template<typename MeshType_>
  void run(const int rank, const int nprocs, const Config& cfg, Control::Domain::DomainControl<MeshType_>& domain)
  {
    // create a time-stamp
    TimeStamp stamp_start;

    // define our mesh type
    typedef MeshType_ MeshType;

    // our dimension
    static constexpr int dim = MeshType::shape_dim;

    // our arch types
    typedef Mem::Main MemType;
    typedef Real DataType;
    typedef Index IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<MeshType_> DomainControlType;

    // define our velocity and pressure system levels
    typedef NavierStokesBlockedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    // define our transfer levels
    typedef Control::StokesBlockedTransferLevel<SystemLevelType> TransferLevelType;

    // define our trafo and FE spaces
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> VeloSpaceType;
    typedef Space::Lagrange1::Element<TrafoType> PresSpaceType;

    // define our assembler level
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef NavierStokesBlockedAssemblerLevel<VeloSpaceType, PresSpaceType> AssemblerLevelType;

    // get our domain level and layer
    typedef typename DomainControlType::LayerType DomainLayerType;
    const DomainLayerType& layer = *domain.get_layers().back();
    const std::deque<DomainLevelType*>& domain_levels = domain.get_levels();

    std::deque<AssemblerLevelType*> asm_levels;
    std::deque<SystemLevelType*> system_levels;
    std::deque<TransferLevelType*> transfer_levels;

    const Index num_levels = Index(domain_levels.size());

    // create a batch of stop-watches
    StopWatch watch_total, watch_asm_rhs, watch_asm_mat, watch_calc_def,
      watch_sol_init, watch_solver_a, watch_solver_s, watch_vtk;

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

    Util::mpi_cout("Creating gates...\n");

    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_gates(layer, *system_levels.at(i));
    }

    /* ***************************************************************************************** */

    Util::mpi_cout("Assembling basic matrices...\n");
    for(Index i(0); i < num_levels; ++i)
    {
      // assemble velocity matrix structure
      asm_levels.at(i)->assemble_velo_struct(*system_levels.at(i));
      // assemble pressure laplace matrix
      asm_levels.at(i)->assemble_pres_laplace(*system_levels.at(i));
    }

    // assemble B/D matrices on finest level
    asm_levels.back()->assemble_grad_div_matrices(*system_levels.back());

    /* ***************************************************************************************** */

    Util::mpi_cout("Assembling system filters...\n");
    for(Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_velo_filter(cfg, *system_levels.at(i));
      asm_levels.at(i)->assemble_pres_filter(cfg, *system_levels.at(i));
    }

    /* ***************************************************************************************** */

    Util::mpi_cout("Assembling transfer matrices...\n");

    for (Index i(0); (i + 1) < num_levels; ++i)
    {
      asm_levels.at(i+1)->assemble_velo_transfer(*transfer_levels.at(i), *asm_levels.at(i));
      asm_levels.at(i+1)->assemble_pres_transfer(*transfer_levels.at(i), *asm_levels.at(i));
    }

    /* ***************************************************************************************** */

    // get our vector types
    typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain_levels.back();
    SystemLevelType& the_system_level = *system_levels.back();
    AssemblerLevelType& the_asm_level = *asm_levels.back();

    // get our fine-level matrices
    typename SystemLevelType::GlobalMatrixBlockA& matrix_a = the_system_level.matrix_a;
    typename SystemLevelType::GlobalMatrixBlockB& matrix_b = the_system_level.matrix_b;
    typename SystemLevelType::GlobalMatrixBlockD& matrix_d = the_system_level.matrix_d;
    typename SystemLevelType::GlobalSchurMatrix& matrix_s = the_system_level.matrix_s;

    // get out fine-level filters
    typename SystemLevelType::GlobalVeloFilter& filter_v = the_system_level.filter_velo;
    typename SystemLevelType::GlobalPresUnitFilter& filter_p = the_system_level.filter_pres_unit;

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    Util::mpi_cout("Setting up Velocity Multigrid...\n");

    // create a multigrid solver for the velocity
    auto multigrid_hierarchy_velo = std::make_shared<Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalMatrixBlockA,
      typename SystemLevelType::GlobalVeloFilter,
      typename TransferLevelType::GlobalVeloTransferMatrix,
      typename TransferLevelType::GlobalVeloTransferMatrix
      >>();

    // use multigrid for A-solver?
    if(cfg.multigrid_a)
    {
      // set up velocity coarse grid solver
      {
        auto& lvl = *system_levels.front();
        auto jac = Solver::new_jacobi_precond(lvl.matrix_a, lvl.filter_velo);
        auto cgs = Solver::new_richardson(lvl.matrix_a, lvl.filter_velo, cfg.smooth_damp_a, jac);
        cgs->set_max_iter(cfg.smooth_steps_a);
        cgs->set_min_iter(cfg.smooth_steps_a);

        multigrid_hierarchy_velo->push_level(lvl.matrix_a, lvl.filter_velo, cgs);
      }

      // loop over all finer levels
      for(Index i(1); i < num_levels; ++i)
      {
        auto& lvl = *system_levels.at(std::size_t(i));
        auto& trs = *transfer_levels.at(std::size_t(i-1));
        auto jac = Solver::new_jacobi_precond(lvl.matrix_a, lvl.filter_velo);
        auto smoother = Solver::new_richardson(lvl.matrix_a, lvl.filter_velo, cfg.smooth_damp_a, jac);
        smoother->set_max_iter(cfg.smooth_steps_a);
        smoother->set_min_iter(cfg.smooth_steps_a);
        multigrid_hierarchy_velo->push_level(
          lvl.matrix_a, lvl.filter_velo,
          trs.prol_velo, trs.rest_velo,
          smoother, smoother, smoother);
      }
    }

    /* ***************************************************************************************** */

    Util::mpi_cout("Setting up Pressure Multigrid...\n");

    // create another multigrid solver for the pressure
    auto multigrid_hierarchy_pres = std::make_shared<Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalSchurMatrix,
      typename SystemLevelType::GlobalPresUnitFilter,
      typename TransferLevelType::GlobalPresTransferMatrix,
      typename TransferLevelType::GlobalPresTransferMatrix
      >>();

    // set up pressure coarse grid solver
    {
      auto& lvl = *system_levels.front();
      auto jac = Solver::new_jacobi_precond(lvl.matrix_s, lvl.filter_pres_unit);
      auto cgs = Solver::new_richardson(lvl.matrix_s, lvl.filter_pres_unit, cfg.smooth_damp_s, jac);
      cgs->set_max_iter(cfg.smooth_steps_s);
      cgs->set_min_iter(cfg.smooth_steps_s);

      multigrid_hierarchy_pres->push_level(lvl.matrix_s, lvl.filter_pres_unit, cgs);
    }

    // loop over all finer levels
    for(Index i(1); i < num_levels; ++i)
    {
      auto& lvl = *system_levels.at(std::size_t(i));
      auto& trs = *transfer_levels.at(std::size_t(i-1));
      auto jac = Solver::new_jacobi_precond(lvl.matrix_s, lvl.filter_pres_unit);
      auto smoother = Solver::new_richardson(lvl.matrix_s, lvl.filter_pres_unit, cfg.smooth_damp_s, jac);
      smoother->set_max_iter(cfg.smooth_steps_s);
      smoother->set_min_iter(cfg.smooth_steps_s);
      multigrid_hierarchy_pres->push_level(
        lvl.matrix_s, lvl.filter_pres_unit,
        trs.prol_pres, trs.rest_pres,
        smoother, smoother, smoother);
    }

    /* ***************************************************************************************** */

    // create a Richardson-MG for the velocity
    std::shared_ptr<Solver::IterativeSolver<GlobalVeloVector>> solver_a;

    if(cfg.multigrid_a)
    {
      // use Multigrid
      auto multigrid_velo = Solver::new_multigrid(multigrid_hierarchy_velo);
      solver_a = Solver::new_richardson(matrix_a, filter_v, 1.0, multigrid_velo);
    }
    else
    {
      // use Richardson-Jacobi
      auto jac = Solver::new_jacobi_precond(the_system_level.matrix_a, the_system_level.filter_velo);
      solver_a = Solver::new_richardson(matrix_a, filter_v, cfg.smooth_damp_a, jac);
    }

    solver_a->set_max_iter(cfg.max_iter_a);
    solver_a->set_tol_rel(cfg.tol_rel_a);

    // for the velocity multigrid/solver, we can only perform symbolic initialisation up to now:
    if(cfg.multigrid_a)
    {
      multigrid_hierarchy_velo->init_branch();
      multigrid_hierarchy_velo->init_symbolic();
    }
    solver_a->init_branch();
    solver_a->init_symbolic();

    /* ***************************************************************************************** */

    std::shared_ptr<Solver::IterativeSolver<GlobalPresVector>> solver_s;

    // create a PCG-MG for the pressure
    auto multigrid_pres = Solver::new_multigrid(multigrid_hierarchy_pres);
    solver_s = Solver::new_pcg(matrix_s, filter_p, multigrid_pres);

    solver_s->set_max_iter(cfg.max_iter_s);
    solver_s->set_tol_rel(cfg.tol_rel_s);

    // for the pressure multigrid, we can perform full initialisation:
    multigrid_hierarchy_pres->init();
    solver_s->init();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    Util::mpi_cout("\n\n");

    // create RHS and SOL vectors
    GlobalVeloVector vec_sol_v = matrix_a.create_vector_l();
    GlobalPresVector vec_sol_p = matrix_s.create_vector_l();
    GlobalVeloVector vec_rhs_v = matrix_a.create_vector_l();
    GlobalPresVector vec_rhs_p = matrix_s.create_vector_l();

    // create defect and correction vectors
    GlobalVeloVector vec_def_v = matrix_a.create_vector_l();
    GlobalPresVector vec_def_p = matrix_s.create_vector_l();
    GlobalVeloVector vec_cor_v = matrix_a.create_vector_l();
    GlobalPresVector vec_cor_p = matrix_s.create_vector_l();
    // create convection vector
    GlobalVeloVector vec_conv = matrix_a.create_vector_l();

    // format the vectors
    vec_sol_v.format();
    vec_sol_p.format();
    vec_rhs_v.format();
    vec_rhs_p.format();

    // apply filter onto solution vector
    filter_v.filter_sol(vec_sol_v);

    // create solution backup vectors; these store vec_sol_v/p of the last two time-steps
    GlobalVeloVector vec_sol_v_1 = vec_sol_v.clone();
    GlobalVeloVector vec_sol_v_2 = vec_sol_v.clone();
    GlobalPresVector vec_sol_p_1 = vec_sol_p.clone();

    // write header line to console
    if(rank == 0)
    {
      const std::size_t nf = stringify_fp_sci(0.0, 3).size();
      String head;
      head += String("Step").pad_front(6) + "  ";
      head += String("Time").pad_back(8) + " ";
      head += String("NL").pad_front(3) + "   ";
      head += String("Def-V").pad_back(nf) + " ";
      head += String("Def-P").pad_back(nf) + "   ";
      head += String("Def-V").pad_back(nf) + " ";
      head += String("Def-P").pad_back(nf) + "   ";
      head += String("IT-A").pad_front(4) + " ";
      head += String("IT-S").pad_front(4) + "   ";
      head += String("Runtime    ");
      std::cout << head << std::endl;
      std::cout << String(head.size(), '-') << std::endl;
    }

    watch_total.start();

    // compute time-step size
    const DataType delta_t = cfg.time_max / DataType(cfg.time_steps);

    // keep track whether something failed miserably...
    bool failure = false;

    // time-step loop
    for(Index time_step(1); time_step <= cfg.max_time_steps; ++time_step)
    {
      // compute current time
      const DataType cur_time = DataType(time_step) * delta_t;

      // assemble RHS vector
      watch_asm_rhs.start();
      vec_rhs_v.format();
      vec_rhs_p.format();
      the_asm_level.assemble_rhs_vector(cfg, delta_t, vec_rhs_v, vec_sol_v);
      // subtract pressure (?)
      //matrix_b.apply(vec_rhs_v, vec_sol_p, vec_rhs_v, -DataType(1));
      watch_asm_rhs.stop();

      // apply RHS filter
      filter_v.filter_rhs(vec_rhs_v);

      // non-linear loop
      for(Index nonlin_step(0); nonlin_step < cfg.nonlin_steps; ++nonlin_step)
      {
        // Phase 1: compute convection vector
        // extrapolate previous time-step solution in first NL step
        if((time_step > Index(1)) && (nonlin_step == Index(0)))
        {
          // linear extrapolation of solution in time
          vec_conv.scale(vec_sol_v_1, DataType(2));
          vec_conv.axpy(vec_sol_v_2, vec_conv, -DataType(1));
        }
        else
        {
          // constant extrapolation of solution in time
          vec_conv.copy(vec_sol_v);
        }

        // Phase 2: loop over all levels and assemble the burgers matrices
        watch_asm_mat.start();
        if(cfg.multigrid_a)
        {
          // assemble burgers matrices on all levels
          for(std::size_t i(0); i < asm_levels.size(); ++i)
          {
            asm_levels.at(i)->assemble_burgers_matrix(cfg, delta_t, system_levels.at(i)->matrix_a, vec_conv);
          }
        }
        else
        {
          // assemble burgers matrices on finest level
          the_asm_level.assemble_burgers_matrix(cfg, delta_t, the_system_level.matrix_a, vec_conv);
        }
        watch_asm_mat.stop();

        // Phase 3: compute non-linear defects
        watch_calc_def.start();
        matrix_a.apply(vec_def_v, vec_sol_v, vec_rhs_v, -DataType(1));
        matrix_b.apply(vec_def_v, vec_sol_p, vec_def_v, -DataType(1));
        matrix_d.apply(vec_def_p, vec_sol_v, vec_rhs_p, -DataType(1));
        filter_v.filter_def(vec_def_v);

        // compute defect norms
        const DataType def_nl1_v = vec_def_v.norm2();
        const DataType def_nl1_p = vec_def_p.norm2();
        watch_calc_def.stop();

        // console output, part 1
        if(rank == 0)
        {
          String line;
          line += stringify(time_step).pad_front(6) + " ";
          line += stringify_fp_fix(cur_time, 5, 8) + " ";
          line += stringify(nonlin_step).pad_front(4);
          line += " : ";
          line += stringify_fp_sci(def_nl1_v, 3) + " ";
          line += stringify_fp_sci(def_nl1_p, 3);
          line += " > ";
          std::cout << line;
        }

        // Phase 4: initialise linear solvers
        watch_sol_init.start();
        if(cfg.multigrid_a)
          multigrid_hierarchy_velo->init_numeric();
        solver_a->init_numeric();
        watch_sol_init.stop();

        // linear solver iterations counts
        Index iter_v(0), iter_p(0);

        // Phase 5: linear DPM loop
        // Note: we need to perform one additional velocity solve,
        // so the break condition of the loop is inside...
        for(Index dpm_step(0); ; ++dpm_step)
        {
          // solve velocity system
          watch_solver_a.start();
          Solver::Status status_a = solver_a->apply(vec_cor_v, vec_def_v);
          watch_solver_a.stop();
          if(!Solver::status_success(status_a))
          {
            Util::mpi_cout("\n\nERROR: velocity solver broke down!\n\n");
            failure = true;
            break;
          }
          iter_v += solver_a->get_num_iter();

          // update velocity solution
          vec_sol_v.axpy(vec_cor_v, vec_sol_v);

          // are we done yet?
          if(dpm_step >= cfg.dpm_steps)
            break;

          // update pressure defect
          watch_calc_def.start();
          matrix_d.apply(vec_def_p, vec_sol_v, vec_rhs_p, -DataType(1));
          filter_p.filter_def(vec_def_p);
          watch_calc_def.stop();

          // solve pressure system
          watch_solver_s.start();
          Solver::Status status_s = solver_s->apply(vec_cor_p, vec_def_p);
          watch_solver_s.stop();
          if(!Solver::status_success(status_s))
          {
            Util::mpi_cout("\n\nERROR: pressure solver broke down!\n\n");
            failure = true;
            break;
          }
          iter_p += solver_s->get_num_iter();

          // update pressure solution
          vec_sol_p.axpy(vec_cor_p, vec_sol_p, -DataType(1) / delta_t);

          // compute new defect
          watch_calc_def.start();
          matrix_a.apply(vec_def_v, vec_sol_v, vec_rhs_v, -DataType(1));
          matrix_b.apply(vec_def_v, vec_sol_p, vec_def_v, -DataType(1));
          filter_v.filter_def(vec_def_v);
          watch_calc_def.stop();
        } // inner Uzawa loop

        // Phase 6: release linear solvers
        solver_a->done_numeric();
        if(cfg.multigrid_a)
          multigrid_hierarchy_velo->done_numeric();

        // epic fail?
        if(failure)
          break;

        // Phase 7: compute final defect and norms (only for console output)
        watch_calc_def.start();
        matrix_a.apply(vec_def_v, vec_sol_v, vec_rhs_v, -DataType(1));
        matrix_b.apply(vec_def_v, vec_sol_p, vec_def_v, -DataType(1));
        matrix_d.apply(vec_def_p, vec_sol_v, vec_rhs_p, -DataType(1));
        filter_v.filter_def(vec_def_v);

        const DataType def_nl2_v = vec_def_v.norm2();
        const DataType def_nl2_p = vec_def_p.norm2();
        watch_calc_def.stop();

        // console output, part 2
        if(rank == 0)
        {
          String line;
          line += stringify_fp_sci(def_nl2_v, 3) + " ";
          line += stringify_fp_sci(def_nl2_p, 3);
          line += " | ";
          line += stringify(iter_v).pad_front(4) + " ";
          line += stringify(iter_p).pad_front(4) + " | ";
          line += stamp_start.elapsed_string_now();
          std::cout << line << std::endl;
        }
      } // non-linear loop

      // epic fail?
      if(failure)
        break;

      // VTK-Export
      if(!cfg.vtk_name.empty() && (cfg.vtk_step > 0) && (time_step % cfg.vtk_step == 0))
      {
        watch_vtk.start();
        String vtk_path = cfg.vtk_name + "." + stringify(nprocs) + "." + stringify(time_step).pad_front(5, '0');

        Geometry::ExportVTK<MeshType> vtk(the_domain_level.get_mesh());

        // write solution
        vtk.add_vertex_vector("v", (*vec_sol_v));
        vtk.add_vertex_scalar("p", (*vec_sol_p).elements());

        // compute and write time-derivatives
        GlobalVeloVector vec_der_v = vec_sol_v.clone();
        GlobalPresVector vec_der_p = vec_sol_p.clone();
        vec_der_v.axpy(vec_sol_v_1, vec_der_v, -DataType(1));
        vec_der_p.axpy(vec_sol_p_1, vec_der_p, -DataType(1));
        vec_der_v.scale(vec_der_v, DataType(1) / delta_t);
        vec_der_p.scale(vec_der_p, DataType(1) / delta_t);

        vtk.add_vertex_vector("v_dt", (*vec_der_v));
        vtk.add_vertex_scalar("p_dt", (*vec_der_p).elements());

        // export
        vtk.write(vtk_path, rank, nprocs);
        watch_vtk.stop();
      }

      // finally, update our solution vector backups
      vec_sol_v_2.copy(vec_sol_v_1);
      vec_sol_v_1.copy(vec_sol_v);
      vec_sol_p_1.copy(vec_sol_p);

      // continue with next time-step
    } // time-step loop

    watch_total.stop();

    // release pressure solvers
    solver_s->done();
    multigrid_hierarchy_pres->done();

    // release velocity solvers
    solver_a->done_symbolic();
    if(cfg.multigrid_a)
      multigrid_hierarchy_velo->done_symbolic();

    // write timings
    if(rank == 0)
    {
      double t_total = watch_total.elapsed();
      double t_asm_mat = watch_asm_mat.elapsed();
      double t_asm_rhs = watch_asm_rhs.elapsed();
      double t_calc_def = watch_calc_def.elapsed();
      double t_sol_init = watch_sol_init.elapsed();
      double t_solver_a = watch_solver_a.elapsed();
      double t_solver_s = watch_solver_s.elapsed();
      double t_vtk = watch_vtk.elapsed();
      double t_sum = t_asm_mat + t_asm_rhs + t_calc_def + t_sol_init + t_solver_a + t_solver_s + t_vtk;
      Util::mpi_cout("\n\n");
      dump_time("Total Solver Time", t_total, t_total);
      dump_time("Matrix Assembly Time", t_asm_mat, t_total);
      dump_time("Vector Assembly Time", t_asm_rhs, t_total);
      dump_time("Defect-Calc Time", t_calc_def, t_total);
      dump_time("Solver-A Init Time", t_sol_init, t_total);
      dump_time("Solver-A Time", t_solver_a, t_total);
      dump_time("Solver-S Time", t_solver_s, t_total);
      dump_time("VTK-Write Time", t_vtk, t_total);
      dump_time("Other Time", t_total-t_sum, t_total);
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

  int main(int argc, char* argv [])
  {
    int rank = 0;
    int nprocs = 0;

    // initialise
    FEAT::Runtime::initialise(argc, argv, rank, nprocs);
#ifdef FEAT_HAVE_MPI
    Util::mpi_cout("NUM-PROCS: " + stringify(nprocs) + "\n");
#endif

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    args.support("help", "\nDisplays this help message.\n");
    args.support("setup", "<config>\nLoads a pre-defined configuration:\n"
      "square    Poiseuille-Flow on Unit-Square\n"
      "nozzle    Jet-Flow through Nozzle domain\n"
      "bench1    Nonsteady Flow Around A Cylinder\n"
      );
    args.support("deformation", "\nUse deformation tensor instead of gradient tensor.\n");
    args.support("nu <nu>", "\nSets the viscosity parameter.\n");
    args.support("time-max", "<T_max>\nSets the maximum simulation time T_max.\n");
    args.support("time-steps", "<N>\nSets the number of time-steps for the time interval.\n");
    args.support("max-time-steps", "<N>\nSets the maximum number of time-steps to perform.\n");
    args.support("part-in", "<name>\nSpecifies the name of the inflow mesh-part.\n");
    args.support("part-out", "<name>\nSpecifies the name of the outflow mesh-part.\n");
    args.support("profile", "<x0> <y0> <x1> <y1>\nSpecifies the line segment coordinates for the inflow profile.\n");
    args.support("level", "<max> [<min>]\nSets the maximum and minimum mesh refinement levels.\n");
    args.support("vtk", "<name> [<step>]\nSets the name for VTK output and the time-stepping for the output (optional).\n");
    args.support("rank-elems", "<n>\nSpecifies the minimum number of elements per rank.\nDefault: 4\n");
    args.support("mesh-file", "<name>\nSpecifies the filename of the input mesh file.\n");
    args.support("mesh-path", "<path>\nSpecifies the path of the directory containing the mesh file.\n");
    args.support("nl-steps", "<N>\nSets the number of non-linear iterations per time-step.\nDefault: 1\n");
    args.support("dpm-steps", "<N>\nSets the number of Discrete-Projection-Method steps per non-linear step.\nDefault: 1\n");
    args.support("no-multigrid-a", "\nUse Jacobi instead of Multigrid as A-Solver.\n");
    args.support("max-iter-a", "<N>\nSets the maximum number of allowed iterations for the A-Solver.\nDefault: 25\n");
    args.support("tol-rel-a", "<eps>\nSets the relative tolerative for the A-Solver.\nDefault: 1E-5\n");
    args.support("smooth-a", "<N>\nSets the number of smoothing steps for the A-Solver.\nDefault: 4\n");
    args.support("damp-a", "<omega>\nSets the smoother daming parameter for the A-Solver.\nDefault: 0.7\n");
    args.support("max-iter-s", "<N>\nSets the maximum number of allowed iterations for the S-Solver.\nDefault: 50\n");
    args.support("tol-rel-s", "<eps>\nSets the relative tolerative for the S-Solver.\nDefault: 1E-5\n");
    args.support("smooth-s", "<N>\nSets the number of smoothing steps for the S-Solver.\nDefault: 4\n");
    args.support("damp-s", "<omega>\nSets the smoother daming parameter for the S-Solver.\nDefault: 0.7\n");

    // no arguments given?
    if((argc <= 1) || (args.check("help") >= 0))
    {
      Util::mpi_cout("\n2D Nonsteady Navier-Stokes CP-Q2/Q1 Toycode Solver (TM)\n\n");
      Util::mpi_cout("The easiest way to make this application do something useful is\n");
      Util::mpi_cout("to load a pre-defined problem configuration by supplying the\n");
      Util::mpi_cout("option '--setup <config>', where <config> may be one of:\n\n");
      Util::mpi_cout("  square    Poiseuille-Flow on Unit-Square\n");
      Util::mpi_cout("  nozzle    Jet-Flow through Nozzle domain\n");
      Util::mpi_cout("  bench1    Nonsteady Flow Around A Cylinder\n");
      Util::mpi_cout("\n");
      Util::mpi_cout("This will pre-configure this application to solve one of the\n");
      Util::mpi_cout("above problems. Note that you can further adjust the configration\n");
      Util::mpi_cout("by specifying additional options to override the default problem\n");
      Util::mpi_cout("configuration.\n");
      if(args.check("help") >= 0)
      {
        Util::mpi_cout("\nSupported Options:\n");
        Util::mpi_cout(args.get_supported_help());
      }
      else
      {
        Util::mpi_cout("\nUse the option '--help' to display a list of all supported options.\n\n");
      }
      return Runtime::finalise();
    }

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      if(rank == 0)
      {
        for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
          std::cerr << "ERROR: unknown option '--" << (*it).second << "'" << std::endl;

        std::cerr << "Supported options are:" << std::endl;
        std::cerr << args.get_supported_help() << std::endl;
      }
      // abort
      FEAT::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;

    // parse our configuration
    Config cfg;
    if(!cfg.parse_args(args))
      FEAT::Runtime::abort();

    // build the mesh file path
    String meshfile = cfg.mesh_path + "/" + cfg.mesh_file;
    int lvl_min = int(cfg.level_min_in);
    int lvl_max = int(cfg.level_max_in);

#ifndef DEBUG
    try
#endif
    {
      TimeStamp stamp1;

      // let's create our domain
      Util::mpi_cout("\nPreparing domain...\n");
      int rank_elems(4);
      args.parse("rank-elems", rank_elems);
      Index min_elems_partitioner = Index(nprocs * rank_elems);

      TimeStamp stamp_partition;
#ifdef FEAT_HAVE_PARMETIS
      Control::Domain::PartitionerDomainControl<Foundation::PExecutorParmetis<Foundation::ParmetisModePartKway>, MeshType> domain(lvl_max, lvl_min, min_elems_partitioner, meshfile);
#elif defined(FEAT_HAVE_MPI)
      Control::Domain::PartitionerDomainControl<Foundation::PExecutorFallback<double, Index>, MeshType> domain(lvl_max, lvl_min, min_elems_partitioner, meshfile);
#else
      Control::Domain::PartitionerDomainControl<Foundation::PExecutorNONE<double, Index>, MeshType> domain(lvl_max, lvl_min, min_elems_partitioner, meshfile);
#endif
      Statistics::toe_partition = stamp_partition.elapsed_now();

      // store levels after partitioning
      cfg.level_max = Index(domain.get_levels().back()->get_level_index());
      cfg.level_min = Index(domain.get_levels().front()->get_level_index());

      // dump our configuration
      cfg.dump();

      // run our application
      run<MeshType>(rank, nprocs, cfg, domain);

      TimeStamp stamp2;

      // get times
      long long time1 = stamp2.elapsed_micros(stamp1);

      // accumulate times over all processes
      long long time2 = time1 * (long long) nprocs;

      // print time
      Util::mpi_cout("Run-Time: " + stringify(TimeStamp::format_micros(time1, TimeFormat::m_s_m)) + " [" +
        stringify(TimeStamp::format_micros(time2, TimeFormat::m_s_m)) + "]\n");
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

    // okay
    return FEAT::Runtime::finalise();
  }
} // namespace NaverStokesCP2D


int main(int argc, char* argv [])
{
  return NaverStokesCP2D::main(argc, argv);
}
