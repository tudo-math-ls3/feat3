// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// The 2D Nonsteady Navier-Stokes CP-Q2/P1dc Toy-Code Solver (TM)
// ------------------------------------------------------------
// This application implements a "simple" parallel non-steady Navier-Stokes solver using
// the "CP" approach with Q2/P1dc space and Crank-Nicolson time discretization.
//
// ---------------
// !!! WARNING !!!
// ---------------
// This application is a "toy code" solver, i.e. it is meant as a playground for the
// HPC guys to tweak their parallel poisson solvers for more interesting scenarios than
// poisson on the unit-square. You can furthermore generate fancy videos of vortex
// streets to impress your girlfriend or to show your parents what you are being paid for.
// But: Do not expect this application to be accurate in time and/or space,
// so do *NOT* use it to perform any serious PDE/FEM analysis work !!!
// </Warning>
//
// This application is a "flow-through-a-domain" solver, i.e. it handles Navier-Stokes
// equations with an inflow and outflow region without any external forces, moving
// boundaries or any other fancy stuff.
//
// This application has four pre-configured benchmark problems, which can be launched
// by specifying the "--setup <config>" command line arguments, where <config> specifies
// one of the following:
//
// --setup square
// Loads the "Poiseuille-Flow-On-Unit-Square" problem.
// This is the most simple of the three pre-configured problems, where the time-dependent
// solution converges to a steady-state Poiseuille-Flow.
//
// --setup nozzle
// Loads the "Jet-Flow-Through-A-Nozzle" problem.
// Similar to the previous problem, but on a slightly more complex domain.
// The solution also converges to a steady-state flow.
//
// --setup fb_c2d_00
// --setup fb_c2d_01
// --setup fb_c2d_02
// --setup fb_c2d_03
// Loads the famous "Flow-Around-A-Cylinder" problem (non-steady version).
// This is the problem which generates the fancy "Von-Karman vortex shedding".
// In contrast to the previous problems, this solution is periodic.
// The three versions differ in the mesh that they use:
// * fb_c2d_00: 130 quads mesh (the old FEAT bench1 mesh)
// * fb_c2d_01: 32 quads mesh
// * fb_c2d_02: 48 quads mesh
// * fb_c2d_03: 64 quads mesh
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
// In the following, the options are categorized:
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
// This application supports only very limited customization of boundary conditions, which
// is limited to specifying
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
// for the whole time interval [0,T]. The "mesh width" of the time discretization is
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
// non-linear iteration. By default, only 1 step is performed.
//
//
// Linear Solver/Preconditioner Options
// ------------------------------------
// This application uses 2 multigrid solvers as linear preconditioners for the DPM:
//  1) a Richardson-Multigrid for the linearized Burgers system in velocity space (A-solver)
//  2) a PCG-Multigrid for the Poisson problem in pressure space (S-solver)
//
// Both multigrids use a damped Jacobi smoother as well as a Jacobi "smoother" as the coarse-grid solver.
//
// These two multigrid solvers are configured by the same set of options with slightly
// different names: options for the A-solver are postfixed by "-a", whereas options for
// the S-solver are postfixed by "-s".
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
// Default: 0.5 for both A- and S-solver
//
// Furthermore, it is possible to use a simple (one-grid) damped Jacobi-Iteration instead
// of multigrid as the A-solver. This can be achieved by supplying the '--no-multigrid-a' option.
//
// Furthermore, it is possible to use a simple (one-grid) damped Jacobi-Iteration instead
// of multigrid as the S-solver. This can be achieved by supplying the '--no-multigrid-s' option.
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
#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/bernstein2/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/global/pmdcdsc_matrix.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/umfpack.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/statistics.hpp>

#include <deque>
#include <numeric>

namespace NavierStokesCP2D
{
  using namespace FEAT;

  // helper functions for padded console output
  static inline void dump_line(const Dist::Comm& comm, String s, String t)
  {
    comm.print(s.pad_back(30, '.') + ": " + t);
  }

  static inline void dump_line(const Dist::Comm& comm, String s, const std::deque<String>& t)
  {
    comm.print(s.pad_back(30, '.') + ": " + stringify_join(t, " "));
  }

  template<typename T_>
  static inline void dump_line(const Dist::Comm& comm, String s, T_ t)
  {
    comm.print(s.pad_back(30, '.') + ": " + stringify(t));
  }

  static inline void dump_time(const Dist::Comm& comm, String s, double t, double total)
  {
    comm.print(s.pad_back(30, '.') + ": " + stringify_fp_fix(t, 3, 10)
      + " (" + stringify_fp_fix(100.0*t/total,3,7) + "%)");
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
    /// path to the mesh directory
    String mesh_path;
    /// names of the mesh files
    std::deque<String> mesh_files;

    /// minimum and maximum levels (as configured)
    Index level_min_in, level_max_in;

    /// minimum and maximum levels (after partitioning)
    Index level_min, level_max;

    /// base-name of VTK files
    String vtk_name;
    /// stepping of VTK output
    Index vtk_step;

    // names of inflow mesh-parts
    std::deque<String> part_names_in;

    // names of outflow mesh-part
    std::deque<String> part_names_out;

    // names of noflow mesh-parts
    std::deque<String> part_names_no;

    /// compute "flow-around-a-cylinder" benchmark quantities?
    bool flowbench_c2d;

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

    // use multigrid for S-solver ?
    bool multigrid_s;

    // use UMFPACK as coarse grid solver?
    bool coarse_umfpack_a;

    // use UMFPACK as coarse grid solver?
    bool coarse_umfpack_s;

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

    // enables verbose statistics output
    bool statistics;

    // specifies whether we run in test mode
    bool test_mode;

  public:
    Config() :
      level_min_in(0),
      level_max_in(0),
      level_min(0),
      level_max(0),
      vtk_step(0),
      flowbench_c2d(false),
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
      multigrid_s(false),
#ifdef FEAT_HAVE_UMFPACK
      coarse_umfpack_a(true),
      coarse_umfpack_s(true),
#else
      coarse_umfpack_a(false),
      coarse_umfpack_s(false),
#endif
      max_iter_a(50),
      tol_rel_a(1E-5),
      smooth_steps_a(4),
      smooth_damp_a(0.5),
      max_iter_s(100),
      tol_rel_s(1E-5),
      smooth_steps_s(4),
      smooth_damp_s(0.5),
      statistics(false),
      test_mode(false)
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
        else if(s.compare_no_case("fb_c2d_00") == 0)
          setup_fb_c2d_00();
        else if(s.compare_no_case("fb_c2d_01") == 0)
          setup_fb_c2d_01();
        else if(s.compare_no_case("fb_c2d_02") == 0)
          setup_fb_c2d_02();
        else if(s.compare_no_case("fb_c2d_03") == 0)
          setup_fb_c2d_03();
        else
        {
          std::cerr << "ERROR: unknown setup '" << s << "'" << "\n";
          return false;
        }
      }

      flowbench_c2d |= (args.check("flowbench") >= 0); // this may already be enabled by a setup above
      deformation = (args.check("deformation") >= 0);
      statistics = (args.check("statistics") >= 0);
      test_mode = (args.check("test-mode") >= 0);

      args.parse("mesh-path", mesh_path);
      if(args.check("mesh-file") > 0)
        mesh_files = args.query("mesh-file")->second;
      if(args.parse("vtk", vtk_name, vtk_step) == 1)
        vtk_step = 1; // vtk-name given, but not vtk-step, so set to 1
      args.parse("nu", nu);
      if(args.check("part-in") > 0)
        part_names_in = args.query("part-in")->second;
      if(args.check("part-out") > 0)
        part_names_out = args.query("part-out")->second;
      if(args.check("part-no") > 0)
        part_names_no = args.query("part-no")->second;
      args.parse("profile", ix0, iy0, ix1, iy1);
      args.parse("vmax", vmax);
      args.parse("time-max", time_max);
      args.parse("time-steps", time_steps);
      if(args.parse("max-time-steps", max_time_steps) < 1)
        max_time_steps = time_steps;
      args.parse("nl-steps", nonlin_steps);
      args.parse("dpm-steps", dpm_steps);
      multigrid_a = (args.check("no-multigrid-a") < 0);
      multigrid_s = (args.check("no-multigrid-s") < 0);
#ifdef FEAT_HAVE_UMFPACK
      coarse_umfpack_a = (args.check("no-umfpack-a") < 0);
      coarse_umfpack_s = (args.check("no-umfpack-s") < 0);
#endif
      args.parse("max-iter-a", max_iter_a);
      args.parse("tol-rel-a", tol_rel_a);
      args.parse("smooth-a", smooth_steps_a);
      args.parse("damp-a", smooth_damp_a);
      args.parse("max-iter-s", max_iter_s);
      args.parse("tol-rel-s", tol_rel_s);
      args.parse("smooth-s", smooth_steps_s);
      args.parse("damp-s", smooth_damp_s);

      // only 5 time-steps in test mode
      if(test_mode)
        max_time_steps = 5;

      return true;
    }

    void dump(const Dist::Comm& comm)
    {
      comm.print("\nConfiguration Summary:");
      dump_line(comm, "Mesh Path", mesh_path);
      dump_line(comm, "Mesh Files", mesh_files);
      dump_line(comm, "Level-Min", stringify(level_min) + " [" + stringify(level_min_in) + "]");
      dump_line(comm, "Level-Max", stringify(level_max) + " [" + stringify(level_max_in) + "]");
      dump_line(comm, "VTK-Name", vtk_name);
      dump_line(comm, "VTK-Step", vtk_step);
      dump_line(comm, "Flow Benchmark", (flowbench_c2d ? "yes" : "no"));
      dump_line(comm, "Inflow-Parts", part_names_in);
      dump_line(comm, "Outflow-Parts", part_names_out);
      dump_line(comm, "Noflow-Parts", part_names_no);
      dump_line(comm, "Inflow-Profile", "( " + stringify(ix0) + " , " + stringify(iy0) + " ) - ( "
        + stringify(ix1) + " , " + stringify(iy1) + " )");
      dump_line(comm, "V-Max", vmax);
      dump_line(comm, "Tensor", (deformation ? "Deformation" : "Gradient"));
      dump_line(comm, "Nu", nu);
      dump_line(comm, "Time-Max", time_max);
      dump_line(comm, "Time-Steps", time_steps);
      dump_line(comm, "Max Time-Steps", max_time_steps);
      dump_line(comm, "Non-Linear Steps", nonlin_steps);
      dump_line(comm, "Linear DPM Steps", dpm_steps);
      dump_line(comm, "A: Solver", (multigrid_a ? "Rich-Multigrid" : "BiCGStab-Jacobi"));
      dump_line(comm, "A: Max-Iter", max_iter_a);
      dump_line(comm, "A: Tol-Rel", tol_rel_a);
      dump_line(comm, "A: Smooth Steps", smooth_steps_a);
      dump_line(comm, "A: Smooth Damp", smooth_damp_a);
      dump_line(comm, "S: Solver", (multigrid_s ? "PCG-Multigrid" : "PCG-Jacobi"));
      dump_line(comm, "S: Max-Iter", max_iter_s);
      dump_line(comm, "S: Tol-Rel", tol_rel_s);
      dump_line(comm, "S: Smooth Steps", smooth_steps_s);
      dump_line(comm, "S: Smooth Damp", smooth_damp_s);
      dump_line(comm, "Test Mode", (test_mode ? "yes" : "no"));
      dump_line(comm, "Statistics", (statistics ? "yes" : "no"));
    }

    // Setup: Poiseuille-Flow on unit-square
    void setup_square()
    {
      mesh_files.push_back("unit-square-quad.xml");
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
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
      mesh_files.push_back("nozzle-2-quad.xml");
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
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

    // auxiliary function: basic flow benchmark setup
    void setup_fb_c2d_aux()
    {
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
      part_names_no.push_back("bnd:c");  // circle
      nu = 1E-3;
      ix0 = 0.0;
      iy0 = 0.0;
      ix1 = 0.0;
      iy1 = 0.41;
      vmax = 1.5;
      time_max = 3.0;
      flowbench_c2d = true;
    }

    // Setup: flow around a cylinder (130 quads)
    void setup_fb_c2d_00()
    {
      setup_fb_c2d_aux();
      mesh_files.push_back("flowbench_c2d_00_quad_130.xml");
      level_min = level_min_in = Index(0);
      level_max = level_max_in = Index(4);
      time_steps = max_time_steps = 1500;
    }

    // Setup: flow around a cylinder (32 quads)
    void setup_fb_c2d_01()
    {
      setup_fb_c2d_aux();
      mesh_files.push_back("flowbench_c2d_01_quad_32.xml");
      level_min = level_min_in = Index(0);
      level_max = level_max_in = Index(5);
      time_steps = max_time_steps = 1500;
    }

    // Setup: flow around a cylinder (48 quads)
    void setup_fb_c2d_02()
    {
      setup_fb_c2d_aux();
      mesh_files.push_back("flowbench_c2d_02_quad_48.xml");
      level_min = level_min_in = Index(0);
      level_max = level_max_in = Index(4);
      time_steps = max_time_steps = 1500;
    }

    // Setup: flow around a cylinder (64 quads)
    void setup_fb_c2d_03()
    {
      setup_fb_c2d_aux();
      mesh_files.push_back("flowbench_c2d_03_quad_64.xml");
      level_min = level_min_in = Index(0);
      level_max = level_max_in = Index(4);
      time_steps = max_time_steps = 1500;
    }
  }; // class Config

  template<typename DataType_>
  class BenchBodyForceAccumulator
  {
  private:
    const bool _defo;
    const DataType_ _nu;
    const DataType_ _v_max;

  public:
    DataType_ drag;
    DataType_ lift;

    explicit BenchBodyForceAccumulator(bool defo, DataType_ nu, DataType_ v_max) :
      _defo(defo), _nu(nu), _v_max(v_max),
      drag(DataType_(0)), lift(DataType_(0))
    {
    }

    /// 2D variant
    template<typename T_>
    void operator()(
      const T_ omega,
      const Tiny::Vector<T_, 2, 2>& /*pt*/,
      const Tiny::Matrix<T_, 2, 1, 2, 1>& jac,
      const Tiny::Vector<T_, 2, 2>& /*val_v*/,
      const Tiny::Matrix<T_, 2, 2, 2, 2>& grad_v,
      const T_ val_p)
    {
      // compute normal and tangential
      const T_ n2 = T_(1) / Math::sqrt(jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0));
      const T_ tx = jac(0,0) * n2;
      const T_ ty = jac(1,0) * n2;
      const T_ nx = -ty;
      const T_ ny =  tx;

      /// \todo adjust this to support the deformation tensor (?)

      Tiny::Matrix<T_, 2, 2, 2, 2> nt;
      nt(0,0) = tx * nx;
      nt(0,1) = tx * ny;
      nt(1,0) = ty * nx;
      nt(1,1) = ty * ny;

      const T_ dut = Tiny::dot(nt, grad_v);
      const T_ dpf1 = _nu;
      const T_ dpf2 = (2.0 / (0.1*Math::sqr(_v_max*(2.0/3.0)))); // = 2 / (rho * U^2 * D)

      drag += DataType_(omega * dpf2 * ( dpf1 * dut * ny - val_p * nx));
      lift += DataType_(omega * dpf2 * (-dpf1 * dut * nx - val_p * ny));
    }

    /// 3D variant
    template<typename T_>
    void operator()(
      const T_ omega,
      const Tiny::Vector<T_, 3, 3>& /*pt*/,
      const Tiny::Matrix<T_, 3, 2, 3, 2>& jac,
      const Tiny::Vector<T_, 3, 3>& /*val_v*/,
      const Tiny::Matrix<T_, 3, 3, 3, 3>& grad_v,
      const T_ val_p)
    {
      // compute normal and tangential
      const T_ n2 = T_(1) / Math::sqrt(jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0));
      const T_ tx = jac(0,0) * n2;
      const T_ ty = jac(1,0) * n2;
      const T_ nx = -ty;
      const T_ ny =  tx;

      /// \todo adjust this to support the deformation tensor (?)

      Tiny::Matrix<T_, 3, 3, 3, 3> nt;
      nt.format();
      nt(0,0) = tx * nx;
      nt(0,1) = tx * ny;
      nt(1,0) = ty * nx;
      nt(1,1) = ty * ny;

      const T_ dut = Tiny::dot(nt, grad_v);
      const T_ dpf1 = _nu;
      const T_ dpf2 = (2.0 / (0.1*Math::sqr(_v_max*(4.0/9.0))* 0.41)); // = 2 / (rho * U^2 * D * H)

      drag += DataType_(omega * dpf2 * ( dpf1 * dut * ny - val_p * nx));
      lift += DataType_(omega * dpf2 * (-dpf1 * dut * nx - val_p * ny));
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(&drag, &drag, std::size_t(1), Dist::op_sum);
      comm.allreduce(&lift, &lift, std::size_t(1), Dist::op_sum);
    }
  }; // class BenchBodyForceAccumulator<...>

  template<typename DataType_>
  class XFluxAccumulator
  {
  public:
    DataType_ flux;

    XFluxAccumulator() :
      flux(DataType_(0))
    {
    }

    template<typename T_, int d_, int d2_>
    void operator()(
      const T_ omega,
      const Tiny::Vector<T_, d_, d_>& /*pt*/,
      const Tiny::Matrix<T_, d_, d2_, d_, d2_>& /*jac*/,
      const Tiny::Vector<T_, d_, d_>& val_v,
      const Tiny::Matrix<T_, d_, d_, d_, d_>& /*grad_v*/,
      const T_ /*val_p*/)
    {
      flux += DataType_(omega * val_v[0]);
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(&flux, &flux, std::size_t(1), Dist::op_sum);
    }
  };

  /**
   * \brief Navier-Stokes System Level class
   *
   * This extends the StokesBlockedSystemLevel by the corresponding filters for
   * the velocity and pressure sub-systems.
   */
  template<
    int dim_,
    typename DataType_ = Real,
    typename IndexType_ = Index,
    typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
    typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
    typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
    typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>>
  class NavierStokesBlockedSystemLevel :
    public Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
  {
  public:
    typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;

    typedef typename BaseClass::GlobalVeloVector GlobalVeloVector;
    typedef typename BaseClass::LocalMatrixBlockA LocalMatrixBlockA;
    typedef typename BaseClass::GlobalMatrixBlockB GlobalMatrixBlockB;
    typedef typename BaseClass::GlobalMatrixBlockD GlobalMatrixBlockD;

    // lumped velocity mass matrix
    GlobalVeloVector lumped_mass_velo;
    // inverse lumped velocity mass matrix
    GlobalVeloVector inverse_lumped_mass_velo;

    typedef Global::PMDCDSCMatrix<GlobalMatrixBlockB, GlobalMatrixBlockD> GlobalSchurMatrix;

    // filtered local matrix A for UMFPACK coarse grid solver
    LocalMatrixBlockA local_matrix_a_filtered;

    // schur matrix
    GlobalSchurMatrix matrix_s;

    NavierStokesBlockedSystemLevel() :
      lumped_mass_velo(&this->gate_velo),
      matrix_s(this->inverse_lumped_mass_velo, this->matrix_b, this->matrix_d)
    {
    }
  }; // class NavierStokesBlockedSystemLevel

  template <typename SystemLevelType,typename DomainLevel_>
  void report_statistics(const Dist::Comm& comm, double t_total,
    std::deque<std::shared_ptr<SystemLevelType>> & system_levels,
      Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    const int nprocs = comm.size();

    /// \todo cover exactly all la op timings (some are not timed yet in the application) and replace t_total by them
    double solver_toe = t_total; //t_solver_a + t_solver_s + t_calc_def;

    comm.print("\nsolver_a:");
    comm.print(FEAT::Statistics::get_formatted_solver_tree("solver_a").trim());
    comm.print("solver_s:");
    comm.print(FEAT::Statistics::get_formatted_solver_tree("solver_s").trim());

    std::size_t la_size(0);
    std::size_t mpi_size(0);
    for(const auto& sl : system_levels)
    {
      la_size += sl->bytes();
      mpi_size += sl->gate_sys.bytes() + sl->coarse_muxer_sys.bytes();
    }
    String op_timings = FEAT::Statistics::get_formatted_times(solver_toe);

    Index cells_coarse_local = domain.back()->get_mesh().get_num_elements();
    Index cells_coarse_max;
    Index cells_coarse_min;
    comm.allreduce(&cells_coarse_local, &cells_coarse_max, std::size_t(1), Dist::op_max);
    comm.allreduce(&cells_coarse_local, &cells_coarse_min, std::size_t(1), Dist::op_min);
    Index cells_fine_local = domain.front()->get_mesh().get_num_elements();
    Index cells_fine_max;
    Index cells_fine_min;
    comm.allreduce(&cells_fine_local, &cells_fine_max, std::size_t(1), Dist::op_max);
    comm.allreduce(&cells_fine_local, &cells_fine_min, std::size_t(1), Dist::op_min);

    Index dofs_coarse_local = system_levels.back()->matrix_a.local().columns() + system_levels.back()->matrix_b.local().columns();
    Index dofs_coarse_max;
    Index dofs_coarse_min;
    comm.allreduce(&dofs_coarse_local, &dofs_coarse_max, std::size_t(1), Dist::op_max);
    comm.allreduce(&dofs_coarse_local, &dofs_coarse_min, std::size_t(1), Dist::op_min);
    Index dofs_fine_local = system_levels.front()->matrix_a.local().columns() + system_levels.front()->matrix_b.local().columns();
    Index dofs_fine_max;
    Index dofs_fine_min;
    comm.allreduce(&dofs_fine_local, &dofs_fine_max, std::size_t(1), Dist::op_max);
    comm.allreduce(&dofs_fine_local, &dofs_fine_min, std::size_t(1), Dist::op_min);

    Index nzes_coarse_local = system_levels.back()->matrix_a.local().used_elements() + system_levels.back()->lumped_mass_velo.local().size();
    Index nzes_coarse_max;
    Index nzes_coarse_min;
    comm.allreduce(&nzes_coarse_local, &nzes_coarse_max, std::size_t(1), Dist::op_max);
    comm.allreduce(&nzes_coarse_local, &nzes_coarse_min, std::size_t(1), Dist::op_min);
    Index nzes_fine_local = system_levels.front()->matrix_a.local().used_elements() + system_levels.back()->lumped_mass_velo.local().size();
    Index nzes_fine_max;
    Index nzes_fine_min;
    comm.allreduce(&nzes_fine_local, &nzes_fine_max, std::size_t(1), Dist::op_max);
    comm.allreduce(&nzes_fine_local, &nzes_fine_min, std::size_t(1), Dist::op_min);

    String flops = FEAT::Statistics::get_formatted_flops(solver_toe, nprocs);
    comm.print(flops + "\n");
    comm.print(op_timings);
    comm.print(FEAT::Statistics::get_formatted_solver_internals("solver_a"));
    comm.print("");
    comm.print(FEAT::Statistics::get_formatted_solver_internals("solver_s"));
    comm.print("");
    comm.print(String("Domain size:").pad_back(20) + stringify(double(domain.bytes())  / (1024. * 1024.))  + " MByte");
    comm.print(String("MPI size:").pad_back(20) + stringify(double(mpi_size) / (1024. * 1024.)) + " MByte");
    comm.print(String("LA size:").pad_back(20) + stringify(double(la_size) / (1024. * 1024.)) + " MByte\n");
    MemoryUsage mu;
    comm.print(mu.get_formatted_memory_usage());
    comm.print(String("#Mesh cells:").pad_back(20) + "coarse " + stringify(cells_coarse_max) + "/" + stringify(cells_coarse_min) + ", fine " +
        stringify(cells_fine_max) + "/" + stringify(cells_fine_min));
    comm.print(String("#DOFs:").pad_back(20) + "coarse " + stringify(dofs_coarse_max) + "/" + stringify(dofs_coarse_min) + ", fine " +
        stringify(dofs_fine_max) + "/" + stringify(dofs_fine_min));
    comm.print(String("#NZEs").pad_back(20) + "coarse " + stringify(nzes_coarse_max) + "/" + stringify(nzes_coarse_min) + ", fine " +
        stringify(nzes_fine_max) + "/" + stringify(nzes_fine_min) + "\n");
  }

  template<typename DomainLevel_>
  void run(Config& cfg, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();
    const int rank = comm.rank();

    // create a time-stamp
    TimeStamp stamp_start;

    // our arch types
    typedef Real DataType;
    typedef Index IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef typename DomainLevelType::TrafoType TrafoType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // define our velocity and pressure system levels
    typedef NavierStokesBlockedSystemLevel<dim, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = Index(domain.size_physical());

    // create a batch of stop-watches
    StopWatch watch_total, watch_asm_rhs, watch_asm_mat, watch_calc_def,
      watch_sol_init, watch_solver_a, watch_solver_s, watch_vtk;

    // create stokes and system levels
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    String cubature_name("auto-degree:7");
    Cubature::DynamicFactory cubature(cubature_name);

    // compute time-step size
    const DataType delta_t = cfg.time_max / DataType(cfg.time_steps);

    /* ***************************************************************************************** */

    comm.print("");
    comm.print("Assembling gates, muxers and transfers...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gates(domain.at(i));
    }

    for (Index i(0); (i < domain.size_physical()) && ((i+1) < domain.size_virtual()); ++i)
    {
      system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
      if((i+1) < domain.size_physical())
        system_levels.at(i)->assemble_transfers(*system_levels.at(i+1), domain.at(i), domain.at(i+1), cubature_name, true);
      else
        system_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature_name, true);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling basic matrices...");

    for(Index i(0); i < num_levels; ++i)
    {
      // assemble velocity matrix structure
      system_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
      // assemble B/D matrices
      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);

      // assemble lumped velocity mass matrix
      {
        /*Assembly::BurgersAssembler<DataType, IndexType, 2> burgers_mat;
        burgers_mat.theta = DataType(1);

        auto& loc_mat_a = system_levels.at(i)->matrix_a.local();
        auto loc_vec_v = loc_mat_a.create_vector_l();
        loc_vec_v.format();
        loc_mat_a.format();
        burgers_mat.assemble_matrix(loc_mat_a, loc_vec_v, domain.at(i)->space_velo, cubature);*/
        auto& loc_mat_a = system_levels.at(i)->matrix_a.local();
        loc_mat_a.format();
        Assembly::Common::IdentityOperatorBlocked<dim> id_op;
        Assembly::BilinearOperatorAssembler::assemble_matrix1(loc_mat_a, id_op, domain.at(i)->space_velo, cubature);
      }

      system_levels.at(i)->lumped_mass_velo = system_levels.at(i)->matrix_a.lump_rows();
      system_levels.at(i)->inverse_lumped_mass_velo.clone(system_levels.at(i)->lumped_mass_velo);
      system_levels.at(i)->inverse_lumped_mass_velo.component_invert(system_levels.at(i)->lumped_mass_velo);

      // coarse grid level?
      if((i+1u) == domain.size_virtual())
      {
        system_levels.at(i)->local_matrix_a_filtered = system_levels.at(i)->matrix_a.convert_to_1();
      }

      // perform symbolic initialization of Schur-complement matrix
      if((i == Index(0)) || cfg.multigrid_s)
      {
        system_levels.at(i)->matrix_s.init_symbolic();
      }
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    // create parabolic inflow profile
    Analytic::Common::ParProfileVector<DataType> inflow(cfg.ix0, cfg.iy0, cfg.ix1, cfg.iy1, cfg.vmax);

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local system filters
      typename SystemLevelType::LocalVeloFilter& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // create unit-filter assemblers
      Assembly::UnitFilterAssembler<MeshType> unit_asm_noflow, unit_asm_inflow;

      // loop over all inflow boundary parts
      for(const auto& name : cfg.part_names_in)
      {
        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);

        // found it?
        XASSERT(mesh_part_node != nullptr);

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        auto* mesh_part = mesh_part_node->get_mesh();
        if(mesh_part == nullptr)
          continue;

        // add to corresponding boundary assembler
        unit_asm_inflow.add_mesh_part(*mesh_part);
      }

      // loop over all noflow boundary parts
      for(const auto& name : cfg.part_names_no)
      {
        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
        XASSERT(mesh_part_node != nullptr);
        auto* mesh_part = mesh_part_node->get_mesh();
        if(mesh_part == nullptr)
          continue;
        unit_asm_noflow.add_mesh_part(*mesh_part);
      }

      // assemble the velocity filter
      unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo);

      // assemble inflow BC
      unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, inflow);

      // apply velocity filter onto inverse lumped mass matrix
      fil_loc_v.filter_cor(system_levels.at(i)->inverse_lumped_mass_velo.local());

      // filter matrix block A if available
      if(system_levels.at(i)->local_matrix_a_filtered.rows() > Index(0))
        fil_loc_v.filter_mat(system_levels.at(i)->local_matrix_a_filtered);

      // perform numeric initialization of Schur-complement matrix
      if((i == Index(0)) || cfg.multigrid_s)
      {
        system_levels.at(i)->matrix_s.init_numeric();
      }
    }

    /* ***************************************************************************************** */

    // get our vector types
    typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // get our fine-level matrices
    typename SystemLevelType::GlobalMatrixBlockA& matrix_a = the_system_level.matrix_a;
    typename SystemLevelType::GlobalMatrixBlockB& matrix_b = the_system_level.matrix_b;
    typename SystemLevelType::GlobalMatrixBlockD& matrix_d = the_system_level.matrix_d;
    typename SystemLevelType::GlobalSchurMatrix& matrix_s = the_system_level.matrix_s;

    // get out fine-level filters
    typename SystemLevelType::GlobalVeloFilter& filter_v = the_system_level.filter_velo;
    typename SystemLevelType::GlobalPresFilter& filter_p = the_system_level.filter_pres;

    /* ***************************************************************************************** */

    // set up body forces and flux assemblers
    Assembly::TraceAssembler<TrafoType> body_force_asm(the_domain_level.trafo);
    Assembly::TraceAssembler<TrafoType> flux_u_asm(the_domain_level.trafo);
    Assembly::TraceAssembler<TrafoType> flux_l_asm(the_domain_level.trafo);
    Trafo::InverseMappingData<DataType, dim> point_pa, point_pe;

    // set up post-processing for flow benchmark?
    if(cfg.flowbench_c2d)
    {
      auto* mesh_part_c = the_domain_level.get_mesh_node()->find_mesh_part("bnd:c");
      auto* mesh_part_u = the_domain_level.get_mesh_node()->find_mesh_part("inner:u");
      auto* mesh_part_l = the_domain_level.get_mesh_node()->find_mesh_part("inner:l");
      if(mesh_part_c != nullptr)
        body_force_asm.add_mesh_part(*mesh_part_c);
      if(mesh_part_u != nullptr)
        flux_u_asm.add_mesh_part(*mesh_part_u);
      if(mesh_part_l != nullptr)
        flux_l_asm.add_mesh_part(*mesh_part_l);

      body_force_asm.compile();
      flux_u_asm.compile();
      flux_l_asm.compile();

      // set up inverse trafo mapping
      typedef Trafo::InverseMapping<TrafoType, DataType> InvMappingType;
      InvMappingType inv_mapping(the_domain_level.trafo);

      // reference pressure points
      typename InvMappingType::ImagePointType v_a, v_e;
      if(dim == 2)
      {
        v_a[0] = 0.15;
        v_e[0] = 0.25;
        v_a[1] = v_e[1] = 0.2;
      }
      else
      {
        v_a[0] = 0.45;
        v_e[0] = 0.55;
        v_a[1] = v_e[1] = 0.2;
        v_a[2] = v_e[2] = 0.205;
      }

      // unmap points
      point_pa = inv_mapping.unmap_point(v_a, true);
      point_pe = inv_mapping.unmap_point(v_e, true);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("Setting up Velocity Burgers Solver...");

    // create a multigrid solver for the velocity
    auto multigrid_hierarchy_velo = std::make_shared<Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalMatrixBlockA,
      typename SystemLevelType::GlobalVeloFilter,
      typename SystemLevelType::GlobalVeloTransfer
      >>(domain.size_virtual());

    // use multigrid for A-solver?
    if(cfg.multigrid_a)
    {
      // loop over all levels
      for (std::size_t i(0); i < system_levels.size(); ++i)
      {
        auto& lvl = *system_levels.at(i);

        if((i+1) < domain.size_virtual())
        {
          // smoother
          auto jac = Solver::new_jacobi_precond(lvl.matrix_a, lvl.filter_velo);
          auto smoother = Solver::new_richardson(lvl.matrix_a, lvl.filter_velo, cfg.smooth_damp_a, jac);
          smoother->set_max_iter(cfg.smooth_steps_a);
          smoother->set_min_iter(cfg.smooth_steps_a);
          multigrid_hierarchy_velo->push_level(lvl.matrix_a, lvl.filter_velo, lvl.transfer_velo,
            smoother, smoother, smoother);
        }
#ifdef FEAT_HAVE_UMFPACK
        else if(cfg.coarse_umfpack_a)
        {
          auto umf = Solver::new_generic_umfpack(lvl.local_matrix_a_filtered);
          auto cgs = Solver::new_schwarz_precond(umf, lvl.filter_velo);
          multigrid_hierarchy_velo->push_level(lvl.matrix_a, lvl.filter_velo, cgs);
        }
#endif // FEAT_HAVE_UMFPACK
        else
        {
          // coarse grid solver
          auto jac = Solver::new_jacobi_precond(lvl.matrix_a, lvl.filter_velo);
          auto cgs = Solver::new_richardson(lvl.matrix_a, lvl.filter_velo, cfg.smooth_damp_a, jac);
          cgs->set_max_iter(cfg.smooth_steps_a);
          cgs->set_min_iter(cfg.smooth_steps_a);
          multigrid_hierarchy_velo->push_level(lvl.matrix_a, lvl.filter_velo, cgs);
        }
      }
    }

    /* ***************************************************************************************** */

    comm.print("Setting up Pressure Laplace Solver...");

    // create another multigrid solver for the pressure
    auto multigrid_hierarchy_pres = std::make_shared<Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalSchurMatrix,
      typename SystemLevelType::GlobalPresFilter,
      typename SystemLevelType::GlobalPresTransfer
      >>(domain.size_virtual());

    if (cfg.multigrid_s)
    {
      // loop over all levels
      for (std::size_t i(0); i < system_levels.size(); ++i)
      {
        auto& lvl = *system_levels.at(i);

        if((i+1) < domain.size_virtual())
        {
          // smoother
          auto jac = Solver::new_jacobi_precond(lvl.matrix_s, lvl.filter_pres);
          auto smoother = Solver::new_richardson(lvl.matrix_s, lvl.filter_pres, cfg.smooth_damp_s, jac);
          smoother->set_max_iter(cfg.smooth_steps_s);
          smoother->set_min_iter(cfg.smooth_steps_s);
          multigrid_hierarchy_pres->push_level(lvl.matrix_s, lvl.filter_pres, lvl.transfer_pres,
            smoother, smoother, smoother);
        }
#ifdef FEAT_HAVE_UMFPACK
        else if(cfg.coarse_umfpack_s)
        {
          auto umf = Solver::new_umfpack(lvl.matrix_s.get_local_schur_matrix());
          auto cgs = Solver::new_schwarz_precond(umf, lvl.filter_pres);
          multigrid_hierarchy_pres->push_level(lvl.matrix_s, lvl.filter_pres, cgs);
        }
#endif // FEAT_HAVE_UMFPACK
        else
        {
          // coarse grid solver
          auto jac = Solver::new_jacobi_precond(lvl.matrix_s, lvl.filter_pres);
          auto cgs = Solver::new_richardson(lvl.matrix_s, lvl.filter_pres, cfg.smooth_damp_s, jac);
          cgs->set_max_iter(cfg.smooth_steps_s);
          cgs->set_min_iter(cfg.smooth_steps_s);
          multigrid_hierarchy_pres->push_level(lvl.matrix_s, lvl.filter_pres, cgs);
        }
      }
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
      // use BiCGStab-Jacobi
      auto jac = Solver::new_jacobi_precond(the_system_level.matrix_a, the_system_level.filter_velo);
      solver_a = Solver::new_bicgstab(matrix_a, filter_v, jac);
    }

    solver_a->set_max_iter(cfg.max_iter_a);
    solver_a->set_tol_rel(cfg.tol_rel_a);

    // for the velocity multigrid/solver, we can only perform symbolic initialization up to now:
    if(cfg.multigrid_a)
    {
      multigrid_hierarchy_velo->init_symbolic();
    }
    solver_a->init_symbolic();

    /* ***************************************************************************************** */

    // create a PCG-MG for the pressure
    std::shared_ptr<Solver::IterativeSolver<GlobalPresVector>> solver_s;

    if(cfg.multigrid_s)
    {
      // use Multigrid
      auto multigrid_pres = Solver::new_multigrid(multigrid_hierarchy_pres);
      solver_s = Solver::new_pcg(matrix_s, filter_p, multigrid_pres);
    }
    else
    {
      // Use PCG-Jacobi
      auto jac = Solver::new_jacobi_precond(matrix_s, filter_p);
      solver_s = Solver::new_pcg(matrix_s, filter_p, jac);
    }

    solver_s->set_max_iter(cfg.max_iter_s);
    solver_s->set_tol_rel(cfg.tol_rel_s);

    // for the pressure multigrid, we can perform full initialization:
    if(cfg.multigrid_s)
    {
      multigrid_hierarchy_pres->init();
    }
    solver_s->init();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("\n");

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
      if(cfg.flowbench_c2d)
      {
        head += String("Drag       ");
        head += String("Lift       ");
        head += String("P-Diff     ");
        head += String("U-Flux     ");
        head += String("L-Flux       ");
      }
      head += String("Runtime    ");
      comm.print(head);
      comm.print(String(head.size(), '-'));
    }

    watch_total.start();

    // keep track whether something failed miserably...
    bool failure = false;

    // set up a burgers assembler for the RHS vector
    Assembly::BurgersAssembler<DataType, IndexType, 2> burgers_rhs;
    burgers_rhs.deformation = cfg.deformation;
    burgers_rhs.nu = -cfg.nu;
    burgers_rhs.beta = -DataType(1);
    burgers_rhs.theta = DataType(2) / delta_t;

    // set up a burgers assembler for the velocity matrix
    Assembly::BurgersAssembler<DataType, IndexType, 2> burgers_mat;
    burgers_mat.deformation = cfg.deformation;
    burgers_mat.nu = cfg.nu;
    burgers_mat.beta = DataType(1);
    burgers_mat.theta = DataType(2) / delta_t;

    Statistics::reset();

    // time-step loop
    for(Index time_step(1); time_step <= cfg.max_time_steps; ++time_step)
    {
      // compute current time
      const DataType cur_time = DataType(time_step) * delta_t;

      // assemble RHS vector
      watch_asm_rhs.start();
      vec_rhs_v.format();
      vec_rhs_p.format();
      burgers_rhs.assemble_vector(
        vec_rhs_v.local(), vec_sol_v.local(), vec_sol_v.local(), the_domain_level.space_velo, cubature);
      matrix_b.local().apply(vec_rhs_v.local(), vec_sol_p.local(), vec_rhs_v.local(), -DataType(1));
      vec_rhs_v.sync_0();
      watch_asm_rhs.stop();

      // apply RHS filter
      filter_v.filter_rhs(vec_rhs_v);

      // non-linear loop
      for(Index nonlin_step(0); nonlin_step < cfg.nonlin_steps; ++nonlin_step)
      {
        if(comm.rank() == 0 && nonlin_step > 0u)
          std::cout << "\n";

        // Phase 1: compute convection vector
        // extrapolate previous time-step solution in first NL step
        if((time_step > Index(2)) && (nonlin_step == Index(0)))
        {
          // linear extrapolation of solution in time
          vec_conv.scale(vec_sol_v_1, DataType(2));
          vec_conv.axpy(vec_sol_v_2,-DataType(1));
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
          // clone convection vector
          auto vec_cv = vec_conv.clone();

          // assemble burgers matrices on all levels
          for(std::size_t i(0); i < system_levels.size(); ++i)
          {
            // assemble burgers matrix on this level
            auto& loc_mat_a = system_levels.at(i)->matrix_a.local();
            loc_mat_a.format();
            burgers_mat.assemble_matrix(loc_mat_a, vec_cv.local(), domain.at(i)->space_velo, cubature);

            // no more virtual levels to restrict to?
            if((i+1) >= domain.size_virtual())
              break;

            // does this process have another system level?
            if((i+1) < system_levels.size())
            {
              // create a coarse mesh velocity vector
              auto vec_crs = system_levels.at(i+1)->matrix_a.create_vector_l();

              // truncate fine mesh velocity vector
              system_levels.at(i)->transfer_velo.trunc(vec_cv, vec_crs);

              // the coarse vector is our next convection vector
              vec_cv = std::move(vec_crs);
            }
            else
            {
              // this process is a child, so send truncation to parent
              system_levels.at(i)->transfer_velo.trunc_send(vec_cv);
            }
          }
        }
        else
        {
          // assemble burgers matrices on finest level only
          the_system_level.matrix_a.local().format();
          burgers_mat.assemble_matrix(
            the_system_level.matrix_a.local(), vec_conv.local(), the_domain_level.space_velo, cubature);
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

        // Phase 4: initialize linear solvers
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
          FEAT::Statistics::expression_target = "solver_a";
          watch_solver_a.start();
          Solver::Status status_a = solver_a->apply(vec_cor_v, vec_def_v);
          watch_solver_a.stop();
          if(!Solver::status_success(status_a))
          {
            comm.print(std::cerr, "\n\nERROR: velocity solver broke down!\n");
            failure = true;
            break;
          }
          iter_v += solver_a->get_num_iter();

          // update velocity solution
          vec_sol_v.axpy(vec_cor_v);

          // are we done yet?
          if(dpm_step >= cfg.dpm_steps)
            break;

          // update pressure defect
          watch_calc_def.start();
          matrix_d.apply(vec_def_p, vec_sol_v, vec_rhs_p, -DataType(1));
          filter_p.filter_def(vec_def_p);
          watch_calc_def.stop();

          // solve pressure system
          FEAT::Statistics::expression_target = "solver_s";
          watch_solver_s.start();
          Solver::Status status_s = solver_s->apply(vec_cor_p, vec_def_p);
          watch_solver_s.stop();
          if(!Solver::status_success(status_s))
          {
            comm.print(std::cerr, "\n\nERROR: pressure solver broke down!\n");
            failure = true;
            break;
          }
          iter_p += solver_s->get_num_iter();

          // update pressure solution
          vec_sol_p.axpy(vec_cor_p, -DataType(1) / delta_t);

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
          line += stringify(iter_p).pad_front(4);
          std::cout << line;// << "\n";
        }
      } // non-linear loop

      // epic fail?
      if(failure)
        break;

      if(cfg.flowbench_c2d)
      {
        // body forces, pressure difference and flux values
        DataType c_drag = DataType(0);
        DataType c_lift = DataType(0);
        DataType p_diff = DataType(0);
        DataType u_flux = DataType(0);
        DataType l_flux = DataType(0);

        // assemble drag and lift forces
        {
          BenchBodyForceAccumulator<DataType> body_force_accum(false, cfg.nu, cfg.vmax);
          body_force_asm.assemble_flow_accum(
            body_force_accum,
            vec_sol_v.local(),
            vec_sol_p.local(),
            the_domain_level.space_velo,
            the_domain_level.space_pres,
            cubature);
          body_force_accum.sync(comm);
          c_drag = body_force_accum.drag;
          c_lift = body_force_accum.lift;
        }

        // compute pressure values and difference
        {
          // evaluate pressure
          auto pval_a = Assembly::DiscreteEvaluator::eval_fe_function(
            point_pa, vec_sol_p.local(), the_domain_level.space_pres);
          auto pval_e = Assembly::DiscreteEvaluator::eval_fe_function(
            point_pe, vec_sol_p.local(), the_domain_level.space_pres);

          // compute pressure mean
          const auto p_a = pval_a.mean_value_dist(comm);
          const auto p_e = pval_e.mean_value_dist(comm);
          p_diff = p_a - p_e;
        }

        // assemble upper flux
        {
          XFluxAccumulator<DataType> flux_accum;
          flux_u_asm.assemble_flow_accum(
            flux_accum,
            vec_sol_v.local(),
            vec_sol_p.local(),
            the_domain_level.space_velo,
            the_domain_level.space_pres,
            cubature);

          flux_accum.sync(comm);
          u_flux = flux_accum.flux / DataType(2);
        }

        // assemble lower flux
        {
          XFluxAccumulator<DataType> flux_accum;
          flux_l_asm.assemble_flow_accum(
            flux_accum,
            vec_sol_v.local(),
            vec_sol_p.local(),
            the_domain_level.space_velo,
            the_domain_level.space_pres,
            cubature);

          flux_accum.sync(comm);
          l_flux = flux_accum.flux / DataType(2);
        }

        // console output, part 3
        if(rank == 0)
        {
          String line;
          line += " | ";
          line += stringify_fp_fix(c_drag, 7).trunc_back(10u).pad_front(10u) + " ";
          line += stringify_fp_fix(c_lift, 7).trunc_back(10u).pad_front(10u) + " ";
          line += stringify_fp_fix(p_diff, 7).trunc_back(10u).pad_front(10u) + " ";
          line += stringify_fp_fix(u_flux, 7).trunc_back(10u).pad_front(10u) + " ";
          line += stringify_fp_fix(l_flux, 7).trunc_back(10u).pad_front(10u);
          std::cout << line;
        }
      }

      if(rank == 0)
      {
        String line(" | ");
        line += stamp_start.elapsed_string_now().pad_front(10);
        std::cout << line << "\n";
      }

      // VTK-Export
      if(!cfg.vtk_name.empty() && (cfg.vtk_step > 0) && (time_step % cfg.vtk_step == 0))
      {
        watch_vtk.start();
        String vtk_path = cfg.vtk_name + "." + stringify(comm.size()) + "." + stringify(time_step).pad_front(5, '0');

        Geometry::ExportVTK<MeshType> vtk(the_domain_level.get_mesh());

        // write solution
        vtk.add_vertex_vector("v", vec_sol_v.local());

        // project pressure
        Cubature::DynamicFactory cub("gauss-legendre:2");
        LAFEM::DenseVector<double, Index> vtx_p, vtx_der_p;
        Assembly::DiscreteCellProjector::project(vtx_p, vec_sol_p.local(), the_domain_level.space_pres, cub);

        // write pressure
        vtk.add_cell_scalar("p", vtx_p.elements());

        // compute and write time-derivatives
        GlobalVeloVector vec_der_v = vec_sol_v.clone();
        GlobalPresVector vec_der_p = vec_sol_p.clone();
        vec_der_v.axpy(vec_sol_v_1, -DataType(1));
        vec_der_p.axpy(vec_sol_p_1, -DataType(1));
        vec_der_v.scale(vec_der_v, DataType(1) / delta_t);
        vec_der_p.scale(vec_der_p, DataType(1) / delta_t);
        Assembly::DiscreteCellProjector::project(vtx_der_p, vec_der_p.local(), the_domain_level.space_pres, cub);

        vtk.add_vertex_vector("v_dt", vec_der_v.local());
        vtk.add_cell_scalar("p_dt", vtx_der_p.elements());

        // export
        vtk.write(vtk_path, comm);
        watch_vtk.stop();
      }

      // finally, update our solution vector backups
      vec_sol_v_2.copy(vec_sol_v_1);
      vec_sol_v_1.copy(vec_sol_v);
      vec_sol_p_1.copy(vec_sol_p);

      // compress all statistics from the current timestep for further analysis after the solution is finished
      FEAT::Statistics::compress_solver_expressions();

      // continue with next time-step
    } // time-step loop

    watch_total.stop();

    // are we in test-mode?
    if(cfg.test_mode)
    {
      if(failure)
        comm.print("\nTest-Mode: FAILED");
      else
        comm.print("\nTest-Mode: PASSED");
    }

    // release pressure solvers
    solver_s->done();
    if(cfg.multigrid_s)
    {
      multigrid_hierarchy_pres->done();
    }

    // release velocity solvers
    solver_a->done_symbolic();
    if(cfg.multigrid_a)
      multigrid_hierarchy_velo->done_symbolic();

    double t_total = watch_total.elapsed();
    double t_asm_mat = watch_asm_mat.elapsed();
    double t_asm_rhs = watch_asm_rhs.elapsed();
    double t_calc_def = watch_calc_def.elapsed();
    double t_sol_init = watch_sol_init.elapsed();
    double t_solver_a = watch_solver_a.elapsed();
    double t_solver_s = watch_solver_s.elapsed();
    double t_vtk = watch_vtk.elapsed();
    double t_sum = t_asm_mat + t_asm_rhs + t_calc_def + t_sol_init + t_solver_a + t_solver_s + t_vtk;

    // write timings
    if(rank == 0)
    {
      comm.print("");
      dump_time(comm, "Total Solver Time", t_total, t_total);
      dump_time(comm, "Matrix Assembly Time", t_asm_mat, t_total);
      dump_time(comm, "Vector Assembly Time", t_asm_rhs, t_total);
      dump_time(comm, "Defect-Calc Time", t_calc_def, t_total);
      dump_time(comm, "Solver-A Init Time", t_sol_init, t_total);
      dump_time(comm, "Solver-A Time", t_solver_a, t_total);
      dump_time(comm, "Solver-S Time", t_solver_s, t_total);
      dump_time(comm, "VTK-Write Time", t_vtk, t_total);
      dump_time(comm, "Other Time", t_total-t_sum, t_total);

      if(cfg.multigrid_a)
      {
        comm.print("\nMultigrid-A Timings:");
        dump_time(comm, "Smoother Time", multigrid_hierarchy_velo->get_time_smooth(), t_total);
        dump_time(comm, "Coarse Grid Time", multigrid_hierarchy_velo->get_time_coarse(), t_total);
        dump_time(comm, "Defect Calc Time", multigrid_hierarchy_velo->get_time_defect(), t_total);
        dump_time(comm, "Grid Transfer Time", multigrid_hierarchy_velo->get_time_transfer(), t_total);
      }

      if(cfg.multigrid_s)
      {
        comm.print("\nMultigrid-S Timings:");
        dump_time(comm, "Smoother Time", multigrid_hierarchy_pres->get_time_smooth(), t_total);
        dump_time(comm, "Coarse Grid Time", multigrid_hierarchy_pres->get_time_coarse(), t_total);
        dump_time(comm, "Defect Calc Time", multigrid_hierarchy_pres->get_time_defect(), t_total);
        dump_time(comm, "Grid Transfer Time", multigrid_hierarchy_pres->get_time_transfer(), t_total);
      }
    }

    if (cfg.statistics)
    {
      report_statistics(comm, t_total, system_levels, domain);
    }
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
    args.support("help", "\nDisplays this help message.\n");
    args.support("setup", "<config>\nLoads a pre-defined configuration:\n"
      "square      Poiseuille-Flow on Unit-Square\n"
      "nozzle      Jet-Flow through Nozzle domain\n"
      "fb_c2d_00   Nonsteady Flow Around A Cylinder (bench1 mesh)\n"
      "fb_c2d_01   Nonsteady Flow Around A Cylinder (32 quad mesh)\n"
      "fb_c2d_02   Nonsteady Flow Around A Cylinder (48 quad mesh)\n"
      "fb_c2d_03   Nonsteady Flow Around A Cylinder (64 quad mesh)\n"
    );
    args.support("flowbench","\nEnables the computation of 'flow around a cylinder' post-processing\nquantities such as drag, lift, etc.\n");
    args.support("deformation", "\nUse deformation tensor instead of gradient tensor.\n");
    args.support("nu", "<nu>\nSets the viscosity parameter.\n");
    args.support("time-max", "<T_max>\nSets the maximum simulation time T_max.\n");
    args.support("time-steps", "<N>\nSets the number of time-steps for the time interval.\n");
    args.support("max-time-steps", "<N>\nSets the maximum number of time-steps to perform.\n");
    args.support("part-in", "<names...>\nSpecifies the names of the inflow mesh-parts.\n");
    args.support("part-out", "<names...>\nSpecifies the names of the outflow mesh-parts.\n");
    args.support("part-no", "<names...>\nSpecifies the names of the noflow mesh-parts.\n");
    args.support("profile", "<x0> <y0> <x1> <y1>\nSpecifies the line segment coordinates for the inflow profile.\n");
    args.support("vmax", "<v-max>\nSpecifies the maximum inflow velocity.\n");
    args.support("level", "<max> [<min>]\nSets the maximum and minimum mesh refinement levels.\n");
    args.support("vtk", "<name> [<step>]\nSets the name for VTK output and the time-stepping for the output (optional).\n");
    args.support("mesh-file", "<name>\nSpecifies the filename of the input mesh file.\n");
    args.support("mesh-path", "<path>\nSpecifies the path of the directory containing the mesh file.\n");
    args.support("nl-steps", "<N>\nSets the number of non-linear iterations per time-step.\nDefault: 1\n");
    args.support("dpm-steps", "<N>\nSets the number of Discrete-Projection-Method steps per non-linear step.\nDefault: 1\n");
    args.support("no-multigrid-a", "\nUse BiCGStab-Jacobi instead of Multigrid as A-Solver.\n");
    args.support("max-iter-a", "<N>\nSets the maximum number of allowed iterations for the A-Solver.\nDefault: 50\n");
    args.support("tol-rel-a", "<eps>\nSets the relative tolerance for the A-Solver.\nDefault: 1E-5\n");
    args.support("smooth-a", "<N>\nSets the number of smoothing steps for the A-Solver.\nDefault: 4\n");
    args.support("damp-a", "<omega>\nSets the smoother damping parameter for the A-Solver.\nDefault: 0.5\n");
    args.support("no-multigrid-s", "\nUse PCG-Jacobi instead of Multigrid as S-Solver.\n");
    args.support("max-iter-s", "<N>\nSets the maximum number of allowed iterations for the S-Solver.\nDefault: 100\n");
    args.support("tol-rel-s", "<eps>\nSets the relative tolerance for the S-Solver.\nDefault: 1E-5\n");
    args.support("smooth-s", "<N>\nSets the number of smoothing steps for the S-Solver.\nDefault: 4\n");
    args.support("damp-s", "<omega>\nSets the smoother damping parameter for the S-Solver.\nDefault: 0.5\n");
    args.support("statistics", "[dump]\nEnables general statistics output.\nAdditional parameter 'dump' enables complete statistics dump.\n");
    args.support("test-mode", "\nRuns the application in regression test mode.\n");

    // no arguments given?
    if((argc <= 1) || (args.check("help") >= 0))
    {
      comm.print("\n2D Nonsteady Navier-Stokes CP-Q2/P1dc Toycode Solver (TM)\n");
      comm.print("The easiest way to make this application do something useful is");
      comm.print("to load a pre-defined problem configuration by supplying the");
      comm.print("option '--setup <config>', where <config> may be one of:\n");
      comm.print("  square      Poiseuille-Flow on Unit-Square");
      comm.print("  nozzle      Jet-Flow through Nozzle domain");
      comm.print("  fb_c2d_00   Nonsteady Flow Around A Cylinder (bench1 mesh)");
      comm.print("  fb_c2d_01   Nonsteady Flow Around A Cylinder (32 quad mesh)");
      comm.print("  fb_c2d_02   Nonsteady Flow Around A Cylinder (48 quad mesh)");
      comm.print("  fb_c2d_03   Nonsteady Flow Around A Cylinder (64 quad mesh)");
      comm.print("This will pre-configure this application to solve one of the");
      comm.print("above problems. Note that you can further adjust the configuration");
      comm.print("by specifying additional options to override the default problem");
      comm.print("configuration.");
      if(args.check("help") >= 0)
      {
        comm.print("\nSupported Options:");
        comm.print(args.get_supported_help());
      }
      else
      {
        comm.print("\nUse the option '--help' to display a list of all supported options.\n");
      }
      return;
    }

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());

      // abort
      FEAT::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // parse our configuration
    Config cfg;
    if(!cfg.parse_args(args))
      FEAT::Runtime::abort();

    // create a time-stamp
    TimeStamp time_stamp;

    // let's create our domain
    comm.print("\nPreparing domain...");

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    // parse arguments and set levels
    domain.parse_args(args);
    if(args.check("level") > 0)
      domain.set_desired_levels(args.query("level")->second);
    else
      domain.set_desired_levels(int(cfg.level_max_in), int(cfg.level_min_in));

    domain.create(cfg.mesh_files, cfg.mesh_path);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // store levels after partitioning
    cfg.level_max_in = Index(domain.get_desired_level_max());
    cfg.level_min_in = Index(domain.get_desired_level_min());
    cfg.level_max = Index(domain.max_level_index());
    cfg.level_min = Index(domain.min_level_index());

    // ensure that we do not use UMFPACK if we have more than 1 coarse mesh processes
    if(domain.back_layer().comm().size() > 1)
      cfg.coarse_umfpack_a = cfg.coarse_umfpack_s = false;

    // dump our configuration
    cfg.dump(comm);

    // run our application
    run(cfg, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }
} // namespace NavierStokesCP2D


int main(int argc, char* argv [])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  try
  {
    NavierStokesCP2D::main(argc, argv);
  }
  catch (const std::exception& exc)
  {
    std::cerr << "ERROR: unhandled exception: " << exc.what() << "\n";
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << "\n";
    FEAT::Runtime::abort();
  }
  return 0;
}
