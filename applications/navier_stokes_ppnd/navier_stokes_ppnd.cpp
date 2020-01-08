// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// global MPSC aka PP2D / PP3D (Projection solution by Projection solver) for the Navier stokes equation
//
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
//
//  Purpose: - solver for the nonstationary incompressible Navier Stokes equations by time adaptive projection methods
//             multigrid for linear problems for p
//           - nonlinear iteration for Burgers equations:
//               - fixed point defect correction as outer iteration
//               - mg as solver for the linear auxiliary problems
//               - nonlinear parameter optimization for the correction from the linear solver step
//
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
//
// The main idea of this method is first to compute a velocity field without taking into account incompressibility,
// and then perform a pressure correction, which is a projection back to the subspace of divergence free vector fields.
//
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
//
// This application is a "flow-through-a-domain" solver, i.e. it handles Navier-Stokes
// equations with an inflow and outflow region in 2D and 3D.
//
// WARNING : 3D have not been tested !!!
//
// This application has a few pre-configured benchmark problems, which can be launched
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
// launch flow around cylinder 2D benchmark:
// bench 1: RE = 20;
// bench 2: RE = 100;
// bench 3: RE in [0,100]; (time dependent)
// --setup fb-c2d-01
// --setup fb-c2d-02
// --setup fb-c2d-03
//
// launch flow around cylinder 3D benchmark:
// --setup fb-c3d-01
// --setup fb-c3d-02
// --setup fb-c3d-03
//
// Moreover, this application can be configured by specifying further options, which can
// be used to define new (non-preconfigured) problems or override the pre-defined settings.
//
// IMPORTANT #1:
// In any case, you will need to specify the path to the mesh directory by using the
// '--mesh-path <directory>' option (see below), as the application will fail to find
// the required mesh-files otherwise.
//
// IMPORTANT #2:
// It is not necessary to demand a high accuracy for the velocity solver.
// Instead of increasing the accuracy of solver_a, it's better to increase
// the number of nonlinear / fixpoint iterations.
// It may be necessary to increase the number of time-steps in some cases,
// as the non-linear solver may run amok otherwise...
//
// IMPORTANT #3:
// You can save states by using '--save savename timestep' and load savestates by using '--load savename'
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
// \author Mirco Arndt
//

// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
// Includes
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>

// FEAT-Assembly includes
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>

// FEAT-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_node.hpp>

// FEAT-Global includes
#include <kernel/global/pmdcdsc_matrix.hpp>

// Misc. FEAT includes
#include <kernel/util/runtime.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/string.hpp>

// FEAT-Space includes
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/lagrange2/element.hpp>

// FEAT-Solver includes
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/umfpack.hpp>

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/statistics.hpp>

#include <deque>
#include <numeric>

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

namespace NavierStokesPP
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

    /// minimum, medium maximum levels (as configured)
    String levels_in;

    /// minimum, medium and maximum levels (after partitioning)
    String levels;

    /// base-name of VTK files
    String vtk_name;
    /// stepping of VTK output
    Index vtk_step;

    /// base-name of savefiles
    String save_file;
    /// stepping of savefile output
    Index save_step;

    /// base-name of the load file (former saved file)
    String load;

    // names of inflow mesh-parts
    std::deque<String> part_names_in;

    // names of outflow mesh-part
    std::deque<String> part_names_out;

    // names of noflow mesh-parts
    std::deque<String> part_names_no;

    /// compute "flow-around-a-cylinder" benchmark quantities?
    bool flowbench_c2d;

    bool flowbench_1;
    bool flowbench_2;
    bool flowbench_3;

    // -------------------------------

    /// use deformation tensor?
    bool deformation;

    /// viscosity (for the velocity)
    Real nu;

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

    // use linear extrapolation (with one fixpoint iteration) of u instead a full fixpoint iteration
    bool no_nonlinear;

    // number of fixpoint iterations (for the velocity)
    Index fixpoint_steps;

    // tolerance for the fixpoint iteration (for the velocity)
    Real tol_fix;

    // use multigrid for A-solver ? (velocity)
    bool multigrid_a;

    // use multigrid for S-solver ? (pressure)
    bool multigrid_s;

    // -------------------------------

    // solver a (velocity)

    // maximum number of iterations for velocity
    Index max_iter_a;

    // relative tolerance
    Real tol_rel_a;

    // smoothing steps for velocity
    Index smooth_steps_a;

    // damping parameter for velocity smoother
    Real smooth_damp_a;


    // -------------------------------

    // solver s (pressure)

    // maximum number of iterations for pressure
    Index max_iter_s;

    // relative tolerance for pressure
    Real tol_abs_s;

    // smoothing steps for pressure
    Index smooth_steps_s;

    // damping parameter for pressure smoother
    Real smooth_damp_s;

    // use UMFPACK as coarse grid solver?
    bool coarse_umfpack_s;

    // -------------------------------

    // alpha for the velocity mass matrix
    //
    // Algorithm:
    // S := alpha M + theta K + nu L (burgers equation)
    // S u + k B p = g
    Real alpha;

    // the alpha's for the pressure update
    //
    // Algorithm:
    // p = p_old + alpha_r * q + alpha_d * matrix_mass_pres^(-1) * f;
    Real alpha_r;
    Real alpha_d;

    // use fractional-step-theta-scheme ?
    bool fractional;

    // theta for the theta-scheme
    //
    // 1.0 = Backward-Euler
    // 0.5 = Crank-Nicolson
    Real theta;

    // specifies whether we run in test mode
    bool test_mode;

    // use a different inflow function for testcase square
    // if inflow_time = false; <- use the flow around cylinder inflow function
    bool inflow_time;

    // use Chorin ?
    bool chorin;

    // residuals (stopping criteria)
    Real residual;
    Real residual_drag;  // drag absolute tolerance
    Real residual_drag2; // drag relative tolerance
    Real residual_lift;  // lift absolute tolerance
    Real residual_lift2; // lift relative tolerance

  public :
    // the default configuration
    Config() :
      levels_in(),
      levels(),
      vtk_step(0),
      save_step(0),
      flowbench_c2d(false),
      flowbench_1(false),
      flowbench_2(false),
      flowbench_3(false),
      deformation(false),
      nu(1.0),
      vmax(1.0),
      time_max(0.0),
      time_steps(0),
      max_time_steps(0),
      no_nonlinear(true),
      fixpoint_steps(100),
      tol_fix(1E-2),
      multigrid_a(true),
      multigrid_s(true),
      max_iter_a(100),
      tol_rel_a(0.1),
      smooth_steps_a(2),
      smooth_damp_a(0.5),
      max_iter_s(100),
      tol_abs_s(1E-14),
      smooth_steps_s(4),
      smooth_damp_s(1.0),
      coarse_umfpack_s(false),
      alpha(1.0),                 // non-stationary: alpha = 1, stationary: alpha = 0
      alpha_r(1.0),               // non-stationary: alpha_r <= 1, stationary: alpha_r = 0; (alpha_r = 1, if the preconditioner is exact)
      alpha_d(-1.0),              // non-stationary: alpha_d <= theta * nu * delta_t, stationary: alpha <= nu
      fractional(false),
      theta(1.0),
      test_mode(false),
      inflow_time(false),
      chorin(false),
      residual(1E-32),            // deactivation of the stopping criteria
      residual_drag(1E-32),       // deactivation of the stopping criteria
      residual_drag2(1E-32),      // deactivation of the stopping criteria
      residual_lift(1E-32),       // deactivation of the stopping criteria
      residual_lift2(1E-32)       // deactivation of the stopping criteria
    {
#ifdef FEAT_HAVE_UMFPACK
      coarse_umfpack_s = true;
#endif
    }

    bool parse_args(SimpleArgParser& args)
    {
      // parsing input arguments

      String s;

      // use a default setup (like flow around the cylinder benchmark)
      if ((args.parse("setup", s) > 0) && load.empty())
      {
        if(s.compare_no_case("square") == 0)
          setup_square();
        else if(s.compare_no_case("square_3d") == 0)
          setup_square_3d();
        else if(s.compare_no_case("nozzle") == 0)
          setup_nozzle();
        else if(s.compare_no_case("fb-c2d-01") == 0)
          setup_fb_c2d_01();
        else if(s.compare_no_case("fb-c2d-02") == 0)
          setup_fb_c2d_02();
        else if(s.compare_no_case("fb-c2d-03") == 0)
          setup_fb_c2d_03();
        else if(s.compare_no_case("fb-c3d-01") == 0)
          setup_fb_c3d_01();
        else if(s.compare_no_case("fb-c3d-02") == 0)
          setup_fb_c3d_02();
        else if(s.compare_no_case("fb-c3d-03") == 0)
          setup_fb_c3d_03();
        else
        {
          std::cerr << "ERROR: unknown setup '" << s << "'" << std::endl;
          return false;
        }
      }

      flowbench_c2d |= (args.check("flowbench") >= 0); // this may already be enabled by a setup above
      flowbench_1 |= (args.check("flowbench-1") >= 0);
      flowbench_2 |= (args.check("flowbench-2") >= 0);
      flowbench_3 |= (args.check("flowbench-3") >= 0);
      deformation = (args.check("deformation") >= 0);
      test_mode = (args.check("test-mode") >= 0);

      if (load.empty())
      {
        args.parse("mesh-path", mesh_path);
        if(args.check("mesh-file") > 0)
          mesh_files = args.query("mesh-file")->second;
      }
      if(args.parse("vtk", vtk_name, vtk_step) == 1)
        vtk_step = 1;  // vtk-name given, but not vtk-step, so set to 1
      if(args.parse("save", save_file, save_step) == 1)
        save_step = 1; // save-name given, but not save-step, so set to 1
      args.parse("nu", nu);
      if (load.empty())
      {
        if(args.check("part-in") > 0)
          part_names_in = args.query("part-in")->second;
        if(args.check("part-out") > 0)
          part_names_out = args.query("part-out")->second;
        if(args.check("part-no") > 0)
          part_names_no = args.query("part-no")->second;
      }
      args.parse("vmax", vmax);
      args.parse("time-max", time_max);
      args.parse("time-steps", time_steps);
      if(args.parse("max-time-steps", max_time_steps) < 1)
        max_time_steps = time_steps;
      args.parse("fix-steps", fixpoint_steps);
      no_nonlinear = (args.check("no-nonlinear") < 0);
      multigrid_a = (args.check("no-multigrid-a") < 0);
      multigrid_s = (args.check("no-multigrid-s") < 0);
#ifdef FEAT_HAVE_UMFPACK
      coarse_umfpack_s = (args.check("no-umfpack-s") < 0);
#endif
      args.parse("max-iter-a", max_iter_a);
      args.parse("tol-rel-a", tol_rel_a);
      args.parse("smooth-a", smooth_steps_a);
      args.parse("damp-a", smooth_damp_a);
      args.parse("max-iter-s", max_iter_s);
      args.parse("tol-abs-s", tol_abs_s);
      args.parse("smooth-s", smooth_steps_s);
      args.parse("damp-s", smooth_damp_s);
      args.parse("alpha", alpha);
      args.parse("alpha-r", alpha_r);
      args.parse("alpha-d", alpha_d);
      args.parse("fractional", fractional);
      args.parse("theta", theta);
      args.parse("chorin", chorin);
      args.parse("tol-fix", tol_fix);
      args.parse("residual", residual);
      args.parse("residual-drag", residual_drag);
      args.parse("residual-drag2", residual_drag2);
      args.parse("residual-lift", residual_lift);
      args.parse("residual-lift2", residual_lift2);

      // only 10 time-steps in test mode
      if(test_mode)
        max_time_steps = 10;

      return true;
    }

    void dump(const Dist::Comm& comm)
    {
      // print the configuration summary
      comm.print("Configuration Summary:");
      // =======================================================================================
      comm.print("\nGeneral Setting:");
      dump_line(comm, "Mesh Path", mesh_path);
      dump_line(comm, "Mesh Files", mesh_files);
      dump_line(comm, "Desired Levels", levels_in);
      dump_line(comm, "Chosen Levels", levels);
      dump_line(comm, "VTK-Name", vtk_name);
      dump_line(comm, "VTK-Step", vtk_step);
      dump_line(comm, "Save-File", save_file);
      dump_line(comm, "Save-Step", save_step);
      dump_line(comm, "Load-File", load);
      // =======================================================================================
      comm.print("\nFlow Setting:");
      dump_line(comm, "Flow Benchmark", (flowbench_c2d ? "yes" : "no"));
      dump_line(comm, "Flow Benchmark 1", (flowbench_1 ? "yes" : "no"));
      dump_line(comm, "Flow Benchmark 2", (flowbench_2 ? "yes" : "no"));
      dump_line(comm, "Flow Benchmark 3", (flowbench_3 ? "yes" : "no"));
      dump_line(comm, "Inflow-Parts", part_names_in);
      dump_line(comm, "Outflow-Parts", part_names_out);
      dump_line(comm, "Noflow-Parts", part_names_no);
      dump_line(comm, "Tensor", (deformation ? "Deformation" : "Gradient"));
      dump_line(comm, "Nu", nu);
      dump_line(comm, "V-Max", vmax);
      // =======================================================================================
      comm.print("\nPP Setting:");
      dump_line(comm, "Alpha", alpha);
      dump_line(comm, "Alpha_r", alpha_r);
      if (alpha_d==-1.0)
      {
      dump_line(comm, "Alpha_d", nu*theta*time_max/time_steps);
      }else{
      dump_line(comm, "Alpha_d", alpha_d);
      }
      if (fractional)
        dump_line(comm, "Fractional-Step-Theta-Scheme", "true");
      else
        dump_line(comm, "Theta", theta);
      dump_line(comm, "Chorin", (chorin ? "yes" : "no"));
      dump_line(comm, "Time-Max", time_max);
      dump_line(comm, "Time-Steps", time_steps);
      dump_line(comm, "Max Time-Steps", max_time_steps);
      dump_line(comm, "Delta-t", time_max / time_steps);
       // =======================================================================================
      comm.print("\nSolver Setting:");
      dump_line(comm, "Iteration", (no_nonlinear ? "Linear Extrapolation" : "Fixpoint"));
      if (!no_nonlinear)
      {
        dump_line(comm, "Fixpoint Iteration", fixpoint_steps);
        dump_line(comm, "Tol-fixpoint", tol_fix);
      }
      dump_line(comm, "A: Solver", (multigrid_a ? "Rich-Multigrid" : "Richardson"));
      dump_line(comm, "A: Max-Iter", max_iter_a);
      dump_line(comm, "A: Tol-Rel", tol_rel_a);
      dump_line(comm, "A: Smooth Steps", smooth_steps_a);
      dump_line(comm, "A: Smooth Damp", smooth_damp_a);
      dump_line(comm, "S: Solver", (multigrid_s ? "PCG-Multigrid" : "PCG-Jacobi"));
      dump_line(comm, "S: Max-Iter", max_iter_s);
      dump_line(comm, "S: Tol-Abs", tol_abs_s);
      dump_line(comm, "S: Smooth Steps", smooth_steps_s);
      dump_line(comm, "S: Smooth Damp", smooth_damp_s);
      dump_line(comm, "S: Coarse Solver", (coarse_umfpack_s ? "UMFPACK" : "Richardson-Jacobi"));

      dump_line(comm, "Residual", residual);
      if (flowbench_c2d)
      {
        dump_line(comm, "Residual Drag", residual_drag);
        dump_line(comm, "Residual Drag2", residual_drag2);
        dump_line(comm, "Residual Lift", residual_lift);
        dump_line(comm, "Residual Lift2", residual_lift2);
      }
      dump_line(comm, "\nTest Mode", (test_mode ? "yes" : "no"));
    }

    // 2D variant
    // Setup: Poiseuille-Flow on unit-square
    void setup_square()
    {
      mesh_files.push_back("unit-square-quad.xml");
      part_names_in.push_back("bnd:t");   // top
      part_names_in.push_back("bnd:b");   // bottom
      part_names_in.push_back("bnd:l");   // left
      part_names_out.push_back("bnd:r");  // right
      levels_in = "4 1";
      nu = 1E-3;
      vmax = 0.3;
      time_max = 2.0;
      time_steps = max_time_steps = 2000; // 20000;
      tol_fix = 0.01;                     // tolerance for the fixpoint iteration (velocity)
      tol_abs_s = 1e-14;                  // absolute tolerance for the pressure solver
      inflow_time = true;                         // use a special inflow function
    }

    // 2D variant
    // Setup: nozzle-jet simulation
    void setup_nozzle()
    {
      mesh_files.push_back("nozzle-2-quad.xml");
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
      levels_in = "6 0";
      nu = 1E-3;
      vmax = 1.0;
      time_max = 7.0;
      time_steps = max_time_steps = 35000;
    }

    // 2D variant
    // auxiliary function: basic flow benchmark setup
    void setup_fb_c2d_aux()
    {
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
      part_names_no.push_back("bnd:c");  // circle
      nu = 1E-3;
      time_max = 8.0;
      smooth_damp_a = smooth_damp_s = 0.5;
      flowbench_c2d = true;              // calculate drag lift ...
    }

    // 2D variant
    // Setup: flow around a cylinder (64 quads)
    // DFG flow around cylinder benchmark 2D-1, laminar case Re=20
    void setup_fb_c2d_01()
    {
      setup_fb_c2d_aux();
      mesh_files.push_back("flowbench_c2d_03_quad_64.xml");
      levels_in = "4 0";
      tol_rel_a = 0.5;
      smooth_steps_a = 1;
      smooth_damp_a = 0.01;
      tol_fix = 1E-1;
      tol_rel_a = 0.1;
      tol_abs_s = 1e-14;
      vmax = 0.3;
//      theta = 0.5;                           // Crank-Nicolson-Scheme
      time_max = 1;                      // delta_t = 1
      time_steps = max_time_steps = 10000;
      flowbench_1 = true;
    }

    // 2D variant
    // Setup: flow around a cylinder (64 quads)
    // DFG flow around cylinder benchmark 2D-2, time-periodic case Re=100
    void setup_fb_c2d_02() // (unsteady)
    {
      setup_fb_c2d_aux();
      mesh_files.push_back("flowbench_c2d_03_quad_64.xml");
      levels_in = "5 0";
      vmax = 1.5;
      theta = 0.5;                           // Crank-Nicolson-Scheme
      time_max = 0.35;
//      time_steps = max_time_steps = 35;      // delta_t = 1/100
//      time_steps = max_time_steps = 70;      // delta_t = 1/200
      time_steps = max_time_steps = 140;     // delta_t = 1/400
      flowbench_2 = true;
    }

    // 2D variant
    // Setup: flow around a cylinder (64 quads)
    // DFG flow around cylinder benchmark 2D-3, fixed time interval (Re=100)
    void setup_fb_c2d_03() // (unsteady)
    {
      setup_fb_c2d_aux();
      mesh_files.push_back("flowbench_c2d_03_quad_64.xml");
      levels_in = "3 0";
      tol_rel_a = 0.1;
      smooth_steps_a = 0;
      smooth_damp_a = 0.01;
      max_iter_a = 1000;
      tol_fix = 1E-2;
      tol_rel_a = 0.1;
      tol_abs_s = 1e-14;
      vmax = 1.5;
      theta = 0.5;                           // Crank-Nicolson-Scheme
      time_max = 8.0;
      time_steps = max_time_steps = 8000;
      flowbench_3 = true;
    }


    /// \todo adjust support for "front" and "verso" side for the 3D model
    // 3D variant
    // Setup: Poiseuille-Flow on unit-square
    void setup_square_3d()
    {
      mesh_files.push_back("unit-cube-hexa.xml");
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
      part_names_no.push_back("bnd:f");  // front
      part_names_no.push_back("bnd:n");  // verso
      levels_in = "5 0";
      nu = 1E-3;
      vmax = 0.45;
      time_max = 3.0;
      time_steps = max_time_steps = 12800;
    }

    // 3D variant
    // Setup: flow around a cylinder (32 quads)
    void setup_fb_c3d_01() // (steady)
    {
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
      part_names_no.push_back("bnd:c");  // circle
      part_names_no.push_back("bnd:f");  // front
      part_names_no.push_back("bnd:n");  // verso
      smooth_damp_a = smooth_damp_s = 0.5;
      nu = 1E-3;
      vmax = 0.3;
      time_max = 10000;                      // delta_t = 1
      time_steps = max_time_steps = 10000;
      flowbench_c2d = true;
      mesh_files.push_back("flowbench_c3d_01_hexa_128.xml");
      levels_in = "3 0";
      flowbench_1 = true;
    }

    // 3D variant
    // Setup: flow around a cylinder (32 quads)
    void setup_fb_c3d_02() // (unsteady)
    {
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
      part_names_no.push_back("bnd:c");  // circle
      part_names_no.push_back("bnd:f");  // front
      part_names_no.push_back("bnd:n");  // verso
      smooth_damp_a = smooth_damp_s = 0.5;
      nu = 1E-3;
      vmax = 1.5;
      theta = 0.5;                           // Crank-Nicolson-Scheme
      time_max = 0.35;
//      time_steps = max_time_steps = 35;      // delta_t = 1/100
//      time_steps = max_time_steps = 70;      // delta_t = 1/200
      time_steps = max_time_steps = 140;     // delta_t = 1/400
      flowbench_c2d = true;
      mesh_files.push_back("flowbench_c3d_01_hexa_128.xml");
      levels_in = "3 0";
      flowbench_2 = true;
    }

    // 3D variant
    // Setup: flow around a cylinder (32 quads)
    void setup_fb_c3d_03() // (unsteady)
    {
      part_names_in.push_back("bnd:l");  // left
      part_names_out.push_back("bnd:r"); // right
      part_names_no.push_back("bnd:t");  // top
      part_names_no.push_back("bnd:b");  // bottom
      part_names_no.push_back("bnd:c");  // circle
      part_names_no.push_back("bnd:f");  // front
      part_names_no.push_back("bnd:n");  // verso
      nu = 1E-3;
      vmax = 1.5;
      theta = 0.5;                           // Crank-Nicolson-Scheme
      time_max = 8.0;
      flowbench_c2d = true;
      mesh_files.push_back("flowbench_c3d_01_hexa_128.xml");
      levels_in = "3 0";
      time_steps = max_time_steps = 12800;
      flowbench_3 = true;
      tol_rel_a = 0.1;
      tol_fix = 1E-2;
      tol_rel_a = 0.1;
      tol_abs_s = 1e-14;
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
  }; // class XFluxAccumulator<...>

  //
  // InflowFunction
  //

  template<int dim_>
  class InflowFunction;

  template<>
  class InflowFunction<2> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 2;
    typedef Analytic::Image::Vector<2> ImageType;
    static constexpr bool can_value = true;

  protected:
    Real _vmax;      // maximum inflow velocity
    Real _cur_time;  // time
    bool _setup;     // setup = true <- use flow around cylinder benchmark inflow function
    bool _bench3;

  public:
    explicit InflowFunction(Real vmax, Real cur_time, bool setup, bool bench3) :
      _vmax(vmax),
      _cur_time(cur_time),
      _setup(setup),
      _bench3(bench3)
    {
    }

    template<typename Traits_>
    class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
    {
    protected:
      typedef typename Traits_::DataType DataType;
      typedef typename Traits_::PointType PointType;
      typedef typename Traits_::ValueType ValueType;

      const DataType _vmax, _d, _den, _cur_time;
      bool _setup, _bench3;

    public:
      explicit Evaluator(const InflowFunction& function) :
        _vmax(function._vmax),
        _d(DataType(0.41)),
        _den(_d*_d),
        _cur_time(function._cur_time),
        _setup(function._setup),
        _bench3(function._bench3)
      {
      }

      ValueType value(const PointType& point)
      {
        ValueType val;
        if (_setup)
        {
          // flow around cylinder benchmark inflow function
          const DataType y = point[1];

          if (_bench3)
            val[0] = (_vmax * Math::sin((Math::template pi<DataType>() * DataType(_cur_time)) / DataType(8.0)) * DataType(4) * y * (_d - y)) / _den;
          else
            val[0] = (_vmax * DataType(4) * y * (_d - y)) / _den;
          val[1] = DataType(0);
        }
        else
        {
          // own inflow function
          const DataType y = point[1];

          val[0] = y * (DataType(1) -y); //DataType(_cur_time) * DataType(_cur_time) * DataType(_cur_time) * Math::sin(Math::template pi<DataType>() * y);
          val[1] = DataType(0);
        }
        return val;
      }
    };
  }; // class InflowFunction<2>

  template<>
  class InflowFunction<3> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 3;
    typedef Analytic::Image::Vector<3> ImageType;
    static constexpr bool can_value = true;

  protected:
    Real _vmax;      // maximum inflow velocity
    Real _cur_time;  // time
    bool _setup;     // setup = true <- use flow around cylinder benchmark inflow function
    bool _bench3;

  public:
    explicit InflowFunction(Real vmax, Real cur_time, bool setup, bool bench3) :
      _vmax(vmax),
      _cur_time(cur_time),
      _setup(setup),
      _bench3(bench3)
    {
    }

    template<typename Traits_>
    class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
    {
    protected:
      typedef typename Traits_::DataType DataType;
      typedef typename Traits_::PointType PointType;
      typedef typename Traits_::ValueType ValueType;

      const DataType _vmax, _d, _den, _cur_time;
      bool _setup, _bench3;

    public:
      explicit Evaluator(const InflowFunction& function) :
        _vmax(function._vmax),
        _d(DataType(0.41)),
        _den(_d*_d*_d*_d),
        _cur_time(function._cur_time),
        _setup(function._setup),
        _bench3(function._bench3)
      {
      }

      ValueType value(const PointType& point)
      {
        ValueType val;
        if (_setup)
        {
          // flow around cylinder benchmark inflow function
          const DataType y = point[1];
          const DataType z = point[2];

          if (_bench3)
            val[0] = (_vmax * Math::sin((Math::template pi<DataType>() * _cur_time) / DataType(8.0)) * DataType(16) * y * (_d - y) * z * (_d - z)) / _den;
          else
            val[0] = (_vmax * DataType(16) * y * (_d - y) * z * (_d - z)) / _den;
        }
        else
        {
          // own inflow function
          const DataType y = point[1];

          val[0] = DataType(_cur_time) * DataType(_cur_time) * DataType(_cur_time) * Math::sin(Math::template pi<DataType>() * y);
          val[1] = DataType(0.0);
        }
        return val;
      }
    };
  }; // class InflowFunction


  /**
  * \brief Navier-Stokes System Level class
  *
  * This extends the StokesBlockedSystemLevel by the corresponding filters for
  * the velocity and pressure sub-systems.
  */
  template
  <
    int dim_,
    typename MemType_ = Mem::Main,
    typename DataType_ = Real,
    typename IndexType_ = Index,
    // MatrixBlockA = standard Matrix Block A
    // Algorithm:
    // S_u := alpha M_u + theta K_u + nu_u L_u (burgers equation for the velocity)
    typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
    // MatrixBlockB = standard Matrix Block B
    // gradient Matrix
    typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
    // MatrixBlockD = standard Matrix Block D
    // divergence Matrix
    typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
    typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>
  >
  class NavierStokesBlockedSystemLevel :
    public Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
  {
  public:
    typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;

    typedef typename BaseClass::GlobalVeloVector GlobalVeloVector;
    typedef typename BaseClass::GlobalPresVector GlobalPresVector;
    typedef typename BaseClass::GlobalMatrixBlockB GlobalMatrixBlockB;
    typedef typename BaseClass::GlobalMatrixBlockD GlobalMatrixBlockD;
    typedef Global::Matrix<ScalarMatrix_, typename BaseClass::PresMirror, typename BaseClass::PresMirror> GlobalPresMatrix;

    // (lumped) pressure mass matrix
    // according to the pressure mass matrix structure, no lumping is required
    GlobalPresMatrix matrix_mass_pres;
    GlobalPresVector inverse_lumped_mass_pres;

    // inverse lumped velocity mass matrix
    GlobalVeloVector inverse_lumped_mass_velo;


    // schur matrix
    //
    // Algorithm:
    // P = B^T M_l^(-1) B
    //
    // FEAT Pseudocode:
    // matrix_s = matrix_d * inv_lumped_mass_velo * matrix_b
    //
    // (FEAT : inv_lumped_mass_velo = matrix_s.inv_lumped_matrix_a)
    typedef Global::PMDCDSCMatrix<GlobalMatrixBlockB, GlobalMatrixBlockD> GlobalSchurMatrix;
    GlobalSchurMatrix matrix_s;

    NavierStokesBlockedSystemLevel() :
      inverse_lumped_mass_pres(&this->gate_pres),
      inverse_lumped_mass_velo(&this->gate_velo),
      matrix_s(this->inverse_lumped_mass_velo, this->matrix_b, this->matrix_d)
    {
    }

  }; // class NavierStokesBlockedSystemLevel

  // load the configuration of a given save state (save file)
  void load_conf(Config& cfg, const Dist::Comm& comm, int &processes)
  {
    PropertyMap savefile;
    savefile.read(comm, cfg.load, true);

    // the number of processes (mpi)
    processes = std::stoi(savefile.get_entry("processes").first);
    cfg.part_names_in = savefile.get_entry("cfg.part_names_in").first.split_by_whitespaces();
    cfg.part_names_out = savefile.get_entry("cfg.part_names_out").first.split_by_whitespaces();
    cfg.part_names_no = savefile.get_entry("cfg.part_names_no").first.split_by_whitespaces();
    cfg.levels_in = savefile.get_entry("cfg.levels_in").first;
    cfg.levels = savefile.get_entry("cfg.levels").first;
    cfg.vtk_step = std::stoul(savefile.get_entry("cfg.vtk_step").first, NULL, 0);
    cfg.vtk_name = savefile.get_entry("cfg.vtk_name").first;
    cfg.save_step = std::stoul(savefile.get_entry("cfg.save_step").first, NULL, 0);
    // add '_restart' to the default savefile name and use it this name for the new save file
    cfg.save_file = cfg.load;
    cfg.save_file = cfg.save_file.substr(0, cfg.save_file.length()-4) + "_restart.ini";
    if (savefile.get_entry("cfg.flowbench_c2d").first == "true")
      cfg.flowbench_c2d = true;
    if (savefile.get_entry("cfg.flowbench_1").first == "true")
      cfg.flowbench_1 = true;
    if (savefile.get_entry("cfg.flowbench_2").first == "true")
      cfg.flowbench_2 = true;
    if (savefile.get_entry("cfg.flowbench_3").first == "true")
      cfg.flowbench_3 = true;
    if (savefile.get_entry("cfg.deformation").first == "true")
      cfg.deformation = true;
    cfg.nu = std::stod(savefile.get_entry("cfg.nu").first);
    cfg.vmax = std::stod(savefile.get_entry("cfg.vmax").first);
    cfg.time_max = std::stod(savefile.get_entry("cfg.time_max").first);
    cfg.time_steps = std::stoul(savefile.get_entry("cfg.time_steps").first, NULL, 0);
    cfg.max_time_steps = std::stoul(savefile.get_entry("cfg.max_time_steps").first, NULL, 0);
    cfg.fixpoint_steps = std::stoul(savefile.get_entry("cfg.fixpoint_steps").first, NULL, 0);
    cfg.tol_fix = std::stod(savefile.get_entry("cfg.tol_fix").first);
    if (savefile.get_entry("cfg.multigrid_a").first == "true")
      cfg.multigrid_a = true;
    if (savefile.get_entry("cfg.multigrid_s").first == "true")
      cfg.multigrid_s = true;
    cfg.max_iter_a = std::stoul(savefile.get_entry("cfg.max_iter_a").first, NULL, 0);
    cfg.tol_rel_a = std::stod(savefile.get_entry("cfg.tol_rel_a").first);
    cfg.smooth_steps_a = std::stoul(savefile.get_entry("cfg.smooth_steps_a").first, NULL, 0);
    cfg.smooth_damp_a = std::stod(savefile.get_entry("cfg.smooth_damp_a").first);
    cfg.max_iter_s = std::stoul(savefile.get_entry("cfg.max_iter_s").first, NULL, 0);
    cfg.tol_abs_s = std::stod(savefile.get_entry("cfg.tol_abs_s").first);
    cfg.smooth_steps_s = std::stoul(savefile.get_entry("cfg.smooth_steps_s").first, NULL, 0);
    cfg.smooth_damp_s = std::stod(savefile.get_entry("cfg.smooth_damp_s").first);
    cfg.alpha = std::stod(savefile.get_entry("cfg.alpha").first);
    cfg.alpha_r = std::stod(savefile.get_entry("cfg.alpha_r").first);
    cfg.alpha_d = std::stod(savefile.get_entry("cfg.alpha_d").first);
    if (savefile.get_entry("cfg.fractional").first == "true")
      cfg.fractional = true;
    cfg.theta = std::stod(savefile.get_entry("cfg.theta").first);
    if (savefile.get_entry("cfg.test_mode").first == "true")
      cfg.test_mode = true;
    if (savefile.get_entry("cfg.inflow_time").first == "true")
      cfg.inflow_time = true;
    if (savefile.get_entry("cfg.chorin").first == "true")
      cfg.chorin = true;
    cfg.residual = std::stod(savefile.get_entry("cfg.residual").first);
    cfg.residual_drag = std::stod(savefile.get_entry("cfg.residual_drag").first);
    cfg.residual_drag2 = std::stod(savefile.get_entry("cfg.residual_drag2").first);
    cfg.residual_lift = std::stod(savefile.get_entry("cfg.residual_lift").first);
    cfg.residual_lift2 = std::stod(savefile.get_entry("cfg.residual_lift2").first);
    cfg.mesh_path = savefile.get_entry("cfg.mesh_path").first;
    cfg.mesh_files = savefile.get_entry("cfg.mesh_files").first.split_by_whitespaces();
  }

  // load the data from the save state
  void load_data(Config& cfg, const Dist::Comm& comm, Index &time_step, Real &c_drag_old, Real &c_lift_old)
  {
    PropertyMap savefile;
    savefile.read(comm, cfg.load, true);

    time_step = std::stoul(savefile.get_entry("time_step").first, NULL, 0);

    c_drag_old = std::stod(savefile.get_entry("c_drag_old").first);
    c_lift_old = std::stod(savefile.get_entry("c_lift_old").first);
  }

  // save the configuration in a savefile
  void save(Config cfg, const Dist::Comm& comm, Index time_step,
             Real c_drag_old, Real c_lift_old)
  {
    PropertyMap savefile;

    savefile.add_entry("processes",stringify_fp_sci(comm.size(),8),true);
    savefile.add_entry("time_step",stringify_fp_sci(time_step,8),true);

    savefile.add_entry("cfg.part_names_in",stringify_join(cfg.part_names_in, " "),true);
    savefile.add_entry("cfg.part_names_out",stringify_join(cfg.part_names_out, " "),true);
    savefile.add_entry("cfg.part_names_no",stringify_join(cfg.part_names_no, " "),true);
    savefile.add_entry("cfg.levels_in",stringify(cfg.levels_in),true);
    savefile.add_entry("cfg.levels",stringify(cfg.levels),true);
    savefile.add_entry("cfg.vtk_step",stringify(cfg.vtk_step),true);
    savefile.add_entry("cfg.vtk_name",cfg.vtk_name,true);
    savefile.add_entry("cfg.save_step",stringify(cfg.save_step),true);
    savefile.add_entry("cfg.flowbench_c2d",stringify(cfg.flowbench_c2d),true);
    savefile.add_entry("cfg.flowbench_1",stringify(cfg.flowbench_1),true);
    savefile.add_entry("cfg.flowbench_2",stringify(cfg.flowbench_2),true);
    savefile.add_entry("cfg.flowbench_3",stringify(cfg.flowbench_3),true);
    savefile.add_entry("cfg.deformation",stringify(cfg.deformation),true);
    savefile.add_entry("cfg.nu",stringify_fp_sci(cfg.nu,8),true);
    savefile.add_entry("cfg.vmax",stringify_fp_sci(cfg.vmax,8),true);
    savefile.add_entry("cfg.time_max",stringify_fp_sci(cfg.time_max,8),true);
    savefile.add_entry("cfg.time_steps",stringify(cfg.time_steps),true);
    savefile.add_entry("cfg.max_time_steps",stringify(cfg.max_time_steps),true);
    savefile.add_entry("cfg.fixpoint_steps",stringify(cfg.fixpoint_steps),true);
    savefile.add_entry("cfg.tol_fix",stringify_fp_sci(cfg.tol_fix,8),true);
    savefile.add_entry("cfg.multigrid_a",stringify(cfg.multigrid_a),true);
    savefile.add_entry("cfg.multigrid_s",stringify(cfg.multigrid_s),true);
    savefile.add_entry("cfg.max_iter_a",stringify(cfg.max_iter_a),true);
    savefile.add_entry("cfg.tol_rel_a",stringify_fp_sci(cfg.tol_rel_a,8),true);
    savefile.add_entry("cfg.smooth_steps_a",stringify(cfg.smooth_steps_a),true);
    savefile.add_entry("cfg.smooth_damp_a",stringify_fp_sci(cfg.smooth_damp_a,8),true);
    savefile.add_entry("cfg.max_iter_s",stringify(cfg.max_iter_s),true);
    savefile.add_entry("cfg.tol_abs_s",stringify_fp_sci(cfg.tol_abs_s,8),true);
    savefile.add_entry("cfg.smooth_steps_s",stringify(cfg.smooth_steps_s),true);
    savefile.add_entry("cfg.smooth_damp_s",stringify_fp_sci(cfg.smooth_damp_s,8),true);
    savefile.add_entry("cfg.alpha",stringify_fp_sci(cfg.alpha,8),true);
    savefile.add_entry("cfg.alpha_r",stringify_fp_sci(cfg.alpha_r,8),true);
    savefile.add_entry("cfg.alpha_d",stringify_fp_sci(cfg.alpha_d,8),true);
    savefile.add_entry("cfg.fractional",stringify(cfg.fractional),true);
    savefile.add_entry("cfg.theta",stringify_fp_sci(cfg.theta,8),true);
    savefile.add_entry("cfg.test_mode",stringify(cfg.test_mode),true);
    savefile.add_entry("cfg.inflow_time",stringify(cfg.inflow_time),true);
    savefile.add_entry("cfg.chorin",stringify(cfg.chorin),true);
    savefile.add_entry("cfg.residual",stringify_fp_sci(cfg.residual,8),true);
    savefile.add_entry("cfg.residual_drag",stringify_fp_sci(cfg.residual_drag,8),true);
    savefile.add_entry("cfg.residual_drag2",stringify_fp_sci(cfg.residual_drag2,8),true);
    savefile.add_entry("cfg.residual_lift",stringify_fp_sci(cfg.residual_lift,8),true);
    savefile.add_entry("cfg.residual_lift2",stringify_fp_sci(cfg.residual_lift2,8),true);
    savefile.add_entry("cfg.mesh_path",cfg.mesh_path,true);
    savefile.add_entry("cfg.mesh_files",stringify_join(cfg.mesh_files, " "),true);

    savefile.add_entry("c_drag_old",stringify_fp_sci(c_drag_old,8),true);
    savefile.add_entry("c_lift_old",stringify_fp_sci(c_lift_old,8),true);

    String filename = cfg.save_file + ".ini";

    // save the configuration in the save file 'filename'
    savefile.write(comm, filename);
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
    typedef Mem::Main MemType;
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
    typedef NavierStokesBlockedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = Index(domain.size_physical());

    // create a batch of stop-watches
    StopWatch watch_total, watch_asm_rhs, watch_asm_mat, watch_calc_def,
      watch_sol_init, watch_solver_a, watch_solver_s, watch_solver_p, watch_vtk, watch_save;

    // create stokes and system levels
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("auto-degree:7");

    // compute time-step size
    DataType delta_t = cfg.time_max / DataType(cfg.time_steps);

    /* ***************************************************************************************** */

    comm.print("");
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

    // assemble velocity truncation operators -- we need those for the assembly of the
    // non-linear burgers operators on the coarser levels
    for (Index i(0); i < num_levels; ++i)
    {
      if(i+1 < num_levels)
        system_levels.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature, system_levels.at(i+1).get());
      else if(i+1 < domain.size_virtual())
        system_levels.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature);
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
        auto& loc_mat_a = system_levels.at(i)->matrix_a.local();
        loc_mat_a.format();
        Assembly::Common::IdentityOperatorBlocked<dim> id_op;
        Assembly::BilinearOperatorAssembler::assemble_block_matrix1(loc_mat_a, id_op, domain.at(i)->space_velo, cubature);

        system_levels.at(i)->inverse_lumped_mass_velo = system_levels.at(i)->matrix_a.lump_rows();
        system_levels.at(i)->inverse_lumped_mass_velo.component_invert(system_levels.at(i)->inverse_lumped_mass_velo);
      }

      // perform symbolic initialisation of Schur-complement matrix
      if((i == Index(0)) || cfg.multigrid_s)
      {
        system_levels.at(i)->matrix_s.init_symbolic();
      }
    }

    // assemble pressure mass matrix on finest level
    {
      // get the local pressure mass matrix
      auto& loc_mat_pres = system_levels.front()->matrix_mass_pres.local();

      // assemble matrix structure
      Assembly::SymbolicAssembler::assemble_matrix_std1(loc_mat_pres, domain.front()->space_pres);

      // assemble matrix
      loc_mat_pres.format();
      Assembly::Common::IdentityOperator id_op;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(loc_mat_pres, id_op, domain.front()->space_pres, cubature, DataType(1));
      system_levels.front()->inverse_lumped_mass_pres = system_levels.front()->matrix_mass_pres.lump_rows();
      system_levels.front()->inverse_lumped_mass_pres.component_invert(system_levels.front()->inverse_lumped_mass_pres);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    {
      // our inflow BC function
      InflowFunction<dim> inflow_func(cfg.vmax,delta_t,!cfg.inflow_time,cfg.flowbench_3);

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

        // assemble the filters
        unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, inflow_func);
        unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo);

        // finally, compile the system filter
        system_levels.at(i)->compile_system_filter();

        // apply velocity filter onto inverse lumped mass matrix
        fil_loc_v.filter_cor(system_levels.at(i)->inverse_lumped_mass_velo.local());

        // perform numeric initialisation of Schur-complement matrix
        if((i == Index(0)) || cfg.multigrid_s)
        {
          system_levels.at(i)->matrix_s.init_numeric();
        }
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

    // set up post-processing for flow benchmark
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

      body_force_asm.compile_facets(true);
      flux_u_asm.compile_facets(false);
      flux_l_asm.compile_facets(false);

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
          multigrid_hierarchy_velo->push_level(lvl.matrix_a, lvl.filter_velo, lvl.transfer_velo,
            smoother, smoother, smoother);
        }
        else
        {
          // coarse grid solver
          auto jac = Solver::new_jacobi_precond(lvl.matrix_a, lvl.filter_velo);
          auto cgs = Solver::new_richardson(lvl.matrix_a, lvl.filter_velo, cfg.smooth_damp_a, jac);
          cgs->set_max_iter(cfg.smooth_steps_a);
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
      solver_a = Solver::new_richardson(matrix_a, filter_v, cfg.smooth_damp_a);
    }

    solver_a->set_max_iter(cfg.max_iter_a);
    solver_a->set_tol_rel(cfg.tol_rel_a);

    // for the velocity multigrid/solver, we can only perform symbolic initialisation up to now:
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
    // use the absolute tolerance instead of the relative tolerance
    solver_s->set_tol_rel(0.99); // 'deactivation' of the relative tolerance
    solver_s->set_tol_abs(cfg.tol_abs_s/ (cfg.time_max / DataType(cfg.time_steps)));

    // for the pressure multigrid, we can perform full initialisation:
    if(cfg.multigrid_s)
    {
      multigrid_hierarchy_pres->init();
    }
    solver_s->init();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("");

    // create RHS and SOL vectors
    GlobalVeloVector vec_sol_v = matrix_a.create_vector_l();
    GlobalPresVector vec_sol_p = matrix_s.create_vector_l();
    GlobalVeloVector vec_rhs = matrix_a.create_vector_l();

    GlobalVeloVector vec_f = matrix_a.create_vector_l();
    GlobalVeloVector vec_f_old = matrix_a.create_vector_l();

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
    vec_rhs.format();

    vec_f.format();
    vec_f_old.format();

    // apply filter onto solution vector
    filter_v.filter_sol(vec_sol_v);

    // create solution backup vectors; these store vec_sol_v/p of the last two time-steps
    GlobalVeloVector vec_sol_v_1 = vec_sol_v.clone();
    GlobalVeloVector vec_sol_v_2 = vec_sol_v.clone();
    GlobalPresVector vec_sol_p_1 = vec_sol_p.clone();

    // write header line to console
    {
      const std::size_t nf = stringify_fp_sci(0.0, 3).size();
      String head;
      head += String("step").pad_front(6) + " | ";
      head += String("time").pad_back(12) + " | ";
      head += String("nl").pad_back(1) + " | ";
      head += String("iter").pad_back(4) + " | ";
      head += String("defect velocity").pad_back(21) + " | ";
      head += String("iter").pad_back(4) + " | ";
      head += String("res-pres").pad_back(nf) + " | ";
      head += String("res-velo").pad_back(nf) + " | ";
      if(cfg.flowbench_c2d)
      {
        head += String("drag       | ");
        head += String("lift       | ");
        head += String("p-diff     | ");
        if (cfg.flowbench_2)
        {
          head += String("drag-min   | ");
          head += String("drag-max   | ");
          head += String("drag-mean  | ");
          head += String("drag-amp   | ");
          head += String("lift-min   | ");
          head += String("lift-max   | ");
          head += String("lift-mean  | ");
          head += String("lift-amp   | ");
        }
        else if (cfg.flowbench_3)
        {
          head += String("t(drag-max)| ");
          head += String("drag-max   | ");
          head += String("t(lift-max)| ");
          head += String("lift-max   | ");
        }
        else
        {
          head += String("drag-abs   | ");
          head += String("drag-rel   | ");
          head += String("lift-abs   | ");
          head += String("lift-rel   | ");
        }
      }
      head += String("runtime    ");
      comm.print(head);
      comm.print(String(head.size(), '-'));
    }

    Statistics::reset();
    watch_total.start();

    // keep track whether something failed miserably...
    bool failure = false;

    Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_rhs;
    Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_mat;

    //
    // theta scheme
    //

    // set up a burgers assembler for the RHS vector
    // Algorithm: g = M u^n + theta k F^(n+1) + (1 - theta) k [F^(n) - (nu L + K)]
    // FEAT: g = M + (1 - theta) delta_t [- (nu L + K)]
    burgers_rhs.deformation = cfg.deformation;
    burgers_rhs.nu = -(DataType(1) - cfg.theta) * delta_t * cfg.nu;
    burgers_rhs.beta = -(DataType(1) - cfg.theta) * delta_t;
    burgers_rhs.theta = cfg.alpha;

    // set up a burgers assembler for the velocity matrix
    // Algorithm: S(u^(n+1)) = M + k[K(u^(n+1)) + nu L]
    // FEAT: S = alpha M + theta delta_t [K + nu L]
    burgers_mat.deformation = cfg.deformation;
    burgers_mat.nu = cfg.theta * delta_t * cfg.nu;
    burgers_mat.beta =  cfg.theta * delta_t;
    burgers_mat.theta = cfg.alpha;

    //*********************************************************************************//

    auto vec_def_v_old_1 = vec_def_v.clone();

    // set alpha_d, if alpha_d was not given by the user
    if (cfg.alpha_d==-1.0)
    {
      cfg.alpha_d = cfg.theta * delta_t * cfg.nu; // alpha_d <= theta * delta_t * nu  (<= k / Reynolds)
    }

    if (cfg.no_nonlinear == true)
    {
      cfg.fixpoint_steps = 1;
    }

    // body forces, pressure difference and flux values
    DataType c_drag(0), c_lift(0), p_diff(0);
    DataType c_drag_time(0), c_lift_time(0);
    DataType c_lift_old(0), c_drag_old(0);
    DataType c_drag_min(1E+32), c_drag_max(-1E+32), c_lift_min(1E+32), c_lift_max(-1E+32);

    Index t_step(0);
    if (!cfg.load.empty())
    {
      load_data(cfg, comm, t_step, c_drag_old, c_lift_old);
      vec_sol_v_2.local().read_from(LAFEM::FileMode::fm_dvb, cfg.load.substr(0, cfg.load.length() - 4) + "_v2_" + stringify(rank));
      vec_sol_v_1.local().read_from(LAFEM::FileMode::fm_dvb, cfg.load.substr(0, cfg.load.length() - 4) + "_v1_" + stringify(rank));
      vec_sol_v.copy(vec_sol_v_1);
      vec_sol_p_1.local().read_from(LAFEM::FileMode::fm_dv, cfg.load.substr(0, cfg.load.length() - 4) + "_p1_" + stringify(rank));
      vec_sol_p.copy(vec_sol_p_1);
      // save given rhs vector (if time dependent)
      vec_f.local().read_from(LAFEM::FileMode::fm_dvb, cfg.load.substr(0, cfg.load.length() - 4) + "_f_" + stringify(rank));
      vec_f_old.local().read_from(LAFEM::FileMode::fm_dvb, cfg.load.substr(0, cfg.load.length() - 4) + "_f_old_" + stringify(rank));
    }

    // fractional_step = 1 means use theta sheme
    Index  fractional_step(1);
    if (cfg.fractional)
      fractional_step = 3;

    // time-step loop
    for(Index time_step(t_step+1); time_step <= cfg.max_time_steps; ++time_step)
    {
      String line;
      // compute current time
      DataType cur_time = DataType(time_step) * delta_t;

      for (Index fractional = 1; fractional <= fractional_step; ++fractional)
      {

        // is inflow function time dependent?
        // -> assemble in inflow function in each time step
        if (cfg.inflow_time || cfg.flowbench_3)
        {
          /// Assembling time dependent system filters...
          for(Index i(0); i < num_levels; ++i)
          {
            Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow;
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
            typename SystemLevelType::LocalVeloFilter& fil_loc_v = system_levels.at(i)->filter_velo.local();
            InflowFunction<dim> inflow_func(cfg.vmax,cur_time,!cfg.inflow_time,cfg.flowbench_3);
            unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, inflow_func);
          }
        } // end if (Assemble inflow function)

        // if vec_f time dependent, update vec_f and vec_f_old (aka F^(n+1) F^n)
        vec_f_old.copy(vec_f);
        vec_f.format();

        if (cfg.fractional)
        {
          // Choose theta such that the scheme is of second order. See Turek book, page 151.
          DataType theta(DataType(1.0) - Math::sqrt(2.0) / DataType(2.0));
          DataType thetas(DataType(1.0) - DataType(2.0) * theta);                                      // thetas stands for \theta'
          DataType alpha_t((DataType(1.0) - DataType(2.0) * theta) / (DataType(1.0) - theta));         // => all coefficient matrices will be the same in all substeps.
          DataType beta_t(DataType(1.0) - alpha_t);
          DataType theta1(0),theta2(0),theta3(0),theta4(0);

          switch (fractional)
          {
            case 1:
              // Fractional-step-theta-scheme
              theta1 = alpha_t * theta * delta_t;                      // K = delta_t*3.0
              theta2 = beta_t * theta * delta_t;
              theta3 = theta * delta_t;                                // rhs
              theta4 = theta * delta_t;                                // pressure
              break;
            case 2:
              // Fractional-step-theta-scheme
              theta1 = beta_t * thetas * delta_t;
              theta2 = alpha_t * thetas * delta_t;
              theta3 = thetas * delta_t;                               // rhs
              theta4 = thetas * delta_t;                               // pressure
              break;
            case 3:
              // Fractional-step-theta-scheme
              theta1 = alpha_t * theta * delta_t;                      // K = delta_t*3.0
              theta2 = beta_t * theta * delta_t;
              theta3 = theta * delta_t;                                // rhs
              theta4 = theta * delta_t;                                // pressure
              break;
          }
          // set up a burgers assembler for the RHS vector
          // Algorithm: g = M u^n - theta2 k (nu L + K) + theta3 k F^(n)
          // F^(n+1) = F^(n) = 0
          // FEAT: g = M + theta2 [- (nu L + K)]
          burgers_rhs.deformation = cfg.deformation;
          burgers_rhs.nu = -theta2 * cfg.nu;
          burgers_rhs.beta = -theta2;
          burgers_rhs.theta = DataType(1);

          // set up a burgers assembler for the velocity matrix
          // Algorithm: S(u^(n+1)) = M + theta1 k[K(u^(n+1)) + nu L]
          // FEAT: S = alpha M + theta1 [K + nu L]
          burgers_mat.deformation = cfg.deformation;
          burgers_mat.nu = theta1 * cfg.nu;
          burgers_mat.beta =  theta1;
          burgers_mat.theta = cfg.alpha;

          cur_time += theta3;

          //
          // assemble RHS vector
          //
          // Algorithm: f := g - k B p^(l-1)
          watch_asm_rhs.start();
          vec_rhs.format();
          // assemble vector vec_rhs aka g
          //
          // Algorithm: g := M u^n + theta2 k [- (nu L + K)] + theta3 k F^(n+1)
          // now : g := M u^n - (1 - theta) k [- (nu L + K)]
          // later : support for the given rhs F^(n+1)
          burgers_rhs.assemble_vector(
            vec_rhs.local(), vec_sol_v.local(), vec_sol_v.local(), the_domain_level.space_velo, cubature);

          // FEAT: rhs_v <- rhs_v + (-theta4) * matrix_b * sol_p
          if (!cfg.chorin)
            matrix_b.local().apply(vec_rhs.local(), vec_sol_p.local(), vec_rhs.local(), -theta4);

          vec_rhs.sync_0();

          // add support for the given righthandside vec_f (aka F^(n+1))
          vec_rhs.axpy(vec_f, vec_rhs,theta3);
        }
        else
        {
          //
          // assemble RHS vector
          //
          // Algorithm: f := g - k B p^(l-1)
          watch_asm_rhs.start();
          vec_rhs.format();
          // assemble vector vec_rhs aka g
          //
          // Algorithm: g := M u^n - (1 - theta) k [K(u^(n+1)) + nu L] + theta k F^n+1 + (1 - theta) k F^n
          // now : g := M u^n - (1 - theta) k [K(u^(n+1)) + nu L]
          // later : support for the given rhs F^(n+1) and F^n
          burgers_rhs.assemble_vector(
            vec_rhs.local(), vec_sol_v.local(), vec_sol_v.local(), the_domain_level.space_velo, cubature);

          // FEAT: rhs <- rhs + (-delta_t) * matrix_b * sol_p
          if (!cfg.chorin)
            matrix_b.local().apply(vec_rhs.local(), vec_sol_p.local(), vec_rhs.local(), -delta_t);

          vec_rhs.sync_0();

          // add support for the given righthandside vec_f and vec_f_old (aka F^(n+1) and F^n)
          vec_rhs.axpy(vec_f, vec_rhs, cfg.theta*delta_t);
          vec_rhs.axpy(vec_f_old, vec_rhs, (DataType(1.0) - cfg.theta)*delta_t);
        }
        watch_asm_rhs.stop();
        // apply RHS filter
        filter_v.filter_rhs(vec_rhs);

        // linear extrapolation of solution in time
        if((time_step > Index(2)) && (cfg.no_nonlinear))
        {
          vec_conv.scale(vec_sol_v_1, DataType(2));
          vec_conv.axpy(vec_sol_v_2, vec_conv, -DataType(1));
        }
        else
        {
          // constant extrapolation of solution in time
          vec_conv.copy(vec_sol_v);
        }

        // Iteration number of solver_a
        Index iter_a(0);

        // non-linear iteration (fixpoint iteration)
        for(Index nonlin_step(0); ; ++nonlin_step)
        {
          //
          // loop over all levels and assemble the burgers matrices
          //
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

          //
          // compute defect d
          //
          // Algorithm : d^(n-1) = - T(tilde_u^(n-1))tilde_u^(n-1) + f
          //                     = - S(tilde_u^(n-1))tilde_u^(n-1) + f
          // -> def_v = - A tilde_u + f_v
          watch_calc_def.start();
          vec_def_v.format();
          matrix_a.apply(vec_def_v, vec_conv);
          vec_def_v.axpy(vec_def_v,vec_rhs,-DataType(1));
          filter_v.filter_def(vec_def_v);
          if (nonlin_step == 0)
            vec_def_v_old_1.copy(vec_def_v);
          watch_calc_def.stop();

          Real norm_old_1 = vec_def_v_old_1.norm2();
          Real norm_v = vec_def_v.norm2();

          // console output
          if(nonlin_step > 0)
          {
            line += stringify(time_step).pad_front(6) + " | ";
            line += stringify_fp_fix(cur_time, 8).pad_front(12) + " | ";
            line += stringify(nonlin_step).pad_front(2) + " | ";
            line += stringify_fp_sci(iter_a, 3).pad_front(4) +  " | ";
            line += stringify_fp_sci(norm_old_1, 3) + " > ";
            line += stringify_fp_sci(norm_v, 3) + " | ";
          }

          if (((Math::abs(norm_v/norm_old_1)  <= cfg.tol_fix) && (nonlin_step > 0)) || (nonlin_step >= cfg.fixpoint_steps))
          {
            break;
          }

          //
          //  initialise linear solvers
          //
          watch_sol_init.start();
          if(cfg.multigrid_a)
            multigrid_hierarchy_velo->init_numeric();
          solver_a->init_numeric();
          watch_sol_init.stop();

          //
          // solve velocity system
          //
          // Algorithm : tilde_T(tilde_u^(n-1)) y^(n-1) = T(tilde_u^(n-1)) y^(n-1) = d^(n-1)
          // FEAT: A * y^(n-1) = d^(n-1)
          //       matrix_a * cor_v = def_v
          vec_cor_v.format();
          watch_solver_a.start();
          // solve for vec_cor_v aka y^(n-1)
          Solver::Status status_a = solver_a->apply(vec_cor_v, vec_def_v);
          watch_solver_a.stop();
          if(!Solver::status_success(status_a))
          {
            comm.print(std::cerr, "\n\nERROR: velocity solver broke down!\n");
            failure = true;
            break;
          }
          // release linear solvers
          solver_a->done_numeric();
          if(cfg.multigrid_a)
            multigrid_hierarchy_velo->done_numeric();
          iter_a = solver_a->get_num_iter();

          vec_conv.axpy(vec_cor_v, vec_conv);
          // apply filter onto solution vector
          filter_v.filter_sol(vec_conv);

          // console output
          if ((nonlin_step <= cfg.fixpoint_steps) && (nonlin_step > 0))
          {
            comm.print(line);
            line = "";
          }

        } // non-linear iteration

        if (failure)
          break;

        //
        // calculate the right hand side for the Pressure-Poisson-Problem
        //
        // Algorithm : f_p = 1/k B^T tilde_u^l
        // FEAT : f = 1/delta_t D tilde_u
        vec_def_p.format();
        matrix_d.apply(vec_def_p, vec_conv);
        vec_def_p.scale(vec_def_p,DataType(1) / delta_t);
        filter_p.filter_def(vec_def_p);

        //
        // Solve the discrete Pressure-Poisson-Problem for q
        //
        // Algorithm : P q = f_p
        //
        //
        auto vec_q = vec_sol_p.clone();
        vec_q.format();
        watch_solver_s.start();
        // solve for vec_q
        Solver::Status status_s = solver_s->apply(vec_q, vec_def_p);
        watch_solver_s.stop();
        if(!Solver::status_success(status_s))
        {
          comm.print(std::cerr, "\n\nERROR: pressure solver broke down!\n");
          failure = true;
          break;
        }
        Index iter_p = solver_s->get_num_iter();

        Real norm_p = vec_def_p.norm2();
        // console output
        {
          line += stringify_fp_sci(iter_p, 3).pad_front(4) +  " | ";
          line += stringify_fp_sci(norm_p, 3) + " | ";
        }

        //
        // Update the new pressure
        //
        // Algorithm : p^l = p^(l-1) + alpha_r q + alpha_d M_p^(-1) f_p

        // update the pressure part I
        //
        // p = p_old + alpha_r q
        //
        // chorin: (means p^(l-1) := 0)
        // p = alpha_r q
        if (cfg.chorin)
          vec_sol_p.scale(vec_q,cfg.alpha_r);
        else
          vec_sol_p.axpy(vec_q,vec_sol_p,cfg.alpha_r);

        // update the pressure part II
        // p = p_old + alpha_d f      (= p_old + alpha_d M_p^(-1) f_p)
        if (cfg.alpha_d > 0)
        {
          vec_def_p.component_product(the_system_level.inverse_lumped_mass_pres, vec_def_p);
          vec_sol_p.axpy(vec_def_p,vec_sol_p,cfg.alpha_d);
        }
        filter_p.filter_sol(vec_sol_p);

        //
        // Update the new velocity
        //
        // Algorithm : u^l = tilde_u^l - k M_l^(-1) B q
        auto vec_tmp = vec_sol_v.clone();
        vec_tmp.format();
        // tmp = B q
        matrix_b.apply(vec_tmp,vec_q);
        // tmp = M_l^(-1) B q
        vec_tmp.component_product(the_system_level.inverse_lumped_mass_velo, vec_tmp);
        // u^l = -k tmp + tilde_u^l = -k M_l^(-1) B q + tilde_u^l
        vec_sol_v.axpy(vec_tmp,vec_conv,-delta_t);
        // apply filter onto solution vector
        filter_v.filter_sol(vec_sol_v);

        // body forces, pressure difference and flux values
        c_drag = DataType(0);
        c_lift = DataType(0);
        p_diff = DataType(0);

        if(cfg.flowbench_c2d)
        {

          // p
          auto vec_sol_p_mid = vec_sol_p.clone();
          // midpoint of p
          //  vec_sol_p_mid.axpy(vec_sol_p, vec_sol_p_1);
          //  vec_sol_p_mid.scale(vec_sol_p_mid, DataType(0.5));
          // extrapolation
          vec_sol_p_mid.axpy(vec_sol_p_1, vec_sol_p, -DataType(1));
          vec_sol_p_mid.axpy(vec_sol_p_mid, vec_sol_p, DataType(0.5));

          {
          // assemble drag and lift forces
            BenchBodyForceAccumulator<DataType> body_force_accum(false, cfg.nu, cfg.vmax);
            body_force_asm.assemble_flow_accum(
              body_force_accum,
              vec_sol_v_1.local(),
              vec_sol_p_mid.local(),
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
              point_pa, vec_sol_p_mid.local(), the_domain_level.space_pres);
            auto pval_e = Assembly::DiscreteEvaluator::eval_fe_function(
              point_pe, vec_sol_p_mid.local(), the_domain_level.space_pres);

            // compute pressure mean
            const auto p_a = pval_a.mean_value_dist(comm);
            const auto p_e = pval_e.mean_value_dist(comm);
            p_diff = p_a - p_e;
          }
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
          LAFEM::DenseVector<Mem::Main, double, Index> vtx_p, vtx_der_p;
          GlobalPresVector vec_sol_p_post = vec_sol_p.clone();
          vec_sol_p_post.axpy(vec_sol_p,vec_sol_p_1);
          vec_sol_p_post.scale(vec_sol_p_post,DataType(0.5));
          Assembly::DiscreteCellProjector::project(vtx_p, vec_sol_p.local(), the_domain_level.space_pres, cub);

          // write pressure
          vtk.add_cell_scalar("p", vtx_p.elements());

          // compute and write time-derivatives
          GlobalVeloVector vec_der_v = vec_sol_v.clone();
          GlobalPresVector vec_der_p = vec_sol_p.clone();
          vec_der_v.axpy(vec_sol_v_1, vec_der_v, -DataType(1));
          vec_der_p.axpy(vec_sol_p_1, vec_der_p, -DataType(1));
          vec_der_v.scale(vec_der_v, DataType(1) / delta_t);
          vec_der_p.scale(vec_der_p, DataType(1) / delta_t);
          Assembly::DiscreteCellProjector::project(vtx_der_p, vec_der_p.local(), the_domain_level.space_pres, cub);

          vtk.add_vertex_vector("v_dt", vec_der_v.local());
          vtk.add_cell_scalar("p_dt", vtx_der_p.elements());

          // export
          vtk.write(vtk_path, comm);
          watch_vtk.stop();
        }

        //
        // compute final Residual (stopping criteria)
        //
        // B^T u^(n+1) = B^T u - delta_t B^T M_l^(-1) B q
        //             = delta_t f_p - delta_t P q
        //             = 0 (residual)
        //
        auto vec_def_res = vec_def_p.clone();
        vec_def_res.format();
        matrix_d.apply(vec_def_res,vec_sol_v);

        DataType norm_res(vec_def_res.norm2());

        // console output
        {
          line += stringify_fp_sci(norm_res, 3) + " | ";
        }


        // console output
        if (cfg.flowbench_c2d)
        {
          line += stringify_fp_sci(c_drag, 3).trunc_back(10u).pad_front(10u) + " | ";
          line += stringify_fp_sci(c_lift, 3).trunc_back(10u).pad_front(10u) + " | ";
          line += stringify_fp_sci(p_diff, 3).trunc_back(10u).pad_front(10u) + " | ";
          if (cfg.flowbench_2)
          {
            // Drag
            c_drag_min = Math::min(c_drag_min, c_drag);
            c_drag_max = Math::max(c_drag_max, c_drag);
            c_lift_min = Math::min(c_lift_min, c_lift);
            c_lift_max = Math::max(c_lift_max, c_lift);
            line += stringify_fp_sci(c_drag_min, 3).trunc_back(10u).pad_front(10u) + " | "; //min
            line += stringify_fp_sci(c_drag_max, 3).trunc_back(10u).pad_front(10u) + " | "; //max
            line += stringify_fp_sci((c_drag_max + c_drag_min)/2.0, 3).trunc_back(10u).pad_front(10u) + " | "; //mean
            line += stringify_fp_sci(c_drag_max - c_drag_min, 3).trunc_back(10u).pad_front(10u) + " | "; //amp
            // Lift
            line += stringify_fp_sci(c_lift_min, 3).trunc_back(10u).pad_front(10u) + " | "; // min
            line += stringify_fp_sci(c_lift_max, 3).trunc_back(10u).pad_front(10u) + " | "; //max
            line += stringify_fp_sci((c_lift_max + c_lift_min)/2.0, 3).trunc_back(10u).pad_front(10u) + " | "; // mean
            line += stringify_fp_sci(c_lift_max - c_lift_min, 3).trunc_back(10u).pad_front(10u) + " | "; //amp
          }
          else if (cfg.flowbench_3)
          {
            if (c_drag > c_drag_max)
            {
              c_drag_max = c_drag;
              c_drag_time = cur_time;
            }
            if (c_lift > c_lift_max)
            {
              c_lift_max = c_lift;
              c_lift_time = cur_time;
            }
            line += stringify_fp_sci(c_drag_time, 3).trunc_back(10u).pad_front(10u) + " | "; //t(drag_max)
            line += stringify_fp_sci(c_drag_max, 3).trunc_back(10u).pad_front(10u) + " | "; //drag_max
            line += stringify_fp_sci(c_lift_time, 3).trunc_back(10u).pad_front(10u) + " | "; //t(lift_max)
            line += stringify_fp_sci(c_lift_max, 3).trunc_back(10u).pad_front(10u) + " | "; //lift_max
          }
          else
          {
            line += stringify_fp_sci(Math::abs(c_drag - c_drag_old), 3).trunc_back(10u).pad_front(10u) + " | ";
            line += stringify_fp_sci(Math::abs((c_drag - c_drag_old)/c_drag), 3).trunc_back(10u).pad_front(10u) + " | ";
            line += stringify_fp_sci(Math::abs(c_lift - c_lift_old), 3).trunc_back(10u).pad_front(10u) + " | ";
            line += stringify_fp_sci(Math::abs((c_lift - c_lift_old)/c_lift), 3).trunc_back(10u).pad_front(10u) + " | ";
          }
        }

        // console output
        {
          line += stamp_start.elapsed_string_now();
          comm.print(line);
          line = "";
        }

        // finally, update our solution vector backups
        vec_sol_v_2.copy(vec_sol_v_1);
        vec_sol_v_1.copy(vec_sol_v);
        vec_sol_p_1.copy(vec_sol_p);

        // compress all statistics from the current timestep for further analysis after the solution is finished
        FEAT::Statistics::compress_solver_expressions();
      } // fractional loop

      // final Residual
      //
      // compute final Residual (stopping criteria)
      //
      // B^T u^(n+1) = B^T u - delta_t B^T M_l^(-1) B q
      //             = delta_t f_p - delta_t P q
      //             = 0 (residual)
      //
      auto vec_def_res = vec_def_p.clone();
      vec_def_res.format();
      matrix_d.apply(vec_def_res,vec_sol_v);
      Real norm_res = vec_def_res.norm2();

      //
      // stopping criteria
      //

      // final Residual
      if (norm_res <= cfg.residual)
      {
        comm.print("\n'Residual'-stopping criteria reached " + stringify_fp_fix(norm_res, 3) + "\n");
        break;
      }
      if ((cfg.flowbench_c2d) && (time_step>2))
      {
        // forces : drag
        if (Math::abs(c_drag - c_drag_old) <= cfg.residual_drag)
        {
          comm.print("\n'Drag'-stopping criteria reached " + stringify_fp_fix(Math::abs((c_drag - c_drag_old)/c_drag), 3) + "\n");
          break;
        }
        // forces : drag-time
        if (Math::abs((c_drag - c_drag_old)/c_drag) <= cfg.residual_drag2)
        {
          comm.print("\n'Drag-time'-stopping criteria reached " + stringify_fp_fix(Math::abs((c_drag - c_drag_old)/delta_t), 3) + "\n");
          break;
        }
        // forces : lift
        if (Math::abs(c_lift - c_lift_old) <= cfg.residual_lift)
        {
          comm.print("\n'Lift'-stopping criteria reached " + stringify_fp_fix(Math::abs((c_lift - c_lift_old)/c_lift), 3) + "\n");
          break;
        }
        // forces : lift-time
        if (Math::abs((c_lift - c_lift_old)/c_lift) <= cfg.residual_lift2)
        {
          comm.print("\n'Lift-time'-stopping criteria reached " + stringify_fp_fix(Math::abs((c_lift - c_lift_old)/delta_t), 3) + "\n");
          break;
        }
      }

      // finally, update our backups vectors

      c_drag_old = c_drag;
      c_lift_old = c_lift;

      //write out save
      if((cfg.save_step > 0) && (time_step % cfg.save_step == 0))
      {
        watch_save.start();

        // get vector data
        vec_sol_v_1.local().write_out(LAFEM::FileMode::fm_dvb, cfg.save_file + "_v1_" + stringify(rank));
        vec_sol_v_2.local().write_out(LAFEM::FileMode::fm_dvb, cfg.save_file + "_v2_" + stringify(rank));
        vec_sol_p_1.local().write_out(LAFEM::FileMode::fm_dv, cfg.save_file + "_p1_" + stringify(rank));
        vec_f.local().write_out(LAFEM::FileMode::fm_dvb, cfg.save_file + "_f_" + stringify(rank));
        vec_f_old.local().write_out(LAFEM::FileMode::fm_dvb, cfg.save_file + "_f_old_" + stringify(rank));
        save(cfg, comm, time_step, c_drag_old, c_lift_old);
        comm.print("Save State: " + cfg.save_file + ".ini");
        watch_save.stop();
      }

      // continue with next time-step
    } // time-step loop

    watch_total.stop();

    // are we in test-mode?
    if(cfg.test_mode)
    {
      // (on rosetyler)
      // LVL 1 : linear extrapolation
      // mpirun -n 12 applications/navier_stokes_ppnd/navier_stokes_ppnd --setup fb-c2d-03 --mesh-path ~/FEAT/feat3/data/meshes/  --level 1 --time-steps 1600  --tol-rel-a 0.99   --test-mode
      // 10 |   0.05000000 |  1 |    1 | 6.388e-05 > 6.292e-05 |   87 | 1.034e-03 | 5.379e-15 |  3.045e-01 | -6.638e-04 |  1.714e-01 |  1.000e-02 |  3.492e-01 |  5.000e-03 | -4.636e-04 | 0.148
      // LVL 1: full fixpoint iteration
      // mpirun -n 12 applications/navier_stokes_ppnd/navier_stokes_ppnd --setup fb-c2d-03 --mesh-path ~/FEAT/feat3/data/meshes/  --level 1 --time-steps 1600  --tol-rel-a 0.99 --test-mode --no-nonlinear --fix-steps 400
      // 10 |   0.05000000 | 311 |    1 | 7.178e-05 > 7.099e-07 |   84 | 4.511e-03 | 9.428e-15 |  1.839e-01 | -2.967e-04 |  1.066e-01 |  5.000e-03 |  2.082e-01 |  1.000e-02 | -2.905e-04 | 1.085

      if (cfg.no_nonlinear)
      {
        if ((c_drag < 3.035e-01) || (c_drag > 3.055e-01)
              || (c_lift < -6.65e-04) || (c_lift > -6.62e-4)
              || (p_diff < 1.70e-01) || (p_diff > 1.73e-01)
              || (c_drag_time != 1.00e-02) || (c_lift_time != 5.00e-03))
                failure = true;
      }
      else
      {
        if ((c_drag < 1.82e-01) || (c_drag > 1.85e-01)
              || (c_lift < -2.97e-04) || (c_lift > -2.96e-4)
              || (p_diff < 1.00e-01) || (p_diff > 1.10e-01)
              || (c_drag_time != 5.00e-03) || (c_lift_time != 1.00e-02))
                failure = true;
      }

      if(failure)
        comm.print(std::cerr, "\nTest-Mode: FAILED");
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

    double t_total = watch_total.elapsed();
    double t_asm_mat = watch_asm_mat.elapsed();
    double t_asm_rhs = watch_asm_rhs.elapsed();
    double t_calc_def = watch_calc_def.elapsed();
    double t_sol_init = watch_sol_init.elapsed();
    double t_solver_a = watch_solver_a.elapsed();
    double t_solver_s = watch_solver_s.elapsed();
    double t_solver_p = watch_solver_p.elapsed();
    double t_vtk = watch_vtk.elapsed();
    double t_save = watch_save.elapsed();
    double t_sum = t_asm_mat + t_asm_rhs + t_calc_def + t_sol_init + t_solver_a + t_solver_s + t_solver_p + t_vtk + t_save;

    // write timings
    {
      comm.print("");
      dump_time(comm, "Total Solver Time", t_total, t_total);
      dump_time(comm, "Matrix Assembly Time", t_asm_mat, t_total);
      dump_time(comm, "Vector Assembly Time", t_asm_rhs, t_total);
      dump_time(comm, "Defect-Calc Time", t_calc_def, t_total);
      dump_time(comm, "Solver-A Init Time", t_sol_init, t_total);
      dump_time(comm, "Solver-A Time", t_solver_a, t_total);
      dump_time(comm, "Solver-S Time", t_solver_s, t_total);
      dump_time(comm, "Solver-P Time", t_solver_p, t_total);
      dump_time(comm, "VTK-Write Time", t_vtk, t_total);
      dump_time(comm, "Save-Write Time", t_save, t_total);
      dump_time(comm, "Other Time", t_total-t_sum, t_total);
    }
  }

  template<typename Shape_>
  void domain_and_run(SimpleArgParser& args, Config& cfg, const Dist::Comm& comm)
  {
    // create a time-stamp
    TimeStamp time_stamp;

    // let's create our domain
    comm.print("\nPreparing domain...");

    typedef Shape_ ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    // parse arguments and set levels
    domain.parse_args(args);
    if ((args.check("level") > 0) && cfg.load.empty())
      domain.set_desired_levels(args.query("level")->second);
    else
      domain.set_desired_levels(cfg.levels_in);

    domain.create(cfg.mesh_files, cfg.mesh_path);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // store levels after partitioning
    if (cfg.load.empty())
    {
      cfg.levels_in = domain.format_desired_levels();
      cfg.levels = domain.format_chosen_levels();
    }

    // ensure that we do not use UMFPACK if we have more than 1 coarse mesh processes
    if(domain.back_layer().comm().size() > 1)
      cfg.coarse_umfpack_s = false;

    // dump our configuration
    cfg.dump(comm);

    // run our application
    run(cfg, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("help", "\nDisplays this help message.\n");
    args.support("setup", "<config>\nLoads a pre-defined configuration:\n"
      "square      Poiseuille-Flow on Unit-Square\n"
      "nozzle      Jet-Flow through Nozzle domain\n"
      "fb-c2d-00   Nonsteady Flow Around A Cylinder (bench1 mesh)\n"
      "fb-c2d-01   Nonsteady Flow Around A Cylinder (32 quad mesh)\n"
      "fb-c2d-02   Nonsteady Flow Around A Cylinder (48 quad mesh)\n"
      "fb-c2d-03   Nonsteady Flow Around A Cylinder (64 quad mesh)\n"
    );
    args.support("alpha", "<alpha>\nSets the scale parameter alpha.\nThis scale the velocity mass matrix.\nDefault = 1\n");
    args.support("alpha-r", "<alpha_r>\nSets the scale parameter alpha-r (used in pressure update).\n");
    args.support("alpha-d", "<alpha_d>\nSets the scale parameter alpha-d (used in pressure update).\nDefault: alpha-d <= theta * nu * delta_t, Stationary: alpha <= nu\n");
    args.support("deformation", "\nUse deformation tensor instead of gradient tensor.\n");
    args.support("flowbench","\nEnables the computation of 'flow around a cylinder' post-processing\nquantities such as drag, lift, etc.\n");
    args.support("nu", "<nu>\nSets the viscosity parameter for the velocity.\n");
    args.support("fractional", "\nUse the fractional-step-theta-scheme instead of the theta-scheme.\nDefault = false\n");
    args.support("theta", "<theta>\nSets the parameter theta.\n1 for Backward Euler, 0.5 for Crank-Nicolson...\nDefault = 1\n");
    args.support("res-fix", "<eps>\nSets the residual for the nonlinear / fixpoint iteration.\n");
    args.support("time-max", "<T_max>\nSets the maximum simulation time T_max.\n");
    args.support("time-steps", "<N>\nSets the number of time-steps for the time interval.\n");
    args.support("max-time-steps", "<N>\nSets the maximum number of time-steps to perform.\n");
    args.support("part-in", "<names...>\nSpecifies the names of the inflow mesh-parts.\n");
    args.support("part-out", "<names...>\nSpecifies the names of the outflow mesh-parts.\n");
    args.support("part-no", "<names...>\nSpecifies the names of the noflow mesh-parts.\n");
    args.support("vmax", "<v_max>\nSpecifies the maximum inflow velocity.\n");
    args.support("level", "<max> [[<med>] <min>]\nSets the maximum, medium and minimum mesh refinement levels.\n");
    args.support("vtk", "<name> [<step>]\nSets the name for VTK output and the time-stepping for the output (optional).\n");
    args.support("save", "<name> [<step>]\nSets the name for savefile and the time-stepping for the output (optional).\n");
    args.support("load", "<name> \nLoad a savefile.\n");
    args.support("mesh-file", "<name>\nSpecifies the filename of the input mesh file.\n");
    args.support("mesh-path", "<path>\nSpecifies the path of the directory containing the mesh file.\n");
    args.support("no-nonlinear", "\n Sets the fixpoint iterations on / off.\n");
    args.support("fix-steps", "<N>\nSets the number of non-linear / fixpoint iterations per time-step.\n");
    args.support("tol-fix", "<eps>\nSets the tolerance of the non-linear / fixpoint iteration.\n");
    args.support("no-multigrid-a", "\nUse no Multigrid-Solver as A-Solver (velocity).\n");
    args.support("no-multigrid-s", "\nUse PCG-Jacobi instead of Multigrid as S-Solver (pressure).\n");
    args.support("no-umfpack-s", "\nUse Richardson-Jacobi instead of UMFPACK as S-Coarse-Grid-Solver.\n");
    args.support("max-iter-a", "<N>\nSets the maximal number of iteration for the A-Solver (velocity).\nDefault: 1E-5\n");
    args.support("tol-rel-a", "<eps>\nSets the relative tolerance for the A-Solver (velocity).\nDefault: 1E-5\n");
    args.support("smooth-a", "<N>\nSets the number of smoothing steps for the A-Solver (velocity).\nDefault: 4\n");
    args.support("damp-a", "<omega>\nSets the smoother damping parameter for the A-Solver (velocity).\nDefault: 0.5\n");
    args.support("max-iter-s", "<N>\nSets the maximal number of iteration for the S-Solver (pressure).\nDefault: 1E-5\n");
    args.support("tol-abs-s", "<eps>\nSets the absolute tolerance for the S-Solver (pressure).\nDefault: 1E-5\n");
    args.support("smooth-s", "<N>\nSets the number of smoothing steps for the S-Solver (pressure).\nDefault: 4\n");
    args.support("damp-s", "<omega>\nSets the smoother damping parameter for the S-Solver (pressure).\nDefault: 0.5\n");
    args.support("test-mode", "\nRuns the application in regression test mode.\n");
    args.support("chorin", "\nUse Chorn mode.\nDefault = false\n");
    args.support("residual", "<eps>\nSet the global residual (stopping criteria).\nDefault: 1E-32 (off)\n");
    args.support("residual_drag", "<eps>\nSet the drag-residual (stopping criteria).\nDefault: 1E-32 (off)\n");
    args.support("residual_drag2", "<eps>\nSet the drag-time-residual (stopping criteria).\nDefault: 1E-32 (off)\n");
    args.support("residual_lift", "<eps>\nSet the lift-residual (stopping criteria).\nDefault: 1E-32 (off)\n");
    args.support("residual_lift2", "<eps>\nSet the lift-time-residual (stopping criteria).\nDefault: 1E-32 (off)\n");

    // no arguments given?
    if((argc <= 1) || (args.check("help") >= 0))
    {
      comm.print("\n2D/3D Navier-Stokes PP-Q2/P1 Solver\n");
      comm.print("The easiest way to make this application do something useful is");
      comm.print("to load a pre-defined problem configuration by supplying the");
      comm.print("option '--setup <config>', where <config> may be one of:\n");
      comm.print("  square      Poiseuille-Flow on Unit-Square");
      comm.print("  nozzle      Jet-Flow through Nozzle domain");
      comm.print("  fb-c2d-00   Nonsteady Flow Around A Cylinder (bench1 mesh)");
      comm.print("  fb-c2d-01   Nonsteady Flow Around A Cylinder (32 quad mesh)");
      comm.print("  fb-c2d-02   Nonsteady Flow Around A Cylinder (48 quad mesh)");
      comm.print("  fb-c2d-03   Nonsteady Flow Around A Cylinder (64 quad mesh)");
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

    // create our domain control
    Config cfg;

    args.parse("load", cfg.load);
    if (!cfg.load.empty())
    {
      int processes(0);
      load_conf(cfg, comm, processes);
      if (comm.size() != processes)
      {
        comm.print(std::cerr, "ERROR: wrong number of processes");

        comm.print(std::cerr, "Required are: " + stringify(processes));
        // abort
        FEAT::Runtime::abort();
      }
    }

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // parse our configuration
    if(!cfg.parse_args(args))
      FEAT::Runtime::abort();

    //
    // create a MeshFileReader to know the dimension of the mesh
    //
    Geometry::MeshFileReader mesh_reader;
    try
    {
      mesh_reader.add_mesh_files(comm, cfg.mesh_files, cfg.mesh_path);
    }
    catch(const std::exception& exc)
    {
      // Something went wrong; probably one of the files could not be opened...
      comm.print(std::cerr, "ERROR: " + stringify(exc.what()));
      Runtime::abort();
    }
    try
    {
      mesh_reader.read_root_markup();
    }
    catch(const std::exception& exc)
    {
      comm.print(std::cerr, "ERROR: " + stringify(exc.what()));
      Runtime::abort();
    }
    const String mesh_type = mesh_reader.get_meshtype_string();
    if(mesh_type.empty())
    {
      comm.print(std::cerr, "ERROR: Mesh file(s) did not provide a mesh-type!\n");
      comm.print(std::cerr, "Did you supply all required mesh-files?\n");
      Runtime::abort();
    }

    if(mesh_type == "conformal:hypercube:2:2") // 2D quadrilateral mesh
      domain_and_run<Shape::Hypercube<2>>(args, cfg, comm);
    else if(mesh_type == "conformal:hypercube:3:3") // 3D hexahedron mesh
      domain_and_run<Shape::Hypercube<3>>(args, cfg, comm);
    else
    {
      // The mesh-type is either invalid or not supported. In fact, there are some
      // other valid mesh-type specifiers (e.g. sub-dimensional meshes), which this
      // (and most other) application does not support by construction.
      comm.print(std::cerr, "ERROR: unsupported mesh type!\n");
      Runtime::abort();
    }
  }
} // namespace NavierStokesPP


int main(int argc, char* argv[])
{
  FEAT::Runtime::initialise(argc, argv);
  try
  {
    NavierStokesPP::main(argc, argv);
  }
  catch(const std::exception& exc)
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
