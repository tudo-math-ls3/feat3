// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// HPC Steady CCnD solver for DFG95 Flow-Around-A-Circle/Square/Cylinder/Cuboid/Sphere Benchmarks
// ------------------------------------------------------------------------------------------------
// This application implements a parallel steady CCND solver, which is pre-configured to solve
// the infamous unsteady "flow-around-an-obstacle" benchmark problem, which is defined in
//
//     M. Schaefer and S. Turek: Benchmark Computations of Laminar Flow Around a Cylinder
//
// The system is discretized using an isoparametric Q2/P1dc finite element discretization.
// The monolithic nonlinear Oseen systems are solved using an adaptive Newton-Multigrid solver
// with an additive matrix-based Vanka smoother ("AmaVanka") and using UMFPACK (if available)
// as a coarse grid solver. This application supports recursive partitioning.
//
// This application is aimed to be the high-performance implementation of the CCND, which makes use
// of all of FEAT's advanced features like multi-precision solver configuration, hybrid parallelism
// by combining MPI, CUDA and OpenMP. The downside is that this application is fairly 'hard-wired'
// to the one benchmark problem.
//
// ------------------------------------
// Basic Setup and Mandatory Parameters
// ------------------------------------
// This application defines default values for most of its parameters, however, two parameters
// are mandatory and always have to be specified explicitly:
//
// --mesh <meshfiles...>
// Specifies the input mesh file(s).
//
// --level <level-max> [levels...] <level-min>
// Specifies the mesh refinement levels in the syntax according to Control::PartiDomainControl.
//
//
// --------------------------------
// Backend and Threading Parameters
// --------------------------------
// This section describes the parameters that can be used to configure the backends and threading
// strategies:
//
// --backend [generic|cuda|mkl]
// Specifies the solver backend to use.
//
// --threads <n>
// Specifies the number of threads to use for multi-threading. Defaults to 1.
//
//
// ----------------------------
// System Definition Parameters
// ----------------------------
// Some of the basic parameters of the underlying system can be configured by the following
// parameters.
//
// --nu <nu>
// Specifies the viscosity parameter nu for the diffusion operator. Defaults to 1E-3.
//
// --deform
// If specified, the deformation tensor formulation is used for the diffusion operator,
// otherwise the gradient tensor formulation is used. Defaults to gradient tensor.
//
// --v-max <vmax>
// Specifies the maximum velocity of the inflow boundary condition.
// Defaults to 0.3 in 2D and 0.45 in 3D.
//
// --upsam <ups>
// Specifies the stabilization parameter <ups> for the streamline diffusion stabilization.
// Defaults to 0, i.e. unstabilized.
//
// --stokes
// If specified, only the steady-state Stokes equations are solved and the nonlinear iteration
// for solving the Navier-Stokes equations is skipped entirely.
//
//
// -------------------------------
// Solver Configuration Parameters
// -------------------------------
// This section describes the parameters that control the non-linear Newton/Picard solver
// as well as its multigrid preconditioner and its smoother component.
//
// --alpine
// If specified, the nonlinear system in each time step will be solved using an alternating
// Picard-Newton Picard iteration instead of the Newton iteration.
//
// --picard
// If specified, the nonlinear system in each time step will be solved using a simple
// Picard iteration instead of the Newton iteration. See \cite
//
// --nl-tol <tol>
// Specifies the normalized tolerance for the nonlinear solver. Defaults to 1E-8.
//
// --nl-max-iter <N>
// Specifies the maximum number of nonlinear (Newton/Picard/Alpine) solver iterations.
// Defaults to 25.
//
// --nl-stag-rate <rho>
// Specifies the stagnation rate criterion of nonlinear (Newton/Picard/Alpine) solver.
// Defaults to 1.05.
//
// --nl-slow-rate <rho>
// Specifies the slowdown rate criterion of nonlinear (Newton/Picard/Alpine) solver.
// Defaults to 0.1.
//
// --nl-damp <omega>
// Specifies the damping parameter for the nonlinear (Newton/Picard/Alpine) solver.
// Defaults to 1.0.
//
// --mg-min-iter <N>
// Specifies the minimum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 1.
//
// --mg-max-iter <N>
// Specifies the maximum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 50.
//
// --mg-max-boot-iter <N>
// Specifies the maximum number of multigrid iterations in the boot-strapping phase of the nonlinear solver.
// Defaults to 3.
//
// --mg-stag-rate <rho>
// Specifies the stagnation rate criterion of multigrid solver.
// Defaults to 1.05.
//
// --mg-damp <omega>
// Specifies the damping parameter for the multigrid solver.
// Defaults to 1.0.
//
// --mg-adapt <fixed|defect|energy>
// Specifies which coarse grid adaption method to use for the multigrid solver.
// Defaults to 'fixed'.
//
// --mg-cycle <V|F|W>
// Specifies which multigrid cycle to use to use for the multigrid solver.
// Defaults to 'V'.
//
// --mg-smooth <pre|post|both>
// Specifies which smoothing strategy to use for the multigrid solver.
// Defaults to 'both'.
//
// --smooth-plot
// If given, specifies that the smoother iterations are to be plotted to the console.
//
// --smooth-steps <N>
// Specifies the number of pre- and post-smoothing steps. Defaults to 8.
//
// --smooth-damp <omega>
// Specifies the damping parameter for the smoother. Defaults to 0.7.
//
// --smooth-gmres
// If given, specifies that the AmaVanka is to be used in a GMRES iteration rather than a Richardson iteration.
//
// --mg-tol-rel <tol>
// If given, specifies the relative tolerance for the multigrid solver.
// If not given, then the tolerance for the multigrid solver is chosen in an adaptive
// manner depending on the two previous nonlinear solver defects, which is the default case.
// The adaptive tolerance is chosen in each nonlinear iteration by analyzing the nonlinear
// defect improvement in the previous nonlinear (Newton/Picard) solver iteration in the
// following manner: Let def_{j} and def_{j-1} denote the two previous nonlinear defect norms,
// then the next nonlinear defect norm def_{j+1} should approximately fulfill the equation
//
//        (def_{j+1} / def_{j}) \approx (def_{j} / def_{j+1})^C
//
// where C \in {1,2} is the convergence speed of the nonlinear solver, i.e. C=2 for Newton
// and C=1 for Picard iteration. Multiplying the above equation by def_{j} gives us an
// estimate for the next nonlinear defect norm def_{j+1}:
//
//        def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})^C
//
// To obtain an absolute tolerance for the next multigrid solver application, we simply
// multiply the above estimate by 0.1.
//
// --flush-mg
// If given, each output line of the multigrid solver is flushed to cout.
//
//
// -----------------------------------------------------
// Initial Solution Read-In and Final Solution Write-Out
// -----------------------------------------------------
//
// --save-sol <filename>
// Specifies that the application should write the final (partitioned) solution (and the
// partitioning) to a single binary output file. The output file can be loaded by the --load-sol
// option if the input mesh, the refinement level as well as the partitioning (and thus the
// process count) are identical.
//
// --load-sol <filename>
// Specifies that the application should read in the initial (partitioned) solution guess
// from a single binary output file, which was written by a --save-sol from a previous run.
// If specified, the solving of the Stokes system to obtain an initial guess is skipped.
//
//
// -----------------------------------------------------
// Post-Processing Parameters
// -----------------------------------------------------
// This section describes parameters that configure the post-processing phase.
//
// --vtk <filename> [<refined-filename>]
// Specifies that the application should write a VTK visualization output file. The second
// optional parameter specifies the filename for a VTK output file on a once refined mesh,
// if given. Note: it is possible to output only the refined VTKs by passing a whitespace
// string as the first filename, e.g.: --vtk " " myvtk
//
// --skip-analysis
// If given, specifies to skip the entire post-processing solution analysis phase.
//
// --p-a <x> <y> [<z>]
// Specifies the coordinates of the pressure difference point p_a.
// Defaults to (0.15, 0.2) in 2D and (0.45, 0.2, 0.205) in 3D.
//
// --p-e <x> <y> [<z>]
// Specifies the coordinates of the pressure difference point p_e.
// Defaults to (0.25, 0.2) in 2D and (0.55, 0.2, 0.205) in 3D.
//
// --ext-stats
// If given, specifies that the application should output extensive statistics at the end of the
// program run, including detailed MPI timings.
//
//
// \author Peter Zajac
//

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/burgers_assembly_job.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/vanka.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/umfpack.hpp>
#include <kernel/solver/direct_stokes_solver.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/checkpoint_control.hpp>
#include <control/statistics.hpp>

#ifdef FEAT_HAVE_OMP
#include <omp.h>
#endif


namespace DFG95
{
  using namespace FEAT;

#ifndef FEAT_DFG95_BENCH1_DIM
#error Macro 'FEAT_DFG95_BENCH1_DIM' must be defined to 2 or 3 by the build system!
#endif

  // our dimension
  static constexpr int dim = FEAT_DFG95_BENCH1_DIM;

  // our one and only index type; use 32 bit for a small performance boost
  typedef std::uint32_t IndexType;

  // select floating point precisions for the system and the solver
  // the solver precision may be lower than the system precision
#ifdef FEAT_DFG95_BENCH1_QUAD
#  define Q_(x) (x##Q)
  typedef __float128 SystemDataType;
  typedef __float128 SolverDataType;
  static constexpr int fp_num_digs = 35;
  static const char* fp_typename = "system: quadruple, solver: quadruple";
#elif defined(FEAT_DFG95_BENCH1_QUAD_DOUBLE)
#  define Q_(x) (x##Q)
  typedef __float128 SystemDataType;
  typedef double SolverDataType;
  static constexpr int fp_num_digs = 35;
  static const char* fp_typename = "system: quadruple, solver: double";
#elif defined(FEAT_DFG95_BENCH1_DOUBLE)
#  define Q_(x) x
  typedef double SystemDataType;
  typedef double SolverDataType;
  static constexpr int fp_num_digs = 17;
  static const char* fp_typename = "system: double, solver: double";
#elif defined(FEAT_DFG95_BENCH1_DOUBLE_SINGLE)
#  define Q_(x) (x##Q)
  typedef double SystemDataType;
  typedef float  SolverDataType;
  static constexpr int fp_num_digs = 17;
  static const char* fp_typename = "system: double, solver: single";
#else
#error Exactly one of the following macros must be defined to specify the floating point formats to use:
#error FEAT_DFG95_BENCH1_DOUBLE: use double for system and solver
#error FEAT_DFG95_BENCH1_DOUBLE_SINGLE: use double for system and single for solver
#error FEAT_DFG95_BENCH1_QUAD: use quadruple for system and solver
#error FEAT_DFG95_BENCH1_QUAD_DOUBLE: use quadruple for system and double for solver
#endif

  struct Counts
  {
    static constexpr std::size_t velo_dofs = 0u;        // number of velocity dofs nodes
    static constexpr std::size_t pres_dofs = 1u;        // number of pressure dofs
    static constexpr std::size_t total_dofs = 2u;       // total dofs = dim*velo_dofs + pres_dofs
    static constexpr std::size_t nonlin_iter = 3u;      // number of nonlinear iterations (newton/picard)
    static constexpr std::size_t linsol_iter = 4u;      // number of linear solver iterations
    static constexpr std::size_t nnze_a = 5u;           // number of non-zero blocks in matrix-block A
    static constexpr std::size_t nnze_b = 6u;           // number of non-zero blocks in matrix-block B
    static constexpr std::size_t nnze_total = 7u;       // total number of non-zero entries in matrix
    static constexpr std::size_t vanka_data = 8u;       // total number of data entries for Vanka
    static constexpr std::size_t fine_elements = 9u;    // total number of fine mesh elements
    static constexpr std::size_t all_elements = 10u;    // total number of all mesh elements
    static constexpr std::size_t count = 11u;
  };

  struct Times
  {
    static constexpr std::size_t total_run = 0u;        // total runtime
    static constexpr std::size_t nonlin_total = 1u;     // nonlinear solver time
    static constexpr std::size_t nonlin_asm_def = 2u;   // nonlinear defect assembly time
    static constexpr std::size_t nonlin_asm_mat = 3u;   // nonlinear matrix assembly time
    static constexpr std::size_t linsol_init = 4u;      // linear solver init time
    static constexpr std::size_t linsol_apply = 5u;     // linear solver apply time
    static constexpr std::size_t mg_defect = 6u;        // multigrid defect compute time
    static constexpr std::size_t mg_smooth = 7u;        // multigrid smooth apply time
    static constexpr std::size_t mg_coarse = 8u;        // multigrid coarse grid solve time
    static constexpr std::size_t mg_transfer = 9u;      // multigrid grid transfer time
    static constexpr std::size_t vanka_init_sym = 10u;  // vanka symbolic factorization time
    static constexpr std::size_t vanka_init_num = 11u;  // vanka numeric factorization time
    static constexpr std::size_t vanka_apply = 12u;     // vanka apply time
    static constexpr std::size_t vtk_write = 13u;       // VTK output time
    static constexpr std::size_t checkpoint = 14u;      // checkpoint write time
    static constexpr std::size_t sol_analysis = 15u;    // solution analysis time
    static constexpr std::size_t count = 16u;
  };

  struct Bytes
  {
    static constexpr std::size_t peak_p = 0u;           // peak physical memory usage
    static constexpr std::size_t peak_v = 1u;           // peak virtual memory usage
    static constexpr std::size_t mesh = 2u;             // mesh node memory usage
    static constexpr std::size_t gate = 3u;             // gate memory usage
    static constexpr std::size_t muxer = 4u;            // muxer memory usage
    static constexpr std::size_t matrix = 5u;           // system matrices memory usage
    static constexpr std::size_t matrix_struct = 6u;    // matrix structures memory usage
    static constexpr std::size_t matrix_values = 7u;    // matrix values memory usage
    static constexpr std::size_t transfer = 8u;         // transfers operators memory usage
    static constexpr std::size_t vanka = 9u;            // vanka memory usage
    static constexpr std::size_t count = 10u;
  };

  static constexpr std::size_t padlen = std::size_t(30);

  // class for collecting benchmark post-processing summary data
  class BenchmarkSummary
  {
  public:
    // raw body forces from surface/volume integration
    Tiny::Vector<SystemDataType, 3> body_forces_raw_surf, body_forces_raw_vol;

    // drag/lift forces by surface integration
    SystemDataType drag_coeff_surf_cylinder, lift_coeff_surf_cylinder, drag_coeff_surf_sphere, lift_coeff_surf_sphere;

    // drag/lift forces by volume integration
    SystemDataType drag_coeff_vol_cylinder, lift_coeff_vol_cylinder, drag_coeff_vol_sphere, lift_coeff_vol_sphere;

    // pressure difference
    SystemDataType pres_diff;

    // pressure drop
    SystemDataType pres_drop_raw, pres_drop_cuboid, pres_drop_pipe;

    // flow through upper/lower region
    SystemDataType flux_upper, flux_lower;

    // velocity field information
    Assembly::VelocityInfo<SystemDataType, dim> velo_info;

    BenchmarkSummary() :
      drag_coeff_surf_cylinder(0.0), lift_coeff_surf_cylinder(0.0),
      drag_coeff_surf_sphere(0.0), lift_coeff_surf_sphere(0.0),
      drag_coeff_vol_cylinder(0.0), lift_coeff_vol_cylinder(0.0),
      drag_coeff_vol_sphere(0.0), lift_coeff_vol_sphere(0.0),
      pres_diff(0.0), pres_drop_raw(0.0), pres_drop_cuboid(0.0), pres_drop_pipe(0.0),
      flux_upper(0.0), flux_lower(0.0)
    {
    }

    String format(int prec = fp_num_digs+5) const
    {
      String s;
      const char pc = '.';
      // append coefficients and velocity info
      s += "Solution Analysis:\n";
      s += String("Drag Body Force [Raw Surface]").pad_back(padlen, pc) + ": " + stringify_fp_sci(body_forces_raw_surf[0], prec) + "\n";
      s += String("Lift Body Force [Raw Surface]").pad_back(padlen, pc) + ": " + stringify_fp_sci(body_forces_raw_surf[1], prec) + "\n";
      s += String("Side Body Force [Raw Surface]").pad_back(padlen, pc) + ": " + stringify_fp_sci(body_forces_raw_surf[2], prec) + "\n";
      s += String("Drag Body Force [Raw Volume]").pad_back(padlen, pc) + ": " + stringify_fp_sci(body_forces_raw_vol[0], prec) + "\n";
      s += String("Lift Body Force [Raw Volume]").pad_back(padlen, pc) + ": " + stringify_fp_sci(body_forces_raw_vol[1], prec) + "\n";
      s += String("Side Body Force [Raw Volume]").pad_back(padlen, pc) + ": " + stringify_fp_sci(body_forces_raw_vol[2], prec) + "\n";
      s += String("Drag Coeff [Cylinder Surface]").pad_back(padlen, pc) + ": " + stringify_fp_fix(drag_coeff_surf_cylinder, prec) + "\n";
      s += String("Lift Coeff [Cylinder Surface]").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff_surf_cylinder, prec) + "\n";
      s += String("Drag Coeff [Cylinder Volume]").pad_back(padlen, pc) + ": " + stringify_fp_fix(drag_coeff_vol_cylinder, prec) + "\n";
      s += String("Lift Coeff [Cylinder Volume]").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff_vol_cylinder, prec) + "\n";
      s += String("Drag Coeff [Sphere Surface]").pad_back(padlen, pc) + ": " + stringify_fp_fix(drag_coeff_surf_sphere, prec) + "\n";
      s += String("Lift Coeff [Sphere Surface]").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff_surf_sphere, prec) + "\n";
      s += String("Drag Coeff [Sphere Volume]").pad_back(padlen, pc) + ": " + stringify_fp_fix(drag_coeff_vol_sphere, prec) + "\n";
      s += String("Lift Coeff [Sphere Volume]").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff_vol_sphere, prec) + "\n";
      s += String("Pressure Difference").pad_back(padlen, pc) + ": " + stringify_fp_fix(pres_diff, prec) + "\n";
      s += String("Pressure Drop [Raw]").pad_back(padlen, pc) + ": " + stringify_fp_sci(pres_drop_raw, prec) + "\n";
      s += String("Pressure Drop [Cuboid]").pad_back(padlen, pc) + ": " + stringify_fp_fix(pres_drop_cuboid, prec) + "\n";
      s += String("Pressure Drop [Pipe]").pad_back(padlen, pc) + ": " + stringify_fp_fix(pres_drop_pipe, prec) + "\n";
      s += String("Upper Flux").pad_back(padlen, pc) + ": " + stringify_fp_fix(flux_upper, prec) + "\n";
      s += String("Lower Flux").pad_back(padlen, pc) + ": " + stringify_fp_fix(flux_lower, prec) + "\n";
      s += String("Kinetic Energy").pad_back(padlen, pc) + ": " + stringify_fp_sci(SystemDataType(0.5)*Math::sqr(velo_info.norm_h0), prec) + "\n";
      s += velo_info.format_string(prec, padlen, pc);

      return s;
    }
  }; // class BenchmarkSummary

  // class for collection benchmark statistics
  class BenchmarkStats
  {
  public:
    // timings (in seconds)
    double times[Times::count];     // of this process
    double times_sum[Times::count]; // sum over all processes
    double times_min[Times::count]; // minimum over all processes
    double times_max[Times::count]; // maximum over all processes

    // counts
    unsigned long long counts[Counts::count];     // of this process
    unsigned long long counts_sum[Counts::count]; // summed up over all processes
    unsigned long long counts_min[Counts::count]; // minimum over all processes
    unsigned long long counts_max[Counts::count]; // maximum over all processes

    // sizes (in bytes)
    unsigned long long bytes[Bytes::count];     // of this process
    unsigned long long bytes_sum[Bytes::count]; // summed up over all processes
    unsigned long long bytes_min[Bytes::count]; // minimum over all processes
    unsigned long long bytes_max[Bytes::count]; // maximum over all processes

    // multigrid timings for each level
    std::vector<double> mg_times_defect, mg_times_smooth, mg_times_transfer, mg_times_coarse;
    std::vector<double> mg_times_defect_sum, mg_times_smooth_sum, mg_times_transfer_sum, mg_times_coarse_sum;
    std::vector<double> mg_times_defect_min, mg_times_smooth_min, mg_times_transfer_min, mg_times_coarse_min;
    std::vector<double> mg_times_defect_max, mg_times_smooth_max, mg_times_transfer_max, mg_times_coarse_max;

    // --------------------------------------------------------------

    explicit BenchmarkStats(std::size_t sv) :
      mg_times_defect(sv+1u, 0.0), mg_times_smooth(sv+1u, 0.0), mg_times_transfer(sv+1u, 0.0), mg_times_coarse(sv+1u, 0.0),
      mg_times_defect_sum(sv+1u, 0.0), mg_times_smooth_sum(sv+1u, 0.0), mg_times_transfer_sum(sv+1u, 0.0), mg_times_coarse_sum(sv+1u, 0.0),
      mg_times_defect_min(sv+1u, 0.0), mg_times_smooth_min(sv+1u, 0.0), mg_times_transfer_min(sv+1u, 0.0), mg_times_coarse_min(sv+1u, 0.0),
      mg_times_defect_max(sv+1u, 0.0), mg_times_smooth_max(sv+1u, 0.0), mg_times_transfer_max(sv+1u, 0.0), mg_times_coarse_max(sv+1u, 0.0)
    {
      for(std::size_t i(0); i < Counts::count; ++i)
        counts[i] = counts_sum[i] = counts_min[i] = counts_max[i] = std::size_t(0);
      for(std::size_t i(0); i < Times::count; ++i)
        times[i] = times_sum[i] = times_min[i] = times_max[i] = std::size_t(0);
      for(std::size_t i(0); i < Bytes::count; ++i)
        bytes[i] = bytes_sum[i] = bytes_min[i] = bytes_max[i] = std::size_t(0);
    }

    String format() const
    {
      String s;
      const char pc = '.';

      // solver statistics
      s += "\nSolver Statistics:\n";
      s += String("Multigrid Iterations").pad_back(padlen, pc) + ": " + stringify(counts[Counts::linsol_iter]) + "\n";
      s += String("Nonlinear Iterations").pad_back(padlen, pc) + ": " + stringify(counts[Counts::nonlin_iter]) + "\n";

      s += String("\nDiscretization Statistics:").pad_back(padlen+15) + String("Sum").pad_back(29)
        + String("Min").pad_back(29) + String("Max").pad_back(22) + "Balance\n";
      s += format_count("Fine Mesh Elements", Counts::fine_elements);
      s += format_count("All Mesh Elements", Counts::all_elements);
      s += format_count("Velocity Dof Nodes", Counts::velo_dofs);
      s += format_count("Pressure Dofs", Counts::pres_dofs);
      s += format_count("Total Dofs", Counts::total_dofs);
      s += format_count("Nonzero Blocks A", Counts::nnze_a);
      s += format_count("Nonzero Blocks B", Counts::nnze_b);
      s += format_count("Total Nonzero Entries", Counts::nnze_total);
      s += format_count("Vanka Nonzero Entries", Counts::vanka_data);

      // append timing info
      s += String("\nTiming Statistics:").pad_back(padlen+15) + String("Sum").pad_back(29)
        + String("Min").pad_back(29) + String("Max").pad_back(22) + "Balance\n";
      s += format_runtime("Total Runtime", Times::total_run);
      s += format_runtime("Nonlinear Solver Total", Times::nonlin_total);
      s += format_runtime("Nonlinear Defect Assembly", Times::nonlin_asm_def);
      s += format_runtime("Nonlinear Matrix Assembly", Times::nonlin_asm_mat);
      s += format_runtime("Multigrid Solver Initialize", Times::linsol_init);
      s += format_runtime("Multigrid Solver Apply", Times::linsol_apply);
      s += format_runtime("Multigrid Defect Compute", Times::mg_defect);
      s += format_runtime("Multigrid Smooth Apply", Times::mg_smooth);
      s += format_runtime("Multigrid Coarse Solve", Times::mg_coarse);
      s += format_runtime("Multigrid Transfer Apply", Times::mg_transfer);
      s += format_runtime("Vanka Symbolic Initialize", Times::vanka_init_sym);
      s += format_runtime("Vanka Numeric Initialize", Times::vanka_init_num);
      s += format_runtime("Vanka Local Apply", Times::vanka_apply);
      s += format_runtime("Checkpointing", Times::checkpoint);
      s += format_runtime("Solution Analysis", Times::sol_analysis);
      s += format_runtime("VTK Write", Times::vtk_write);

      // append memory info
      s += String("\nMemory Statistics:").pad_back(padlen+15) + String("Sum").pad_back(29)
        + String("Min").pad_back(29) + String("Max").pad_back(22) + "Balance\n";;
      s += format_memuse("Peak Physical Memory", Bytes::peak_p);
      s += format_memuse("Peak Virtual Memory", Bytes::peak_v);
      s += format_memuse("Mesh-Node Size", Bytes::mesh);
      s += format_memuse("Gate Size", Bytes::gate);
      s += format_memuse("Muxer Size", Bytes::muxer);
      s += format_memuse("Matrix Total Size", Bytes::matrix);
      s += format_memuse("Matrix Struct Size", Bytes::matrix_struct);
      s += format_memuse("Matrix Values Size", Bytes::matrix_values);
      s += format_memuse("Transfer Size", Bytes::transfer);
      s += format_memuse("Vanka Size", Bytes::vanka);
      s += format_memuse("Vanka Size", Bytes::vanka);

      // append multigrid timings
      s += String("\nMultigrid Timings #1:").pad_back(40)
        + String("Sum").pad_back(12) + String("Min").pad_back(12) + String("Max").pad_back(6) + String("Balance").pad_back(22)
        + String("Sum").pad_back(12) + String("Min").pad_back(12) + String("Max").pad_back(6) + String("Balance") + "\n";
        for(std::size_t i(0); i < mg_times_defect.size(); ++i)
      {
        if(i == 0u)
          s += "Defect/Smooth   Level Sum:";
        else
          s += "Defect/Smooth   Level" + stringify(i-1u).pad_front(4) + ":";
        s += stringify_fp_fix(mg_times_defect_sum[i], 3, 16) + " /";
        s += stringify_fp_fix(mg_times_defect_min[i], 3, 10) + " /";
        s += stringify_fp_fix(mg_times_defect_max[i], 3, 10) + " /";
        s += stringify_fp_fix(mg_times_defect_max[i] > 0.0 ? mg_times_defect_min[i] / mg_times_defect_max[i] : 0.0, 5, 8) + " |";

        s += stringify_fp_fix(mg_times_smooth_sum[i], 3, 16) + " /";
        s += stringify_fp_fix(mg_times_smooth_min[i], 3, 10) + " /";
        s += stringify_fp_fix(mg_times_smooth_max[i], 3, 10) + " /";
        s += stringify_fp_fix(mg_times_smooth_max[i] > 0.0 ? mg_times_smooth_min[i] / mg_times_smooth_max[i] : 0.0, 5, 8) + "\n";
      }

      s += String("\nMultigrid Timings #2:").pad_back(40)
        + String("Sum").pad_back(12) + String("Min").pad_back(12) + String("Max").pad_back(6) + String("Balance").pad_back(22)
        + String("Sum").pad_back(12) + String("Min").pad_back(12) + String("Max").pad_back(6) + String("Balance") + "\n";
      for(std::size_t i(0); i < mg_times_defect.size(); ++i)
      {
        if(i == 0u)
          s += "Transfer/Coarse Level Sum:";
        else
          s += "Transfer/Coarse Level" + stringify(i-1u).pad_front(4) + ":";
        s += stringify_fp_fix(mg_times_transfer_sum[i], 3, 16) + " /";
        s += stringify_fp_fix(mg_times_transfer_min[i], 3, 10) + " /";
        s += stringify_fp_fix(mg_times_transfer_max[i], 3, 10) + " /";
        s += stringify_fp_fix(mg_times_transfer_max[i] > 0.0 ? mg_times_transfer_min[i] / mg_times_transfer_max[i] : 0.0, 5, 8) + " |";

        s += stringify_fp_fix(mg_times_coarse_sum[i], 3, 16) + " /";
        s += stringify_fp_fix(mg_times_coarse_min[i], 3, 10) + " /";
        s += stringify_fp_fix(mg_times_coarse_max[i], 3, 10) + " /";
        s += stringify_fp_fix(mg_times_coarse_max[i] > 0.0 ? mg_times_coarse_min[i] / mg_times_coarse_max[i] : 0.0, 5, 8) + "\n";
      }

      return s;
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(times, times_sum, Times::count, Dist::op_sum);
      comm.allreduce(times, times_min, Times::count, Dist::op_min);
      comm.allreduce(times, times_max, Times::count, Dist::op_max);
      comm.allreduce(counts, counts_sum, Counts::count, Dist::op_sum);
      comm.allreduce(counts, counts_min, Counts::count, Dist::op_min);
      comm.allreduce(counts, counts_max, Counts::count, Dist::op_max);
      comm.allreduce(bytes, bytes_sum, Bytes::count, Dist::op_sum);
      comm.allreduce(bytes, bytes_min, Bytes::count, Dist::op_min);
      comm.allreduce(bytes, bytes_max, Bytes::count, Dist::op_max);
      comm.allreduce(mg_times_defect.data(), mg_times_defect_sum.data(), mg_times_defect.size(), Dist::op_sum);
      comm.allreduce(mg_times_defect.data(), mg_times_defect_min.data(), mg_times_defect.size(), Dist::op_min);
      comm.allreduce(mg_times_defect.data(), mg_times_defect_max.data(), mg_times_defect.size(), Dist::op_max);
      comm.allreduce(mg_times_smooth.data(), mg_times_smooth_sum.data(), mg_times_smooth.size(), Dist::op_sum);
      comm.allreduce(mg_times_smooth.data(), mg_times_smooth_min.data(), mg_times_smooth.size(), Dist::op_min);
      comm.allreduce(mg_times_smooth.data(), mg_times_smooth_max.data(), mg_times_smooth.size(), Dist::op_max);
      comm.allreduce(mg_times_transfer.data(), mg_times_transfer_sum.data(), mg_times_transfer.size(), Dist::op_sum);
      comm.allreduce(mg_times_transfer.data(), mg_times_transfer_min.data(), mg_times_transfer.size(), Dist::op_min);
      comm.allreduce(mg_times_transfer.data(), mg_times_transfer_max.data(), mg_times_transfer.size(), Dist::op_max);
      comm.allreduce(mg_times_coarse.data(), mg_times_coarse_sum.data(), mg_times_coarse.size(), Dist::op_sum);
      comm.allreduce(mg_times_coarse.data(), mg_times_coarse_min.data(), mg_times_coarse.size(), Dist::op_min);
      comm.allreduce(mg_times_coarse.data(), mg_times_coarse_max.data(), mg_times_coarse.size(), Dist::op_max);
    }

    String format_count(String s, std::size_t cc) const
    {
      return s.pad_back(padlen, '.') + ":"
        + stringify_fp_fix(counts_sum[cc], 3, 16).pad_back(35)
        + stringify_fp_fix(counts_min[cc], 3, 10).pad_back(29)
        + stringify_fp_fix(counts_max[cc], 3, 10).pad_back(29)
        + stringify_fp_fix(counts_max[cc] > 0u ? double(counts_min[cc]) / double(counts_max[cc]) : 0.0, 5, 7) + "\n";
    }

    String format_runtime(String s, std::size_t tt) const
    {
      return s.pad_back(padlen, '.') + ":"
        + stringify_fp_fix(times_sum[tt], 3, 16) + " s  [" + stringify_fp_fix(100.0*times_sum[tt]/times_sum[Times::total_run], 3, 8) + "% ] / "
        + stringify_fp_fix(times_min[tt], 3, 10) + " s  [" + stringify_fp_fix(100.0*times_min[tt]/times_max[Times::total_run], 3, 8) + "% ] / "
        + stringify_fp_fix(times_max[tt], 3, 10) + " s  [" + stringify_fp_fix(100.0*times_max[tt]/times_max[Times::total_run], 3, 8) + "% ] / "
        + stringify_fp_fix(times_max[tt] > 1e-9 ? times_min[tt] / times_max[tt] : 0.0, 5, 7) + "\n";
    }

    String format_memuse(String s, std::size_t mm) const
    {
      return s.pad_back(padlen, '.') + ":"
        + stringify_fp_fix(double(bytes_sum[mm]) / 1073741824.0, 3, 16) + " GB [" + stringify_fp_fix(100.0*double(bytes_sum[mm])/double(bytes_sum[Bytes::peak_p]), 3, 8) + "% ] / "
        + stringify_fp_fix(double(bytes_min[mm]) / 1073741824.0, 3, 10) + " GB [" + stringify_fp_fix(100.0*double(bytes_min[mm])/double(bytes_max[Bytes::peak_p]), 3, 8) + "% ] / "
        + stringify_fp_fix(double(bytes_max[mm]) / 1073741824.0, 3, 10) + " GB [" + stringify_fp_fix(100.0*double(bytes_max[mm])/double(bytes_max[Bytes::peak_p]), 3, 8) + "% ] / "
        + stringify_fp_fix(bytes_max[mm] > 0u ? double(bytes_min[mm]) / double(bytes_max[mm]) : 0.0, 5, 7) + "\n";
    }

  }; // class BenchmarkSummary

  // accumulator for benchmark body forces (i.e. drag and lift)
  // this class is used for the computation of the 'surface integration' variants of drag and lift
  class BenchBodyForceAccumulator
  {
  public:
    const SystemDataType _nu;

    // old: computation method from DFG95 paper, works only for 2.5D (cylinder), gives wrong results for sphere obstacle
    SystemDataType drag_raw_old, lift_raw_old;

    // new: computation method from Giles paper, always gives correct results, but may be slightly less accurate on cylinder
    SystemDataType drag_raw_new, lift_raw_new, side_raw_new;

    explicit BenchBodyForceAccumulator(SystemDataType nu) :
      _nu(nu),
      drag_raw_old(SystemDataType(0)), lift_raw_old(SystemDataType(0)),
      drag_raw_new(SystemDataType(0)), lift_raw_new(SystemDataType(0)), side_raw_new(SystemDataType(0))
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

      Tiny::Matrix<T_, 2, 2, 2, 2> nt;
      nt(0,0) = tx * nx;
      nt(0,1) = tx * ny;
      nt(1,0) = ty * nx;
      nt(1,1) = ty * ny;

      const T_ dut = Tiny::dot(nt, grad_v);

      drag_raw_old += SystemDataType(omega * ( _nu * dut * ny - val_p * nx));
      lift_raw_old += SystemDataType(omega * (-_nu * dut * nx - val_p * ny));

      const Tiny::Vector<T_, 2> eta = Tiny::orthogonal(jac).normalize();
      drag_raw_new -= SystemDataType(omega * (_nu*(T_(2) * grad_v(0,0) * eta[0] + (grad_v(0, 1) + grad_v(1, 0)) * eta[1]) - val_p * eta[0]));
      lift_raw_new -= SystemDataType(omega * (_nu*(T_(2) * grad_v(1,1) * eta[1] + (grad_v(1, 0) + grad_v(0, 1)) * eta[0]) - val_p * eta[1]));
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

      Tiny::Matrix<T_, 3, 3, 3, 3> nt;
      nt.format();
      nt(0,0) = tx * nx;
      nt(0,1) = tx * ny;
      nt(1,0) = ty * nx;
      nt(1,1) = ty * ny;

      const T_ dut = Tiny::dot(nt, grad_v);

      drag_raw_old += SystemDataType(omega * ( _nu * dut * ny - val_p * nx));
      lift_raw_old += SystemDataType(omega * (-_nu * dut * nx - val_p * ny));

      const Tiny::Vector<T_, 3> eta = Tiny::orthogonal(jac).normalize();
      drag_raw_new -= SystemDataType(omega * (_nu*(T_(2) * grad_v(0,0) * eta[0] + (grad_v(0, 1) + grad_v(1, 0)) * eta[1] + (grad_v(0, 2) + grad_v(2, 0)) * eta[2]) - val_p * eta[0]));
      lift_raw_new -= SystemDataType(omega * (_nu*(T_(2) * grad_v(1,1) * eta[1] + (grad_v(1, 2) + grad_v(2, 1)) * eta[2] + (grad_v(1, 0) + grad_v(0, 1)) * eta[0]) - val_p * eta[1]));
      side_raw_new -= SystemDataType(omega * (_nu*(T_(2) * grad_v(2,2) * eta[2] + (grad_v(2, 0) + grad_v(0, 2)) * eta[0] + (grad_v(2, 1) + grad_v(1, 2)) * eta[1]) - val_p * eta[2]));
    }

    void sync(const Dist::Comm& comm)
    {
      SystemDataType v[] =
      {
        drag_raw_old, lift_raw_old, drag_raw_new, lift_raw_new, side_raw_new
      };
      comm.allreduce(v, v, std::size_t(5), Dist::op_sum);
      drag_raw_old = v[0];
      lift_raw_old = v[1];
      drag_raw_new = v[2];
      lift_raw_new = v[3];
      side_raw_new = v[4];
    }
  }; // class BenchBodyForces<...,2>

  // accumulator for computation of X-flux
  class XFluxAccumulator
  {
  public:
    SystemDataType flux;

    XFluxAccumulator() :
      flux(SystemDataType(0))
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
      flux += SystemDataType(omega * val_v[0]);
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(&flux, &flux, std::size_t(1), Dist::op_sum);
    }
  }; // class XFluxAccumulator

  // computes the body forces by the volumetric 'defect vector' approach
  template<typename DT_, typename IT_, int dim_>
  void assemble_bdforces_vol(
    Tiny::Vector<DT_, 3>& forces,
    const LAFEM::DenseVectorBlocked<DT_, IT_, dim_>& vec_def_v,
    const LAFEM::DenseVector<DT_, IT_>& vec_char)
  {
    static_assert(dim_ <= 3, "invalid forces size");
    Tiny::Vector<DT_, dim_, 3> frc(DT_(0));

    XASSERT(vec_def_v.size() == vec_char.size());

    // get the vector arrays
    const Index n = vec_def_v.size();
    const auto* vdef = vec_def_v.elements();
    const auto* vchr = vec_char.elements();
    for(Index i(0); i < n; ++i)
    {
      frc.axpy(vchr[i], vdef[i]);
    }
    forces.template copy_n<dim_>(frc);
  }

  template<typename Vector_>
  void vec2buf(std::vector<char>& buffer, const Vector_& vector)
  {
    const auto& loc_v = vector.template at<0>();
    const auto& loc_p = vector.template at<1>();
    const std::size_t nv = loc_v.template size<LAFEM::Perspective::pod>() * sizeof(SystemDataType);
    const std::size_t np = loc_p.template size<LAFEM::Perspective::pod>() * sizeof(SystemDataType);
    std::uint64_t header[4u] =
    {
      std::uint64_t(dim),
      sizeof(SystemDataType),
      loc_v.template size<LAFEM::Perspective::native>(),
      loc_p.template size<LAFEM::Perspective::native>()
    };
    buffer.resize(nv + np + 32ull);
    char* buf_h = buffer.data();
    char* buf_v = &buf_h[32ull];
    char* buf_p = &buf_v[nv];
    memcpy(buf_h, header, 32ull);
    memcpy(buf_v, loc_v.elements(), nv);
    memcpy(buf_p, loc_p.elements(), np);
  }

  template<typename Vector_>
  void buf2vec(Vector_& vector, const std::vector<char>& buffer)
  {
    auto& loc_v = vector.template at<0>();
    auto& loc_p = vector.template at<1>();
    const std::size_t nv = loc_v.template size<LAFEM::Perspective::pod>() * sizeof(SystemDataType);
    const std::size_t np = loc_p.template size<LAFEM::Perspective::pod>() * sizeof(SystemDataType);

    XASSERTM(buffer.size() == (nv+np+32ull), "invalid vector buffer size!");
    const char* buf_h = buffer.data();
    const char* buf_v = &buf_h[32ull];
    const char* buf_p = &buf_v[nv];

    std::uint64_t header[4u] = {0u, 0u, 0u, 0u};
    memcpy(&header, buf_h, 32ull);
    XASSERTM(int(header[0]) == dim, "invalid vector buffer dimension!");
    XASSERTM(header[1] == sizeof(SystemDataType), "invalid vector buffer data type!");
    XASSERTM(header[2] == loc_v.template size<LAFEM::Perspective::native>(), "invalid velocity vector buffer size!");
    XASSERTM(header[3] == loc_p.template size<LAFEM::Perspective::native>(), "invalid pressure vector buffer size!");

    memcpy(loc_v.elements(), buf_v, nv);
    memcpy(loc_p.elements(), buf_p, np);
  }
} // namespace DFG95

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  using namespace FEAT;
  using namespace DFG95;

  // create runtime scope guard
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // create world communicator
  Dist::Comm comm(Dist::Comm::world());

  // print number of processes
  comm.print(String("Number of MPI Processes").pad_back(padlen, '.') + ": " + stringify(comm.size()));

  // create command line argument parser
  SimpleArgParser args(argc, argv);

  // check command line arguments
  Control::Domain::add_supported_pdc_args(args);
  args.support("threads");
  args.support("backend");
  args.support("level");
  args.support("vtk");
  args.support("mesh");
  args.support("nu");
  args.support("upsam");
  args.support("nl-tol");
  args.support("nl-damp");
  args.support("nl-stag-rate");
  args.support("nl-slow-rate");
  args.support("nl-max-iter");
  args.support("mg-min-iter");
  args.support("mg-boot-iter");
  args.support("mg-max-iter");
  args.support("mg-stag-rate");
  args.support("mg-adapt");
  args.support("mg-damp");
  args.support("mg-cycle");
  args.support("mg-smooth");
  args.support("smooth-plot");
  args.support("smooth-steps");
  args.support("smooth-damp");
  args.support("smooth-gmres");
  args.support("alpine");
  args.support("picard");
  args.support("v-max");
  args.support("deform");
  args.support("stokes");
  args.support("ext-stats");
  args.support("load-sol");
  args.support("save-sol");
  args.support("p-a");
  args.support("p-e");
  args.support("flush-mg");
  args.support("skip-analysis");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if (!unsupported.empty())
  {
    // print all unsupported options to cerr
    for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
      comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");
    comm.print_flush();

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

  int num_threads = 1;
  if((args.parse("threads", num_threads) < 0) || (num_threads < 1))
  {
    comm.print(std::cerr, "ERROR: Failed to parse '--threads <n>' option!");
    FEAT::Runtime::abort();
  }
  comm.print(String("Threads per Process").pad_back(padlen, '.') + ": " + stringify(num_threads));
#ifdef FEAT_HAVE_OMP
  omp_set_num_threads(num_threads);
  comm.print(String("Maximum OpenMP Threads").pad_back(padlen, '.') + ": " + stringify(omp_get_max_threads()));
#endif

  comm.print(String("Total CPU Cores").pad_back(padlen, '.') + ": " + stringify(std::size_t(num_threads) * comm.size()));

  if(args.check("backend") > 0)
  {
    String sbackend;
    args.parse("backend", sbackend);
    if(sbackend.compare_no_case("cuda") == 0)
      Backend::set_preferred_backend(PreferredBackend::cuda);
    else if(sbackend.compare_no_case("mkl") == 0)
      Backend::set_preferred_backend(PreferredBackend::mkl);
    else if(sbackend.compare_no_case("generic") == 0)
      Backend::set_preferred_backend(PreferredBackend::generic);
    else
    {
      comm.print(std::cerr, "ERROR: unknown backend specified: '" + sbackend + "'!");
      FEAT::Runtime::abort();
    }
  }
  comm.print(String("Floating Point Formats").pad_back(padlen, '.') + ": " + String(fp_typename));
  comm.print(String("Chosen Backend").pad_back(padlen, '.') + ": " + stringify(Backend::get_preferred_backend()));

  // create our mesh file reader
  std::unique_ptr<Geometry::MeshFileReader> mesh_reader(new Geometry::MeshFileReader);

  // read in the mesh files
  mesh_reader->add_mesh_files(comm, args.query("mesh")->second);

  // read the mesh file root markups
  mesh_reader->read_root_markup();
  String mesh_type = mesh_reader->get_meshtype_string();

  // is this the correct mesh type?
  if constexpr (dim == 2)
  {
    if(mesh_type != "conformal:hypercube:2:2")
    {
      comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'; expected 'conformal:hypercube:2:2'");
      FEAT::Runtime::abort();
    }
  }
  if constexpr (dim == 3)
  {
    if(mesh_type != "conformal:hypercube:3:3")
    {
      comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'; expected 'conformal:hypercube:3:3'");
      FEAT::Runtime::abort();
    }
  }

  // define our mesh type
  typedef Shape::Hypercube<dim> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, SystemDataType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
  //typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
  typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

  // create a time-stamp
  TimeStamp time_stamp;

  // create our domain control
  typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
  Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

  domain.parse_args(args);
  domain.set_desired_levels(args.query("level")->second);

  // create domain
  domain.create(*mesh_reader);
  domain.add_trafo_mesh_part_charts();

  // delete mesh reader
  mesh_reader.reset();

  // print partitioning info
  comm.print(String("Partition Info").pad_back(padlen, '.') + ": " + domain.get_chosen_parti_info());

  // plot our levels
  comm.print(String("Desired Levels").pad_back(padlen, '.') + ": " + domain.format_desired_levels());
  comm.print(String("Chosen  Levels").pad_back(padlen, '.') + ": " + domain.format_chosen_levels());
  comm.print(String("Mesh File Type").pad_back(padlen, '.') + ": " + mesh_type);
  comm.print(String("Transformation").pad_back(padlen, '.') + ": Isoparametric:2");
  comm.print(String("Velocity Space").pad_back(padlen, '.') + ": " + SpaceVeloType::name());
  comm.print(String("Pressure Space").pad_back(padlen, '.') + ": " + SpacePresType::name());

  // define our system level with (higher) system precision
  typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim, SystemDataType, IndexType> SystemLevelType;

  // define our solver level with (lower) solver precision
  typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim, SolverDataType, IndexType> SolverLevelType;

  BenchmarkSummary summary;
  BenchmarkStats statistics(domain.size_virtual());
  StopWatch watch_total_run;
  watch_total_run.start();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // solver Navier-Stokes or just Stokes?
  const bool navier = (args.check("stokes") < 0);
  // use nonlinear AlPiNe solver?
  const bool alpine = (args.check("alpine") >= 0);
  // use nonlinear Newton solver? (alpine overrides Newton)
  const bool newton = !alpine && (args.check("picard") < 0);
  // use deformation tensor or gradient tensor?
  const bool deform = (args.check("deform") >= 0);
  // unit-test mode?
  const bool testmode = (args.check("test-mode") >= 0);
  // collect/print extended solver statistics?
  const bool ext_stats = (args.check("ext-stats") >= 0);
  // print smoother iterations?
  const bool smooth_plot = (args.check("smooth-plot") >= 0);
  // use GMRES smoother?
  const bool smooth_gmres = (args.check("smooth-gmres") >= 0);
  // flush print after each MG iteration?
  const bool flush_mg = (args.check("flush-mg") >= 0);
  // skip solution analysis?
  const bool skip_analysis = (args.check("skip-analysis") >= 0);

  // viscosity parameter
  const SystemDataType nu = args.parse_default("nu", SystemDataType(1e-3));
  // maximum velocity, default: 2D: 0.3, 3D: 0.45
  const SystemDataType v_max = args.parse_default("v-max", SystemDataType(dim) * SystemDataType(0.15));
  // streamline diffusion parameter
  const SystemDataType upsam = args.parse_default("upsam", SystemDataType(0));
  // max. nonlinear solver iterations
  const Index nl_max_iter = args.parse_default("nl-max-iter", Index(25));
  // min. multigrid iterations
  const Index mg_min_iter = args.parse_default("mg-min-iter", Index(1));
  // max. multigrid boot iterations
  const Index mg_boot_iter = args.parse_default("mg-boot-iter", Index(3));
  // max. multigrid iterations
  const Index mg_max_iter = args.parse_default("mg-max-iter", Index(50));
  // number of smoothing steps
  const Index smooth_steps = args.parse_default("smooth-steps", Index(8));
  // smoother damping parameter
  const SolverDataType smooth_damp = args.parse_default("smooth-damp", SolverDataType(0.7));
  // multigrid damping parameter
  const SolverDataType mg_damp = args.parse_default("mg-damp", SolverDataType(1.0));
  // multigrid stagnation rate
  const SolverDataType mg_stag_rate = args.parse_default("mg-stag-rate", SolverDataType(1.05));
  // nonlinear solver damping parameter
  const SolverDataType nl_damp = args.parse_default("nl-damp", SolverDataType(1.0));
  // nonlinear stagnation rate
  const SystemDataType nl_stag_rate = args.parse_default("nl-stag-rate", SystemDataType(1.05));
  // nonlinear slowdown rate
  const SystemDataType nl_slow_rate = args.parse_default("nl-slow-rate", SystemDataType(0.1));
  // normalized tolerance for nonlinear solver
  const SystemDataType nl_tol_norm = args.parse_default("nl-tol", SystemDataType(1E-8));
  // multigrid cycle
  const Solver::MultiGridCycle mg_cycle = args.parse_default("mg-cycle", Solver::MultiGridCycle::V);
  // multigrid adaptive coarse grid correction
  Solver::MultiGridAdaptCGC mg_adapt = Solver::MultiGridAdaptCGC::Fixed;
  if(args.check("mg-adapt") > 0)
  {
    String s;
    args.parse("mg-adapt", s);
    if(s.compare_no_case("fixed") == 0)
      mg_adapt = Solver::MultiGridAdaptCGC::Fixed;
    else if(s.compare_no_case("defect") == 0)
      mg_adapt = Solver::MultiGridAdaptCGC::MinDefect;
    else if(s.compare_no_case("energy") == 0)
      mg_adapt = Solver::MultiGridAdaptCGC::MinEnergy;
    else
    {
      comm.print("ERROR: invalid parameter for --mg-adapt: '" + s + "'");
      Runtime::abort();
    }
  }
  /// multigrid smoothing type
  enum class MGSmooth
  {
    pre_only,
    post_only,
    both
  };
  MGSmooth mg_smooth = MGSmooth::both;
  if(args.check("mg-smooth") > 0)
  {
    String s;
    args.parse("mg-smooth", s);
    if(s.compare_no_case("pre") == 0)
      mg_smooth = MGSmooth::pre_only;
    else if(s.compare_no_case("post") == 0)
      mg_smooth = MGSmooth::post_only;
    else if(s.compare_no_case("both") == 0)
      mg_smooth = MGSmooth::both;
    else
    {
      comm.print("ERROR: invalid parameter for --mg-adapt: '" + s + "'");
      Runtime::abort();
    }
  }

  // can we use a direct coarse grid solver?
  const bool direct_coarse_solver = (domain.back_layer().comm().size() == 1) && Solver::direct_stokes_solver_available;

  // is our system actually non-linear?
  const bool non_linear = navier; // currently, only Navier-Stokes is non-linear

  {
    static constexpr std::size_t pl = padlen;
    static constexpr char pc = '.';
    comm.print("\nProblem Parameters:");
    comm.print(String("Nu").pad_back(pl, pc) + ": " + stringify(nu));
    comm.print(String("V-Max").pad_back(pl, pc) + ": " + stringify(v_max));
    comm.print(String("Upsam").pad_back(pl, pc) + ": " + stringify(upsam));
    comm.print(String("System").pad_back(pl, pc) + ": " + (navier ? "Navier-Stokes" : "Stokes"));
    comm.print(String("System Linearity").pad_back(pl, pc) + ": " + (non_linear ? "non-linear" : "linear"));
    comm.print(String("Diffusion Tensor").pad_back(pl, pc) + ": " + (deform ? "Deformation" : "Gradient"));
    comm.print(String("Nonlinear Solver").pad_back(pl, pc) + ": " + (alpine ? "AlPiNe" : newton ? "Newton" : "Picard"));
    //comm.print(String("Nonlinear Absolute Tolerance").pad_back(pl, pc) + ": " + stringify_fp_sci(nl_tol_abs));
    comm.print(String("Nonlinear Normalized Tolerance").pad_back(pl, pc) + ": " + stringify_fp_sci(nl_tol_norm));
    comm.print(String("Nonlinear Max Iterations").pad_back(pl, pc) + ": " + stringify(nl_max_iter));
    comm.print(String("Nonlinear Damping").pad_back(pl, pc) + ": " + stringify(nl_damp));
    comm.print(String("Nonlinear Stagnation Rate").pad_back(pl, pc) + ": " + stringify(nl_stag_rate));
    comm.print(String("Nonlinear Slowdown Rate").pad_back(pl, pc) + ": " + stringify(nl_slow_rate));
    comm.print(String("Multigrid Min  Iterations").pad_back(pl, pc) + ": " + stringify(mg_min_iter));
    comm.print(String("Multigrid Boot Iterations").pad_back(pl, pc) + ": " + stringify(mg_boot_iter));
    comm.print(String("Multigrid Max  Iterations").pad_back(pl, pc) + ": " + stringify(mg_max_iter));
    comm.print(String("Multigrid Stagnation Rate").pad_back(pl, pc) + ": " + stringify(mg_stag_rate));
    comm.print(String("Multigrid Damping").pad_back(pl, pc) + ": " + stringify(mg_damp));
    comm.print(String("Multigrid Cycle").pad_back(pl, pc) + ": " + stringify(mg_cycle));
    switch(mg_smooth)
    {
    case MGSmooth::pre_only:
      comm.print(String("Multigrid Smoothing").pad_back(pl, pc) + ": pre-smoothing only");
      break;
    case MGSmooth::post_only:
      comm.print(String("Multigrid Smoothing").pad_back(pl, pc) + ": post-smoothing only");
      break;
    case MGSmooth::both:
      comm.print(String("Multigrid Smoothing").pad_back(pl, pc) + ": both pre- and post-smoothing");
      break;
    }
    switch(mg_adapt)
    {
    case Solver::MultiGridAdaptCGC::Fixed:
      comm.print(String("Multigrid Adaptivity").pad_back(pl, pc) + ": disabled");
      break;
    case Solver::MultiGridAdaptCGC::MinEnergy:
      comm.print(String("Multigrid Adaptivity").pad_back(pl, pc) + ": min energy");
      break;
    case Solver::MultiGridAdaptCGC::MinDefect:
      comm.print(String("Multigrid Adaptivity").pad_back(pl, pc) + ": min defect");
      break;
    }
    comm.print(String("Smoother Type").pad_back(pl, pc) + ": " + String(smooth_gmres ? "GMRES-AmaVanka" : "Richardson-AmaVanka"));
    comm.print(String("Smoother Steps").pad_back(pl, pc) + ": " + stringify(smooth_steps));
    comm.print(String("Smoother Damping").pad_back(pl, pc) + ": " + stringify(smooth_damp));
    if(direct_coarse_solver)
      comm.print(String("Coarse Grid Solver").pad_back(pl, pc) + ": direct solver");
    else
      comm.print(String("Coarse Grid Solver").pad_back(pl, pc) + ": GMRES[16]-AmaVanka");
  }

  // enable solver expressions if extended statistics are desired
  if(ext_stats)
    Statistics::enable_solver_expressions = true;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::deque<std::shared_ptr<SolverLevelType>> solver_levels;

  const Index num_levels = Index(domain.size_physical());

  // create system levels
  for (Index i(0); i < num_levels; ++i)
  {
    solver_levels.push_back(std::make_shared<SolverLevelType>());
  }

  // get the finest solver level
  SolverLevelType& the_solver_level = *solver_levels.front();

  // create a system level for the finest domain level
  SystemLevelType the_system_level;

  // cubature for assembly
  const String cubature("gauss-legendre:3");

  // cubature for post-processing
  Cubature::DynamicFactory cubature_postproc("gauss-legendre:5");

  StopWatch watch_asm;
  watch_asm.start();

  // fetch our finest level
  DomainLevelType& the_domain_level = *domain.front();

  for (Index i(0); i < num_levels; ++i)
  {
    if(num_threads > 1)
    {
      domain.at(i)->domain_asm.set_threading_strategy(Assembly::ThreadingStrategy::layered);
      domain.at(i)->domain_asm.set_max_worker_threads(num_threads);
    }
    domain.at(i)->domain_asm.compile_all_elements();
  }

  // assemble gates
  the_system_level.assemble_gates(domain.front());
  the_solver_level.convert(the_system_level);
  for (Index i(1); i < num_levels; ++i)
  {
    solver_levels.at(i)->assemble_gates(domain.at(i));
  }

  // assemble muxers and transfers
  for (Index i(0); (i < domain.size_physical()) && ((i+1) < domain.size_virtual()); ++i)
  {
    solver_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
    if((i+1) < domain.size_physical())
      solver_levels.at(i)->assemble_transfers(*solver_levels.at(i+1), domain.at(i), domain.at(i+1), cubature, true);
    else
      solver_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature, true);
  }

  // assemble base splitter on finest level if required
  if((args.check("save-joined-sol") >= 0) || (args.check("load-joined-sol") >= 0))
    the_system_level.assemble_base_splitters(domain.front());

  // collect some finest-level statistics
  {
    auto tv = the_system_level.gate_velo._freqs.clone(LAFEM::CloneMode::Deep);
    tv.format(1.0);
    Index velo_dofs = Index(the_system_level.gate_velo.dot(tv, tv));
    Index locp_dofs = the_system_level.gate_pres._freqs.size();
    Index pres_dofs = Index(the_system_level.gate_pres.sum(SystemDataType(locp_dofs)));
    statistics.counts[Counts::velo_dofs] = velo_dofs/dim;
    statistics.counts[Counts::pres_dofs] = pres_dofs;
    statistics.counts[Counts::total_dofs] = velo_dofs+pres_dofs;
    statistics.counts[Counts::fine_elements] = domain.front()->get_mesh().get_num_elements();
  }

  // collect some all-level statistics
  for (Index i(1); i < num_levels; ++i)
  {
    statistics.counts[Counts::all_elements] += domain.at(i)->get_mesh().get_num_elements();
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // assemble basic matrices
  the_system_level.assemble_grad_div_matrices(domain.front()->domain_asm, domain.front()->space_velo, domain.front()->space_pres, cubature);
  for(Index i(0); i < num_levels; ++i)
  {
    if(i == Index(0))
    {
      solver_levels.at(i)->matrix_b.local().convert(the_system_level.matrix_b.local());
      solver_levels.at(i)->matrix_d.local().convert(the_system_level.matrix_d.local());
    }
    else
    {
      solver_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->domain_asm, domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
    }

    // assemble matrix structure for A and format local A matrix to 1, because
    // some solvers/smoothers may not perform symbolic factorization otherwise
    solver_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
    solver_levels.at(i)->matrix_a.local().format(SolverDataType(1));

    solver_levels.at(i)->compile_system_matrix();
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // our inflow BC function for bench1
  Analytic::Common::DFG95SteadyInflowFunction<dim, SystemDataType> inflow_func(v_max);

  // our inflow BC function for bench7
  Tiny::Vector<SystemDataType, dim> pipe_origin, pipe_axis;
  pipe_origin = SystemDataType(0.205);
  pipe_axis = SystemDataType(0);
  pipe_axis[0] = SystemDataType(1);
  Analytic::Common::PoiseuillePipeFlow<SystemDataType, dim> pipe_inflow_func(pipe_origin, pipe_axis, SystemDataType(0.205), v_max);

  // the names of the mesh parts on which to assemble
  std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

  for(Index i(0); i < num_levels; ++i)
  {
    // create unit-filter assembler
    Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow_1, unit_asm_inflow_2, unit_asm_noflow;

    // loop over all boundary parts except for the right one, which is outflow
    for(const auto& name : part_names)
    {
      // skip non-boundary mesh-parts
      if(!name.starts_with("bnd:"))
        continue;

      // try to fetch the corresponding mesh part node
      auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
      XASSERT(mesh_part_node != nullptr);

      auto* mesh_part = mesh_part_node->get_mesh();
      if (mesh_part != nullptr)
      {
        if(name == "bnd:l")
        {
          // inflow (bench1)
          unit_asm_inflow_1.add_mesh_part(*mesh_part);
        }
        else if(name == "bnd:in")
        {
          // inflow (bench7)
          unit_asm_inflow_2.add_mesh_part(*mesh_part);
        }
        else if((name != "bnd:r") && (name != "bnd:out"))
        {
          // outflow
          unit_asm_noflow.add_mesh_part(*mesh_part);
        }
      }
    }

    // assemble system filters on finest level
    if(i == Index(0))
    {
      auto& fil_loc_v1 = the_system_level.filter_velo.local();
      unit_asm_inflow_1.assemble(fil_loc_v1, domain.at(i)->space_velo, inflow_func);
      unit_asm_inflow_2.assemble(fil_loc_v1, domain.at(i)->space_velo, pipe_inflow_func);
      unit_asm_noflow.assemble(fil_loc_v1, domain.at(i)->space_velo);
      the_system_level.compile_system_filter();
    }

    // assemble in SystmemDataType and then convert to SolverDataType
    SystemLevelType::LocalVeloFilter fil_loc_v;
    unit_asm_inflow_1.assemble(fil_loc_v, domain.at(i)->space_velo, inflow_func);
    unit_asm_inflow_2.assemble(fil_loc_v, domain.at(i)->space_velo, pipe_inflow_func);
    unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo);
    solver_levels.at(i)->filter_velo.local().convert(fil_loc_v);
    solver_levels.at(i)->compile_system_filter();
  }

  // finally, compile the local type-1 matrices
  for(Index i(0); i < num_levels; ++i)
  {
    solver_levels.at(i)->compile_local_matrix_sys_type1();
  }

  watch_asm.stop();
  Statistics::toe_assembly = watch_asm.elapsed();

  // accumulate sizes
  for(Index i(0); i < num_levels; ++i)
  {
    statistics.bytes[Bytes::mesh] += domain.at(i)->get_mesh_node()->bytes();
    statistics.bytes[Bytes::gate] += solver_levels.at(i)->gate_sys.bytes();
    statistics.bytes[Bytes::muxer] += solver_levels.at(i)->coarse_muxer_sys.bytes();
    statistics.bytes[Bytes::matrix] += solver_levels.at(i)->matrix_sys.local().bytes();
    statistics.bytes[Bytes::transfer] += solver_levels.at(i)->transfer_sys.bytes();
    const auto& loc_a = solver_levels.at(i)->matrix_sys.local().block_a();
    const auto& loc_b = solver_levels.at(i)->matrix_sys.local().block_b();
    const auto& loc_d = solver_levels.at(i)->matrix_sys.local().block_d();
    statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_a.used_elements() + loc_a.rows() + Index(1));
    statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_b.used_elements() + loc_b.rows() + Index(1));
    statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_d.used_elements() + loc_d.rows() + Index(1));
    statistics.bytes[Bytes::matrix_values] += sizeof(SolverDataType) * std::size_t(loc_a.template used_elements<LAFEM::Perspective::pod>());
    statistics.bytes[Bytes::matrix_values] += sizeof(SolverDataType) * std::size_t(loc_b.template used_elements<LAFEM::Perspective::pod>());
    statistics.bytes[Bytes::matrix_values] += sizeof(SolverDataType) * std::size_t(loc_d.template used_elements<LAFEM::Perspective::pod>());
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // get our global solver system types
  typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
  typedef typename SolverLevelType::GlobalSystemVector GlobalSolverVector;

  // create solver vectors
  GlobalSolverVector vec_mg_sol = the_solver_level.create_global_vector_sys();
  GlobalSolverVector vec_mg_rhs = the_solver_level.create_global_vector_sys();
  GlobalSolverVector vec_mg_def = the_solver_level.create_global_vector_sys();
  GlobalSolverVector vec_mg_cor = the_solver_level.create_global_vector_sys();

  // create new vectors
  GlobalSystemVector vec_sol = the_system_level.create_global_vector_sys();
  GlobalSystemVector vec_def = the_system_level.create_global_vector_sys();
  GlobalSystemVector vec_cor = the_system_level.create_global_vector_sys();
  GlobalSystemVector vec_def_unsynced = the_system_level.create_global_vector_sys();

  // format the vectors
  vec_sol.format();
  vec_def.format();

  // assemble a temporary vector representing the constant 1-functional and compute its norm;
  // this is used to normalize the defect norms, which depend on the mesh resolution
  {
    Analytic::Common::ConstantVectorFunction<dim, SystemDataType> one_v(SystemDataType(1));
    Analytic::Common::ConstantFunction<dim, SystemDataType> one_p(SystemDataType(1));
    vec_cor.format();
    Assembly::assemble_force_function_vector(the_domain_level.domain_asm, vec_cor.local().template at<0>(), one_v, the_domain_level.space_velo, cubature);
    Assembly::assemble_force_function_vector(the_domain_level.domain_asm, vec_cor.local().template at<1>(), one_p, the_domain_level.space_pres, cubature);
    vec_cor.sync_0();
  }
  const SystemDataType norm1_v = Math::sqrt(the_system_level.gate_velo.dot(vec_cor.local().template at<0>(), vec_cor.local().template at<0>()));
  const SystemDataType norm1_p = Math::sqrt(the_system_level.gate_pres.dot(vec_cor.local().template at<1>(), vec_cor.local().template at<1>()));
  const SystemDataType norm1 = SystemDataType(1);//Math::sqrt(norm1_v*norm1_v + norm1_p*norm1_p);

  comm.print(String("\nNorm of 1-Functional V/P").pad_back(padlen, '.') + ".: " + stringify_fp_sci(norm1_v) + " / " + stringify_fp_sci(norm1_p));
  comm.print(String("Defect Normalization Factor").pad_back(padlen, '.') + ": " + stringify_fp_sci(SystemDataType(1) / norm1));

  // and filter them
  the_system_level.filter_sys.filter_sol(vec_sol);
  the_system_level.filter_sys.filter_rhs(vec_def);

  {
    // count non-zeros in a and b
    statistics.counts[Counts::nnze_a] = the_solver_level.matrix_sys.local().block_a().used_elements();
    statistics.counts[Counts::nnze_b] = the_solver_level.matrix_sys.local().block_b().used_elements();
    statistics.counts[Counts::nnze_total] = the_solver_level.matrix_sys.local().template used_elements<LAFEM::Perspective::pod>();
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const MeshPartType* mesh_part_bnd_c = the_domain_level.get_mesh_node()->find_mesh_part("bnd:c");
  const MeshPartType* mesh_part_bnd_s = the_domain_level.get_mesh_node()->find_mesh_part("bnd:s");
  const MeshPartType* mesh_part_bnd_sphere = the_domain_level.get_mesh_node()->find_mesh_part("bnd:sphere");
  const MeshPartType* mesh_part_bnd_l = the_domain_level.get_mesh_node()->find_mesh_part("bnd:l");
  const MeshPartType* mesh_part_bnd_r = the_domain_level.get_mesh_node()->find_mesh_part("bnd:r");
  const MeshPartType* mesh_part_bnd_in = the_domain_level.get_mesh_node()->find_mesh_part("bnd:in");
  const MeshPartType* mesh_part_bnd_out = the_domain_level.get_mesh_node()->find_mesh_part("bnd:out");
  const MeshPartType* mesh_part_inner_u = the_domain_level.get_mesh_node()->find_mesh_part("inner:u");
  const MeshPartType* mesh_part_inner_l = the_domain_level.get_mesh_node()->find_mesh_part("inner:l");

  // create trace assembler for body force assembly
  Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> body_force_asm(the_domain_level.trafo);
  if(mesh_part_bnd_c != nullptr)
    body_force_asm.add_mesh_part(*mesh_part_bnd_c);
  if(mesh_part_bnd_s != nullptr)
    body_force_asm.add_mesh_part(*mesh_part_bnd_s);
  if(mesh_part_bnd_sphere != nullptr)
    body_force_asm.add_mesh_part(*mesh_part_bnd_sphere);
  body_force_asm.compile();

  // create trace assembler for upper flux
  Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> flux_asm_u(the_domain_level.trafo);
  if(mesh_part_inner_u != nullptr)
    flux_asm_u.add_mesh_part(*mesh_part_inner_u);
  flux_asm_u.compile();

  // create trace assembler for lower flux
  Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> flux_asm_l(the_domain_level.trafo);
  if(mesh_part_inner_l != nullptr)
    flux_asm_l.add_mesh_part(*mesh_part_inner_l);
  flux_asm_l.compile();

  // create trace assembler for in flux
  Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> flux_asm_in(the_domain_level.trafo);
  if(mesh_part_bnd_l != nullptr)
    flux_asm_in.add_mesh_part(*mesh_part_bnd_l);
  if(mesh_part_bnd_in != nullptr)
    flux_asm_in.add_mesh_part(*mesh_part_bnd_in);
  flux_asm_in.compile();

  // create trace assembler for in flux
  Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> flux_asm_out(the_domain_level.trafo);
  if(mesh_part_bnd_r != nullptr)
    flux_asm_out.add_mesh_part(*mesh_part_bnd_r);
  if(mesh_part_bnd_out != nullptr)
    flux_asm_out.add_mesh_part(*mesh_part_bnd_out);
  flux_asm_out.compile();

  // unmap pressure evaluation points p_a and p_e
  Trafo::InverseMappingData<SystemDataType, dim> point_iv_a, point_iv_e;
  {
    typedef Trafo::InverseMapping<typename SpacePresType::TrafoType, SystemDataType> InvMappingType;
    InvMappingType inv_mapping(the_domain_level.trafo);

    // reference pressure points
    typename InvMappingType::ImagePointType v_a, v_e;
    if constexpr(dim == 2)
    {
      v_a[0] = SystemDataType(0.15);
      v_e[0] = SystemDataType(0.25);
      v_a[1] = v_e[1] = SystemDataType(0.2);
      args.parse("p-a", v_a[0], v_a[1]);
      args.parse("p-e", v_e[0], v_e[1]);
    }
    else
    {
      v_a[0] = SystemDataType(0.45);
      v_e[0] = SystemDataType(0.55);
      v_a[1] = v_e[1] = SystemDataType(0.2);
      v_a[2] = v_e[2] = SystemDataType(0.205);
      args.parse("p-a", v_a[0], v_a[1], v_a[2]);
      args.parse("p-e", v_e[0], v_e[1], v_e[2]);
    }

    // unmap points
    point_iv_a = inv_mapping.unmap_point(v_a, true);
    point_iv_e = inv_mapping.unmap_point(v_e, true);
  }

  // create characteristic function vector for circle/sphere boundary
  // this is needed for the volumetric drag/lift computation
  LAFEM::DenseVector<SystemDataType, IndexType> vec_char(vec_sol.local().template at<0>().size());
  vec_char.format(SystemDataType(0));
  {
    LAFEM::UnitFilter<SystemDataType, IndexType> filter_char;
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    if(mesh_part_bnd_c != nullptr)
      unit_asm.add_mesh_part(*mesh_part_bnd_c);
    if(mesh_part_bnd_s != nullptr)
      unit_asm.add_mesh_part(*mesh_part_bnd_s);
    if(mesh_part_bnd_sphere != nullptr)
      unit_asm.add_mesh_part(*mesh_part_bnd_sphere);
    unit_asm.assemble(filter_char, the_domain_level.space_velo);
    filter_char.get_filter_vector().format(SystemDataType(1));
    filter_char.filter_sol(vec_char);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // create a multigrid solver
  auto multigrid_hierarchy = std::make_shared<
    Solver::MultiGridHierarchy<
      typename SolverLevelType::GlobalSystemMatrix,
      typename SolverLevelType::GlobalSystemFilter,
      typename SolverLevelType::GlobalSystemTransfer>>(domain.size_virtual());

  // array of AmaVanka pointers - this is only required to collect the memory usage
  // statistics of the Vankas, as we need to query the memory usage after factorization
  std::deque<
    std::shared_ptr<
      Solver::AmaVanka<
        typename SolverLevelType::LocalSystemMatrix,
        typename SolverLevelType::LocalSystemFilter>>> ama_vankas;

  // push levels into multigrid
  for(std::size_t i(0); i < solver_levels.size(); ++i)
  {
    SolverLevelType& lvl = *solver_levels.at(i);

    // not the coarse level?
    if((i+1) < domain.size_virtual())
    {
      // create an GMRES-Schwarz-AmaVanka smoother
      auto vanka = Solver::new_amavanka(lvl.local_matrix_sys_type1, lvl.filter_sys.local());
      vanka->set_skip_singular(true);
      ama_vankas.push_back(vanka);
      auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
      schwarz->set_ignore_status(true);
      std::shared_ptr<Solver::IterativeSolver<GlobalSolverVector>> smoother;
      if(smooth_gmres)
        smoother = Solver::new_gmres(lvl.matrix_sys, lvl.filter_sys, smooth_steps, 0.0, schwarz);
      else
        smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);
      //smoother = Solver::new_descent(lvl.matrix_sys, lvl.filter_sys, Solver::DescentVariant::defect, schwarz);
      smoother->set_min_iter(smooth_steps);
      smoother->set_max_iter(smooth_steps);
      smoother->set_plot_name("Smoother[" + stringify(i) + "]");
      smoother->set_plot_mode(smooth_plot ? Solver::PlotMode::iter : Solver::PlotMode::none);
      // pre-/post-smoother?
      auto pre_smoother = smoother;
      auto post_smoother = smoother;
      if(mg_smooth == MGSmooth::pre_only)
        post_smoother.reset();
      if(mg_smooth == MGSmooth::post_only)
        pre_smoother.reset();
      multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, pre_smoother, post_smoother, smoother);
    }
    else if(direct_coarse_solver)
    {
      // create a direct coarse grid solver
      auto coarse_solver = Solver::new_direct_stokes_solver(lvl.matrix_sys, lvl.filter_sys);
      multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, coarse_solver);
    }
    else
    {
      // create FGMRES-AmaVanka coarse grid solver
      auto vanka = Solver::new_amavanka(lvl.local_matrix_sys_type1, lvl.filter_sys.local());
      vanka->set_skip_singular(true);
      ama_vankas.push_back(vanka);
      auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
      schwarz->set_ignore_status(true);
      auto coarse_solver = Solver::new_gmres(lvl.matrix_sys, lvl.filter_sys, 16, SolverDataType(0.9), schwarz);
      //auto coarse_solver = Solver::new_bicgstab(lvl.matrix_sys, lvl.filter_sys, schwarz);
      coarse_solver->set_max_iter(500);
      coarse_solver->set_tol_rel(SolverDataType(1e-3));
      coarse_solver->set_plot_name("CoarseSolver");
      multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, coarse_solver);
    }
  }

  // create our multigrid preconditioner object
  auto multigrid = Solver::new_multigrid(multigrid_hierarchy, mg_cycle);
  multigrid->set_adapt_cgc(mg_adapt);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // create a few watches
  StopWatch watch_nonlin_loop;
  StopWatch watch_nonlin_def_asm;
  StopWatch watch_nonlin_mat_asm;
  StopWatch watch_nonlin_solver_init;
  StopWatch watch_nonlin_solver_apply;

  // initialize solver
  watch_nonlin_solver_init.start();
  multigrid_hierarchy->init_symbolic();
  multigrid->init_symbolic();
  watch_nonlin_solver_init.stop();

  // accumulate vanka sizes
  if(!ama_vankas.empty())
    statistics.counts[Counts::vanka_data] = Index(ama_vankas.front()->data_size());
  statistics.bytes[Bytes::vanka] = 0ull;
  for(auto& v : ama_vankas)
    statistics.bytes[Bytes::vanka] += v->bytes();

  // load distributed solution?
  bool loaded_solution = false;
  if(args.check("load-sol") > 0)
  {
    StopWatch watch_checkpoint;
    watch_checkpoint.start();

    String save_name;
    args.parse("load-sol", save_name);
    comm.print("\nReading (partitioned) initial solution from '" + save_name + "'...");

    // read file into streams
    BinaryStream bs_com, bs_sol;
    DistFileIO::read_combined(bs_com, bs_sol, save_name, comm);

    // parse vector from stream
    buf2vec(vec_sol.local(), bs_sol.container());

    // apply our solution filter in case the BCs have changed
    the_system_level.filter_sys.filter_sol(vec_sol);

    watch_checkpoint.stop();
    statistics.times[Times::checkpoint] += watch_checkpoint.elapsed();

    loaded_solution = true;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // vector of all non-linear absolute defect norms (front: latest defect, back: first defect)
  std::deque<SystemDataType> nl_defs;
  //std::deque<SolverDataType> mg_defs;

  watch_nonlin_loop.start();

  // print header
  comm.print(String("\n") + String(100, '=') + "\n");
  //         "Newton...:   1 : 1.887920e-04 / 9.286588e-08 / 1.887920e-04 / 1.650402e-02 / 1.650e-02 | 2.7238e-05"
  comm.print("Solver      it : Momentum     / Continuity   / Absolute     / Normalized   / Improve   | MG Norm Tol");
  comm.print_flush();

  // keep track whether solvers are initialized numerically
  bool solvers_numeric_init = false;

  // nonlinear loop
  for(Index nl_iter(0); nl_iter <= nl_max_iter; ++nl_iter)
  {
    statistics.counts[Counts::nonlin_iter] = nl_iter;
    comm.print(String(100, '-'));

    // --------------------------------------------------------------------------------------------
    // PHASE 1: assemble non-linear defect vector
    // --------------------------------------------------------------------------------------------

    // set up Burgers assembly job for our defect vector
    Assembly::BurgersBlockedVectorAssemblyJob<typename SystemLevelType::LocalVeloVector, SpaceVeloType> burgers_def_job(
      vec_def.local().template at<0>(), vec_sol.local().template at<0>(),
      vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature);
    burgers_def_job.deformation = deform;
    burgers_def_job.nu = -nu;
    // enable convection if we're solving Navier-Stokes
    burgers_def_job.beta = -SystemDataType(navier ? 1 : 0);

    // assemble nonlinear defect vector
    watch_nonlin_def_asm.start();
    vec_def.format();
    // assemble burgers operator defect
    the_domain_level.domain_asm.assemble(burgers_def_job);
    // compute remainder of defect vector
    the_system_level.matrix_b.local().apply(
      vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
    the_system_level.matrix_d.local().apply(
      vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);
    // store the unsynchronized and unfiltered defect for later body forces computation
    vec_def_unsynced.copy(vec_def);
    // synchronize and filter the defect
    vec_def.sync_0();
    the_system_level.filter_sys.filter_def(vec_def);
    watch_nonlin_def_asm.stop();

    // --------------------------------------------------------------------------------------------
    // PHASE 2: analyze defect and check stopping criterions
    // --------------------------------------------------------------------------------------------

    // compute defect norm for momentum and continuity equations
    const SystemDataType nl_def_v = Math::sqrt(the_system_level.gate_velo.dot(vec_def.local().template at<0>(), vec_def.local().template at<0>()));
    const SystemDataType nl_def_p = Math::sqrt(the_system_level.gate_pres.dot(vec_def.local().template at<1>(), vec_def.local().template at<1>()));
    const SystemDataType nl_def_abs = Math::sqrt(nl_def_v*nl_def_v + nl_def_p*nl_def_p);

    // compute defect norm
    const SystemDataType nl_def_norm = nl_def_abs / norm1;
    const SystemDataType nl_def_prev = (nl_defs.empty() ? SystemDataType(1) : nl_defs.front());
    const SystemDataType nl_def_improve = (nl_defs.empty() ? SystemDataType(0) : nl_def_abs / nl_def_prev);
    nl_defs.push_front(nl_def_abs);

    // build our non-linear iteration plot line
    String line = !navier ? "Stokes...:" : (newton || (alpine && ((nl_iter & 1u) != 0u))) ? "Newton...:" : "Picard...:";
    line += stringify(nl_iter).pad_front(4) + " : ";
    line += stringify_fp_sci(nl_def_v, 6) + " / ";
    line += stringify_fp_sci(nl_def_p, 6) + " / ";
    line += stringify_fp_sci(nl_def_abs, 6) + " / ";
    line += stringify_fp_sci(nl_def_norm) + " / ";
    line += stringify_fp_sci(nl_def_improve, 3);

    if(nl_def_norm > SystemDataType(1E+3))
    {
      comm.print(line);
      comm.print("\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
      if(testmode)
        comm.print("Test-Mode: FAILED");
      comm.print_flush();
      return 1;
    }
    /*else if(nl_def_abs < nl_tol_abs)
    {
      comm.print(line);
      comm.print("\nNonlinear solver converged!");
      comm.print_flush();
      break;
    }*/
    else if(nl_def_norm < nl_tol_norm)
    {
      comm.print(line);
      comm.print("\nNonlinear solver converged!");
      comm.print_flush();
      break;
    }
    else if(nl_iter >= nl_max_iter)
    {
      comm.print(line);
      comm.print("\nMaximum iterations reached!");
      comm.print_flush();
      break;
    }
    // the defect should always decrease after a slowdown phase
    else if((nl_iter >= 3) && (nl_def_prev < nl_stag_rate*nl_def_abs))
    {
      comm.print(line);
      comm.print("\nNonlinear solver stagnated!");
      comm.print_flush();
      break;
    }
    // in the case of Newton, the convergence speed should accelerate after a short warmup phase, i.e. the defect
    // improvement of this iteration should be better than the defect improvement of the previous iteration, so:
    // defs[1]/defs[2] < defs[0]/defs[1] <==> defs[1]^2 < defs[0]*defs[2]
    // if this criterion fails, then this usually means that we have reached the end of our solver precision
    // and there is no point in iterating even further, as the quality of our solution cannot improve anymore
    // we damp the right-hand-side by a user-chosen slowdown rate for better control in case of problems
    else if(newton && (nl_iter >= 5) && (nl_defs[1u]*nl_defs[1u] < nl_slow_rate*nl_defs[0u]*nl_defs[2u]))
    {
      comm.print(line);
      comm.print("\nNonlinear solver convergence slowdown!");
      comm.print_flush();
      break;
    }
    // we also check for slowdown in AlPiNe, but only in Newton iterations and we have to compare to the previous
    // Newton defect, so in the AlPiNe case we have to check
    // defs[2]/defs[4] < defs[0]/defs[2] <==> defs[2]^2 < defs[0]*defs[4]
    // note that we have to check (nl_iter % 2u == 0u) to check whether the last iteration was a Newton iteration
    else if(alpine && (nl_iter >= 5) && (nl_iter % 2u == 0u) && (nl_defs[2u]*nl_defs[2u] < nl_slow_rate*nl_defs[0u]*nl_defs[4u]))
    {
      comm.print(line);
      comm.print("\nNonlinear solver convergence slowdown!");
      comm.print_flush();
      break;
    }

    // --------------------------------------------------------------------------------------------
    // PHASE 3: compute multigrid solver tolerance and print line
    // --------------------------------------------------------------------------------------------

    // compute tolerance for multigrid; note that we will normalize the RHS vector for the MG,
    // so therefore absolute and relative stopping criterions are identical for the MG here
    SolverDataType mg_tol = SolverDataType(1E-4);
    if(!non_linear)
    {
      // linear system: just try to reach the nonlinear tolerance and add 5% as a buffer
      mg_tol = SolverDataType(0.95) * SolverDataType(nl_tol_norm);
    }
    else if(nl_iter > Index(0))
    {
      // in the case of AlPiNe, things are a bit more complicated
      if(alpine)
      {
        // the first three iterations are a warm-up phase consisting of "Picard-damped-Newton-Picard" iteration
        if(nl_iter > Index(2))
        {
          // let's also compute the previous improvement of the defect
          const SystemDataType nl_def_improve_prev = nl_defs.at(1u) / nl_defs.at(2u);

          // we expect quadratic convergence for the Newton steps, but only linear in Picard steps
          if((nl_iter % 2u) != 0u)
          {
            // Newton iteration: expect quadratic convergence w.r.t. last Newton iteration
            mg_tol = SolverDataType(nl_def_norm * nl_def_improve_prev * nl_def_improve_prev) * SolverDataType(0.1);
          }
          else
          {
            // Picard iteration: expect linear convergence w.r.t. last Picard iteration
            mg_tol = SolverDataType(nl_def_norm * nl_def_improve_prev) * SolverDataType(0.1);
          }
        }
        else
        {
          // warm-up phase: only expect linear convergence anyways
          mg_tol = SolverDataType(nl_def_norm * nl_def_improve) * SolverDataType(0.1);
        }
      }
      // in the case of Newton, we expect quadratic convergence, so try to square the defect
      // improvement of the last step and try to gain one more digit as a buffer
      else if(newton)
        mg_tol = SolverDataType(nl_def_norm * nl_def_improve * nl_def_improve) * SolverDataType(0.1);
      // and in the case of Picard, just linear convergence, so try to keep the improvement
      else
        mg_tol = SolverDataType(nl_def_norm * nl_def_improve) * SolverDataType(0.1);

      // don't try go below 0.01 * nl_tol_norm, as this is usually just a waste of time
      mg_tol = Math::max(mg_tol, SolverDataType(nl_tol_norm) * SolverDataType(0.01));
    }

    // make sure that we gain at least 2 digits
    mg_tol = Math::min(mg_tol, SolverDataType(nl_def_norm * 0.01));

    // in any case, do not try to solve for more than eps^(0.9), since this would likely
    // stagnate or even diverge depending on the chosen smoother
    mg_tol = Math::max(mg_tol, SolverDataType(nl_def_norm)*Math::pow(Math::eps<SolverDataType>(), SolverDataType(0.9)));

    // print our current nonlinear iteration line
    line += " | ";
    line += stringify_fp_sci(mg_tol, 4);
    comm.print(line);
    comm.print_flush();

    // --------------------------------------------------------------------------------------------
    // PHASE 4: assemble Burgers matrices on all levels
    // --------------------------------------------------------------------------------------------

    // we only have to assemble the system matrix in the non-linear case or in the first step of
    // a linear case
    watch_nonlin_mat_asm.start();
    if(non_linear || (nl_iter == 0u))
    {
      // get a clone of the global velocity vector
      auto vec_sol_vc = vec_sol.local().template at<0>().clone();
      typename SolverLevelType::GlobalVeloVector vec_conv(&the_solver_level.gate_velo);
      vec_conv.local().convert(vec_sol_vc);

      // initialize velocity norm for streamline diffusion (if enabled)
      SolverDataType sd_v_norm = SolverDataType(0);

      // loop over all system levels
      for(std::size_t i(0); i < solver_levels.size(); ++i)
      {
        auto& loc_mat_a = solver_levels.at(i)->matrix_sys.local().block_a();
        loc_mat_a.format();

        // set up Burgers assembly job
        Assembly::BurgersBlockedMatrixAssemblyJob<
          typename SolverLevelType::LocalMatrixBlockA, SpaceVeloType, typename SolverLevelType::LocalVeloVector>
          burgers_mat_job(loc_mat_a, vec_conv.local(), domain.at(i)->space_velo, cubature);

        burgers_mat_job.deformation = deform;
        burgers_mat_job.nu = SolverDataType(nu);
        burgers_mat_job.beta = SolverDataType(navier && (loaded_solution || (nl_iter > 0)) ? 1 : 0);
        burgers_mat_job.frechet_beta = SolverDataType(0);
        if(navier)
        {
          if(alpine)
          {
            // If we start from scratch, then damp the first Newton iteration by 0.5
            if(!loaded_solution && (nl_iter == 1u))
              burgers_mat_job.frechet_beta = SolverDataType(0.5);
            // In all other cases, use full alternating Picard-Newton scheme
            else
              burgers_mat_job.frechet_beta = SolverDataType(nl_iter & 1u);
          }
          else if(newton)
          {
            // If we loaded a solution from file, then go directly to full Newton
            if(loaded_solution)
              burgers_mat_job.frechet_beta = SolverDataType(1);
            // Otherwise make a single Picard iteration in the first nonlinear step
            else if(nl_iter == 0)
              burgers_mat_job.frechet_beta = SolverDataType(0);
            // Followed by a single damped Newton iteration in the second nonlinear step
            else if(nl_iter == 1)
              burgers_mat_job.frechet_beta = SolverDataType(0.5);
            // And switch to full Newton from the third nonlinear step on
            else
              burgers_mat_job.frechet_beta = SolverDataType(1);
          }
        }

        // set convection vector in case streamline diffusion is enabled
        burgers_mat_job.sd_delta = SolverDataType(upsam);
        burgers_mat_job.sd_nu = SolverDataType(nu);
        if(i == size_t(0))
        {
          burgers_mat_job.set_sd_v_norm(vec_conv);
          sd_v_norm = burgers_mat_job.sd_v_norm;
        }
        else
        {
          // use fine mesh norm
          burgers_mat_job.sd_v_norm = sd_v_norm;
        }

        // assemble our system matrix
        domain.at(i)->domain_asm.assemble(burgers_mat_job);
        solver_levels.at(i)->compile_local_matrix_sys_type1();

        // restrict our convection vector
        if((i+1) >= domain.size_virtual())
          break;

        // does this process have another system level?
        if((i+1) < solver_levels.size())
        {
          // create a coarse mesh velocity vector
          auto vec_crs = solver_levels.at(i+1)->matrix_a.create_vector_l();

          // truncate fine mesh velocity vector
          solver_levels.at(i)->transfer_velo.trunc(vec_conv, vec_crs);

          // the coarse vector is our next convection vector
          vec_conv = std::move(vec_crs);
        }
        else
        {
          // this process is a child, so send truncation to parent
          solver_levels.at(i)->transfer_velo.trunc_send(vec_conv);
        }
      }
    }
    watch_nonlin_mat_asm.stop();

    // --------------------------------------------------------------------------------------------
    // PHASE 4: solve linear system by (lower-precision) multigrid
    // --------------------------------------------------------------------------------------------

    // initialize linear solver
    if(non_linear || (nl_iter == 0u))
    {
      watch_nonlin_solver_init.start();
      multigrid_hierarchy->init_numeric();
      multigrid->init_numeric();
      watch_nonlin_solver_init.stop();
      solvers_numeric_init = true;
    }

    // solve linear system
    watch_nonlin_solver_apply.start();

    // we normalize our RHS vector here to avoid small vector entries in case of mixed-precision
    vec_mg_rhs.local().convert(vec_def.local());
    vec_mg_rhs.scale(vec_mg_rhs, SolverDataType(1) / SolverDataType(nl_def_abs));
    vec_mg_def.copy(vec_mg_rhs);
    vec_mg_sol.format();
    vec_mg_cor.format();

    // multigrid loop
    Index mg_iter = 1u;
    bool mg_failed = false;
    bool mg_finished = false;
    SolverDataType mg_def_prev = SolverDataType(1.0);
    for(; mg_iter <= mg_max_iter; ++mg_iter)
    {
      // apply multigrid preconditioner
      Solver::Status status = multigrid->apply(vec_mg_cor, vec_mg_def);

      // update solution vector
      vec_mg_sol.axpy(vec_mg_cor, mg_damp);
      the_solver_level.filter_sys.filter_cor(vec_mg_sol);

      // compute new defect
      the_solver_level.matrix_sys.apply(vec_mg_def, vec_mg_sol, vec_mg_rhs, -SolverDataType(1));
      the_solver_level.filter_sys.filter_def(vec_mg_def);

      // compute defect norm (and keep in mind that we have normalized the RHS)
      const SolverDataType mg_def_v = Math::sqrt(the_solver_level.gate_velo.dot(vec_mg_def.local().template at<0>(), vec_mg_def.local().template at<0>()));
      const SolverDataType mg_def_p = Math::sqrt(the_solver_level.gate_pres.dot(vec_mg_def.local().template at<1>(), vec_mg_def.local().template at<1>()));
      const SolverDataType mg_def = Math::sqrt(mg_def_v*mg_def_v + mg_def_p*mg_def_p);

      // normalize defect
      const SolverDataType mg_def_norm = mg_def * SolverDataType(nl_def_norm);

      // print defect norm, but scale by non-linear defect to avoid confusion
      String mg_line = "Multigrid:" + stringify(mg_iter).pad_front(4) + " : "
        + stringify_fp_sci(mg_def_v * SolverDataType(nl_def_abs)) + " / "
        + stringify_fp_sci(mg_def_p * SolverDataType(nl_def_abs)) + " / "
        + stringify_fp_sci(mg_def * SolverDataType(nl_def_abs)) + " / "
        + stringify_fp_sci(mg_def_norm) + " / "
        + stringify_fp_sci(mg_def / mg_def_prev, 3);

      // multigrid failure?
      if(!Solver::status_success(status))
      {
        mg_line += " > solver breakdown!";
        mg_failed = true;
      }
      // NaNs or infinity?
      else if(!Math::isfinite(mg_def))
      {
        mg_line += " > invalid!";
        mg_failed = true;
      }
      // convergence?
      else if(mg_def_norm < mg_tol)
      {
        mg_line += " > converged!";
        mg_finished = true;
      }
      // divergence?
      else if(SolverDataType(1E+3) < mg_def_norm)
      {
        mg_line += " > diverged!";
        mg_failed = true;
      }
      // stagnation? (don't check for stagnation in first two non-linear iterations)
      else if((!navier || (nl_iter > 1u)) && (mg_iter > mg_min_iter) && (mg_def_prev < mg_stag_rate*mg_def))
      {
        mg_line += " > stagnated!";
        mg_finished = true;
      }
      // in the first two non-linear iterations, we typically stop after less max iterations
      else if(!loaded_solution && navier && (nl_iter <= 1u) && (mg_iter >= mg_boot_iter) && (mg_iter > mg_min_iter))
      {
        mg_line += " > max boot iters!";
        mg_finished = true;
      }
      // maximum number of iterations reached?
      else if(mg_iter >= mg_max_iter)
      {
        mg_line += " > max iters!";
      }

      // print current multigrid iteration line
      comm.print(mg_line);
      if(flush_mg)
        comm.print_flush();
      mg_def_prev = mg_def;

      // and stop iterating if we're finished or if we've failed
      if(mg_finished || mg_failed)
        break;

    } // multigrid loop

    // --------------------------------------------------------------------------------------------
    // PHASE 5: release linear solver and update non-linear solution
    // --------------------------------------------------------------------------------------------

    // re-scale solution vector due to previous RHS normalization
    vec_mg_sol.scale(vec_mg_sol, SolverDataType(nl_def_abs));
    vec_cor.local().convert(vec_mg_sol.local());

    watch_nonlin_solver_apply.stop();

    statistics.counts[Counts::linsol_iter] += mg_iter;

    // release linear solver if solving a non-linear system
    if(non_linear)
    {
      multigrid->done_numeric();
      multigrid_hierarchy->done_numeric();
      solvers_numeric_init = false;
    }

    // multigrid failure?
    if(mg_failed)
    {
      comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
      if(testmode)
      {
        comm.print("Test-Mode: FAILED");
        return 1;
      }
      break;
    }

    // update solution
    vec_sol.axpy(vec_cor, nl_damp);

    // compress solver expressions for statistics
    FEAT::Statistics::compress_solver_expressions();

    // next non-linear iteration
  }

  watch_nonlin_loop.stop();

  // release linear solver if not already released in non-linear case
  if(solvers_numeric_init)
  {
    multigrid->done_numeric();
    multigrid_hierarchy->done_numeric();
    solvers_numeric_init = false;
  }

  // release multigrid
  multigrid->done_symbolic();
  multigrid_hierarchy->done_symbolic();

  // save timings
  statistics.times[Times::nonlin_total] = watch_nonlin_loop.elapsed();
  statistics.times[Times::nonlin_asm_def] = watch_nonlin_def_asm.elapsed();
  statistics.times[Times::nonlin_asm_mat] = watch_nonlin_mat_asm.elapsed();
  statistics.times[Times::linsol_init] = watch_nonlin_solver_init.elapsed();
  statistics.times[Times::linsol_apply] = watch_nonlin_solver_apply.elapsed();

  // get multigrid timings
  statistics.times[Times::mg_defect] = multigrid_hierarchy->get_time_defect();
  statistics.times[Times::mg_smooth] = multigrid_hierarchy->get_time_smooth();
  statistics.times[Times::mg_coarse] = multigrid_hierarchy->get_time_coarse();
  statistics.times[Times::mg_transfer] = multigrid_hierarchy->get_time_transfer();

  // accumulate vanka timings
  for(auto& v : ama_vankas)
  {
    statistics.times[Times::vanka_init_sym] += v->time_init_symbolic();
    statistics.times[Times::vanka_init_num] += v->time_init_numeric();
    statistics.times[Times::vanka_apply] += v->time_apply();
  }

  // get multigrid timings
  for(Index i(0); i < multigrid_hierarchy->size_physical(); ++i)
  {
    statistics.mg_times_defect.front()   += statistics.mg_times_defect[i+1u]   = multigrid_hierarchy->get_time_defect(int(i));
    statistics.mg_times_smooth.front()   += statistics.mg_times_smooth[i+1u]   = multigrid_hierarchy->get_time_smooth(int(i));
    statistics.mg_times_transfer.front() += statistics.mg_times_transfer[i+1u] = multigrid_hierarchy->get_time_transfer(int(i));
    statistics.mg_times_coarse.front()   += statistics.mg_times_coarse[i+1u]   = multigrid_hierarchy->get_time_coarse(int(i));
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  StopWatch watch_sol_analysis;
  watch_sol_analysis.start();

  // compute drag/lift scaling factors for cylinder and sphere (identical in 2D)
  const SystemDataType dpf2_cylinder = SystemDataType(2) / (dim == 2 ?
    SystemDataType(0.100)*Math::sqr(v_max*(SystemDataType(2)/SystemDataType(3))) : // = 2 / (rho * U^2 * D)
    SystemDataType(0.041)*Math::sqr(v_max*(SystemDataType(4)/SystemDataType(9)))); // = 2 / (rho * U^2 * D * H)
  const SystemDataType dpf2_sphere = SystemDataType(2) / (dim == 2 ?
    SystemDataType(0.100)*Math::sqr(v_max*(SystemDataType(2)/SystemDataType(3))) : // = 2 / (rho * U^2 * D)
    SystemDataType(0.01)*Math::sqr(v_max*(SystemDataType(0.5)))); // = 2 / (rho * U^2 * D^2)

  // compute drag & lift coefficients by line integration using the trace assembler
  if(!skip_analysis)
  {
    BenchBodyForceAccumulator body_force_accum(nu);
    body_force_asm.assemble_flow_accum(
      body_force_accum,
      vec_sol.local().template at<0>(),
      vec_sol.local().template at<1>(),
      the_domain_level.space_velo,
      the_domain_level.space_pres,
      cubature_postproc);
    body_force_accum.sync(comm);

    summary.body_forces_raw_surf[0] = body_force_accum.drag_raw_new;
    summary.body_forces_raw_surf[1] = body_force_accum.lift_raw_new;
    summary.body_forces_raw_surf[2] = body_force_accum.side_raw_new;
    summary.drag_coeff_surf_cylinder = body_force_accum.drag_raw_old * dpf2_cylinder;
    summary.lift_coeff_surf_cylinder = body_force_accum.lift_raw_old * dpf2_cylinder;
    summary.drag_coeff_surf_sphere = body_force_accum.drag_raw_new * dpf2_sphere;
    summary.lift_coeff_surf_sphere = body_force_accum.lift_raw_new * dpf2_sphere;
  }

  // compute drag & lift coefficients via volume integration from unsynchronized final defect
  if(!skip_analysis)
  {
    assemble_bdforces_vol(summary.body_forces_raw_vol, vec_def_unsynced.local().first(), vec_char);
    comm.allreduce(summary.body_forces_raw_vol.v, summary.body_forces_raw_vol.v, 3u, Dist::op_sum);

    // compute bench1 coefficients from raw forces
    summary.drag_coeff_vol_cylinder = dpf2_cylinder * summary.body_forces_raw_vol[0];
    summary.lift_coeff_vol_cylinder = dpf2_cylinder * summary.body_forces_raw_vol[1];

    // compute bench7 coefficients from raw forces
    summary.drag_coeff_vol_sphere = dpf2_sphere * summary.body_forces_raw_vol[0];
    summary.lift_coeff_vol_sphere = dpf2_sphere * summary.body_forces_raw_vol[1];
  }

  // compute pressure difference around obstacle
  if(!skip_analysis)
  {
    // evaluate pressure
    auto pval_a = Assembly::DiscreteEvaluator::eval_fe_function(
      point_iv_a, vec_sol.local().template at<1>(), the_domain_level.space_pres);
    auto pval_e = Assembly::DiscreteEvaluator::eval_fe_function(
      point_iv_e, vec_sol.local().template at<1>(), the_domain_level.space_pres);

    // compute pressure mean
    const SystemDataType p_a = pval_a.mean_value_dist(comm);
    const SystemDataType p_e = pval_e.mean_value_dist(comm);
    const SystemDataType d_p = p_a - p_e;

    // compute error to reference values
    summary.pres_diff = d_p;
  }

  // compute pressure drop between inflow and outflow
  if(!skip_analysis)
  {
    // compute pressure integrals
    SystemDataType pv[2] =
    {
      flux_asm_in.assemble_discrete_integral(vec_sol.local().template at<1>(), the_domain_level.space_pres, cubature_postproc),
      flux_asm_out.assemble_discrete_integral(vec_sol.local().template at<1>(), the_domain_level.space_pres, cubature_postproc)
    };

    comm.allreduce(pv, pv, 2u, Dist::op_sum);

    // compute pressure drop average over inflow area
    summary.pres_drop_raw = pv[0] - pv[1];
    summary.pres_drop_pipe = (pv[0] - pv[1]) / (dim == 2 ? SystemDataType(0.41) : SystemDataType(0.205*0.205)*Math::pi<SystemDataType>());
    summary.pres_drop_cuboid = (pv[0] - pv[1]) / (dim == 2 ? SystemDataType(0.41) : SystemDataType(0.41*0.41));
  }

  // compute flux through region above cylinder
  if(!skip_analysis)
  {
    XFluxAccumulator flux_accum_u;
    flux_asm_u.assemble_flow_accum(
      flux_accum_u,
      vec_sol.local().template at<0>(),
      vec_sol.local().template at<1>(),
      the_domain_level.space_velo,
      the_domain_level.space_pres,
      cubature_postproc);
    flux_accum_u.sync(comm);

    summary.flux_upper = flux_accum_u.flux / SystemDataType(2);
  }

  // compute flux through region below cylinder
  if(!skip_analysis)
  {
    XFluxAccumulator flux_accum_l;
    flux_asm_l.assemble_flow_accum(
      flux_accum_l,
      vec_sol.local().template at<0>(),
      vec_sol.local().template at<1>(),
      the_domain_level.space_velo,
      the_domain_level.space_pres,
      cubature_postproc);
    flux_accum_l.sync(comm);

    summary.flux_lower = flux_accum_l.flux / SystemDataType(2);
  }

  // perform analysis of velocity field
  if(!skip_analysis)
  {
    Assembly::VelocityInfo<SystemDataType, dim> vi = Assembly::VelocityAnalyser::compute(
      vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature_postproc);
    vi.synchronize(comm);

    summary.velo_info = vi;
  }

  watch_sol_analysis.stop();
  statistics.times[Times::sol_analysis] = watch_sol_analysis.elapsed();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(args.check("save-sol") >= 0)
  {
    StopWatch watch_checkpoint;
    watch_checkpoint.start();

    String save_name;
    if(args.parse("save-sol", save_name) < 1)
    {
      save_name = String("dfg95-cc") + stringify(dim) + "d-bench1";
      save_name += "-lvl" + stringify(the_domain_level.get_level_index());
      save_name += "-n" + stringify(comm.size()) + ".sol";
    }

    comm.print("\nWriting (partitioned) solution to '" + save_name + "'");

    // serialize the partitioning
    std::vector<char> buf_pdc = domain.serialize_partitioning();

    // serialize solution vector
    //BinaryStream bs_sol;
    //vec_sol.local().write_out(LAFEM::FileMode::fm_binary, bs_sol);
    std::vector<char> buf_sol;
    vec2buf(buf_sol, vec_sol.local());

    // write to combined output file
    DistFileIO::write_combined(buf_pdc, buf_sol, save_name, comm);

    watch_checkpoint.stop();
    statistics.times[Times::checkpoint] += watch_checkpoint.elapsed();
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(args.check("vtk") >= 0)
  {
    StopWatch watch_vtk;
    watch_vtk.start();

    // build VTK name
    String vtk_name, vtk_name2;
    int npars = args.parse("vtk", vtk_name, vtk_name2);
    if(npars < 1)
    {
      vtk_name = String("dfg95-cc") + stringify(dim) + "d-bench1";
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(comm.size());
    }

    comm.print("\nWriting VTK output to '" + vtk_name + ".pvtu'");

    if(!vtk_name.empty())
    {
      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

      // project velocity and pressure
      exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());

      // project pressure
      Cubature::DynamicFactory cub("gauss-legendre:2");
      LAFEM::DenseVector<SystemDataType, Index> vtx_p;
      Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), the_domain_level.space_pres, cub);

      // write pressure
      exporter.add_cell_scalar("pressure", vtx_p.elements());

      // finally, write the VTK file
      exporter.write(vtk_name, comm);
    }

    // write refined VTK?
    if(!vtk_name2.empty())
    {
      comm.print("Writing refined VTK output to '" + vtk_name2 + ".pvtu'");

      // refine mesh
      Geometry::StandardRefinery<MeshType> refinery(the_domain_level.get_mesh());
      MeshType mesh(refinery);

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(mesh);

      // project velocity and pressure
      exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());

      // finally, write the VTK file
      exporter.write(vtk_name2, comm);
    }

    watch_vtk.stop();
    statistics.times[Times::vtk_write] = watch_vtk.elapsed();
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  watch_total_run.stop();
  statistics.times[Times::total_run] = watch_total_run.elapsed();

  {
    MemoryUsage mi;
    statistics.bytes[Bytes::peak_p] = mi.get_peak_physical();
    statistics.bytes[Bytes::peak_v] = mi.get_peak_virtual();
  }

  statistics.sync(comm);

  comm.print(String("\n") + String(100, '=') + "\n");
  if(!skip_analysis)
    comm.print(summary.format());
  comm.print(statistics.format());

  // print extended statistics if desired
  if(ext_stats)
  {
    comm.print("\n");
    comm.print(FEAT::Statistics::get_formatted_flops(statistics.times[Times::linsol_apply], comm.size()));
    comm.print(FEAT::Statistics::get_formatted_times(statistics.times[Times::linsol_apply]));
    comm.print(FEAT::Statistics::get_formatted_solver_internals("default"));
    comm.print("\n");
    comm.print(FEAT::Statistics::get_formatted_solver_tree("default").trim());
  }

  if(testmode)
    comm.print("\nTest-Mode: PASSED");

  // print elapsed runtime
  comm.print("\nTotal Wallclock Runtime: " + time_stamp.elapsed_string_now(TimeFormat::s_m) + " seconds");
  return 0;
}
