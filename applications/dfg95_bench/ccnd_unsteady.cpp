// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Unsteady CCnD solver for DFG95 Flow-Around-A-Cylinder Benchmarks
// ------------------------------------------------------------------------------------------------
// This application implements a parallel unsteady CCND solver, which is pre-configured to solve
// the infamous unsteady "flow-around-a-cylinder" benchmark problems, which are defined in
//
//     M. Schaefer and S. Turek: Benchmark Computations of Laminar Flow Around a Cylinder
//
// The system is discretized using an isoparametric Q2/P1dc finite element discretization with a
// BDF(2) time stepping scheme. The monolithic nonlinear Oseen systems are solved using an
// adaptive Newton-Multigrid solver with an additive matrix-based Vanka smoother ("AmaVanka")
// and using UMFPACK (if available) as a coarse grid solver.
// This application supports recursive partitioning as well as a checkpoint/restart system.
//
//
// ------------------------------------
// Basic Setup and Mandatory Parameters
// ------------------------------------
// This application defines default values for most of its parameters, however, three parameters
// are mandatory and always have to be specified explicitly:
//
// --bench <2|3>
// Specifies which of the 2 unsteady benchmarks is to be solved:
//   --bench 2 corresponds to the test cases 2D-2 and 3D-2Z, which is defined by a
//             steady inflow boundary condition function and is usually solved on
//             the time interval [0, 30]
//   --bench 3 corresponds to the test cases 2D-3 and 3D-3Z, which is defined by an
//             unsteady inflow boundary condition function and is solved on the
//             time interval [0, 8]
//
// --mesh <meshfiles...>
// Specifies the input mesh file(s).
//
// --level <level-max> [levels...] <level-min>
// Specifies the mesh refinement levels in the syntax according to Control::PartiDomainControl.
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
// --defo
// If specified, the deformation tensor formulation is used for the diffusion operator,
// otherwise the gradient tensor formulation is used. Defaults to gradient tensor.
//
// --v-max <vmax>
// Specifies the maximum velocity of the inflow boundary condition.
// Defaults to 1.5 in 2D and 2.25 in 3D.
//
// --upsam <ups>
// Specifies the stabilization parameter <ups> for the streamline diffusion stabilization.
// Defaults to 0, i.e. unstabilized.
//
//
// ------------------------
// Time Stepping Parameters
// ------------------------
// This section describes the parameters that control the time stepping scheme.
//
// --t-max <T>
// Specifies the maximum simulation time T. Defaults to 30.0 for bench 2 and to 8.0 for bench 3.
//
// --delta-t <dt>
// Specifies the time step size. Defaults to 1E-2.
//
// --t-expo <0|1|2>
// Specifies the time step extrapolation order, which is used to define an initial guess u_k for
// the non-linear solver in each time step k, where:
//   --t-expo 0 use constant extrapolation from previous time step, i.e.
//              u_k := u_{k-1}
//   --t-expo 1 use linear extrapolation from two previous time steps, i.e.
//              u_k := 2*u_{k-1} - u_{k-2}
//   --t-expo 2 use quadratic extrapolation from three previous time steps, i.e.
//              u_k := 3*u_{k-1} - 3*u_{k-2}) + u_{k-3}
// Defaults to 2.
//
// --steady-tol <tol>
// Specifies the velocity time-derivative tolerance for steady-state detection.
// The time-stepping loop will be terminated before reaching T-max if the L2-norm
// of the time-derivative of the velocity field drops below this tolerance, which
// indicates that the flow has reached a steady state prematurely. Defaults to 1E-3.
//
//
// -------------------------------
// Solver Configuration Parameters
// -------------------------------
// This section describes the parameters that control the non-linear Newton/Picard solver
// as well as its multigrid preconditioner and its smoother component.
//
// --picard
// If specified, the nonlinear system in each time step will be solved using a simple
// Picard iteration instead of the Newton iteration.
//
// --min-nl-iter <N>
// Specifies the minimum number of nonlinear (Newton/Picard) solver iterations per time step.
// Defaults to 1.
//
// --max-nl-iter <N>
// Specifies the maximum number of nonlinear (Newton/Picard) solver iterations per time step.
// Defaults to 10.
//
// --min-mg-iter <N>
// Specifies the minimum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 1.
//
// --max-mg-iter <N>
// Specifies the maximum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 5.
//
// --smooth-steps <N>
// Specifies the number of pre- and post-smoothing AmaVanka steps. Defaults to 8.
//
// --smooth-damp <omega>
// Specifies the damping parameter for the AmaVanka smoother. Defaults to 0.7.
//
// --no-umfpack
// If specified, the multigrid solver will use a BiCGStab-AmaVanka solver as the coarse grid solver
// instead of the UMFPACK direct solver. Note that UMFPACK is only used if it is included in the
// build id and if the coarse system is solved on a single process.
//
// --nl-tol-abs <tol>
// Specifies the absolute tolerance for the nonlinear solver. Defaults to 1E-8.
//
// --mg-tol-rel <tol>
// If given, specifies the relative tolerance for the multigrid solver.
// If not given, then the tolerance for the multigrid solver is chosen in an adaptive
// manner depending on the two previous nonlinear solver defects, which is the default case.
// The adaptive tolerance is chosen in each nonlinear iteration by analysing the nonlinear
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
//
// ---------------------------------------------
// Checkpoint & Restart Configuration Parameters
// ---------------------------------------------
// This application implements a checkpoint-&-restart system based on the CheckpointControl
// class, which allows to restart a simulation that has been aborted at a later time.
// Note that a checkpoint file contains three consecutive solution vectors u_{k}, u_{k-1} and
// u_{k-1} as well as the time-step index, the current simulation time and the time step size,
// so these three quantities are restored from a checkpoint in a restart and do not need to be
// specified explicitly.
//
// Important Note:
// It is very important that the application, which is restarted from a previously saved
// checkpoint, is launched with the same number of MPI processes, the same mesh and the same
// maximum refinement level as the application run that was used to create the checkpoints.
// However, many other parameters (e.g. viscosity, gradient/deformation tensor, linear solver
// configuration) may be changed upon restart. Note that it is also possible to change the
// time step size 'delta-t' upon restart, in which case the three saved solution vectors will
// be inter-/extrapolated to the new time stepping size.
//
//
// --checkpoint <filename> [<stepping>] [<modulus>]
// Specifies that the application should write out a checkpoint every <stepping> time steps
// and that it should overwrite previously written checkpoint files every <modulus> checkpoints.
// The <modulus> parameter is optional and it defaults to 2. Although it is possible to set the
// modulus to 1, thus ensuring that only a single checkpoint file is created and is continuously
// overwritten, it is strongly discouraged to do so, because this might lead to a corrupted
// checkpoint file if the application crashes while trying to write the checkpoint file.
// Keep in mind that I/O errors are one of the most common reasons for application crashes
// on clusters, thus setting the modulus to 1 is a recipe for disaster.
//
// Example: --checkpoint savegame 5 3
// The parameters in this example will cause the application to write a checkpoint file
// whose filename begins with 'savegame' every 5-th time step and it only will keep the last
// 3 checkpoints on disk by overwriting previous ones.
// So it the application will write checkpoints in the following time-steps:
// Time-Step  5: write to file  "savegame.0.cp"
// Time-Step 10: write to file  "savegame.1.cp"
// Time-Step 15: write to file  "savegame.2.cp"
// Time-Step 20: overwrite file "savegame.0.cp"
// Time-Step 25: overwrite file "savegame.1.cp"
// Time-Step 30: overwrite file "savegame.2.cp"
// Time-Step 35: overwrite file "savegame.0.cp"
// ...
//
//
// --restart <filename> [<curtime>] [<timestep>]
// Specifies that the application should restart from a previously saved checkpoint and
// continue from that checkpoint rather than starting from scratch. The second and third
// optional parameters <curtime> and <timestep> can be used to change the current simulation
// time as well as the current time step. If these parameters are not given, the simulation
// will continue from the simulation time and time step from when the checkpoint was written.
//
//
// ------------------------
// Miscellaneous Parameters
// ------------------------
// This section describes miscellaneous parameters that do not fit into any other section and
// which do not deserve a custom section of their own.
//
// --no-stokes
// If given, specifies that the application should start with a null solution instead of a steady-
// state Stokes solution time step 0 for the bench 2 benchmark. Has no effect for bench 3.
//
// --vtk <filename> [<stepping>]
// Specifies that the application should write a VTK visualization output file every <stepping>
// time steps. The stepping parameter is optional and defaults to 1.
//
// --ext-stats
// If given, specifies that the application should output extensive statistics at the end of the
// program run, including detailed MPI timinigs.
//
// --test-mode
// If given, specifies that the application should run in test-mode rather than its normal mode.
// In test mode, the application only perform 3 time steps and writes some additional output which
// is parsed by the test system.
//
//
// Happy vortex shedding! (^_^)
//
// \author Peter Zajac
//

#include "ccnd_common.hpp"

namespace DFG95
{
  template<typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainLevelType::SpaceVeloType SpaceVeloType;
    typedef typename DomainLevelType::SpacePresType SpacePresType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // define our system level
    typedef NavierStokesBlockedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    BenchmarkSummary<DataType, dim> summary;
    BenchmarkStats statistics;
    StopWatch watch_total_run;
    watch_total_run.start();

    /* ****************************************************************************************** */

    // solve stokes for initial solution? (bench2 only)
    const bool stokes = (args.check("no-stokes") < 0);
    const bool newton = (args.check("picard") < 0);
    const bool defo = (args.check("defo") >= 0);
    const bool adapt_tol = (args.check("mg-tol-rel") < 0);
    const bool testmode = (args.check("test-mode") >= 0);
    const bool ext_stats = (args.check("ext-stats") >= 0);

    // which benchmark?
    const int bench = parse(args, "bench", 0);
    if((bench < 2) || (bench > 3))
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--bench <2|3>' is missing!");
      FEAT::Runtime::abort();
    }

    // use UMFPACK as coarse grid solver for multigrid?
#ifdef FEAT_HAVE_UMFPACK
    const bool umf_cgs = (domain.back_layer().comm().size() == 1) && (args.check("no-umfpack") < 0);
#else
    const bool umf_cgs = false;
#endif

    // viscosity parameter
    const DataType nu = parse(args, "nu", DataType(1e-3));
    // maximum velocity, default: 2D: 1.5, 3D: 2.25
    const DataType v_max = parse(args, "v-max", DataType(dim) * DataType(0.75));
    // streamline diffusion parameter
    const DataType upsam = parse(args, "upsam", DataType(0));
    // maximum simulation time
    const DataType t_max = parse(args, "t-max", DataType(bench == 3 ? 8 : 30));
    // time step size
    const DataType delta_t = parse(args, "delta-t", DataType(0.01));
    // tolerance for steady-state detection
    const DataType steady_tol = parse(args, "steady-tol", DataType(1E-3));
    // extrapolation order: 0, 1 or 2
    const Index t_expo = parse(args, "t-expo", Index(2));
    // min. nonlinear solver iterations
    const Index min_nl_iter = parse(args, "min-nl-iter", Index(1));
    // max. nonlinear solver iterations
    const Index max_nl_iter = parse(args, "max-nl-iter", Index(10));
    // min. multigrid iterations
    const Index min_mg_iter = parse(args, "min-mg-iter", Index(1));
    // max. multigrid iterations
    const Index max_mg_iter = parse(args, "max-mg-iter", Index(5));
    // number of smoothing steps
    const Index smooth_steps = parse(args, "smooth-steps", Index(8));
    // damping parameter for smoother
    const DataType smooth_damp = parse(args, "smooth-damp", DataType(0.7));
    // rel. tolerance for linear solver
    const DataType mg_tol_rel = parse(args, "mg-tol-rel", DataType(1E-5));
    // absolute tolerance for nonlinear solver
    const DataType nl_tol_abs = parse(args, "nl-tol-abs", DataType(1E-8));

    // name of checkpoint input file to restart from
    String restart_name;
    // time of restart checkpoint; 0 => use checkpoint data
    DataType restart_time = DataType(0);
    // timestep of restart checkpoint; 0 => use checkpoint data
    Index restart_step = Index(0);
    const bool restart = (args.parse("restart", restart_name, restart_time, restart_step) > 0);

    // name of checkpoint output files
    String check_name = String("dfg95-cc") + stringify(dim) + "d-bench" + stringify(bench);
    check_name += "-lvl" + stringify(domain.front()->get_level_index());
    check_name += "-n" + stringify(comm.size());
    Index check_step = Index(args.check("checkpoint") >= 0 ? 1 : 0);
    Index check_mod = 2;
    args.parse("checkpoint", check_name, check_step, check_mod);

    // name of VTK output files
    String vtk_name = String("dfg95-cc") + stringify(dim) + "d-bench" + stringify(bench);
    vtk_name += "-lvl" + stringify(domain.front()->get_level_index());
    vtk_name += "-n" + stringify(comm.size());
    Index vtk_step = Index(args.check("vtk") >= 0 ? 1 : 0);
    args.parse("vtk", vtk_name, vtk_step);

    // dump some information to the console
    {
      static constexpr std::size_t pl = 30u;
      static constexpr char pc = '.';
      comm.print("\nProblem Parameters:");
      if(bench == 2)
        comm.print(String("Benchmark").pad_back(pl, pc) + ": bench2 (steady inflow BCs)");
      else
        comm.print(String("Benchmark").pad_back(pl, pc) + ": bench3 (unsteady inflow BCs)");
      comm.print(String("T-max").pad_back(pl, pc) + ": " + stringify(t_max));
      comm.print(String("Time Step Length").pad_back(pl, pc) + ": " + stringify(delta_t));
      comm.print(String("Steady State Derivative Tol").pad_back(pl, pc) + ": " + stringify_fp_sci(steady_tol));
      comm.print(String("Time Extrapolation Order").pad_back(pl, pc) + ": " + stringify(t_expo));
      comm.print(String("Viscosity Parameter Nu").pad_back(pl, pc) + ": " + stringify(nu));
      comm.print(String("V-Max").pad_back(pl, pc) + ": " + stringify(v_max));
      comm.print(String("Upsam").pad_back(pl, pc) + ": " + stringify(upsam));
      comm.print(String("Tensor").pad_back(pl, pc) + ": " + (defo ? "Deformation" : "Gradient"));
      comm.print(String("Start with Stokes").pad_back(pl, pc) + ": " + (stokes ? "yes" : "no"));
      comm.print(String("Nonlinear Solver").pad_back(pl, pc) + ": " + (newton ? "Newton" : "Picard"));
      comm.print(String("AmaVanka Smoother Steps").pad_back(pl, pc) + ": " + stringify(smooth_steps));
      comm.print(String("AmaVanka Smoother Damping").pad_back(pl, pc) + ": " + stringify(smooth_damp));
      comm.print(String("Nonlinear Absolute Tol").pad_back(pl, pc) + ": " + stringify_fp_sci(nl_tol_abs));
      comm.print(String("Multigrid Relative Tol").pad_back(pl, pc) + ": " + (adapt_tol ? String("adaptive") : stringify_fp_sci(mg_tol_rel)));
      comm.print(String("Min Multigrid Iterations").pad_back(pl, pc) + ": " + stringify(min_mg_iter));
      comm.print(String("Max Multigrid Iterations").pad_back(pl, pc) + ": " + stringify(max_mg_iter));
      comm.print(String("Min Nonlinear Iterations").pad_back(pl, pc) + ": " + stringify(min_nl_iter));
      comm.print(String("Max Nonlinear Iterations").pad_back(pl, pc) + ": " + stringify(max_nl_iter));
      if(umf_cgs)
        comm.print(String("Coarse Solver").pad_back(pl, pc) + ": UMFPACK");
      else
        comm.print(String("Coarse Solver").pad_back(pl, pc) + ": BiCGStab-AmaVanka");
      if(vtk_step > Index(0))
      {
        comm.print(String("VTK Output Filename").pad_back(pl, pc) + ": '" + vtk_name + "'");
        comm.print(String("VTK Output Stepping").pad_back(pl, pc) + ": " + stringify(vtk_step));
      }
      else
      {
        comm.print(String("VTK Output Filename").pad_back(pl, pc) + ": N/A");
        comm.print(String("VTK Output Stepping").pad_back(pl, pc) + ": N/A");
      }
      if(check_step > Index(0))
      {
        comm.print(String("Checkpoint Filename").pad_back(pl, pc) + ": '" + check_name + "'");
        comm.print(String("Checkpoint Stepping").pad_back(pl, pc) + ": " + stringify(check_step));
        comm.print(String("Checkpoint Overwrite").pad_back(pl, pc) + ": " + stringify(check_mod));
      }
      else
      {
        comm.print(String("Checkpoint Filename").pad_back(pl, pc) + ": N/A");
        comm.print(String("Checkpoint Stepping").pad_back(pl, pc) + ": N/A");
        comm.print(String("Checkpoint Overwrite").pad_back(pl, pc) + ": N/A");
      }
      if(restart)
      {
        comm.print(String("Restart Filename").pad_back(pl, pc) + ": '" + restart_name + "'");
        comm.print(String("Restart Time").pad_back(pl, pc) + ": " + stringify(restart_time));
        comm.print(String("Restart Timestep").pad_back(pl, pc) + ": " + stringify(restart_step));
      }
      else
      {
        comm.print(String("Restart Filename").pad_back(pl, pc) + ": N/A");
        comm.print(String("Restart Time").pad_back(pl, pc) + ": N/A");
        comm.print(String("Restart Timestep").pad_back(pl, pc) + ": N/A");
      }
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = Index(domain.size_physical());

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    // cubature for assembly
    Cubature::DynamicFactory cubature("gauss-legendre:3");

    // cubature for post-processing
    Cubature::DynamicFactory cubature_postproc("gauss-legendre:5");

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    // assemble gates, muxers and transfers
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

    // collect some finest-level statistics
    {
      auto tv = system_levels.front()->gate_velo._freqs.clone(LAFEM::CloneMode::Deep);
      tv.format(1.0);
      Index velo_dofs = Index(system_levels.front()->gate_velo.dot(tv, tv));
      Index locp_dofs = system_levels.front()->gate_pres._freqs.size();
      Index pres_dofs = Index(system_levels.front()->gate_pres.sum(DataType(locp_dofs)));
      statistics.counts[Counts::velo_dofs] = velo_dofs/dim;
      statistics.counts[Counts::pres_dofs] = pres_dofs;
      statistics.counts[Counts::total_dofs] = velo_dofs+pres_dofs;
      statistics.counts[Counts::elements] = domain.front()->get_mesh().get_num_elements();
    }

    /* ***************************************************************************************** */

    // assemble basic matrices
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
      system_levels.at(i)->assemble_pres_struct(domain.at(i)->space_pres);
      system_levels.at(i)->assemble_velocity_mass_matrix(domain.at(i)->space_velo, cubature);
      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
      system_levels.at(i)->compile_system_matrix();
    }

    /* ***************************************************************************************** */

    // our inflow BC function
    SteadyInflowFunction<dim> steady_inflow_func(v_max);

    // the names of the mesh parts on which to assemble
    std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local velocity filter
      auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow, unit_asm_noflow;

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
          // bench3 starts with homogeneous inflow
          if((name == "bnd:l") && (bench == 2))
          {
            // inflow
            unit_asm_inflow.add_mesh_part(*mesh_part);
          }
          else if(name != "bnd:r")
          {
            // outflow
            unit_asm_noflow.add_mesh_part(*mesh_part);
          }
        }
      }

      // assemble the no-flow filter first
      unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo);

      // depending on whether we solve bench2 or bench3, either assemble the
      // inflow bc or clone the noflow filter for later use
      if(bench == 2)
        unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, steady_inflow_func);
      else
        system_levels.at(i)->local_velo_filter_noflow = fil_loc_v.clone();

      // finally, compile the system filter
      system_levels.at(i)->compile_system_filter();
    }

    // finally, compile the local type-1 matrices
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->compile_local_matrix();
    }

    Statistics::toe_assembly = stamp_ass.elapsed_now();

    // accumulate sizes
    for(Index i(0); i < num_levels; ++i)
    {
      statistics.bytes[Bytes::mesh] += domain.at(i)->get_mesh_node()->bytes();
      statistics.bytes[Bytes::gate] += system_levels.at(i)->gate_sys.bytes();
      statistics.bytes[Bytes::muxer] += system_levels.at(i)->coarse_muxer_sys.bytes();
      statistics.bytes[Bytes::matrix] += system_levels.at(i)->matrix_sys.local().bytes();
      statistics.bytes[Bytes::transfer] += system_levels.at(i)->transfer_sys.bytes();
      const auto& loc_a = system_levels.at(i)->matrix_sys.local().block_a();
      const auto& loc_b = system_levels.at(i)->matrix_sys.local().block_b();
      const auto& loc_d = system_levels.at(i)->matrix_sys.local().block_d();
      statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_a.used_elements() + loc_a.rows() + Index(1));
      statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_b.used_elements() + loc_b.rows() + Index(1));
      statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_d.used_elements() + loc_d.rows() + Index(1));
      statistics.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_a.template used_elements<LAFEM::Perspective::pod>());
      statistics.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_b.template used_elements<LAFEM::Perspective::pod>());
      statistics.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_d.template used_elements<LAFEM::Perspective::pod>());
    }

    /* ***************************************************************************************** */

    // get our global solver system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalSystemTransfer GlobalSystemTransfer;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // get our global solve matrix and filter
    GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
    GlobalSystemFilter& filter = the_system_level.filter_sys;

    // create new vectors
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_def = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_cor = the_system_level.matrix_sys.create_vector_r();

    // format the vectors
    vec_sol.format();
    vec_rhs.format();

    // and filter them
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_rhs);

    // clone solution vector twice to store previous timestep solutions u_{k-1} and u_{k-2}
    GlobalSystemVector vec_sol_1 = vec_sol.clone(); // u_{k-1}
    GlobalSystemVector vec_sol_2 = vec_sol.clone(); // u_{k-2}

    {
      // count non-zeros in a and b
      statistics.counts[Counts::nnze_a] = the_system_level.matrix_sys.local().block_a().used_elements();
      statistics.counts[Counts::nnze_b] = the_system_level.matrix_sys.local().block_b().used_elements();
      statistics.counts[Counts::nnze_total] = the_system_level.matrix_sys.local().template used_elements<LAFEM::Perspective::pod>();
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    const auto* mesh_part_bnd_c = the_domain_level.get_mesh_node()->find_mesh_part("bnd:c");
    const auto* mesh_part_inner_u = the_domain_level.get_mesh_node()->find_mesh_part("inner:u");
    const auto* mesh_part_inner_l = the_domain_level.get_mesh_node()->find_mesh_part("inner:l");

    // create trace assembler for body force assembly
    Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> body_force_asm(the_domain_level.trafo);
    if(mesh_part_bnd_c != nullptr)
      body_force_asm.add_mesh_part(*mesh_part_bnd_c);
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

    // unmap pressure evaluation points p_a and p_e
    Trafo::InverseMappingData<DataType, dim> point_iv_a, point_iv_e;
    {
      typedef Trafo::InverseMapping<typename SpacePresType::TrafoType, DataType> InvMappingType;
      InvMappingType inv_mapping(the_domain_level.trafo);

      // reference pressure points
      typename InvMappingType::ImagePointType v_a, v_e;
      if(dim == 2)
      {
        // pressure evaluations points in 2D: p_a = (0.15, 0.2); p_e = (0.25, 0.2)
        v_a[0] = DataType(0.15);
        v_e[0] = DataType(0.25);
        v_a[1] = v_e[1] = DataType(0.2);
      }
      else
      {
        // pressure evaluations points in 3D: p_a = (0.45, 0.2); p_e = (0.55, 0.2)
        v_a[0] = DataType(0.45);
        v_e[0] = DataType(0.55);
        v_a[1] = v_e[1] = DataType(0.2);
        v_a[2] = v_e[2] = DataType(0.205);
      }

      // unmap points
      point_iv_a = inv_mapping.unmap_point(v_a, true);
      point_iv_e = inv_mapping.unmap_point(v_e, true);
    }

    // create characteristic function vector for circle boundary
    // this is needed for the volumetric drag/lift computation
    LAFEM::DenseVector<Mem::Main, DataType, IndexType> vec_char(vec_sol.local().template at<0>().size());
    vec_char.format(DataType(0));
    if(mesh_part_bnd_c != nullptr)
    {
      LAFEM::UnitFilter<Mem::Main, DataType, IndexType> filter_char;
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      unit_asm.add_mesh_part(*mesh_part_bnd_c);
      unit_asm.assemble(filter_char, the_domain_level.space_velo);
      filter_char.get_filter_vector().format(DataType(1));
      filter_char.filter_sol(vec_char);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a multigrid solver
    auto multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>>(domain.size_virtual());

    // array of Vanka pointers - this is only required to collect the memory usage
    // statistics of the Vankas, as we need the memory usage after factorization
    std::deque<
      std::shared_ptr<
        Solver::AmaVanka<
          typename SystemLevelType::LocalSystemMatrix,
          typename SystemLevelType::LocalSystemFilter>>> ama_vankas;

    // push levels into multigrid
    for(std::size_t i(0); i < system_levels.size(); ++i)
    {
      SystemLevelType& lvl = *system_levels.at(i);

      if((i+1) < domain.size_virtual())
      {
        auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
        ama_vankas.push_back(vanka);
        auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);
        smoother->set_min_iter(smooth_steps);
        smoother->set_max_iter(smooth_steps);
        smoother->skip_defect_calc(true); // skip defect calculation
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
#ifdef FEAT_HAVE_UMFPACK
      else if(umf_cgs)
      {
        // create UMFPACK coarse grid solver
        auto umfpack = Solver::new_generic_umfpack(lvl.local_matrix_sys);
        auto cgsolver = Solver::new_schwarz_precond(umfpack, lvl.filter_sys);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
#endif //  FEAT_HAVE_UMFPACK
      else
      {
        // create BiCGStab-AmaVanka coarse grid solver
        auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
        ama_vankas.push_back(vanka);
        auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
        auto cgsolver = Solver::new_bicgstab(lvl.matrix_sys, lvl.filter_sys, schwarz);
        cgsolver->set_max_iter(500);
        cgsolver->set_tol_rel(1e-3);
        //cgsolver->set_plot_mode(Solver::PlotMode::summary);

        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
    }

    // create our multigrid solver
    auto multigrid = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

    // create our solver
    auto solver = Solver::new_richardson(matrix, filter, 1.0, multigrid);

    solver->set_plot_name("Multigrid");

    solver->set_min_iter(min_mg_iter);
    solver->set_max_iter(max_mg_iter);
    solver->set_tol_rel(mg_tol_rel);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a few watches
    StopWatch watch_stokes_solve;
    StopWatch watch_nonlin_loop;
    StopWatch watch_nonlin_def_asm;
    StopWatch watch_nonlin_mat_asm;
    StopWatch watch_nonlin_solver_init;
    StopWatch watch_nonlin_solver_apply;
    StopWatch watch_vtk_write;
    StopWatch watch_checkpoint;
    StopWatch watch_sol_analysis;

    // initialize solver
    watch_nonlin_solver_init.start();
    multigrid_hierarchy->init_symbolic();
    solver->init_symbolic();
    watch_nonlin_solver_init.stop();

    // accumulate vanka sizes
    statistics.counts[Counts::vanka_data] = Index(ama_vankas.front()->data_size());
    statistics.bytes[Bytes::vanka] = 0ull;
    for(auto& v : ama_vankas)
      statistics.bytes[Bytes::vanka] += v->bytes();

    // solve Stokes for bench 2
    if(stokes && (bench == 2) && !restart)
    {
      comm.print("\nSolving Stokes system...");

      // assemble stokes matrices
      for(Index i(0); i < num_levels; ++i)
      {
        system_levels.at(i)->assemble_velocity_laplace_matrix(domain.at(i)->space_velo, cubature, nu, defo);
        system_levels.at(i)->compile_system_matrix();
        system_levels.at(i)->compile_local_matrix();
      }

      // initialize solver
      watch_nonlin_solver_init.start();
      multigrid_hierarchy->init_numeric();
      solver->init_numeric();
      watch_nonlin_solver_init.stop();

      // solve Stokes system
      solver->set_plot_mode(Solver::PlotMode::iter);
      watch_stokes_solve.start();
      Solver::Status stokes_status = Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);
      watch_stokes_solve.stop();
      solver->set_plot_mode(Solver::PlotMode::none);

      statistics.counts[Counts::linsol_iter] = solver->get_num_iter();

      // release solvers
      solver->done_numeric();
      multigrid_hierarchy->done_numeric();

      FEAT::Statistics::compress_solver_expressions();

      if(!Solver::status_success(stokes_status))
      {
        comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
        if(testmode)
          comm.print("Test-Mode: FAILED");
        return;
      }

      // write VTK files
      if(vtk_step > 0)
      {
        watch_vtk_write.start();
        String now_vtk_name = vtk_name + "." + stringify(0).pad_front(5, '0');
        String line = "Writing VTK file: '" + now_vtk_name + "'";
        comm.print(line);
        {
          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

          // project velocity and pressure
          exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());

          // project pressure
          Cubature::DynamicFactory cub("gauss-legendre:2");
          LAFEM::DenseVector<Mem::Main, DataType, Index> vtx_p;
          Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), the_domain_level.space_pres, cub);

          // write pressure
          exporter.add_cell_scalar("pressure", vtx_p.elements());

          // finally, write the VTK file
          exporter.write(now_vtk_name, comm);
        }
        watch_vtk_write.stop();
      }
    }

    comm.print("\nSolving nonsteady Navier-Stokes system...");

    // setup burgers assembler for matrix
    Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_mat;
    burgers_mat.deformation = defo;
    burgers_mat.nu = nu;
    burgers_mat.beta = DataType(1);
    burgers_mat.frechet_beta = DataType(newton ? 1 : 0);
    burgers_mat.sd_delta = upsam;
    burgers_mat.sd_nu = nu;
    burgers_mat.theta = DataType(1) / delta_t; // implicit Euler in first time step

    // setup burgers assembler for defect vector
    Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_def;
    burgers_def.deformation = defo;
    burgers_def.nu = nu;
    burgers_def.beta = DataType(1);
    burgers_def.theta = DataType(1) / delta_t; // implicit Euler in first time step

    // vector of all non-linear defect norms
    std::vector<DataType> nl_defs;

    // setup current time
    Index time_step(0u);
    DataType cur_time(0.0);

    // ============================================================================================
    // RESTART FROM CHECKPOINT (IF DESIRED)
    // ============================================================================================
    if(restart)
    {
      watch_checkpoint.start();

      comm.print("\n" + String(100, '=') + "\nReading restart checkpoint '" + restart_name + "'...");

      // create checkpoint control
      LAFEM::SerialConfig serial_config(false, false);
      Control::CheckpointControl check_ctrl(comm, serial_config);

      // create checkpoint info vector
      LAFEM::DenseVector<Mem::Main, DataType, IndexType> vcp_info;

      // load checkpoint
      check_ctrl.load(restart_name);

      // restore all three solution vectors from checkpoint
      check_ctrl.restore_object("u[k-2]", vec_sol_2.local());
      check_ctrl.restore_object("u[k-1]", vec_sol_1.local());
      check_ctrl.restore_object("u[k-0]", vec_sol.local());
      check_ctrl.restore_object("info", vcp_info);

      // extract old checkpoint info
      XASSERTM(vcp_info.size() == Index(3), "invalid info vector size");
      DataType cp_cur_time = vcp_info(Index(0));
      DataType cp_delta_t  = vcp_info(Index(1));
      Index cp_time_step   = Index(vcp_info(Index(2)));

      // choose the restart time based on parameters
      time_step = (restart_step <= Index(0) ? cp_time_step : restart_step - Index(1));
      cur_time = (restart_time <= DataType(0) ? cp_cur_time : restart_time - delta_t);

      // same timestep size?
      const bool same_dt = Math::abs(cp_delta_t - delta_t) < delta_t*DataType(1E-10);

      // compute delta t change factor
      // q < 1 ==> new time step size is smaller
      // q > 1 ==> new time step size is bigger
      const DataType q = delta_t / cp_delta_t;

      // print checkpoint and restart info
      comm.print("Checkpoint Delta-T.: " + stringify_fp_fix(cp_delta_t));
      comm.print("Checkpoint Time....: " + stringify_fp_fix(cp_cur_time));
      comm.print("Checkpoint Timestep: " + stringify(cp_time_step));
      comm.print("Restart Delta-T....: " + stringify_fp_fix(delta_t)
        + (same_dt ? " (same as before)" : " (changed)"));
      comm.print("Restart Time.......: " + stringify_fp_fix(cur_time)
        + (restart_time <= DataType(0) ? " (restored)" : " (override)"));
      comm.print("Restart Timestep...: " + stringify(time_step)
        + (restart_step <= Index(0) ? " (restored)" : " (override)"));


      // do we need to reconstruct the vectors due to changed time step size?
      if(!same_dt && (t_expo >= Index(2)))
      {
        // We have the following situation:
        // The three vectors {u[k-2],u[k-1],u[k-0]} that we have restored from the checkpoint
        // correspond to the old timestep size 'd' and where u[k] represents the solution at time T.
        // Since the timestep size has changed to 'b', we have to perform an interpolation and/or
        // extrapolation to obtain the three vectors {v[j-2],v[j-1],v[j]} as shown below:
        //
        //  u[k-2]              u[k-1]              u[k]
        //  |-------------------|-------------------|
        //  T-2*d               T-d                 T
        //              T-2*b         T-b           T
        //              |-------------|-------------|
        //              v[j-2]        v[j-1]        v[j]
        //
        // As a first step, we define the timestep size ratio q:= b/d and transform our old
        // tome interval [T-2*d,T] to the reference interval [-2,0]. Now, all we have to do
        // is to use the Lagrange polynomials to define the interpolation polynomial for the
        // three old vectors {u[k-2],u[k-1],u[k]} and evaluate this polynomial in the points
        // -2*q and -q to obtain the transformation coefficients.
        //
        //  u[k-2]              u[k-1]              u[k]
        //  |-------------------|-------------------|
        // -2         -2*q     -1    -q             0
        //              |-------------|-------------|
        //              v[j-2]        v[j-1]        v[j]
        //
        // The three Lagrange polynomials for the reference points {-2, -1, 0} are given as:
        // L_0(t) := (t-1)*(t-2)/2
        // L_1(t) := t*(2-t)
        // L_2(t) := t*(t-1)/2

        comm.print("Performing quadratic transformation of restored solution vectors");

        // backup old vectors
        const auto old_sol_1 = vec_sol_1.clone(); // u[k-1]
        const auto old_sol_2 = vec_sol_2.clone(); // u[k-2]

        // transform v[k-1] :=       L_0(-q)*u[k] + L_1(-q)*u[k-1] +   L_2(-q)*u[k-2]
        //                   = (q-1)*(q-2)/2*u[k] + q*(2-q)*u[k-1] + q*(q-1)/2*u[k-2]
        vec_sol_1.scale(vec_sol, (q - DataType(1))*(q - DataType(2)) / DataType(2));
        vec_sol_1.axpy(old_sol_1, vec_sol_1, q * (DataType(2) - q));
        vec_sol_1.axpy(old_sol_2, vec_sol_1, q * (q - DataType(1)) / DataType(2));

        // transform v[k-2] :=      L_0(-2q)*u[k] +  L_1(-2q)*u[k-1] +  L_2(-2q)*u[k-2]
        //                   = (q-1)*(2*q-1)*u[k] + 4*q*(1-q)*u[k-1] + q*(2*q-1)*u[k-2]
        vec_sol_2.scale(vec_sol, (q - DataType(1))*(DataType(2)*q - DataType(1)));
        vec_sol_2.axpy(old_sol_1, vec_sol_2, DataType(4) * q * (DataType(1) - q));
        vec_sol_2.axpy(old_sol_2, vec_sol_2, q * (DataType(2)*q - DataType(1)));
      }
      else if(!same_dt && (t_expo >= Index(1)))
      {
        // same approach as the previous case, but only interpolate {v[k-1],v[k]} linearly
        // from {u[k-1],u[k]} by using the first order Lagrange polynomials for {-1, 0}:
        // L_0(t) := x + 1
        // L_1(t) := -x

        comm.print("Performing linear transformation of restored solution vectors");

        // transform v[k-1] = (1-q)*u[k] + q*v[k-1]
        vec_sol_1.scale(vec_sol_1, q);
        vec_sol_1.axpy(vec_sol, vec_sol_1, DataType(1) - q);
      }
      // If t_expo is < 1, i.e. constant extrapolation, then no transformation is
      // required even if the timestep size has changed.

      comm.print(String(100, '='));
      watch_checkpoint.stop();
    }

    // ============================================================================================
    // ============================================================================================
    // START OF TIME STEPPING LOOP
    // ============================================================================================
    // ============================================================================================

    // time-step loop
    while((cur_time += delta_t) <= t_max)
    {
      // next time step
      ++time_step;

      // stop after 3 timesteps in test-mode
      if(testmode && (time_step > 3))
        break;

      comm.print(String(100, '-'));

      // okay, another time step
      String t_line = stringify(time_step).pad_front(6) + ":" + stringify_fp_fix(cur_time, 8, 12);

      // print current runtime
      comm.print(t_line + " > Runtime: " + watch_total_run.elapsed_string(TimeFormat::h_m_s_m));

      // extrapolate initial u[k] from u[k-1] and u[k-2] for this time-step
      if(time_step > Index(1))
      {
        auto& vec_sol_3 = vec_rhs; // abuse RHS vector, which is formatted below anyways

        // shift solution vectors
        vec_sol_3.copy(vec_sol_2); // u_{k-3}
        vec_sol_2.copy(vec_sol_1); // u_{k-2}
        vec_sol_1.copy(vec_sol);   // u_{k-1}

        // extrapolate to current time step
        if((t_expo >= Index(2)) && (time_step > Index(2)))
        {
          // perform quadratic extrapolation of solution to current timestep
          // u_{k} := 3*u_{k-1} - 3*u_{k-2}) + u_{k-3}
          vec_sol.axpy(vec_sol_1, vec_sol_3, DataType(3));
          vec_sol.axpy(vec_sol_2, vec_sol, -DataType(3));
        }
        else if((t_expo >= Index(1)) && (time_step > Index(1)))
        {
          // perform linear extrapolation of solution to current timestep
          // u_{k} := 2*u_{k-1} - u_{k-2}
          vec_sol.scale(vec_sol_1, DataType(2));
          vec_sol.axpy(vec_sol_2, vec_sol, -DataType(1));
        }
        // else-case is constant extrapolation, which is done already by copy
      }

      // in case of bench3, we have to reassemble the inflow BCs now
      if(bench == 3)
      {
        UnsteadyInflowFunction<dim> unsteady_inflow_func(v_max, cur_time);

        for(Index i(0); i < num_levels; ++i)
        {
          // get our local velocity filter
          auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

          // first, clone the noflow filter
          fil_loc_v.clone(system_levels.at(i)->local_velo_filter_noflow);

          // create unit-filter assembler
          Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow;

          // try to fetch the inflow mesh part node
          auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node("bnd:l");
          XASSERT(mesh_part_node != nullptr);

          auto* mesh_part = mesh_part_node->get_mesh();
          if (mesh_part != nullptr)
          {
            unit_asm_inflow.add_mesh_part(*mesh_part);
            unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, unsteady_inflow_func);
          }

          // finally, compile the system filter
          system_levels.at(i)->compile_system_filter();
        }
      }

      // apply the solution filter now
      filter.filter_sol(vec_sol);

      // compute right hand side for this time-step
      // note: we keep vec_rhs as a type-0 vector, because it will be
      // converted to type-1 only as a defect vector in the nonlinear loop below
      vec_rhs.format();
      if(time_step == Index(1))
      {
        // first time step ==> implicit Euler
        // f_k := 1/dt * M * u_{k-1}
        the_system_level.local_velo_mass_matrix.apply(
          vec_rhs.local().template at<0>(),
          vec_sol_1.local().template at<0>());
        vec_rhs.scale(vec_rhs, DataType(1) / delta_t);
      }
      else
      {
        // we're beyond the first time step ==> BDF(2)
        // First, adjust the mass matrix parameter in the burgers assemblers
        burgers_mat.theta = DataType(1.5) / delta_t;
        burgers_def.theta = DataType(1.5) / delta_t;

        // f_k := 3/(2*dt) * (4/3 * M * u_{k-1} - 1/3 M * u_{k-2}
        //      = -1/(2*dt) * M * (u_{k-2} - 4*u_{k-1})
        vec_def.axpy(vec_sol_1, vec_sol_2, -DataType(4));
        the_system_level.local_velo_mass_matrix.apply(
          vec_rhs.local().template at<0>(),
          vec_def.local().template at<0>());
        vec_rhs.scale(vec_rhs, -DataType(0.5) / delta_t);
      }

      // ------------------------------------------------------------------------------------------
      // START OF NONLINEAR SOLVER
      // ------------------------------------------------------------------------------------------

      watch_nonlin_loop.start();
      nl_defs.clear();

      // nonlinear solver successful?
      bool nl_success = false;

      // nonlinear solver loop
      for(Index nl_step(0); nl_step <= max_nl_iter; ++nl_step)
      {
        // yet another nonlinear iteration
        ++statistics.counts[Counts::nonlin_iter];

        // assemble nonlinear defect vector
        watch_nonlin_def_asm.start();
        vec_def.copy(vec_rhs);
        // assemble burgers operator defect
        burgers_def.assemble_vector(vec_def.local().template at<0>(), vec_sol.local().template at<0>(),
          vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature, -1.0);
        // compute remainder of defect vector
        the_system_level.matrix_sys.local().block_b().apply(
          vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
        the_system_level.matrix_sys.local().block_d().apply(
          vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);
        // sync and filter
        vec_def.sync_0();
        filter.filter_def(vec_def);
        watch_nonlin_def_asm.stop();

        // compute defect norm
        const DataType def_prev = (nl_defs.empty() ? DataType(1) : nl_defs.back());
        const DataType def_nl = vec_def.norm2();
        const DataType def_improve = def_nl / def_prev;
        nl_defs.push_back(def_nl);

        // pre-format output line
        String line = t_line + " | ";
        line += stringify(nl_step).pad_front(2) + ": ";
        line += stringify_fp_sci(def_nl, 6) + " / ";
        line += stringify_fp_sci(def_nl/nl_defs.front()) + " / ";
        line += stringify_fp_sci(def_nl/def_prev, 3);

        // diverged, converged, maximum iterations reached?
        if(def_nl > nl_defs.front() * DataType(1E+3))
        {
          comm.print(line);
          comm.print("\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
          if(testmode)
            comm.print("Test-Mode: FAILED");
          return;
        }
        else if(nl_step < min_nl_iter)
        {
          // nothing to do here; this else-case exists just to ensure
          // that none of the following cases fires and breaks the loop
        }
        else if(def_nl < nl_tol_abs)
        {
          line += " > converged!";
          comm.print(line);
          nl_success = true;
          break;
        }
        else if(nl_step >= max_nl_iter)
        {
          line += " > maximum iterations reached!";
          comm.print(line);
          //nl_success = true;
          break;
        }
        else if((nl_step >= 3) && (DataType(0.95)*def_prev < def_nl))
        {
          line += " > stagnated!";
          comm.print(line);
          nl_success = true;
          break;
        }

        // assemble burgers matrices on all levels
        watch_nonlin_mat_asm.start();
        {
          // get a clone of the global velocity vector
          typename SystemLevelType::GlobalVeloVector vec_conv(
            &the_system_level.gate_velo, vec_sol.local().template at<0>().clone());

          // set velocity norm for streamline diffusion (if enabled)
          if(Math::abs(upsam) > DataType(0))
          {
            // set norm by convection vector; we can use this for all levels
            burgers_mat.set_sd_v_norm(vec_conv);
          }

          // loop over all system levels
          for(std::size_t i(0); i < system_levels.size(); ++i)
          {
            // assemble our system matrix
            auto& loc_mat_a = system_levels.at(i)->matrix_sys.local().block_a();
            loc_mat_a.format();
            burgers_mat.assemble_matrix(loc_mat_a, vec_conv.local(), domain.at(i)->space_velo, cubature);
            system_levels.at(i)->compile_local_matrix();

            // restrict our convection vector
            if((i+1) >= domain.size_virtual())
              break;

            // does this process have another system level?
            if((i+1) < system_levels.size())
            {
              // create a coarse mesh velocity vector
              auto vec_crs = system_levels.at(i+1)->matrix_a.create_vector_l();

              // truncate fine mesh velocity vector
              system_levels.at(i)->transfer_velo.trunc(vec_conv, vec_crs);

              // the coarse vector is our next convection vector
              vec_conv = std::move(vec_crs);
            }
            else
            {
              // this process is a child, so send truncation to parent
              system_levels.at(i)->transfer_velo.trunc_send(vec_conv);
            }
          }
        }
        watch_nonlin_mat_asm.stop();

        // initialize linear solver
        watch_nonlin_solver_init.start();
        multigrid_hierarchy->init_numeric();
        solver->init_numeric();
        watch_nonlin_solver_init.stop();

        // specify adaptive tolerance?
        if(adapt_tol && (nl_step > Index(0)))
        {
          if(newton)
          {
            // We're using Newton as the nonlinear solver, which optimally should
            // result in quadratic convergence, i.e. let def_{j} and def_{j-1}
            // denote the two previous nonlinear defect norms, then the next
            // defect norm def_{j+1} should fulfill
            //
            //     (def_{j+1} / def_{j}) \approx (def_{j} / def_{j+1})^2
            //
            // If the multiply the above equation by def_{j}, we can therefore
            // estimate the next def_{j+1} based on the two previous defects:
            //
            //      def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})^2
            //
            // We now multiply the approximation by 0.1, which gives us an absolute
            // tolerance for the multigrid solver for this nonlinear iteration.
            // (Note that def_improve := def_{j} / def_{j+1})
            DataType abs_tol = def_nl * def_improve * def_improve * DataType(0.1);
            // We furthermore limit this absolute tolerance to ensure that we do not
            // overshoot the mark by overoptimistic quadratic convergence expectations.
            solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
            // Also make sure that we gain at least 2 digits.
            solver->set_tol_rel(1E-2);
          }
          else
          {
            // In the case if Picard itertion, we only expect linear convergence,
            // which (in analogy to Newton) leads us to the following estimate:
            //
            //      def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})
            //
            DataType abs_tol = def_nl * def_improve * DataType(0.1);
            solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
            solver->set_tol_rel(1E-2);
          }
        }
        else
        {
          // disable absolute tolerance
          solver->set_tol_abs(1E+10);
          solver->set_tol_rel(mg_tol_rel);
        }

        // solve linear system
        watch_nonlin_solver_apply.start();
        Solver::Status status = solver->apply(vec_cor, vec_def);
        watch_nonlin_solver_apply.stop();

        statistics.counts[Counts::linsol_iter] += solver->get_num_iter();

        line += String(" | ") + stringify(solver->get_num_iter()).pad_front(3) + ": "
          + stringify_fp_sci(solver->get_def_final(), 4) + " / "
          + stringify_fp_sci(solver->get_def_final() / solver->get_def_initial(), 4);
        if(adapt_tol && (nl_step > Index(0)))
          line += String(" [") + stringify_fp_sci(solver->get_tol_abs(), 4) + "]";
        comm.print(line);

        // release linear solver
        solver->done_numeric();
        multigrid_hierarchy->done_numeric();

        if(!Solver::status_success(status))
        {
          comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
          if(testmode)
            comm.print("Test-Mode: FAILED");
          return;
        }

        // update solution
        vec_sol.axpy(vec_cor, vec_sol, DataType(1));

        FEAT::Statistics::compress_solver_expressions();
        // next non-linear iteration
      }

      watch_nonlin_loop.stop();

      // non-linear solver failure?
      if(!nl_success)
        break;

      // ------------------------------------------------------------------------------------------
      // END OF NONLINEAR SOLVER
      // ------------------------------------------------------------------------------------------

      // write checkpoint if desired
      if((check_step > Index(0)) && (time_step % check_step == 0))
      {
        // it's checkpointing time!
        watch_checkpoint.start();

        // generate checkpoint name
        String check_fn = check_name + "." + stringify(((time_step-1) / check_step) % check_mod) + ".cp";
        comm.print(t_line + " > Writing Checkpoint '" + check_fn + "'");

        // create checkpoint control
        LAFEM::SerialConfig serial_config(false, false);
        Control::CheckpointControl check_ctrl(comm, serial_config);

        // add all three known solution vectors to checkpoint
        check_ctrl.add_object("u[k-2]", vec_sol_2.local());
        check_ctrl.add_object("u[k-1]", vec_sol_1.local());
        check_ctrl.add_object("u[k-0]", vec_sol.local());

        // set up info vector and write that one, too
        LAFEM::DenseVector<Mem::Main, DataType, IndexType> vcp_info(Index(3));
        vcp_info(Index(0), cur_time);
        vcp_info(Index(1), delta_t);
        vcp_info(Index(2), DataType(time_step));
        check_ctrl.add_object("info", vcp_info);

        // save checkpoint
        check_ctrl.save(check_fn);

        watch_checkpoint.stop();
      }

      // perform post-processing
      watch_sol_analysis.start();
      // compute drag & lift coefficients by line integration using the trace assembler
      {
        BenchBodyForceAccumulator<DataType> body_force_accum(defo, nu, v_max);
        body_force_asm.assemble_flow_accum(
          body_force_accum,
          vec_sol.local().template at<0>(),
          vec_sol.local().template at<1>(),
          the_domain_level.space_velo,
          the_domain_level.space_pres,
          cubature_postproc);
        body_force_accum.sync(comm);

        summary.drag_coeff_line = body_force_accum.drag;
        summary.lift_coeff_line = body_force_accum.lift;
      }

      // compute drag & lift coefficients via volume integration
      {
        const auto& vec_sol_v = vec_sol.local().template at<0>();
        const auto& vec_sol_p = vec_sol.local().template at<1>();

        Tiny::Matrix<DataType, 2, dim> bdf;
        assemble_bdforces_vol<DataType, dim>(bdf, vec_sol_v, vec_sol_p, vec_char,
          the_domain_level.space_velo, the_domain_level.space_pres, cubature_postproc);

        const DataType dpf1 = nu;
        const DataType dpf2 = DataType(2) / (dim == 2 ?
          DataType(0.100)*Math::sqr(v_max*(DataType(2)/DataType(3))) : // = 2 / (rho * U^2 * D)
          DataType(0.041)*Math::sqr(v_max*(DataType(4)/DataType(9)))); // = 2 / (rho * U^2 * D * H)

        summary.drag_coeff_vol = the_system_level.gate_sys.sum(dpf2 * (dpf1 * bdf[0][0] + bdf[1][0]));
        summary.lift_coeff_vol = the_system_level.gate_sys.sum(dpf2 * (dpf1 * bdf[0][1] + bdf[1][1]));
      }

      // compute pressure difference
      {
        // evaluate pressure
        auto pval_a = Assembly::DiscreteEvaluator::eval_fe_function(
          point_iv_a, vec_sol.local().template at<1>(), the_domain_level.space_pres);
        auto pval_e = Assembly::DiscreteEvaluator::eval_fe_function(
          point_iv_e, vec_sol.local().template at<1>(), the_domain_level.space_pres);

        // compute pressure mean
        const DataType p_a = pval_a.mean_value_dist(comm);
        const DataType p_e = pval_e.mean_value_dist(comm);
        const DataType d_p = p_a - p_e;

        // compute error to reference values
        summary.pres_diff = d_p;
      }

      // compute flow through upper region
      {
        XFluxAccumulator<DataType> flux_accum_u;
        flux_asm_u.assemble_flow_accum(
          flux_accum_u,
          vec_sol.local().template at<0>(),
          vec_sol.local().template at<1>(),
          the_domain_level.space_velo,
          the_domain_level.space_pres,
          cubature_postproc);
        flux_accum_u.sync(comm);

        summary.flux_upper = flux_accum_u.flux / DataType(2);
      }

      // compute flow through lower region
      {
        XFluxAccumulator<DataType> flux_accum_l;
        flux_asm_l.assemble_flow_accum(
          flux_accum_l,
          vec_sol.local().template at<0>(),
          vec_sol.local().template at<1>(),
          the_domain_level.space_velo,
          the_domain_level.space_pres,
          cubature_postproc);
        flux_accum_l.sync(comm);

        summary.flux_lower = flux_accum_l.flux / DataType(2);
      }

      // perform analysis of velocity field
      {
        Assembly::VelocityInfo<DataType, dim> vi = Assembly::VelocityAnalyser::compute(
          vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature_postproc);
        vi.synchronize(comm);

        summary.velo_info = vi;
      }

      // print summary line
      comm.print(summary.format_compact(t_line + " > "));

      // compute time derivative vector
      vec_def.axpy(vec_sol_1, vec_sol, -DataType(1));
      vec_def.scale(vec_def, DataType(1) / delta_t);

      // compute L2-norm of velocity time-derivative by utilising the velocity mass matrix
      // |dt v|_L2 = sqrt( (dt v)^T * M * (dt v) )
      the_system_level.local_velo_mass_matrix.apply(
        vec_cor.local().template at<0>(), vec_def.local().template at<0>());
      const DataType norm_der_v = Math::sqrt(the_system_level.gate_sys.sum(
        vec_cor.local().template at<0>().dot(vec_def.local().template at<0>())));

      // print time derivative
      {
        String line = t_line + " > velocity time derivative norm: ";
        line += stringify_fp_sci(norm_der_v);
        comm.print(line);
      }
      watch_sol_analysis.stop();

      // write VTK files
      if((vtk_step > 0) && (time_step % vtk_step == 0))
      {
        watch_vtk_write.start();
        String now_vtk_name = vtk_name + "." + stringify(time_step).pad_front(5, '0');
        String line = t_line +  " > Writing VTK file: '" + now_vtk_name + "'";
        comm.print(line);
        {
          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

          // project velocity and pressure
          exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());
          exporter.add_vertex_vector("velo_der", vec_def.local().template at<0>());

          // project pressure
          Cubature::DynamicFactory cub("gauss-legendre:2");
          LAFEM::DenseVector<Mem::Main, DataType, Index> vtx_p, der_p;
          Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), the_domain_level.space_pres, cub);
          Assembly::DiscreteCellProjector::project(der_p, vec_def.local().template at<1>(), the_domain_level.space_pres, cub);

          // write pressure
          exporter.add_cell_scalar("pressure", vtx_p.elements());
          exporter.add_cell_scalar("pres_der", der_p.elements());

          // finally, write the VTK file
          exporter.write(now_vtk_name, comm);
        }
        watch_vtk_write.stop();
      }

      // steady state reached?
      if(norm_der_v <= steady_tol)
      {
        comm.print("\nSteady State reached!");
        break;
      }

      // next time step
    }

    // ============================================================================================
    // ============================================================================================
    // END OF TIME STEPPING LOOP
    // ============================================================================================
    // ============================================================================================

    // release solver
    solver->done_symbolic();
    multigrid_hierarchy->done_symbolic();

    // save timings
    statistics.times[Times::stokes_solve] = watch_stokes_solve.elapsed();
    statistics.times[Times::nonlin_total] = watch_nonlin_loop.elapsed();
    statistics.times[Times::nonlin_asm_def] = watch_nonlin_def_asm.elapsed();
    statistics.times[Times::nonlin_asm_mat] = watch_nonlin_mat_asm.elapsed();
    statistics.times[Times::linsol_init] = watch_nonlin_solver_init.elapsed();
    statistics.times[Times::linsol_apply] = watch_nonlin_solver_apply.elapsed();
    statistics.times[Times::checkpoint] = watch_checkpoint.elapsed();
    statistics.times[Times::vtk_write] = watch_vtk_write.elapsed();
    statistics.times[Times::sol_analysis] = watch_sol_analysis.elapsed();

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

    watch_total_run.stop();
    statistics.times[Times::total_run] = watch_total_run.elapsed();

    {
      MemoryUsage mi;
      statistics.bytes[Bytes::peak_p] = mi.get_peak_physical();
      statistics.bytes[Bytes::peak_v] = mi.get_peak_virtual();
    }

    statistics.sync(comm);

    comm.print(String("\n") + String(80, '=') + "\n");
    comm.print(statistics.format());

    // print multigrid timings
    if(comm.rank() == 0)
    {
      comm.print("Multigrid Timings:");
      comm.print("              Defect /   Smoother /   Transfer /     Coarse");
      comm.print("Overall : " +
          stringify_fp_fix(multigrid_hierarchy->get_time_defect(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_smooth(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_transfer(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_coarse(), 3, 10));
      for(int i(0); i < int(multigrid_hierarchy->size_physical()); ++i)
      {
        comm.print("Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(multigrid_hierarchy->get_time_defect(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_smooth(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_transfer(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_coarse(i), 3, 10));
      }
    }

    // print extended statistics if desired
    if(ext_stats)
    {
      comm.print("\n");
      comm.print(FEAT::Statistics::get_formatted_flops(statistics.times[Times::stokes_solve] + statistics.times[Times::linsol_apply], comm.size()));
      comm.print(FEAT::Statistics::get_formatted_times(statistics.times[Times::stokes_solve] + statistics.times[Times::linsol_apply]));
      comm.print(FEAT::Statistics::get_formatted_solver_internals("default"));
      comm.print("\n");
      comm.print(FEAT::Statistics::get_formatted_solver_tree("default").trim());
    }

    if(testmode)
      comm.print("\nTest-Mode: PASSED");
  }

  template<int dim_>
  void run_dim(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
    // define our mesh type
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // create a time-stamp
    TimeStamp time_stamp;

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(mesh_reader);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());
    comm.print("Transformation: Standard");

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }

  template<int dim_>
  void run_dim_iso(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
    // define our mesh type
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, DataType> MeshType;
    typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // create a time-stamp
    TimeStamp time_stamp;

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(mesh_reader);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());
    comm.print("Transformation: Isoparametric:2");

    // get circle chart
    auto chart = domain.get_atlas().find_mesh_chart("circle");

    // loop over all physical levels and add meshpart charts
    for(Index ilvl(0); ilvl < domain.size_physical(); ++ilvl)
    {
      auto mesh_node = domain.at(ilvl)->get_mesh_node();
      const auto mpart = mesh_node->find_mesh_part("bnd:c");
      if((chart != nullptr) && (mpart != nullptr))
        domain.at(ilvl)->trafo.add_meshpart_chart(*mpart, *chart);
    }

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }

  void main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));
    comm.print("Floating Point Type: " + String(fp_typename) + " precision");

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("level");
    args.support("t-max");
    args.support("delta-t");
    args.support("t-expo");
    args.support("vtk");
    args.support("mesh");
    args.support("nu");
    args.support("upsam");
    args.support("bench");
    args.support("min-nl-iter");
    args.support("max-nl-iter");
    args.support("min-mg-iter");
    args.support("max-mg-iter");
    args.support("smooth-steps");
    args.support("smooth-damp");
    args.support("mg-tol-rel");
    args.support("nl-tol-abs");
    //args.support("no-adapt");
    args.support("steady-tol");
    args.support("picard");
    args.support("v-max");
    args.support("defo");
    args.support("no-stokes");
    args.support("test-mode");
    args.support("ext-stats");
    args.support("no-umfpack");
    //args.support("isoparam");
    args.support("checkpoint");
    args.support("restart");

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

    if(args.check("bench") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--bench <2|3>' is missing!");
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

    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run 2D or 3D ?
    //if(args.check("isoparam") >= 0)
    {
      if(mesh_type == "conformal:hypercube:2:2")
        run_dim_iso<2>(args, comm, mesh_reader);
      else if(mesh_type == "conformal:hypercube:3:3")
        run_dim_iso<3>(args, comm, mesh_reader);
      else
      {
        comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
        FEAT::Runtime::abort();
      }
    }
    /*else
    {
      if(mesh_type == "conformal:hypercube:2:2")
        run_dim<2>(args, comm, mesh_reader);
      //else if(mesh_type == "conformal:hypercube:3:3")
        //run_dim<3>(args, comm, mesh_reader);
      else
      {
        comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
        FEAT::Runtime::abort();
      }
    }*/
  }
} // namespace DFG95

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialize(argc, argv);
  try
  {
    DFG95::main(argc, argv);
  }
  catch(std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return FEAT::Runtime::finalize();
}
