// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
// Specifies which of the 3 unsteady benchmarks is to be solved:
//   --bench 2 corresponds to the test cases 2D-2 and 3D-2Z, which is defined by a steady inflow
//             boundary condition function and is usually solved on the time interval [0, 30]
//   --bench 3 corresponds to the test cases 2D-3 and 3D-3Z, which is defined by an unsteady inflow
//             boundary condition function and is solved on the time interval [0, 8]
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
// -----------------------------------------------------
// Initial Solution Read-In and Final Solution Write-Out
// -----------------------------------------------------
//
// --save-joined-sol <filename>
// Specifies that the application should write the final solution into a joined binary output
// file by utilizing the base splitter. The output file can be loaded by the --load-joined-sol
// option if the input mesh and refinement level are identical, however, the process count
// and/or partitioning may differ. This feature should only be used with at most one parallel
// domain layer and moderate process counts.
//
// --load-joined-sol <filename>
// Specifies that the application should read in the initial joined solution guess from a
// single binary output file, which was written by a --save-joined-sol from a previous run.
// The second option argument <scale> is given, then the loaded solution is scaled by that factor,
// which can be used if the loaded solution was computed with a different inflow velocity.
// If specified, the solving of the Stokes system to obtain an initial guess is skipped.
//
// --save-joined-velo <base-filename>
// Specifies that the application should write the velocities in each time step into a joined
// binary output file by utilizing the base splitter.
// The output file can be loaded by the --load-joined-sol
// option if the input mesh and refinement level are identical, however, the process count
// and/or partitioning may differ. This feature should only be used with at most one parallel
// domain layer and moderate process counts.
//
// --save-joined-pres <base-filename>
// Specifies that the application should write the velocities in each time step into a joined
// binary output file by utilizing the base splitter.
// The output file can be loaded by the --load-joined-sol
// option if the input mesh and refinement level are identical, however, the process count
// and/or partitioning may differ. This feature should only be used with at most one parallel
// domain layer and moderate process counts.
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
// --vtk <filename> [<stepping>] [<refined-filename>]
// Specifies that the application should write a VTK visualization output file every <stepping>
// time steps. The stepping parameter is optional and defaults to 1. The third optional parameter
// specifies the filename for a VTK output file on a once refined mesh, if given.
// Note: it is possible to output only the refined VTKs by passing a whitespace string as the
// first filename, e.g.: --vtk " " 1 myvtk
//
// --ext-stats
// If given, specifies that the application should output extensive statistics at the end of the
// program run, including detailed MPI timings.
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

// include common CCND headers
#include "ccnd_common.hpp"
#include "ccnd_common_dfg95.hpp"

// open the namespace and define the DomainLevel and SystemLevel classes
namespace CCND
{
  /// our domain level
  typedef CCND::DomainLevelBase DomainLevel;

  /// define our system level
  typedef CCND::SystemLevelBase SystemLevel;

} //  namespace CCND

  // include unsteady application base header
#include "ccnd_unsteady_appbase.hpp"

namespace CCND
{
  /// our actual Application class
  class Application :
    public UnsteadyAppBase
  {
  public:
    /// our base class
    typedef UnsteadyAppBase BaseClass;

    /// which benchmark to solve? 2 or 3?
    int bench = 2;

    /// maximum inflow velocity: 1.5 for 2D bench 2; 2.25 for 3D bench 2
    DataType v_max = DataType(dim) * DataType(0.75);

    /// DFG95 benchmark summary
    DFG95::BenchmarkAnalysis bench_analysis;

    explicit Application(const Dist::Comm& comm_, SimpleArgParser& args_) :
      BaseClass(comm_, args_)
    {
      args.support("bench");
      args.support("v-max");
      args.support("fbm");

      // this application has homogeneous RHS
      homogeneous_unsteady_rhs = true;
    }

    virtual void parse_args() override
    {
      args.parse("bench", bench);
      if((bench != 2) && (bench != 3))
      {
        comm.print("ERROR: invalid benchmark specified: " + stringify(bench));
        Runtime::abort();
      }

      // in case of bench 3, set the time interval to 8 and delta-t to 1/400
      if(bench == 3)
      {
        time_max = DataType(8);
        max_time_steps = 3200;
      }

      // enable FBM if desired
      enable_fbm = (args.check("fbm") >= 0);

      args.parse("v-max", v_max);

      // set default filename
      default_filename = String("cc") + stringify(dim) + "d-unsteady-dfg95-bench" + stringify(bench);
      if(enable_fbm)
        default_filename += "-fbm";

      // pre-define default parameters
      nu = DataType(0.001);

      // parse remaining arguments
      BaseClass::parse_args();
    }

    virtual void print_problem() override
    {
      BaseClass::print_problem();
      comm.print("\nDFG95 Benchmark Parameters:");
      comm.print(String("Benchmark").pad_back(pad_len, pad_char) + ": bench " + stringify(bench));
      comm.print(String("V-Max").pad_back(pad_len, pad_char) + ": " + stringify(v_max));
    }

    virtual void create_domain() override
    {
      BaseClass::create_domain();

      // create FBM meshpart?
      if(enable_fbm)
      {
        // obstacle midpoint: 2D: (0.2,0.2)
        Tiny::Vector<DataType, dim> mp(DataType(0.2));
        if constexpr (dim == 3)
          mp[0] = DataType(0.5);

        Geometry::SphereHitTestFunction<DataType, dim> fbm_hit_func(mp, 0.05);

        // create meshparts on all levels
        for(std::size_t i(0); i < domain.size_physical(); ++i)
        {
          auto* mesh_node = domain.at(i)->get_mesh_node();
          Geometry::HitTestFactory<decltype(fbm_hit_func), MeshType> hit_factory(fbm_hit_func, *mesh_node->get_mesh());
          mesh_node->add_mesh_part("fbm", hit_factory.make_unique());

          // create FBM assembler
          domain.at(i)->create_fbm_assembler(domain.at(i).layer().comm(), "fbm");
        }
      }

      // add all isoparametric mesh-parts
      BaseClass::add_all_isoparam_parts();
    }

    virtual void create_benchmark_analysis()
    {
      // let the analysis get its mesh parts and charts
      bench_analysis.create(domain.get_atlas(), *domain.front()->get_mesh_node(), domain.front()->space_velo);
    }

    virtual void assemble_filters()
    {
      // our inflow BC function
      DFG95::SteadyInflowFunction<dim> inflow_func(v_max);

      // the names of the mesh parts on which to assemble
      std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

      for(std::size_t i(0); i < domain.size_physical(); ++i)
      {
        // get our local velocity filters
        auto& filter_v_noflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("noflow");
        auto& filter_v_inflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("inflow");

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
            if(name == "bnd:l")
            {
              // inflow
              unit_asm_inflow.add_mesh_part(*mesh_part);
            }
            else if((name != "bnd:r") && (name != "bnd:out"))
            {
              // outflow
              unit_asm_noflow.add_mesh_part(*mesh_part);
            }
          }
        }

        // assemble the filters
        unit_asm_inflow.assemble(filter_v_inflow, domain.at(i)->space_velo, inflow_func);
        unit_asm_noflow.assemble(filter_v_noflow, domain.at(i)->space_velo);

        // assemble FBM filters?
        if(enable_fbm)
        {
          auto& fbm_asm = *domain.at(i)->fbm_assembler;
          auto& filter_fbm_p = system.at(i)->get_local_pres_unit_filter();
          auto& filter_fbm_v = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("fbm");
          auto& filter_fbm_int_v = system.at(i)->filter_interface_fbm;

          // assemble velocity unit filter
          fbm_asm.assemble_inside_filter(filter_fbm_v, domain.at(i)->space_velo);
          fbm_asm.assemble_inside_filter(filter_fbm_p, domain.at(i)->space_pres);

          // assemble interface filter
          fbm_asm.assemble_interface_filter(filter_fbm_int_v, domain.at(i)->space_velo, system.at(i)->matrix_a, system.at(i)->velo_mass_matrix);

          // assemble mask vectors on finest level
          if(i == 0u)
          {
            auto& mask_v = system.at(i)->fbm_mask_velo;
            mask_v.reserve(domain.at(i)->space_velo.get_num_dofs());
            for(int d(0); d <= dim; ++d)
            {
              for(auto k : fbm_asm.get_fbm_mask_vector(d))
                mask_v.push_back(k);
            }
            system.at(i)->fbm_mask_pres = fbm_asm.get_fbm_mask_vector(dim);
          }
        }

        // compile system filter
        system.at(i)->compile_system_filter();
      }
    }

    virtual void update_filters()
    {
      // reassemble inflow BCs for bench 3
      if(bench == 3)
      {
        DFG95::UnsteadyInflowFunction<dim> unsteady_inflow_func(v_max, cur_time);

        for(std::size_t i(0); i < domain.size_physical(); ++i)
        {
          // get our local velocity filter
          auto& filter_v_inflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("inflow");

          // create unit-filter assembler
          Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow;

          // try to fetch the inflow mesh part node
          auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node("bnd:l");
          XASSERT(mesh_part_node != nullptr);

          auto* mesh_part = mesh_part_node->get_mesh();
          if (mesh_part != nullptr)
          {
            unit_asm_inflow.add_mesh_part(*mesh_part);
            unit_asm_inflow.assemble(filter_v_inflow, domain.at(i)->space_velo, unsteady_inflow_func);
          }

          // finally, compile the system filter
          system.at(i)->compile_system_filter();
        }
      }
    }

    virtual void create_vectors() override
    {
      BaseClass::create_vectors();

      // format the vectors
      vec_sol.format();
      vec_rhs.format();

      // apply filters
      system.front()->filter_sys.filter_sol(vec_sol);
      system.front()->filter_sys.filter_rhs(vec_rhs);
    }

    virtual void perform_solution_analysis()
    {
      // only perform the analysis in main steps
      if((main_step == 0u) || (cur_step  % main_step > 0u))
        return;

      watch_sol_analysis.start();

      // compute drag & lift by line integration
      bench_analysis.compute_body_forces_line(comm, vec_sol.local().at<0>(), vec_sol.local().at<1>(),
        domain.front()->space_velo, domain.front()->space_pres, nu, v_max);

      // compute drag & lift coefficients via volume integration from unsynchronized final defect
      bench_analysis.compute_body_forces_vol(comm, vec_def_unsynced.local().first(), nu, v_max, bench);

      // compute pressure difference
      bench_analysis.compute_pressure_difference(comm, vec_sol.local().at<1>(), domain.front()->space_pres);

      // compute fluxes
      bench_analysis.compute_fluxes(comm, vec_sol.local().first(), domain.front()->space_velo);

      // perform analysis of velocity field
      bench_analysis.compute_velocity_info(comm, vec_sol.local().first(), domain.front()->space_velo);

      // print the analysis
      comm.print(bench_analysis.format_compact(line_prefix() + " > "));

      watch_sol_analysis.stop();
      stats.times[Times::sol_analysis] = watch_sol_analysis.elapsed();
    }
  }; // class Application


  void main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));
    comm.print("Floating Point Type: " + String(fp_typename) + " precision");

    // create arg parser
    SimpleArgParser args(argc, argv);

    // create application
    Application app(comm, args);

    // check for unsupported arguments
    app.check_args();

    // parse arguments
    app.parse_args();

    // read mesh and create domain
    app.create_domain();

    // print problem parameters
    app.print_problem();

    // create benchmark analysis stuff
    app.create_benchmark_analysis();

    // create system
    app.create_system(true);

    // assemble filters
    app.assemble_filters();

    // compile systems
    app.compile_system_matrices();

    // create vectors
    app.create_vectors();

    // collect system statistics
    app.collect_system_statistics();

    // create solver
    app.create_solver();

    // initialize solver
    app.init_solver_symbolic();

    // todo: restart from checkpoint?
    if(app.load_checkpoint())
    {
      // ok, application restarted from checkpoint
    }
    else if(app.load_initial_solution())
    {
      // ok, initial solution read in from file
    }
    else if(app.bench == 2)
    {
      // solve stokes to obtain initial guess
      app.solve_stokes();
    }

    comm.print("\nStarting Time-Loop...\n");

    // write initial VTK if desired
    app.write_vtk();

    // time step loop
    while(app.next_time_step())
    {
      // update filters for this time-step
      app.update_filters();

      // initialize solution for this time-step
      app.initialize_step_solution();

      // compute RHS for time step
      app.assemble_step_rhs();

      // solve Navier-Stokes for this step
      if(!app.solve_time_step())
      {
        break;
      }

      // perform solution analysis
      app.perform_solution_analysis();

      // analyze time derivative
      app.analyze_time_derivative();

      // write checkpoint
      app.save_checkpoint();

      // write VTK files if desired
      app.write_vtk();
    } // end of time step loop

    // collect statistics
    app.collect_solver_statistics();

    // release solver
    app.done_solver_symbolic();

    // write solutions if desired
    app.save_final_solution();

    // collect final statistics
    app.collect_final_statistics();

    // print statistics
    app.print_statistics();
  }
} // namespace CCND

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  try
  {
    CCND::main(argc, argv);
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
  return 0;
}
