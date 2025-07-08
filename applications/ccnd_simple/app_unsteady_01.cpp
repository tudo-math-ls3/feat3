// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// ================================================================================================
// Simple Unsteady CCND incompressible Navier-Stokes Solver Application
// ------------------------------------------------------------------------------------------------
// This application implements a simple unsteady monolithic incompressible Navier-Stokes solver
// based on a Q2/P1dc spatial finite element discretization and an implicit BDF(2) time stepping
// scheme, which is  configured to solve a small set of predefined benchmark problems.
//
// To run this application, you have to specify at least the following three mandatory command line
// arguments:
//
// --problem <problem>
// Specifies which benchmark problem to solve; must be one of the following:
// * "simple"     Solves a simple parabolic Poiseuille flow on a 2D/3D unit-square/-cube domain.
// * "conexp"     Simulates a contraction-expansion flow on a 2D domain.
// * "dricav"     Solves a lid-driven cavity problem on a 2D/3D unit-square/-cube domain.
// * "bench2"     Solves a DFG95 flow-around-a-cylinder problem on a flowbench_c2d/c3d domain.
// * "bench8"     Solves a flow-around-a-sphere-inside-a-cylinder problem on a flowbench_p3d domain.
// * "taygre"     Solves a 2D Taylor-Green vortex problem on a unit-square domain.
//                This is the only benchmark which also computes H^k velocity and pressure errors
//                against an analytic reference solution, which can be used to verify the correct
//                order of convergence.
//
// --mesh <meshfile>
// Specifies the mesh file to read; this file must of course match the domain of the problem.
//
// --level <level-max> [[<levels-med...>] <level-min>]
// Specifies the refinement levels for the multigrid hierarchy. The mandatory first parameter
// <level-max> specifies the maximum refinement level on which the system is actually solved.
// The optional final parameter <level-min> specifies the minimum refinement level, which
// represents the coarse mesh level for the multigrid solver. The optional intermediary parameters
// <level-mid...> specify a sequence of partition change levels of the form 'level:ranks' where
// 'level' is the refinement level index and 'ranks' specifies the number of MPI ranks that should
// be used from this level to the next coarser partitioning change level or the coarse level, if
// no further partition change levels were given.
//
// To obtain a list of all other command line arguments, which are currently supported by this
// application, simply compile and run this application without any command line arguments.
//
// To obtain a list of all predefined parameters, simply run this application with a valid problem
// type, mesh and level, and utilize the visual sensors in your skull to parse the program output.
//
// See 'readme.txt' for further information.
//
// \author Peter Zajac
// ================================================================================================

#include "stokes_solver.hpp"
#include "time_stepping.hpp"
#include "check_point.hpp"
#include "vtk_writer.hpp"

int main(int argc, char* argv[])
{
  using namespace CCNDSimple;

  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // create world communicator
  Dist::Comm comm(Dist::Comm::world());

  // print number of processes
  print_pad(comm, "Number of MPI Processes", stringify(comm.size()));

  // create argument parser
  SimpleArgParser args(argc, argv);

  // check command line arguments
  DomainControl::add_supported_args(args);
  StokesSolver::add_supported_args(args);
  TimeStepping::add_supported_args(args);
  CheckPoint::add_supported_args(args);
  VtkWriter::add_supported_args(args);

  args.support("problem",
    "<type>\nSpecifies the problem to solve; valid options are:\n"
    "simple   Simulate a simple Poiseuille flow.\n"
    "conexp   Simulate a contraction-expansion flow.\n"
    "bench1   Simulate a Flow Around A Cylinder."
  );

  // no arguments given?
  if(args.num_args() <= 1)
  {
    comm.print("\nNo arguments given; supported arguments:\n");
    comm.print(args.get_supported_help());
    return 0;
  }

  // check for unsupported options
  {
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      // abort
      return 1;
    }
  }

  // query and verify our problem type
  String problem("simple");
  args.parse("problem", problem);
  if(!problem.is_one_of("simple bench1 conexp"))
  {
    comm.print("ERROR: invalid problem type: '" + problem + "'");
    return 1;
  }

  // a stop watch for our total runtime measurement
  StopWatch watch_total;
  watch_total.start();

  // create domain control object and print info
  DomainControl domain_control(comm);
  domain_control.create_domain(args);
  domain_control.print_info();

  // create Stokes solver, set default parameters, parse arguments and print configuration
  StokesSolver stokes_solver(domain_control);
  stokes_solver.nu = 1E-3;
  stokes_solver.parse_args(args);
  stokes_solver.print_config();

  // create time stepping, set default parameters, parse arguments and print configuration
  TimeStepping time_stepping(comm);
  time_stepping.t_max = 3.0;
  time_stepping.delta_t = 0.01;
  time_stepping.parse_args(args);
  time_stepping.print_config();

  // create checkpoint, parse arguments and print configuration
  CheckPoint check_point(comm);
  check_point.parse_args(args);
  check_point.print_config();

  // create VTK writer, set default parameters, parse arguments and print configuration
  VtkWriter vtk_writer(domain_control);
  vtk_writer.name_prefix = String("ccnd-simple-unsteady-01-") + problem + "-" + stringify(dim) + "d";
  vtk_writer.parse_args(args);
  vtk_writer.print_config();

  // create levels
  stokes_solver.create_levels();

  // get a reference to the finest stokes level
  StokesLevel& the_stokes_level = *stokes_solver.stokes_levels.front();

  // inflow function for parabolic 2D inflow ("simple" problem)
  Analytic::Common::ParProfileVector<DataType> parabolic_inflow_2d(DataType(0), DataType(0), DataType(0), DataType(1));

  // inflow function for parabolic 3D inflow ("simple" problem)
  auto parabolic_inflow_3d = Analytic::create_lambda_function_vector_3d(
    [] (DataType, DataType y, DataType z) {return DataType(16)*y*(DataType(1)-y)*z*(DataType(1)-z);},
    [] (DataType, DataType, DataType) {return DataType(0);},
    [] (DataType, DataType, DataType) {return DataType(0);}
  );

  // 2D inflow function for DFG95 benchmark
  auto dfg95_inflow_2d = Analytic::create_lambda_function_vector_2d(
    [] (DataType, DataType y) {return DataType(6)/DataType(0.1681)*y*(DataType(0.41)-y);},
    [] (DataType, DataType  ) {return DataType(0);});

  // 3D inflow function for DFG95 benchmark
  auto dfg95_inflow_3d = Analytic::create_lambda_function_vector_3d(
    [] (DataType, DataType y, DataType z) {return DataType(36)/DataType(0.02825761)*y*(DataType(0.41)-y)*z*(DataType(0.41)-z);},
    [] (DataType, DataType  , DataType  ) {return DataType(0);},
    [] (DataType, DataType  , DataType  ) {return DataType(0);});

  // inflow function for contraction-expansion flow
  Analytic::Common::ParProfileVector<DataType> conexp_inflow_2d(DataType(0), DataType(0), DataType(0), DataType(1));

  // assemble corresponding boundary conditions depending on which dimension we're in
  if(problem == "simple")
  {
    // Poiseuille-Flow: inflow, no-flow and outflow
    if constexpr (dim == 3)
      stokes_solver.assemble_boundary_conditions("bnd:t bnd:b bnd:n bnd:f", "", "bnd:l", parabolic_inflow_3d, false);
    else
      stokes_solver.assemble_boundary_conditions("bnd:t bnd:b", "", "bnd:l", parabolic_inflow_2d, false);
  }
  else if(problem == "bench1")
  {
    // Flow Around A Cylinder: inflow, no-flow and outflow
    if constexpr (dim == 3)
      stokes_solver.assemble_boundary_conditions("bnd:t bnd:b bnd:n bnd:f bnd:c", "", "bnd:l", dfg95_inflow_3d, false);
    else
      stokes_solver.assemble_boundary_conditions("bnd:t bnd:b bnd:c", "", "bnd:l", dfg95_inflow_2d, false);
  }
  else if(problem == "conexp")
  {
    // contraction-expansion: inflow and slip boundary conditions
    if constexpr(dim == 2)
    {
      stokes_solver.assemble_boundary_conditions(
        "bnd:bl bnd:br bnd:bcl bnd:bcr bnd:bct bnd:tl bnd:tr bnd:tcl bnd:tcr bnd:tcb", "", "bnd:l", conexp_inflow_2d, false);
    }
  }

  // compile local systems for Vanka smoother
  stokes_solver.compile_local_systems();

  // create multigrid stokes_solver
  stokes_solver.create_multigrid_solver();

  // enable or disable full plot mode
  stokes_solver.plot_nonlinear = stokes_solver.plot_nonlinear_header = time_stepping.full_plot;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // create RHS vector
  GlobalStokesVector vec_stokes_rhs = stokes_solver.create_rhs_vector();

  // create initial solution vector
  GlobalStokesVector vec_stokes_sol  = stokes_solver.create_sol_vector();

  // we need up to four additional stokes vectors for our time stepping scheme
  GlobalStokesVector vec_stokes_sol_1 = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);
  GlobalStokesVector vec_stokes_sol_2 = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);
  GlobalStokesVector vec_stokes_sol_3 = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);
  GlobalStokesVector vec_stokes_tmp = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);

  // register Stokes solution vector with VTK writer
  vtk_writer.register_stokes_vector("stokes", vec_stokes_sol);

  // register all our Stokes vectors to the checkpoint
  check_point.register_stokes_vector("stokes_sol_1", vec_stokes_sol_1);
  check_point.register_stokes_vector("stokes_sol_2", vec_stokes_sol_2);
  check_point.register_stokes_vector("stokes_sol_3", vec_stokes_sol_3);

  // restart from checkpoint or load solution vector?
  if(check_point.read_registered(time_stepping))
  {
    vec_stokes_sol.copy(vec_stokes_sol_1);
  }
  else if(!stokes_solver.load_joined_sol_vector(vec_stokes_sol))
  {
    // in any other case, we start by solving a Stokes system
    comm.print("\nSolving steady-state Stokes system...\n");

    // solve Stokes system
    if(!stokes_solver.solve_linear(vec_stokes_sol, vec_stokes_rhs))
      return 1;

    vec_stokes_sol_1.copy(vec_stokes_sol);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Let's go for it
  comm.print(String("\n") + String(100u, '=') + "\n");
  comm.print("Solving unsteady Navier-Stokes system...\n");

  //         "       1      0.1000000000 [  1.00%]       0.000   5: 1.060e-01 > 3.824e-11 :   15;"
  comm.print("    Step       Time           Done       Runtime  NL  Def-Init    Def-Final     MG");
  comm.print("----------------------------------------------------------------------------------");

  time_stepping.begin_loop();

  // write VTK file if desired
  vtk_writer.write_registered(time_stepping);

  // time-stepping loop
  while(time_stepping.advance_step())
  {
    // ============================================================================================
    // Step 1: extrapolate previous Navier-Stokes solution to obtain an initial guess
    // --------------------------------------------------------------------------------------------

    // initialize solution for this time-step by extrapolating the previous ones
    time_stepping.extrapolate_sol(vec_stokes_sol, vec_stokes_sol_1, vec_stokes_sol_2, vec_stokes_sol_3);

    // apply filter
    the_stokes_level.filter_sys.filter_sol(vec_stokes_sol);

    // ============================================================================================
    // Step 2: assemble RHS vector for Navier-Stokes
    // --------------------------------------------------------------------------------------------

    // format RHS
    vec_stokes_rhs.format();

    // combine terms from time stepping scheme into a temporary vector;
    // this one has to be multiplied by the velocity mass matrix and added onto the rhs vector
    time_stepping.combine_mass_rhs(vec_stokes_tmp, vec_stokes_sol_1, vec_stokes_sol_2);

    // add terms from the time stepping scheme to rhs vector
    the_stokes_level.update_unsteady_rhs(vec_stokes_rhs, vec_stokes_tmp);

    // synchronize to obtain a type-1 RHS vector
    vec_stokes_rhs.sync_0();

    // apply filter
    the_stokes_level.filter_sys.filter_rhs(vec_stokes_rhs);

    // ============================================================================================
    // Step 3: solve Navier-Stokes equations
    // --------------------------------------------------------------------------------------------

    // set theta for Burgers matrix assembly
    stokes_solver.theta = time_stepping.theta[0];

    // solve the actual Navier-Stokes system for this time-step
    if(!stokes_solver.solve_nonlinear(vec_stokes_sol, vec_stokes_rhs))
      return 1;

    // ============================================================================================
    // Step 4: write VTK file if desired
    // --------------------------------------------------------------------------------------------

    // write VTK file if desired
    vtk_writer.write_registered(time_stepping);

    // ============================================================================================
    // Step 5: update Stokes vectors
    // --------------------------------------------------------------------------------------------

    // shift the Stokes solution vectors
    vec_stokes_sol_3.copy(vec_stokes_sol_2); // u_{k-3}
    vec_stokes_sol_2.copy(vec_stokes_sol_1); // u_{k-2}
    vec_stokes_sol_1.copy(vec_stokes_sol);   // u_{k-1}

    // ============================================================================================
    // Step 6: write checkpoint if desired
    // --------------------------------------------------------------------------------------------

    // save checkpoint if desired
    check_point.write_registerd(time_stepping);

    // ============================================================================================
    // Step 7: print plot line
    // --------------------------------------------------------------------------------------------

    // print plot line if full plot mode is disabled
    if(!time_stepping.full_plot)
    {
      String plot_line = time_stepping.plot_line;
      plot_line += stokes_solver.plot_line;
      plot_line += check_point.plot_line;
      plot_line += vtk_writer.plot_line;
      comm.print(plot_line);
    }

    // flush cout buffer
    comm.print_flush();

  } // time stepping loop

  time_stepping.finish_loop();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // write joined solution vector if desired
  stokes_solver.save_joined_sol_vector(vec_stokes_sol);

  // release the stokes stokes_solver
  stokes_solver.release();

  // make sure everyone is finished
  comm.barrier();

  // stop watch
  watch_total.stop();

  // print final statistics
  comm.print(String("\n") + String(100u, '=') + "\n");
  print_pad(comm, "Total Runtime", watch_total.elapsed_string().pad_front(10) + " sec {Balance} [Fraction]");
  domain_control.print_runtime(watch_total.elapsed());
  time_stepping.print_runtime(watch_total.elapsed());
  stokes_solver.print_runtime(watch_total.elapsed());
  check_point.print_runtime(watch_total.elapsed());
  vtk_writer.print_runtime(watch_total.elapsed());
  print_memory_usage(comm);

  // okay, we're done
  return 0;
}
