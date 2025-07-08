// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// ================================================================================================
// Simple Steady-State CCND incompressible Navier-Stokes Solver Application
// ------------------------------------------------------------------------------------------------
// This application implements a simple steady-state monolithic incompressible Navier-Stokes solver
// based on a Q2/P1dc spatial finite element discretization, which is  configured to solve a small
// set of predefined benchmark problems.
//
// To run this application, you have to specify at least the following three mandatory command line
// arguments:
//
// --problem <problem>
// Specifies which benchmark problem to solve; must be one of the following:
// * "simple"     Solves a simple parabolic Poiseuille flow on a 2D/3D unit-square/-cube domain.
// * "conexp"     Simulates a contraction-expansion flow on a 2D domain.
// * "dricav"     Solves a lid-driven cavity problem on a 2D/3D unit-square/-cube domain.
// * "bench1"     Solves a DFG95 flow-around-a-cylinder problem on a flowbench_c2d/c3d domain.
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
#include "vtk_writer.hpp"

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  using namespace CCNDSimple;

  // create world communicator
  Dist::Comm comm(Dist::Comm::world());

  // print number of processes
  print_pad(comm, "Number of MPI Processes", stringify(comm.size()));

  // create argument parser
  SimpleArgParser args(argc, argv);

  // check command line arguments
  DomainControl::add_supported_args(args);
  StokesSolver::add_supported_args(args);
  VtkWriter::add_supported_args(args);

  args.support("problem",
    "<type>\nSpecifies the problem to solve; valid options are:\n"
    "simple   Simulate a simple Poiseuille flow.\n"
    "conexp   Simulate a contraction-expansion flow.\n"
    "dricav   Simulate a Driven Cavity flow.\n"
    "bench1   Simulate a Flow Around A Cylinder.\n"
    "taygre   Simulate a Taylor-Green vortex (2D only)."
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
  if(!problem.is_one_of("simple dricav bench1 taygre conexp"))
  {
    comm.print("ERROR: invalid problem type: '" + problem + "'");
    return 1;
  }
  if((problem == "taygre") && (dim != 2))
  {
    comm.print("ERROR: Taylor-Green simulation is only available in 2D");
    return 1;
  }

  StopWatch watch_total;
  watch_total.start();

  // create domain control object and print info
  DomainControl domain_control(comm);
  domain_control.create_domain(args);
  domain_control.print_info();

  // create Stokes solver object, set default parameter, parse arguments and print info
  StokesSolver stokes_solver(domain_control);
  stokes_solver.parse_args(args);
  stokes_solver.print_config();

  // create VTK writer
  VtkWriter vtk_writer(domain_control);
  vtk_writer.name_prefix = String("ccnd-simple-steady-01-") + problem + "-" + stringify(dim) + "d";
  vtk_writer.parse_args(args);
  vtk_writer.print_config();

  // create levels
  stokes_solver.create_levels();

  // null functions
  Analytic::Common::ConstantVectorFunction<dim, DataType> null_function;

  // analytic solution for 2D Taylor-Green problem
  Analytic::Common::TaylorGreenVortexVelo2D<DataType> taylor_green_velo_2d;
  Analytic::Common::TaylorGreenVortexPres2D<DataType> taylor_green_pres_2d;

  // inflow function for parabolic 2D inflow ("simple" problem)
  Analytic::Common::ParProfileVector<DataType> parabolic_inflow_2d(DataType(0), DataType(0), DataType(0), DataType(1));

  // inflow function for parabolic 3D inflow ("simple" problem)
  auto parabolic_inflow_3d = Analytic::create_lambda_function_vector_3d(
    [] (DataType, DataType y, DataType z) {return DataType(16)*y*(DataType(1)-y)*z*(DataType(1)-z);},
    [] (DataType, DataType, DataType) {return DataType(0);},
    [] (DataType, DataType, DataType) {return DataType(0);}
  );

  // inflow function for parabolic 2D inflow for DFG95 benchmark
  auto dfg95_inflow_2d = Analytic::create_lambda_function_vector_2d(
    [] (DataType, DataType y) {return DataType(12)/DataType(1.681)*y*(DataType(0.41)-y);},
    [] (DataType, DataType) {return DataType(0);}
  );

  // inflow function for parabolic 3D inflow for DFG95 benchmark
  auto dfg95_inflow_3d = Analytic::create_lambda_function_vector_3d(
    [] (DataType, DataType y, DataType z) {return DataType(72)/DataType(0.2825761)*y*(DataType(0.41)-y)*z*(DataType(0.41)-z);},
    [] (DataType, DataType, DataType) {return DataType(0);},
    [] (DataType, DataType, DataType) {return DataType(0);}
  );

  // boundary function for 2D driven cavity
  auto dricav_flow_2d = Analytic::create_lambda_function_vector_2d(
    [] (DataType x, DataType) {return (x > DataType(0.0001)) && (x < DataType(0.9999)) ? DataType(1) : DataType(0);},
    [] (DataType, DataType) {return DataType(0);}
  );

  // boundary function for 3D driven cavity
  auto dricav_flow_3d = Analytic::create_lambda_function_vector_3d(
    [] (DataType x, DataType, DataType) {return (x > DataType(0.0001)) && (x < DataType(0.9999)) ? DataType(1) : DataType(0);},
    [] (DataType, DataType, DataType) {return DataType(0);},
    [] (DataType, DataType, DataType) {return DataType(0);}
  );

  // inflow function for contraction-expansion flow
  Analytic::Common::ParProfileVector<DataType> conexp_inflow_2d(DataType(0), DataType(0), DataType(0), DataType(1));

  // assemble boundary conditions
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
    if constexpr(dim == 3)
      stokes_solver.assemble_boundary_conditions("bnd:t bnd:b bnd:n bnd:f bnd:c", "", "bnd:l", dfg95_inflow_3d, false);
    else
      stokes_solver.assemble_boundary_conditions("bnd:t bnd:b bnd:c", "", "bnd:l", dfg95_inflow_2d, false);
  }
  else if(problem == "dricav")
  {
    // Driven Cavity: all Dirichlet boundary conditions
    if constexpr(dim == 3)
      stokes_solver.assemble_boundary_conditions("bnd:l bnd:r bnd:b bnd:n bnd:f", "", "bnd:t", dricav_flow_3d, true);
    else
      stokes_solver.assemble_boundary_conditions("bnd:l bnd:r bnd:b", "", "bnd:t", dricav_flow_2d, true);
  }
  else if(problem == "taygre")
  {
    // Taylor Green: all slip boundary conditions
    stokes_solver.assemble_boundary_conditions("", "bnd:l | bnd:r | bnd:b | bnd:t", "", null_function, true);
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

  // create multigrid solver
  stokes_solver.create_multigrid_solver();

  // create initial solution vector
  GlobalStokesVector vec_sol = stokes_solver.create_sol_vector();

  // create null RHS vector
  GlobalStokesVector vec_rhs = stokes_solver.create_rhs_vector();

  // in case of Taylor-Green, assemble the corresponding RHS vector
  if(problem == "taygre")
  {
    if constexpr(dim == 2)
    {
      // declare constant pi
      const DataType pi = Math::pi<DataType>();
      const DataType nu = stokes_solver.nu;

      // RHS function for 2D Taylor-Green problem
      auto taylor_green_rhs_2d = Analytic::create_lambda_function_vector_2d(
        [pi,nu] (DataType x, DataType y) {return  DataType(2)*pi*pi*nu*Math::sin(pi*x)*Math::cos(pi*y);},
        [pi,nu] (DataType x, DataType y) {return -DataType(2)*pi*pi*nu*Math::cos(pi*x)*Math::sin(pi*y);}
      );

      // create RHS vector from force function
      vec_rhs = stokes_solver.create_rhs_vector(taylor_green_rhs_2d);
    }
  }

  // load joined solution if desired
  if(!stokes_solver.load_joined_sol_vector(vec_sol))
  {
    // no solution loaded; solve Stokes instead
    comm.print("\nSolving Stokes system...\n");

    // solve Stokes system
    if(!stokes_solver.solve_linear(vec_sol, vec_rhs))
      return 1;
  }

  // now let's go for Navier-Stokes
  comm.print("\nSolving Navier-Stokes system...\n");

  // solve the actual non-linear system
  if(!stokes_solver.solve_nonlinear(vec_sol, vec_rhs))
    return 1;

  // post-process solution in dependence of the chosen problem
  if(problem == "taygre")
  {
    if constexpr(dim == 2)
    {
      // compute errors against reference solution
      stokes_solver.compute_errors(vec_sol, taylor_green_velo_2d, taylor_green_pres_2d);
    }
  }
  if(problem == "bench1")
  {
    // compute body forces on circle/cylinder
    stokes_solver.compute_body_forces("bnd:c", vec_sol, vec_rhs,
      dim == 2 ? DataType(500) : DataType(50000) / DataType(41));
  }

  // write joined solution vector if desired
  stokes_solver.save_joined_sol_vector(vec_sol);

  // release the stokes stokes_solver
  stokes_solver.release();

  // write VTK file if desired
  if(vtk_writer.prepare_write())
  {
    comm.print("\nWriting VTK output to '" + vtk_writer.vtk_name + ".pvtu'");
    vtk_writer.add_stokes_vector(vec_sol);
    vtk_writer.write();
  }

  // make sure everyone is finished
  comm.barrier();

  // stop watch and print statistics
  watch_total.stop();
  comm.print(String("\n") + String(100u, '=') + "\n");
  print_pad(comm, "Total Runtime", watch_total.elapsed_string().pad_front(10) + " sec {Balance} [Fraction]", 40);
  domain_control.print_runtime(watch_total.elapsed());
  stokes_solver.print_runtime(watch_total.elapsed());
  vtk_writer.print_runtime(watch_total.elapsed());
  print_memory_usage(comm);

  // okay, we're done
  return 0;
}
