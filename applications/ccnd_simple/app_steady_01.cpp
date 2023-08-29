// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
// * "dricav"     Solves a lid-driven cavity problem on a 2D/3D unit-square/-cube domain.
// * "bench1"     Solves a DFG95 flow-around-a-cylinder problem on a flowbench_c2d/c3d domain.
// * "bench7"     Solves a flow-around-a-sphere-inside-a-cylinder problem on a flowbench_p3d domain.
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

#include "steady_solver.hpp"
#include "vtk_writer.hpp"

namespace CCNDSimple
{
  int main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    print_pad(comm, "Number of MPI Processes", stringify(comm.size()));

    // create argument parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    DomainControl::add_supported_args(args);
    SteadySolver::add_supported_args(args);
    VtkWriter::add_supported_args(args);

    args.support("problem",
      "<type>\nSpecifies the problem to solve; valid options are:\n"
      "simple   Simulate a simple Poiseuille flow.\n"
      "dricav   Simulate a Driven Cavity flow.\n"
      "bench1   Simulate a Flow Around A Cylinder.\n"
      "bench7   Simulate a Flow Inside A Cylinder (3D only).\n"
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
    if(!problem.is_one_of("simple dricav bench1 bench7 taygre"))
    {
      comm.print("ERROR: invalid problem type: '" + problem + "'");
      return 1;
    }
    if((problem == "bench7") && (dim != 3))
    {
      comm.print("ERROR: bench7 simulation is only available in 3D");
      return 1;
    }
    if((problem == "taygre") && (dim != 2))
    {
      comm.print("ERROR: Taylor-Green simulation is only available in 2D");
      return 1;
    }

    StopWatch watch_total;
    watch_total.start();

    // our domain control object
    DomainControl domain(comm);

    // create the domain
    domain.create(args);

    // print domain info
    domain.print_info();

    // our Stokes solver
    SteadySolver solver(domain);

    // parse arguments and print parsed configuration
    solver.parse_args(args);
    solver.print_config();

    // create VTK writer
    VtkWriter writer(domain, String("ccnd-simple-steady-01-") + problem + "-" + stringify(dim) + "d");
    writer.parse_args(args);
    writer.print_config();

    // create levels
    solver.create_levels();

    // analytic solution for Taylor-Green
    Analytic::Common::TaylorGreenVortexVelo2D<DataType> taylor_green_velo;
    Analytic::Common::TaylorGreenVortexPres2D<DataType> taylor_green_pres;

    // assemble boundary conditions
    if(problem == "simple")
    {
      // Poiseuille-Flow: inflow, no-flow and outflow
      if(dim == 3)
        solver.assemble_boundary_conditions("bnd:t bnd:b bnd:n bnd:f", "", "bnd:l", "16*y*(1-y)*z*(1-z) ' 0 ' 0", false);
      else
        solver.assemble_boundary_conditions("bnd:t bnd:b", "", "bnd:l", "4*y*(1-y) ' 0", false);
    }
    else if(problem == "bench1")
    {
      // Flow Around A Cylinder: inflow, no-flow and outflow
      if(dim == 3)
        solver.assemble_boundary_conditions("bnd:t bnd:b bnd:n bnd:f bnd:c", "", "bnd:l", "72/0.2825761*y*(0.41-y)*z*(0.41-z) ' 0 ' 0", false);
      else
        solver.assemble_boundary_conditions("bnd:t bnd:b bnd:c", "", "bnd:l", "12/1.681*y*(0.41-y) ' 0", false);
    }
    else if(problem == "bench7")
    {
      // Flow Inside A Cylinder: inflow, no-flow and outflow
      solver.assemble_boundary_conditions("bnd:pipe bnd:sphere", "", "bnd:in", "40000/1681*(0.205^2-(y-0.2)^2-(z-0.205)^2) ' 0 ' 0", false);
    }
    else if(problem == "dricav")
    {
      // Driven Cavity: all Dirichlet boundary conditions
      if(dim == 3)
        solver.assemble_boundary_conditions("bnd:l bnd:r bnd:b bnd:n bnd:f", "", "bnd:t", "if((x>0.0001) & (x<0.9999),1,0) ' 0 ' 0", true);
      else
        solver.assemble_boundary_conditions("bnd:l bnd:r bnd:b", "", "bnd:t", "1 ' 0", true);
    }
    else if(problem == "taygre")
    {
      // Taylor Green: all slip boundary conditions
      solver.assemble_boundary_conditions("", "bnd:l | bnd:r | bnd:b | bnd:t", "", "", true);
    }

    // compile local systems for Vanka smoother
    solver.compile_local_systems();

    // create multigrid solver
    solver.create_multigrid_solver();

    // create initial solution vector
    GlobalStokesVector vec_sol = solver.create_sol_vector();

    // create RHS vector
    GlobalStokesVector vec_rhs = solver.create_rhs_vector();

    // in case of Taylor-Green, assemble the corresponding RHS vector
    if(problem == "taygre")
    {
      String formula_vx = " 2*(pi^2)*(" + stringify(solver.nu) + ")*sin(pi*x)*cos(pi*y)";
      String formula_vy = "-2*(pi^2)*(" + stringify(solver.nu) + ")*cos(pi*x)*sin(pi*y)";
      if(!solver.navier_stokes)
      {
        formula_vx += "-pi*cos(pi*x)*sin(pi*x)";
        formula_vy += "-pi*sin(pi*y)*cos(pi*y)";
      }
      vec_rhs = solver.create_rhs_vector(formula_vx + " ' " + formula_vy);
    }

    // load joined solution if desired
    if(!solver.load_joined_sol_vector(vec_sol))
    {
      // no solution loaded; solve Stokes instead
      comm.print("\nSolving Stokes system...\n");

      // solve Stokes system
      if(!solver.solve_stokes(vec_sol, vec_rhs))
        return 1;
    }

    // now let's go for Navier-Stokes
    comm.print("\nSolving Navier-Stokes system...\n");

    if(!solver.solve_navier_stokes(vec_sol, vec_rhs))
      return 1;

    // post-process solution in dependence of the chosen problem
    if(problem == "taygre")
    {
#ifndef FEAT_CCND_SIMPLE_3D
      // compute errors against reference solution
      solver.compute_errors(vec_sol, taylor_green_velo, taylor_green_pres);
#endif
    }
    if(problem == "bench1")
    {
      // compute body forces on circle/cylinder
      solver.compute_body_forces("bnd:c", vec_sol, vec_rhs);
    }
    if(problem == "bench7")
    {
      // compute body forces on sphere
      solver.compute_body_forces("bnd:sphere", vec_sol, vec_rhs);
    }

    // write joined solution vector if desired
    solver.save_joined_sol_vector(vec_sol);

    // release the stokes solver
    solver.release();

    // write VTK file if desired
    if(writer.prepare_write())
    {
      comm.print("\nWriting VTK output to '" + writer.vtk_name + ".pvtu'");
      writer.add_stokes_vector(vec_sol);
      writer.write();
    }

    // stop watch and print statistics
    watch_total.stop();
    print_final_stats(comm, watch_total);

    // okay, we're done
    return 0;
  }
} // namespace CCNDSimple

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  int rtn = 1;
  try
  {
    rtn = CCNDSimple::main(argc, argv);
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
  return rtn;
}
