// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include "base.hpp"

namespace Turb
{
  /**
   * \brief Parameter settings for the turbulence/laminar solver
   *
   * This class gathers all model and solver parameters that are typically
   * provided via the command line.  It offers three convenience functions
   * that follow the same convention as the TimeStepping class:
   *   - add_supported_args() registers all supported command‑line flags.
   *   - parse_args() reads the values and performs basic consistency checks.
   *   - print_config() writes a nicely formatted summary to the console.
   *
   * All parameters are initialized with reasonable default values so that
   * the simulation can be launched without specifying any optional flag.
   */
  class ParameterSettings
  {
  public:
    /// communicator
    const Dist::Comm& comm;

    // Reynolds numbers and mean velocity
    DataType Re_laminar   = DataType(1000);     // laminar Reynolds number
    DataType Re_turbulent = DataType(47625);    // turbulent Reynolds number
    DataType u_mean       = DataType(1.0);      //  mean velocity

    // turbulence model constants
    DataType C_1      = DataType(1.44);
    DataType C_2      = DataType(1.92);
    DataType C_mu     = DataType(0.09);
    DataType c_bc     = DataType(0.007);        // blending constant c_bc in [0.003,0.01]
    DataType l_0      = DataType(0.08);         // l in [l_min, l_max]
    DataType sigma_e  = DataType(1.3);
    DataType sigma_k  = DataType(1.0);

    // viscosity / length limits and flow regime switch
    DataType nu_min   = DataType(1.0e-7);       // lower bound for viscosity

    // global time management (separate limits for laminar vs. turbulent)
    bool     do_laminar      = true;            // run laminar simulation instead of k‑e
    DataType t_max_turb      = DataType(50.0);
    DataType delta_t_turb    = DataType(0.001);
    DataType t_max_laminar   = DataType(0.002);
    DataType delta_t_laminar = DataType(0.001);
    IndexType vtk_stepping   = IndexType(100);

    // initial conditions (computed in constructor / parse_args)
    DataType nu_0  = DataType(0.09);
    DataType k_0   = DataType(0.0);             // will be computed from c_bc & u_mean
    DataType e_0   = DataType(0.0);             // will be computed from k_0 & l_0

    // scaling parameters for convection, diffusion, rhs
    DataType scal_conv = DataType(1.0);
    DataType scal_diff = DataType(1.0);
    DataType scal_reac = DataType(1.0);
    DataType scal_rhs  = DataType(1.0);

    // cubature rule for numerical integration
    String cubature_name = String("gauss-legendre:3");

    // use velocity solution from previous iteration?
    bool use_previous_velo = false;

    // number of coupling steps for outer k loop
    IndexType coupling_steps = IndexType(1);

    // Recycle outflow as inflow in iterative boundary treatment
    bool use_iter_inflow = false;

    // Type of boundary conditions (Dirichlet or Neumann)
    String problem = "dirichlet";

  public:
    // ctor – computes derived ICs from defaults
    explicit ParameterSettings(const Dist::Comm& comm_) : comm(comm_)
    {
      compute_initial_conditions();
    }

    // dtor
    virtual ~ParameterSettings() = default;

    // non‑copyable
    ParameterSettings(const ParameterSettings&) = delete;
    ParameterSettings& operator=(const ParameterSettings&) = delete;

    // registers all command‑line flags that this class understands
    static void add_supported_args(SimpleArgParser& args)
    {
      // time stepping parameter
      args.support("vtk-stepping",    "<n>\nWrite VTK output every n time steps.");

      // problem parameter
      args.support("problem",    "<n>\nType of boundary condition (dirichlet or neumann)");
    }

    // parses the supplied command‑line arguments and updates the members
    bool parse_args(SimpleArgParser& args)
    {
      args.parse("vtk-stepping", vtk_stepping);
      if(args.check("problem") < 1)
      {
        comm.print(std::cerr, "ERROR: Mandatory option '--problem <problem>' is missing!");
        FEAT::Runtime::abort();
      }
      args.parse("problem", problem);
      if(problem != "dirichlet" && problem != "neumann")
      {
        comm.print(std::cerr, "ERROR: Selected option '--problem <problem>' is not possible only dirichlet or neumann!");
      }
      compute_initial_conditions();
      return true;
    }

    // prints a nicely formatted summary of the parameter configuration
    void print_config() const
    {
      using namespace std; // print_pad is usually in global ns
      print_pad(comm, "Re_laminar",   stringify_fp_fix(Re_laminar));
      print_pad(comm, "Re_turbulent", stringify_fp_fix(Re_turbulent));
      print_pad(comm, "u_mean",       stringify_fp_fix(u_mean));
      print_pad(comm, "C_1",     stringify_fp_fix(C_1));
      print_pad(comm, "C_2",     stringify_fp_fix(C_2));
      print_pad(comm, "C_mu",    stringify_fp_fix(C_mu));
      print_pad(comm, "c_bc",    stringify_fp_fix(c_bc));
      print_pad(comm, "l_0",     stringify_fp_fix(l_0));
      print_pad(comm, "sigma_e", stringify_fp_fix(sigma_e));
      print_pad(comm, "sigma_k", stringify_fp_fix(sigma_k));
      print_pad(comm, "nu_min",  stringify_fp_fix(nu_min));
      print_pad(comm, "do_laminar", (do_laminar ? "true" : "false"));
      print_pad(comm, "t_max_turb",      stringify_fp_fix(t_max_turb));
      print_pad(comm, "delta_t_turb",    stringify_fp_fix(delta_t_turb));
      print_pad(comm, "t_max_laminar",   stringify_fp_fix(t_max_laminar));
      print_pad(comm, "delta_t_laminar", stringify_fp_fix(delta_t_laminar));
      print_pad(comm, "vtk_stepping",    stringify(vtk_stepping));
      print_pad(comm, "nu_0", stringify_fp_fix(nu_0));
      print_pad(comm, "k_0",  stringify_fp_fix(k_0));
      print_pad(comm, "e_0",  stringify_fp_fix(e_0));
      print_pad(comm, "scal_conv", stringify_fp_fix(scal_conv));
      print_pad(comm, "scal_diff", stringify_fp_fix(scal_diff));
      print_pad(comm, "scal_reac", stringify_fp_fix(scal_reac));
      print_pad(comm, "scal_rhs",  stringify_fp_fix(scal_rhs));
      print_pad(comm, "cubature", cubature_name);
      print_pad(comm, "use_previous_velo", (use_previous_velo ? "true" : "false"));
      print_pad(comm, "coupling_steps", stringify(coupling_steps));
      print_pad(comm, "problem", stringify(problem));
    }

  private:
    /// (re)computes k_0 and e_0 from the current parameters
    void compute_initial_conditions()
    {
      k_0 = c_bc * u_mean * u_mean;
      e_0 = C_mu * Math::sqrt(k_0) * k_0 / l_0;
    }
  }; // class ParameterSettings
} // namespace Turb
