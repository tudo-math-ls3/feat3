// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Unsteady CCnD solver application base class header
// ------------------------------------------------------------------------------------------------
// This header file defines the CCND::UnsteadyAppBase class, which is itself derived from the
// CCND::SteadyAppBase class and which can be used as a base class for CCnD type applications to
// outsource commonly used functionality. This comment block describes the command line parameters
// that can be used by all applications derived from this class that are not already offered by the
// base class CCND::SteadyAppBase, unless the application overrides the parameters explicitly.
//
// The system is discretized using an isoparametric Q2/P1dc finite element discretization with
// a BDF(2) time stepping scheme with constant time step length. The monolithic nonlinear Oseen
// systems are solved using an adaptive Newton-Multigrid solver with an additive matrix-based
// Vanka smoother ("AmaVanka") and using UMFPACK (if available) as a coarse grid solver. The
// smoother iteration as well as the outer multigrid iteration is usually performed by a Richardson
// iteration scheme, however, it is also possible to switch to a FGMRES(k) solver instead.
// This application supports recursive partitioning.
//
//
// ----------------------------------------------
// Mandatory Time Discretization Setup Parameters
// ----------------------------------------------
// This application class defines default values for most of its parameters, however, the following
// parameters are mandatory and always have to be specified explicitly:
//
// --time-max <T-max>
// Specifies the maximum simulation time, i.e. the end T of the simulation time interval (0,T].
// Not to be confused with runtime!
//
// --time-steps <N>
// Specifies the total number of time-steps to perform over the entire time interval.
//
//
// -----------------------------------------------
// Additional Time Discretization Setup Parameters
// -----------------------------------------------
// The following parameters allow to further configure the time discretization.
//
// --time-expo <0|1|2>
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
// ---------------------------------------------
// Checkpoint & Restart Configuration Parameters
// ---------------------------------------------
// This application base implements a checkpoint-&-restart system based on the CheckpointControl
// class, which allows to restart a simulation that has been aborted at a later time.
// Note that a checkpoint file contains three consecutive solution vectors u_{k}, u_{k-1} and
// u_{k-2} as well as the time-step index, the current simulation time and the time step size,
// so these three quantities are restored from a checkpoint in a restart and do not need to be
// specified explicitly.
//
// Important Note:
// It is very important that the application, which is restarted from a previously saved
// checkpoint, is launched with the same number of MPI processes, the same mesh and the same
// maximum refinement level as the application run that was used to create the checkpoints.
// However, many other parameters (e.g. viscosity, gradient/deformation tensor, linear solver
// configuration) may be changed upon restart.
//
// --checkpoint <filename> [<stepping> [<modulus>]]
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
// --vtk-step <step>
// Specifies that the VTK files should only be written out every <step> time steps rather than in
// every single time step.
//
// --main-step <step>
// Specifies that the post-processing of the solution is only to be performed every <step> time
// steps rather than in every single time step.
//
// \author Peter Zajac
//
#pragma once

#include "ccnd_steady_appbase.hpp"

namespace CCND
{
  class UnsteadyAppBase :
    public SteadyAppBase
  {
  public:
    /// our base class
    typedef SteadyAppBase BaseClass;

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Important Note: All the default values below may be overridden by the constructors of derived classes!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    /// maximum simulation time
    DataType time_max = DataType(1);

    /// number of time steps
    Index max_time_steps = Index(1000u);

    /// current time step
    Index cur_step = Index(0);

    /// current time
    DataType cur_time = DataType(0);

    /// time step size delta_t = t_max / max_time_steps
    DataType delta_t = 0.01;

    /// steady-state detection tolerance
    DataType steady_tol = DataType(1E-3);

    /// time extrapolation order: 0, 1 or 2
    IndexType time_expo = Index(2);

    /// current velocity derivative norm
    DataType velo_deriv_norm = DataType(0);

    /// checkpoint output filename
    String checkpoint_filename;
    /// checkpoint stepping; set to 0 to disable checkpointing
    IndexType checkpoint_step = Index(1);
    /// checkpoint stepping modulo, i.e. number of checkpoints to store; should be at least 2
    IndexType checkpoint_mod = Index(2);

    /// restart input filename
    String restart_filename;
    /// time of restart checkpoint; 0 => use checkpoint data
    DataType restart_time = DataType(0);
    /// time-step of restart checkpoint; 0 => use checkpoint data
    Index restart_step = Index(0);

    /// main stepping; all post-processing will only be performed in each main step
    Index main_step = 1;

    /// VTK output stepping; set to 0 to disable VTK output
    Index vtk_step = 0;

    /// specifies whether the simulation has a non-zero right-hand-side for non-steady simulations
    bool homogeneous_unsteady_rhs = false;

    // more vectors
    GlobalSystemVector vec_sol_1, vec_sol_2, vec_der, vec_tmp;

    // time step time stamp
    TimeStamp stamp_cur_step;

    // more stop watches
    StopWatch watch_checkpoint_save, watch_checkpoint_load;

  public:
    explicit UnsteadyAppBase(const Dist::Comm& comm_, SimpleArgParser& args_) :
      BaseClass(comm_, args_)
    {
      // override default values of base class
      max_mg_iter = 5;
      mg_tol_rel = 1E-5;

      args.support("time-max");
      args.support("time-steps");
      args.support("time-expo");
      args.support("steady-tol");
      args.support("checkpoint");
      args.support("restart");
      args.support("vtk-step");
      args.support("main-step");
      //args.support("save-joined-velo");
      //args.support("save-joined-pres");
    }

    virtual String line_prefix() const
    {
      return stringify(cur_step).pad_front(6) + ":" + stringify_fp_fix(cur_time, 8, 12);
    }

    virtual void parse_args() override
    {
      BaseClass::parse_args();

      args.parse("time-max", time_max);
      args.parse("time-steps", max_time_steps);
      args.parse("time-expo", time_expo);
      args.parse("steady-tol", steady_tol);

      delta_t = time_max / DataType(max_time_steps);

      args.parse("checkpoint", checkpoint_filename, checkpoint_step, checkpoint_mod);
      args.parse("restart", restart_filename, restart_time, restart_step);

      vtk_step = Index(args.check("vtk") >= 0 ? 1 : 0);
      args.parse("vtk", vtk_filename, vtk_step);
      args.parse("vtk-step", vtk_step);
      args.parse("main-step", main_step);

      want_vtk = want_vtk || (args.check("vtk-step") >= 0);
    }

    virtual void print_problem() override
    {
      BaseClass::print_problem();

      comm.print(String("T-max").pad_back(pad_len, pad_char) + ": " + stringify(time_max));
      comm.print(String("Number of Time Steps").pad_back(pad_len, pad_char) + ": " + stringify(max_time_steps));
      comm.print(String("Main Stepping").pad_back(pad_len, pad_char) + ": " + stringify(main_step));
      comm.print(String("Time Step Length").pad_back(pad_len, pad_char) + ": " + stringify(delta_t));
      comm.print(String("Time Extrapolation Order").pad_back(pad_len, pad_char) + ": " + stringify(time_expo));
      comm.print(String("Steady State Derivative Tol").pad_back(pad_len, pad_char) + ": " + stringify_fp_sci(steady_tol));

      if(!checkpoint_filename.empty())
      {
        comm.print(String("Checkpoint Filename").pad_back(pad_len, pad_char) + ": " + checkpoint_filename);
        comm.print(String("Checkpoint Stepping").pad_back(pad_len, pad_char) + ": " + stringify(checkpoint_step));
        comm.print(String("Checkpoint Modulus").pad_back(pad_len, pad_char) + ": " + stringify(checkpoint_mod));
      }
      else
        comm.print(String("Checkpoint Filename").pad_back(pad_len, pad_char) + ": -N/A-");

      if(!restart_filename.empty())
      {
        comm.print(String("Restart Filename").pad_back(pad_len, pad_char) + ": " + restart_filename);
        comm.print(String("Restart Time").pad_back(pad_len, pad_char) + ": " + stringify(restart_time));
        comm.print(String("Restart Step").pad_back(pad_len, pad_char) + ": " + stringify(restart_step));
      }
      else
        comm.print(String("Restart Filename").pad_back(pad_len, pad_char) + ": -N/A-");

    }

    virtual void create_vectors() override
    {
      BaseClass::create_vectors();

      watch_create_system.start();
      vec_sol_1 = vec_sol.clone(LAFEM::CloneMode::Weak);
      vec_sol_2 = vec_sol.clone(LAFEM::CloneMode::Weak);
      vec_der = vec_sol.clone(LAFEM::CloneMode::Weak);
      vec_tmp = vec_sol.clone(LAFEM::CloneMode::Weak);
      watch_create_system.stop();
    }

    virtual bool load_checkpoint()
    {
      // no restart?
      if(restart_filename.empty())
        return false;

      watch_checkpoint_load.start();

      //if(log_step > Index(0))
        comm.print("\n" + String(100, '=') + "\nReading restart checkpoint '" + restart_filename + "'...");

      // create checkpoint control
      LAFEM::SerialConfig serial_config(false, false);
      Control::CheckpointControl check_ctrl(comm, serial_config);

      // create checkpoint info vector
      LAFEM::DenseVector<DataType, IndexType> vcp_info;

      // load checkpoint
      check_ctrl.load(restart_filename);

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
      cur_step = (restart_step <= Index(0) ? cp_time_step : restart_step - Index(1));
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
      comm.print("Restart Timestep...: " + stringify(cur_step)
        + (restart_step <= Index(0) ? " (restored)" : " (override)"));

      // do we need to reconstruct the vectors due to changed time step size?
      if(!same_dt && (time_expo >= Index(2)))
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
        vec_sol_1.axpy(old_sol_1, q * (DataType(2) - q));
        vec_sol_1.axpy(old_sol_2, q * (q - DataType(1)) / DataType(2));

        // transform v[k-2] :=      L_0(-2q)*u[k] +  L_1(-2q)*u[k-1] +  L_2(-2q)*u[k-2]
        //                   = (q-1)*(2*q-1)*u[k] + 4*q*(1-q)*u[k-1] + q*(2*q-1)*u[k-2]
        vec_sol_2.scale(vec_sol, (q - DataType(1))*(DataType(2)*q - DataType(1)));
        vec_sol_2.axpy(old_sol_1, DataType(4) * q * (DataType(1) - q));
        vec_sol_2.axpy(old_sol_2, q * (DataType(2)*q - DataType(1)));
      }
      else if(!same_dt && (time_expo >= Index(1)))
      {
        // same approach as the previous case, but only interpolate {v[k-1],v[k]} linearly
        // from {u[k-1],u[k]} by using the first order Lagrange polynomials for {-1, 0}:
        // L_0(t) := x + 1
        // L_1(t) := -x

        comm.print("Performing linear transformation of restored solution vectors");

        // transform v[k-1] = (1-q)*u[k] + q*v[k-1]
        vec_sol_1.scale(vec_sol_1, q);
        vec_sol_1.axpy(vec_sol, DataType(1) - q);
      }
      // If t_expo is < 1, i.e. constant extrapolation, then no transformation is
      // required even if the timestep size has changed.

      comm.print(String(100, '='));
      watch_checkpoint_load.stop();

      return true;
    }

    virtual void save_checkpoint()
    {
      if(checkpoint_filename.empty())
        return;

      if((checkpoint_step > Index(1)) && (cur_step % checkpoint_step != 0))
        return;

      // it's checkpointing time!
      watch_checkpoint_save.start();

      // generate checkpoint name
      String check_fn = checkpoint_filename + "." + stringify(((cur_step-1) / checkpoint_step) % checkpoint_mod) + ".cp";
      comm.print(line_prefix() + " > Writing Checkpoint '" + check_fn + "'");

      // create checkpoint control
      LAFEM::SerialConfig serial_config(false, false);
      Control::CheckpointControl check_ctrl(comm, serial_config);

      // add all three known solution vectors to checkpoint
      check_ctrl.add_object("u[k-2]", vec_sol_2.local());
      check_ctrl.add_object("u[k-1]", vec_sol_1.local());
      check_ctrl.add_object("u[k-0]", vec_sol.local());

      // set up info vector and write that one, too
      LAFEM::DenseVector<DataType, IndexType> vcp_info(Index(3));
      vcp_info(Index(0), cur_time);
      vcp_info(Index(1), delta_t);
      vcp_info(Index(2), DataType(cur_step));
      check_ctrl.add_object("info", vcp_info);

      // save checkpoint
      check_ctrl.save(check_fn);

      watch_checkpoint_save.stop();
    }

    virtual bool next_time_step()
    {
      // print runtime of previous time step
      TimeStamp stamp_next;
      if(cur_step > Index(0))
      {
        comm.print(line_prefix()
          + " > Step Runtime: " + stamp_next.elapsed_string(stamp_cur_step, TimeFormat::h_m_s_m)
          + " [ Total: " + watch_total_run.elapsed_string(TimeFormat::h_m_s_m) + " ]");
      }

      comm.print(String(120, '-'));

      stamp_cur_step = stamp_next;

      // are we done?
      if(cur_step >= max_time_steps)
        return false;

      // steady state reached?
      if((main_step > 0u) && (cur_step > 1u) && (cur_step  % main_step == 0u) && (velo_deriv_norm <= steady_tol))
      {
        comm.print("\nSteady State reached!\ndt Velo = " + stringify_fp_sci(velo_deriv_norm) + " <= " + stringify_fp_sci(steady_tol) + "\n");
        return false;
      }

      // advance time step
      ++cur_step;

      // compute simulation time
      cur_time = time_max * DataType(cur_step) / DataType(max_time_steps);

      // print header line
      //         "     1:  0.10000000 |  1: 4.032408e-03 / 2.160869e-01 / 2.161e-01 |   4: 8.38977e-06 / 2.08058e-03 [1.88287e-05]"
      comm.print("  Step   Time       |  #  Defect (abs)   Defect (rel)   Improve   |  MG  fin abs Def   fin rel Def  abs Tol");

      // keep going
      return true;
    }

    virtual void initialize_step_solution()
    {
      // get a temporary vector
      GlobalSystemVector& vec_sol_3 = vec_tmp;

      // extrapolate initial u_{k} from u_{k-1}, u_{k-2} and u_{k-3} for this time-step
      if(cur_step > Index(1))
      {
        // shift solution vectors
        vec_sol_3.copy(vec_sol_2); // u_{k-3}
        vec_sol_2.copy(vec_sol_1); // u_{k-2}
        vec_sol_1.copy(vec_sol);   // u_{k-1}

        // extrapolate to current time step
        if((time_expo >= Index(2)) && (cur_step > Index(2)))
        {
          // perform quadratic extrapolation of solution to current time-step
          // u_{k} := 3*u_{k-1} - 3*u_{k-2}) + u_{k-3}
          vec_sol.copy(vec_sol_3);
          vec_sol.axpy(vec_sol_1, DataType(3));
          vec_sol.axpy(vec_sol_2, -DataType(3));
        }
        else if((time_expo >= Index(1)) && (cur_step > Index(1)))
        {
          // perform linear extrapolation of solution to current time-step
          // u_{k} := 2*u_{k-1} - u_{k-2}
          vec_sol.scale(vec_sol_1, DataType(2));
          vec_sol.axpy(vec_sol_2, -DataType(1));
        }
        // else-case is constant extrapolation, which is done already by copy
      }
      else
      {
        vec_sol_1.copy(vec_sol);
      }

      // apply current solution filter
      system.front()->filter_sys.filter_sol(vec_sol);
    }

    virtual void assemble_step_rhs()
    {
      /// \todo: figure out what to do for inhom. RHS: type 0 vs type 1
      if(homogeneous_rhs)
        vec_rhs.format();

      // compute right hand side for this time-step
      // note: we keep vec_rhs as a type-0 vector, because it will be
      // converted to type-1 only as a defect vector in the nonlinear loop below
      if(cur_step == Index(1))
      {
        // first time step ==> implicit Euler
        // f_k := 1/dt * M * u_{k-1}
        system.front()->velo_mass_matrix.local().apply(
        //the_system_level.local_velo_mass_matrix.apply(
          vec_rhs.local().template at<0>(),
          vec_sol_1.local().template at<0>());
        vec_rhs.scale(vec_rhs, DataType(1) / delta_t);
      }
      else
      {
        // we're beyond the first time step ==> BDF(2)
        // f_k := 3/(2*dt) * (4/3 * M * u_{k-1} - 1/3 M * u_{k-2}
        //      = -1/(2*dt) * M * (u_{k-2} - 4*u_{k-1})
        vec_tmp.copy(vec_sol_2);
        vec_tmp.axpy(vec_sol_1, -DataType(4));
        system.front()->velo_mass_matrix.local().apply(
          //the_system_level.local_velo_mass_matrix.apply(
          vec_rhs.local().template at<0>(),
          vec_tmp.local().template at<0>());
        vec_rhs.scale(vec_rhs, -DataType(0.5) / delta_t);
      }
      //vec_rhs.sync_0();
    }

    virtual void initialize_nonlinear_defect() override
    {
      /// \todo: figure out what to do for inhom. RHS: type 0 vs type 1
      vec_def.copy(vec_rhs);
    }

    virtual void setup_burgers_defect_job(Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType>& burgers_def_job) override
    {
      BaseClass::setup_burgers_defect_job(burgers_def_job);
      if(cur_step == Index(1))
        burgers_def_job.theta = -DataType(1) / delta_t; // implicit Euler in first time step
      else
        burgers_def_job.theta = -DataType(1.5) / delta_t; // BDF(2) in all further time steps
    }

    virtual void setup_burgers_matrix_job(Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>& burgers_mat_job) override
    {
      BaseClass::setup_burgers_matrix_job(burgers_mat_job);
      if(cur_step == Index(1))
        burgers_mat_job.theta = DataType(1) / delta_t; // implicit Euler in first time step
      else
        burgers_mat_job.theta = DataType(1.5) / delta_t; // BDF(2) in all further time steps
    }

    virtual bool solve_time_step()
    {
      return solve_nonlinear_system(line_prefix() + " | ");
    }

    virtual void analyze_time_derivative()
    {
      if((main_step == 0u) || (cur_step  % main_step > 0u))
        return;

      watch_sol_analysis.start();

      // compute time derivative vector
      vec_der.copy(vec_sol);
      vec_der.axpy(vec_sol_1, -DataType(1));
      vec_der.scale(vec_der, DataType(1) / delta_t);

      // compute time derivative norm
      system.front()->velo_mass_matrix.local().apply(vec_tmp.local().first(), vec_der.local().first());
      velo_deriv_norm = Math::sqrt(system.front()->gate_sys.sum(vec_tmp.local().first().dot(vec_der.local().first())));
      comm.print(line_prefix() + " > Velocity Time-Der Norm: " + stringify_fp_sci(velo_deriv_norm));

      watch_sol_analysis.stop();
    }

    virtual void write_vtk() override
    {
      if(!want_vtk)
        return;

      if((vtk_step > 1u) && (cur_step % vtk_step != 0))
        return;

      watch_vtk.start();

      // build VTK name
      String vtk_name = vtk_filename;
      if(vtk_filename.empty())
      {
        vtk_name = default_filename;
        vtk_name += "-lvl" + stringify(domain.front()->get_level_index());
        vtk_name += "-n" + stringify(comm.size());
      }
      // append time step
      vtk_name += "." + stringify(cur_step).pad_front(std::size_t(Math::max(5, int(Math::ilog10(max_time_steps)))), '0');

      comm.print(line_prefix() + " > Writing VTK output to '" + vtk_name + ".pvtu'");

      // get the mesh for the VTK output
      const MeshType* mesh = &domain.front()->get_mesh();
      if(refined_mesh_vtk)
        mesh = refined_mesh_vtk.get();

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(*mesh);

      // write velocity
      exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());
      exporter.add_vertex_vector("velo_der", vec_der.local().template at<0>());

      // project pressure
      LAFEM::DenseVector<DataType, Index> vtx_p, vtx_dp;
      if(refined_mesh_vtk)
      {
        Assembly::DiscreteCellProjector::project_refined(vtx_p, vec_sol.local().template at<1>(), domain.front()->space_pres);
        Assembly::DiscreteCellProjector::project_refined(vtx_dp, vec_der.local().template at<1>(), domain.front()->space_pres);
      }
      else
      {
        Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), domain.front()->space_pres, "gauss-legendre:2");
        Assembly::DiscreteCellProjector::project(vtx_dp, vec_der.local().template at<1>(), domain.front()->space_pres, "gauss-legendre:2");
      }

      // write pressure
      exporter.add_cell_scalar("pressure", vtx_p.elements());
      exporter.add_cell_scalar("pres_der", vtx_dp.elements());

      // finally, write the VTK file
      exporter.write(vtk_name, comm);

      watch_vtk.stop();
    }

    virtual void collect_final_statistics() override
    {
      stats.times[Times::checkpoint_save] = watch_checkpoint_save.elapsed();
      stats.times[Times::checkpoint_load] = watch_checkpoint_load.elapsed();
      BaseClass::collect_final_statistics();
    }
  }; // class UnsteadyAppBase<...>
} // namespace CCND
