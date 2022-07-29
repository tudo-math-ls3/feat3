// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "unsteady_solver.hpp"
#include <control/checkpoint_control.hpp>

#include <filesystem>

namespace CCNDSimple
{
  UnsteadySolver::UnsteadySolver(DomainControl& domain_) :
    BaseClass(domain_)
  {
    // disable navier-stokes plot header
    plot_navier_header = false;
  }

  UnsteadySolver::~UnsteadySolver()
  {
  }

  void UnsteadySolver::add_supported_args(SimpleArgParser& args)
  {
    BaseClass::add_supported_args(args);

    args.support("full-plot", "Enables the full unsteady solver plot.");
    args.support("t-max", "<T>\nSets the maximum simulation time T.");
    args.support("delta-t", "<dt>\nSets the time step size delta_t.");
    args.support("bdf-type", "<k>\nSets the type of the BDF(k) time stepping scheme; must be 1 or 2.");
    args.support("sol-expo", "<k>\nSets the solution vector time extrapolation order; must be 0, 1 or 2.");
    args.support("checkpoint", "<filename> [<step>] [<modulo>]\nWrites a checkpoint every <step> timesteps; up to <modulo> total files.");
    args.support("restart", "<filename>\nRestart from a previously saved checkpoint file.");
  }

  bool UnsteadySolver::parse_args(SimpleArgParser& args)
  {
    if(!BaseClass::parse_args(args))
      return false;

    // full plot?
    full_plot = (args.check("full-plot") >= 0);
    if(!full_plot)
      this->plot_navier = false;

    // viscosity parameter nu
    args.parse("t-max", t_max);
    args.parse("delta-t", delta_t);
    args.parse("bdf-type", bdf_type);
    args.parse("sol-expo", sol_expo);
    args.parse("checkpoint", check_name, check_step, check_mod);
    args.parse("restart", restart_name);

    return true;
  }

  void UnsteadySolver::print_config() const
  {
    print_pad(comm, "T-Max", stringify_fp_fix(t_max));
    print_pad(comm, "delta-t", stringify_fp_fix(delta_t));
    print_pad(comm, "BDF-Type", stringify(bdf_type));
    print_pad(comm, "Solution Extrapol Order", stringify(sol_expo));
    if(restart_name.empty())
      print_pad(comm, "Checkpoint Restart", "- N/A -");
    else
      print_pad(comm, "Checkpoint Restart", restart_name);
    if(check_name.empty())
      print_pad(comm, "Checkpoint Name", "- N/A -");
    else
      print_pad(comm, "Checkpoint Name", check_name);
    print_pad(comm, "Checkpoint Steps", stringify(check_step));
    print_pad(comm, "Checkpoint Modulo", stringify(check_mod));
    BaseClass::print_config();
  }

  void UnsteadySolver::create_levels()
  {
    BaseClass::create_levels();

    // assemble velocity mass matrix on finest level
    this->stokes_levels.front()->assemble_velocity_mass_matrix(cubature);
  }

  void UnsteadySolver::write_checkpoint(const String& filename)
  {
    // create checkpoint control
    LAFEM::SerialConfig serial_config(false, false);
    Control::CheckpointControl check_ctrl(comm, serial_config);

    // add all three known solution vectors to checkpoint
    check_ctrl.add_object("u[k-2]", vec_sol_3.local());
    check_ctrl.add_object("u[k-1]", vec_sol_2.local());
    check_ctrl.add_object("u[k-0]", vec_sol_1.local());

    // set up info vector and write that one, too
    LAFEM::DenseVector<DataType, IndexType> vcp_info(Index(3));
    vcp_info(Index(0), cur_time);
    vcp_info(Index(1), delta_t);
    vcp_info(Index(2), DataType(time_step));
    check_ctrl.add_object("info", vcp_info);

    // save checkpoint
    check_ctrl.save(filename);
  }

  bool UnsteadySolver::checkpoint_restart(GlobalStokesVector& vec_sol)
  {
    // no restart file?
    if(restart_name.empty())
      return false;

    comm.print("\nReading restart checkpoint '" + restart_name + "'...");

    // create checkpoint control
    LAFEM::SerialConfig serial_config(false, false);
    Control::CheckpointControl check_ctrl(comm, serial_config);

    // create checkpoint info vector
    LAFEM::DenseVector<DataType, IndexType> vcp_info;

    // load checkpoint
    check_ctrl.load(restart_name);

    // restore all three solution vectors from checkpoint
    check_ctrl.restore_object("u[k-2]", vec_sol_3.local());
    check_ctrl.restore_object("u[k-1]", vec_sol_2.local());
    check_ctrl.restore_object("u[k-0]", vec_sol_1.local());
    check_ctrl.restore_object("info", vcp_info);

    // ensure that the dimensions match
    XASSERTM(vec_sol_1.local().at<0>().size() == vec_sol.local().at<0>().size(), "invalid vector size");
    XASSERTM(vec_sol_1.local().at<1>().size() == vec_sol.local().at<1>().size(), "invalid vector size");
    vec_sol.copy(vec_sol_1);

    // extract old checkpoint info
    XASSERTM(vcp_info.size() == Index(3), "invalid info vector size");
    DataType cp_cur_time = vcp_info(Index(0));
    DataType cp_delta_t  = vcp_info(Index(1));
    Index cp_time_step   = Index(vcp_info(Index(2)));

    // choose the restart time based on parameters
    delta_t = cp_delta_t;
    time_step = cp_time_step;
    cur_time = cp_cur_time;

    // print checkpoint and restart info
    comm.print("Checkpoint Delta-T.: " + stringify_fp_fix(cp_delta_t));
    comm.print("Checkpoint Time....: " + stringify_fp_fix(cp_cur_time));
    comm.print("Checkpoint Timestep: " + stringify(cp_time_step));

    // we're restarting
    return true;
  }

  bool UnsteadySolver::start_time_loop(GlobalStokesVector& vec_sol)
  {
    // check for 'STOP' file
    if(std::filesystem::exists("./STOP"))
    {
      //  STOP file found
      comm.print("\nSTOP file found; delete the file 'STOP' and restart program");
      return false;
    }

    // make sure the vectors are allocated
    if(vec_sol_3.local().first().empty())
      vec_sol_3 = vec_sol.clone(LAFEM::CloneMode::Deep);
    if(vec_sol_2.local().first().empty())
      vec_sol_2 = vec_sol.clone(LAFEM::CloneMode::Deep);
    if(vec_sol_1.local().first().empty())
      vec_sol_1 = vec_sol.clone(LAFEM::CloneMode::Deep);
    if(vec_tmp.local().first().empty())
      vec_tmp = vec_sol.clone(LAFEM::CloneMode::Layout);

    // start the loop watch
    watch_loop.start();

    // print time-loop header
    if(!full_plot)
    {
      comm.print("   Step     Time             Done   NL  Def-Init    Def-Final     MG       Runtime" );
      comm.print("-----------------------------------------------------------------------------------");
    }

    // let's go
    return true;
  }

  bool UnsteadySolver::begin_step()
  {
    // update time step and simulation time
    ++time_step;
    cur_time += delta_t;

    // clear statistics
    this->mg_iters.clear();
    this->nl_defs.clear();

    // write output line
    if(full_plot)
    {
      comm.print(String(100u, '-') +
        "Step " + stringify(time_step).pad_front(6) +
        "   Time " + stringify_fp_fix(cur_time, 10, 15) +
        "   [" + stringify_fp_fix(100.0*cur_time/t_max, 2, 6) + "%]" +
        "   Runtime " + watch_loop.elapsed_string().pad_front(12) +
        "\n");
    }

    // check for stop
    return (cur_time <= t_max + 1e-10); // add some tolerance for rounding errors
  }

  bool UnsteadySolver::finish_step(const GlobalStokesVector& vec_sol)
  {
    // save the current vector
    vec_sol_3.copy(vec_sol_2); // u_{k-3}
    vec_sol_2.copy(vec_sol_1); // u_{k-2}
    vec_sol_1.copy(vec_sol);   // u_{k-1}

    String line;
    if(!full_plot)
    {
      line += stringify(time_step).pad_front(7);
      line += stringify_fp_fix(cur_time, 10, 17);
      line += " [" + stringify_fp_fix(100.0*cur_time/t_max, 2, 6) + "%]:";
      line += stringify(Math::max(int(this->nl_defs.size()),1)-1).pad_front(3) + ": ";
      line += stringify_fp_sci(this->nl_defs.front(), 3) + " > ";
      line += stringify_fp_sci(this->nl_defs.back(), 3) + " : ";
      Index mgi(0u);
      for(auto& i : this->mg_iters)
        mgi += i;
      line += stringify(mgi).pad_front(4) + " | ";
      line += "[" + watch_loop.elapsed_string().pad_front(10) + "]";
    }

    // write checkpoint?
    if((check_step > 0u) && (time_step % check_step == 0u))
    {
      // generate checkpoint name
      String check_fn = check_name + "." + stringify(((time_step-1) / check_step) % check_mod) + ".cp";
      if(full_plot)
        comm.print("\nWriting Checkpoint '" + check_fn + "'...");
      else
        line += "; CP > " + check_name;

      write_checkpoint(check_fn);
    }

    // print short summary?
    if(!full_plot)
      comm.print(line);

    // was this the last time step?
    if(t_max < cur_time + 1E-10)
    {
      comm.print("\nMaximum simulation time " + stringify_fp_fix(t_max) + " reached; finishing time step loop");
      return false;
    }

    // check for 'STOP' file
    if(!std::filesystem::exists("./STOP"))
      return true; // continue iteration

    //  STOP file found
    comm.print("\nSTOP file found; stopping time step loop; delete the file 'STOP' and restart program");

    // read the first line of the file to obtain a checkpoint filename
    String check_fn;
    {
      std::stringstream sstr;
      DistFileIO::read_common(sstr, "./STOP", comm);
      std::getline(sstr, check_fn);
    }
    if(check_fn.trim_me().empty())
      return false;
    if(!check_fn.ends_with(".cp") && !check_fn.ends_with(".zcp"))
      check_fn += ".cp";

    // okay, write a checkpoint with that filename
    comm.print("Writing Checkpoint '" + check_fn + "' as requested by the STOP file");
    write_checkpoint(check_fn);

    // abort time step loop
    return false;
  }

  void UnsteadySolver::initialize_current_sol(GlobalStokesVector& vec_sol)
  {
    if((time_step > Index(2)) && (sol_expo >= Index(2)))
    {
      // perform quadratic extrapolation of solution to current time-step
      // u_{k} := 3*u_{k-1} - 3*u_{k-2}) + u_{k-3}
      vec_sol.axpy(vec_sol_1, vec_sol_3, DataType(3));
      vec_sol.axpy(vec_sol_2, vec_sol, -DataType(3));
    }
    else if((time_step > Index(1)) && (sol_expo >= Index(1)))
    {
      // perform linear extrapolation of solution to current time-step
      // u_{k} := 2*u_{k-1} - u_{k-2}
      vec_sol.scale(vec_sol_1, DataType(2));
      vec_sol.axpy(vec_sol_2, vec_sol, -DataType(1));
    }
    // else: leave vec_sol untouched, which results in constant extrapolation

    // apply filter
    stokes_levels.front()->filter_sys.filter_sol(vec_sol);
  }

  void UnsteadySolver::assemble_current_rhs(GlobalStokesVector& vec_rhs)
  {
    // get the finest system level
    auto& the_system_level = *this->stokes_levels.front();

    // convert RHS from type-1 to type-0
    vec_rhs.from_1_to_0();

    if((time_step == Index(1)) || (bdf_type < Index(2)))
    {
      // first time step ==> implicit Euler
      // f_k := 1/dt * M * u_{k-1}
      the_system_level.local_velo_mass_matrix.apply(
        vec_rhs.local().template at<0>(),
        vec_sol_1.local().template at<0>(),
        vec_rhs.local().template at<0>(),
        DataType(1) / delta_t);
    }
    else
    {
      // we're beyond the first time step ==> BDF(2)
      // f_k := 3/(2*dt) * (4/3 * M * u_{k-1} - 1/3 M * u_{k-2}
      //      = -1/(2*dt) * M * (u_{k-2} - 4*u_{k-1})
      vec_tmp.axpy(vec_sol_1, vec_sol_2, -DataType(4));
      the_system_level.local_velo_mass_matrix.apply(
        vec_rhs.local().template at<0>(),
        vec_tmp.local().template at<0>(),
        vec_rhs.local().template at<0>(),
        -DataType(0.5) / delta_t);
    }

    // sync to obtain a type-1 vector
    vec_rhs.sync_0();

    // apply filter
    stokes_levels.front()->filter_sys.filter_rhs(vec_rhs);
  }

  void UnsteadySolver::setup_defect_burgers_job(Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType>& burgers_def_job)
  {
    burgers_def_job.deformation = deform_tensor;
    burgers_def_job.nu = -nu;
    burgers_def_job.beta = -DataType(1);
    if((time_step == Index(1)) || (bdf_type < Index(2)))
      burgers_def_job.theta = -DataType(1) / delta_t; // implicit Euler in first time step
    else
      burgers_def_job.theta = -DataType(1.5) / delta_t; // BDF(2) in all further time steps
  }

  void UnsteadySolver::setup_matrix_burgers_job(Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>& burgers_mat_job)
  {
    burgers_mat_job.deformation = deform_tensor;
    burgers_mat_job.nu = nu;
    burgers_mat_job.beta = DataType(1);
    burgers_mat_job.frechet_beta = DataType(newton_solver ? 1 : 0);
    burgers_mat_job.sd_delta = upsam;
    burgers_mat_job.sd_nu = nu;
    if((time_step == Index(1)) || (bdf_type < Index(2)))
      burgers_mat_job.theta = DataType(1) / delta_t; // implicit Euler in first time step
    else
      burgers_mat_job.theta = DataType(1.5) / delta_t; // BDF(2) in all further time steps
  }
} // namespace CCNDSimple
