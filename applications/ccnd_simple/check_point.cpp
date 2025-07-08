// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "check_point.hpp"

#include <control/checkpoint_control.hpp>

namespace CCNDSimple
{
  CheckPoint::CheckPoint(const Dist::Comm& comm_) :
    comm(comm_)
  {
  }

  CheckPoint::~CheckPoint()
  {
  }

  void CheckPoint::add_supported_args(SimpleArgParser& args)
  {
    args.support("checkpoint", "<filename> [<step>] [<modulo>]\nWrites a checkpoint every <step> time steps; up to <modulo> total files.");
    args.support("continue", "<filename>\nContinue from a previously saved checkpoint file.");
  }

  bool CheckPoint::parse_args(SimpleArgParser& args)
  {
    args.parse("checkpoint", checkpoint_name, check_step, check_mod);
    args.parse("continue", continue_name);

    return true;
  }

  void CheckPoint::print_config() const
  {
    if(continue_name.empty())
      print_pad(comm, "Checkpoint Restart", "- N/A -");
    else
      print_pad(comm, "Checkpoint Restart", continue_name);
    if(checkpoint_name.empty())
      print_pad(comm, "Checkpoint Name", "- N/A -");
    else
      print_pad(comm, "Checkpoint Name", checkpoint_name);
    print_pad(comm, "Checkpoint Steps", stringify(check_step));
    print_pad(comm, "Checkpoint Modulo", stringify(check_mod));
  }

  bool CheckPoint::register_stokes_vector(const String& name, GlobalStokesVector& vector)
  {
    return stokes_vectors.emplace(name, &vector).second;
  }

  void CheckPoint::unregister_stokes_vectors()
  {
    stokes_vectors.clear();
  }

  bool CheckPoint::write_registerd(const TimeStepping& time_stepping)
  {
    plot_line.clear();

    // save checkpoint this time step?
    if((check_step == 0u) || (time_stepping.time_step % check_step != 0u))
      return false;

    watch_save.start();

    // generate checkpoint name
    String check_fn = checkpoint_name + "." + stringify(((time_stepping.time_step-1) / check_step) % check_mod) + ".cp";
    if(time_stepping.full_plot)
      comm.print("\nSaving Checkpoint '" + check_fn + "'...");

    plot_line = " ; > '" + check_fn + "'";

    // create checkpoint control
    LAFEM::SerialConfig serial_config(false, false);
    Control::CheckpointControl check_ctrl(comm, serial_config);

    // set up and write info vector
    LAFEM::DenseVector<DataType, IndexType> vcp_info(Index(3));
    vcp_info(Index(0), time_stepping.cur_time);
    vcp_info(Index(1), time_stepping.delta_t);
    vcp_info(Index(2), DataType(time_stepping.time_step));
    check_ctrl.add_object("{{info}}", vcp_info);

    // add Stokes vectors
    for(const auto& v : stokes_vectors)
      check_ctrl.add_object(v.first, v.second->local());

    // save checkpoint
    check_ctrl.save(check_fn);
    watch_save.stop();
    return true;
  }

  bool CheckPoint::read_registered(TimeStepping& time_stepping)
  {
    // no restart file?
    if(continue_name.empty())
      return false;

    watch_load.start();

    comm.print("\nLoading checkpoint '" + continue_name + "'...");

    // create checkpoint control
    LAFEM::SerialConfig serial_config(false, false);
    Control::CheckpointControl check_ctrl(comm, serial_config);

    // create checkpoint info vector
    LAFEM::DenseVector<DataType, IndexType> vcp_info;

    // load checkpoint
    check_ctrl.load(continue_name);

    // restore info object
    check_ctrl.restore_object("{{info}}", vcp_info);

    // restore all three solution vectors from checkpoint
    for(auto& v : stokes_vectors)
      check_ctrl.restore_object(v.first, v.second->local());

    // ensure that the dimensions match
    //XASSERTM(vec_sol_1.local().at<0>().size() == vec_sol.local().at<0>().size(), "invalid vector size");
    //XASSERTM(vec_sol_1.local().at<1>().size() == vec_sol.local().at<1>().size(), "invalid vector size");
    //vec_sol.copy(vec_sol_1);

    // extract old checkpoint info
    XASSERTM(vcp_info.size() == Index(3), "invalid info vector size");
    DataType cp_cur_time = vcp_info(Index(0));
    DataType cp_delta_t  = vcp_info(Index(1));
    Index cp_time_step   = Index(vcp_info(Index(2)));

    // choose the restart time based on parameters
    time_stepping.delta_t = cp_delta_t;
    time_stepping.time_step = cp_time_step;
    time_stepping.cur_time = cp_cur_time;

    // print checkpoint and restart info
    print_pad(comm, "Checkpoint Delta-T", stringify_fp_fix(cp_delta_t));
    print_pad(comm, "Checkpoint Time", stringify_fp_fix(cp_cur_time));
    print_pad(comm, "Checkpoint Step", stringify(cp_time_step));

    watch_load.stop();

    // we're restarting
    return true;
  }

  void CheckPoint::print_runtime(double total_time)
  {
    print_time(comm, "CheckPoint Load Time", watch_load.elapsed(), total_time);
    print_time(comm, "CheckPoint Save Time", watch_save.elapsed(), total_time);
  }
} // namespace CCNDSimple
