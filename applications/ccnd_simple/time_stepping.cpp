// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "time_stepping.hpp"
#include <control/checkpoint_control.hpp>

namespace CCNDSimple
{
  TimeStepping::TimeStepping(const Dist::Comm& comm_) :
    comm(comm_)
  {
  }

  TimeStepping::~TimeStepping()
  {
  }

  void TimeStepping::add_supported_args(SimpleArgParser& args)
  {
    args.support("full-plot", "\nEnables the full unsteady solver plot.");
    args.support("t-max", "<T>\nSets the maximum simulation time T.");
    args.support("delta-t", "<dt>\nSets the time step size delta_t.");
    args.support("bdf-type", "<k>\nSets the type of the BDF(k) time stepping scheme; must be 1 or 2.");
    args.support("sol-expo", "<k>\nSets the solution vector time extrapolation order; must be 0, 1 or 2.");
  }

  bool TimeStepping::parse_args(SimpleArgParser& args)
  {
    // full plot?
    full_plot = (args.check("full-plot") >= 0);

    // viscosity parameter nu
    args.parse("t-max", t_max);
    args.parse("delta-t", delta_t);
    args.parse("bdf-type", bdf_type);
    args.parse("sol-expo", sol_expo);

    return true;
  }

  void TimeStepping::print_config() const
  {
    print_pad(comm, "T-Max", stringify_fp_fix(t_max));
    print_pad(comm, "delta-t", stringify_fp_fix(delta_t));
    print_pad(comm, "BDF-Type", stringify(bdf_type));
    print_pad(comm, "Solution Extrapolation Order", stringify(sol_expo));
  }

  void TimeStepping::begin_loop()
  {
    watch_loop.start();
  }

  void TimeStepping::finish_loop()
  {
    watch_loop.stop();
  }

  bool TimeStepping::advance_step()
  {
    // update time step and simulation time
    ++time_step;
    cur_time += delta_t;

    // are we beyond the maximum simulation time?
    if(t_max + 1E+10 < cur_time)
    {
      comm.print("\nMaximum simulation time " + stringify_fp_fix(t_max) + " reached; finishing time step loop");
      comm.print_flush();
      return false;
    }

    // compute theta
    theta.format();
    if((time_step == Index(1)) || (bdf_type < Index(2)))
    {
      // first time step or BDF(1) ==> implicit Euler
      theta[0] = DataType(1) / delta_t; // LHS: 1/dt * M * u_k
      theta[1] = DataType(1) / delta_t; // RHS: 1/dt * M * u_{k-1}
    }
    else
    {
      // we're beyond the first time step ==> BDF(2)
      // f_k := 3/(2*dt) * (4/3 * M * u_{k-1} - 1/3 M * u_{k-2}
      //      = -1/(2*dt) * M * (u_{k-2} - 4*u_{k-1})
      theta.format();
      theta[0] =  DataType(1.5) / delta_t; // LHS:  3/(2*dt) * M * u_k
      theta[1] =  DataType(2.0) / delta_t; // RHS:  2/dt     * M * u_{k-1}
      theta[2] = -DataType(0.5) / delta_t; // RHS: -1/(2*dt) * M * u_{k-2}
    }

    // compute extrapolation coefficients
    expolc.format();
    if((time_step > Index(2)) && (sol_expo >= Index(2)))
    {
      // perform quadratic extrapolation of solution to current time-step
      // u_{k} := 3*u_{k-1} - 3*u_{k-2}) + u_{k-3}
      expolc[1] =  DataType(3);
      expolc[2] = -DataType(3);
      expolc[3] =  DataType(1);
    }
    else if((time_step > Index(1)) && (sol_expo >= Index(1)))
    {
      // perform linear extrapolation of solution to current time-step
      // u_{k} := 2*u_{k-1} - u_{k-2}
      expolc[1] =  DataType(2);
      expolc[2] = -DataType(1);
    }
    else
    {
      // perform constant extrapolation
      expolc[1] = DataType(1);
    }

    // write output line
    if(full_plot)
    {
      comm.print(String(120u, '#'));
      comm.print("Step " + stringify(time_step).pad_front(6) +
        "   Time " + stringify_fp_fix(cur_time, 10, 15) +
        "   [" + stringify_fp_fix(100.0*cur_time/t_max, 2, 6) + "%]" +
        "   Runtime " + watch_loop.elapsed_string().pad_front(12));
      comm.print(String(120u, '-') + "\n");
    }
    else
    {
      plot_line.clear();
      plot_line += stringify(time_step).pad_front(8);
      plot_line += stringify_fp_fix(cur_time, 10, 18);
      plot_line += " [" + stringify_fp_fix(100.0*cur_time/t_max, 2, 6) + "%]";
      plot_line += watch_loop.elapsed_string().pad_front(12);
    }

    // check for stop
    return (cur_time <= t_max + 1e-10); // add some tolerance for rounding errors
  }

  void TimeStepping::print_runtime(double total_time)
  {
    print_time(comm, "Total Time Loop Time", watch_loop.elapsed(), total_time);
  }
} // namespace CCNDSimple
