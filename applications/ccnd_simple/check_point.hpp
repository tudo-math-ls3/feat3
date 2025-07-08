// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This header/source file pair defines the CheckPoint class, which can be used by long running applications to save
// the current simulation state to a checkpoint file (aka "savegame") that can be loaded at another time to continue
// the simulation from the state that was stoked in the checkpoint.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "base.hpp"
#include "stokes_level.hpp"
#include "time_stepping.hpp"

namespace CCNDSimple
{
  /**
   * \brief Checkpoint class for CCND apps
   *
   * This class implements a checkpointing system that can be used to write checkpoints periodically during a long
   * running simulation and it enables the user to continue a previous simulation from a checkpoint in the case that
   * the simulation aborts prematurely for whatever reason.
   *
   * \author Peter Zajac
   */
  class CheckPoint
  {
  public:
    /// our communicator
    const Dist::Comm& comm;

    // ----------------
    // input attributes
    // ----------------

    /// desired checkpointing time steps; set to > 0 to enable checkpointing
    IndexType check_step = IndexType(0);

    /// desired number of checkpoints to keep; should be > 1
    IndexType check_mod = IndexType(2);

    /// checkpoint name
    String checkpoint_name;

    /// restart filename
    String continue_name;

    // -----------------
    // output attributes
    // -----------------

    /// plot line for short plot
    String plot_line;

    /// watches for checkpoint save/load
    StopWatch watch_save, watch_load;

    // ----------------
    // state attributes
    // ----------------

    /// a map of Stokes vectors
    std::map<String, GlobalStokesVector*> stokes_vectors;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  public:
    /// constructor
    explicit CheckPoint(const Dist::Comm& comm_);

    /// destructor
    virtual ~CheckPoint();

    /// non-copyable
    CheckPoint(const CheckPoint&) = delete;

    /// non-copyable
    CheckPoint& operator=(const CheckPoint&) = delete;

    /// adds all supported arguments to the argument parser
    static void add_supported_args(SimpleArgParser& args);

    /// parses arguments from command line
    virtual bool parse_args(SimpleArgParser& args);

    /// prints configuration to console
    virtual void print_config() const;

    /// adds a stokes vector reference to the checkpoint system
    virtual bool register_stokes_vector(const String& name, GlobalStokesVector& vector);

    /// removes all previously registered Stokes vectors from the checkpoint system
    virtual void unregister_stokes_vectors();

    /**
     * \brief Writes the current simulation state to a checkpoint file
     *
     * This functions writes all registered vectors to a checkpoint file if the user has requested
     * that checkpoints are to be saved and if the current time step is a checkpoint time step.
     *
     * \param[in] time_stepping
     * The time stepping object that represents the current time step.
     */
    virtual bool write_registerd(const TimeStepping& time_stepping);

    /**
     * \brief Restarts the simulation from a previously saved checkpoint; if desired
     *
     * If the user has requested a restart by supplying the --restart <filename> command line option,
     * then this function will read in that checkpoint and configure the solver to continue the
     * simulation from that checkpoint instead of starting the simulation from scratch.
     *
     * \attention
     * This function will override the time step, current simulation time and time step size to
     * the values stored in the checkpoint file.
     *
     * \param[inout] time_stepping
     * The time stepping object that represents the next time step.
     */
    virtual bool read_registered(TimeStepping& time_stepping);

    /// prints runtime information for this object
    virtual void print_runtime(double total_time);

  }; // class CheckPoint
} // namespace CCNDSimple
