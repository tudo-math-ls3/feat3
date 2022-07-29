// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef FEAT_APPLICATIONS_CCND_SIMPLE_UNSTEADY_SOLVER_HPP
#define FEAT_APPLICATIONS_CCND_SIMPLE_UNSTEADY_SOLVER_HPP 1

#include "steady_solver.hpp"

namespace CCNDSimple
{
  /**
   * \brief Unsteady Navier-Stokes solver class
   *
   * This class provides a simple monolithic unsteady Navier-Stokes solver.
   *
   * This class derives from the SteadySolver class, thus inheriting its functionality to solve
   * the nonlinear system in each time step.
   *
   * \author Peter Zajac
   */
  class UnsteadySolver :
    public SteadySolver
  {
  public:
    /// our base-class typedef
    typedef SteadySolver BaseClass;

    // -------------
    // configuration
    // -------------

    /// the maximum simulation time
    DataType t_max = DataType(1);
    /// the time-step size
    DataType delta_t = DataType(0.1);
    /// desired BDF time stepping scheme; either 1 or 2
    IndexType bdf_type = IndexType(2);
    /// desired solution time-extrapolation for each time-step; either 0, 1 or 2
    IndexType sol_expo = IndexType(2);
    /// desired checkpointing time steps; set to > 0 to enable checkpointing
    IndexType check_step = IndexType(0);
    /// desired number of checkpoints to keep; should be > 1
    IndexType check_mod = IndexType(2);
    /// checkpoint name
    String check_name;
    /// restart filename
    String restart_name;

    /// full plot or short
    bool full_plot = false;

    // -------------------
    // time-stepping state
    // -------------------

    /// the current time step
    IndexType time_step = IndexType(0);
    /// the current simulation time = time_step * delta_t
    DataType cur_time = DataType(0);

    /// three temporary vectors for the previous time step solutions
    GlobalStokesVector vec_tmp, vec_sol_1, vec_sol_2, vec_sol_3;

    /// watch for current run-time
    StopWatch watch_loop;

  public:
    /// constructor
    explicit UnsteadySolver(DomainControl& domain_);

    /// destructor
    virtual ~UnsteadySolver();

    /// non-copyable
    UnsteadySolver(const UnsteadySolver&) = delete;

    /// non-copyable
    UnsteadySolver& operator=(const UnsteadySolver&) = delete;

    /// adds all supported arguments to the argument parser
    static void add_supported_args(SimpleArgParser& args);

    /// parses arguments from command line
    virtual bool parse_args(SimpleArgParser& args) override;

    /// prints configuration to console
    virtual void print_config() const override;

    /// creates the actual levels along with gates, muxers and transfers
    virtual void create_levels() override;

    /**
     * \brief Writes the current simulation state to a checkpoint file
     *
     * \param[in] filename
     * The checkpoint filename
     */
    virtual void write_checkpoint(const String& filename);

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
     * \param[inout] vec_sol
     * A \transient reference to the solution vector that is to be restored from a checkpoint.
     *
     * \returns \c true, if the user has asked to restart from a checkpoint or \c false, if the
     * simulation is to be started from scratch.
     */
    virtual bool checkpoint_restart(GlobalStokesVector& vec_sol);

    /**
     * \brief Starts the time stepping loop
     *
     * \param[in] vec_sol
     * A \transient reference to the initial solution vector.
     */
    virtual bool start_time_loop(GlobalStokesVector& vec_sol);

    /**
     * \brief Advances the time stepping to the next iteration
     *
     * This function should be called in the while-statement of the time stepping loop and its
     * return value should be used as a stopping criterion for the while loop.
     *
     * \returns \c true, if the next time step is to be performed, or \c false, if the time step
     * loop has reached its end.
     */
    virtual bool begin_step();

    /**
     * \brief Finishes the current time step
     *
     * This function saves the current time step solution internally and performs any other
     * finalization tasks for the current time step, which may include writing a checkpoint,
     * if the current time step is a checkpoint step.
     *
     * This function also checks whether a file named \b STOP exists in the working directory of
     * the application and, if so, stops the simulation. Furthermore, if the STOP file is not
     * empty, this function will write a checkpoint to a file, whose name is given by the first
     * line of the STOP file before the time loop is exited.
     *
     * \param[in] vec_sol
     * A \transient reference to the final solution vector of the current time step.
     *
     * \returns \c true, if the next time step is to be performed, or \c false, if the time step
     * loop has reached its end.
     */
    virtual bool finish_step(const GlobalStokesVector& vec_sol);

    /**
     * \brief Initializes the solution vector for the current time step
     *
     * This function is used to store the solution vector of the previous time step and perform an
     * appropriate extrapolation to obtain a initial solution for the current time step.
     *
     * \param[inout] vec_sol
     * A \transient reference to the solution vector, which contains the solution of the previous time step.
     * Its contents are overwritten by the initial solution for the current time step.
     */
    virtual void initialize_current_sol(GlobalStokesVector& vec_sol);

    /**
     * \brief Assembles the right-hand-side vector for the current time step
     *
     * This function is used to assemble the additional RHS terms that arise from the time discretization
     * onto the RHS vector, which may either be a null vector or which may already contain the terms from
     * the 'actual' right-hand-side for this time step.
     *
     * \param[inout] vec_rhs
     * A \transient reference to the RHS vector, which contains the RHS for the current time step.
     * The additional terms that arise from the time discretization are added onto this vector.
     */
    virtual void assemble_current_rhs(GlobalStokesVector& vec_rhs);

    /**
     * \brief Sets up the Burgers assembly job for the nonlinear defect assembly.
     *
     * This function is called by assemble_nonlinear_defect().
     */
    virtual void setup_defect_burgers_job(Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType>& burgers_def_job) override;

    /**
     * \brief Sets up the Burgers assembly job for the nonlinear matrix assembly.
     *
     * This function is called by assemble_burgers_matrices().
     */
    virtual void setup_matrix_burgers_job(Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>& burgers_mat_job) override;

  }; // class UnsteadySolver
} // namespace CCNDSimple

#endif // FEAT_APPLICATIONS_CCND_SIMPLE_UNSTEADY_SOLVER_HPP
