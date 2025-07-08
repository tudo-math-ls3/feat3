// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This header/source file pair defines the TimeStepping class, which implements BDF(1) and BDF(2) time stepping
// schemes. This class controls the internal time stepping loop state and it provides the necessary coefficients
// that have to be assembled to solve the implicit time step system.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "base.hpp"

namespace CCNDSimple
{
  /**
   * \brief Time-Stepping class
   *
   * This class is responsible for computing the coefficients of a simple BDF(1) (aka "implicit Euler") or BDF(2)
   * time stepping scheme and it also manages the current time step state.
   *
   * \author Peter Zajac
   */
  class TimeStepping
  {
  public:
    /// our communicator
    const Dist::Comm& comm;

    // ----------------
    // input attributes
    // ----------------

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

    // -----------------
    // output attributes
    // -----------------

    /// plot line for short plot
    String plot_line;

    /// watch for current run-time
    StopWatch watch_loop;

    // ----------------
    // state attributes
    // ----------------

    /// the current time step
    IndexType time_step = IndexType(0);

    /// the current simulation time = time_step * delta_t
    DataType cur_time = DataType(0);

    /// thetas for the matrix and the right hand side vector components
    /// theta[0] = coefficient for system matrix
    /// theta[i] = coefficient for RHS vector of u_{k-i}
    Tiny::Vector<DataType, 3> theta;

    /// extrapolation coefficients
    /// expolc[0] = 0 (ignored)
    /// expolc[i] = coefficient for sol vector of u_{k-i}
    Tiny::Vector<DataType, 4> expolc;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  public:
    /// constructor
    explicit TimeStepping(const Dist::Comm& comm_);

    /// destructor
    virtual ~TimeStepping();

    /// non-copyable
    TimeStepping(const TimeStepping&) = delete;

    /// non-copyable
    TimeStepping& operator=(const TimeStepping&) = delete;

    /// adds all supported arguments to the argument parser
    static void add_supported_args(SimpleArgParser& args);

    /// parses arguments from command line
    virtual bool parse_args(SimpleArgParser& args);

    /// prints configuration to console
    virtual void print_config() const;

    /**
     * \brief Extrapolates a solution vector from the previous time steps
     *
     * \param[out] vec_sol
     * The vector that receives the extrapolated solution
     *
     * \param[in] vec_sol_1, vec_sol_2, vec_sol_3
     * The solution vectors of the three previous time steps
     */
    template<typename Vector_>
    void extrapolate_sol(Vector_& vec_sol, const Vector_& vec_sol_1, const Vector_& vec_sol_2, const Vector_& vec_sol_3)
    {
      // expolc[0] is ignored
      vec_sol.scale(vec_sol_1, expolc[1]);
      if(Math::abs(expolc[2]) > DataType(0))
        vec_sol.axpy(vec_sol_2, expolc[2]);
      if(Math::abs(expolc[3]) > DataType(0))
        vec_sol.axpy(vec_sol_3, expolc[3]);
    }

    /**
     * \brief Combines all terms from the previous time steps for the RHS that have to be multiplied by a mass matrix
     *
     * \param[out] vec_tmp
     * A temporary vector that receives all combined terms; this vector has to be multiplied by the mass matrix
     *
     * \param[in] vec_sol_1, vec_sol_2
     * The solution vectors of the three previous time steps
     */
    template<typename Vector_>
    void combine_mass_rhs(Vector_& vec_tmp, const Vector_& vec_sol_1, const Vector_& vec_sol_2)
    {
      // theta[0] is ignored
      vec_tmp.scale(vec_sol_1, theta[1]);
      if(Math::abs(theta[2]) > 0.0)
        vec_tmp.axpy(vec_sol_2, theta[2]);
    }

    /**
     * \brief Begins the time stepping loop
     *
     * This function should be called directly before the time-stepping loop in the application.
     */
    virtual void begin_loop();

    /**
     * \brief Finishes the time stepping loop
     *
     * This function should be called directly after the time-stepping loop in the application.
     */
    virtual void finish_loop();

    /**
     * \brief Advances to the next time step
     *
     * This function should be called at the beginning of each time step and the return value of
     * this function can be used as the stopping criterion for the time step while-loop.
     *
     * \returns
     * \c true, if a next time step is to be performed or \c false, if the time loop has reached
     * the end of the time interval.
     */
    virtual bool advance_step();

    /// prints runtime information for this object
    virtual void print_runtime(double total_time);

  }; // class TimeStepping
} // namespace CCNDSimple
