#pragma once
#ifndef KERNEL_PROCESS_HPP
#define KERNEL_PROCESS_HPP 1

// includes, system
#include <stdlib.h>
#include <mpi.h>

// includes, Feast
#include <kernel/base_header.hpp>

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief struct encapsulating an MPI process
  *
  * For each MPI process spawned at program start, one static Process object will be created on this MPI process. It
  * will live throughout the program.
  *
  * \author Hilmar Wobker
  * \author Dominik Goeddeke
  */
  struct Process
  {

  public:

    /* *****************
    * member variables *
    *******************/
    /// rank of the process within MPI_COMM_WORLD, set once via constructor and never changed again
    static int rank;

    /**
    * \brief rank of the master process within MPI_COMM_WORLD, set once via constructor and never changed again
    *
    * Every process has to know the rank of the master process in order to trigger screen output (although this will
    * be mainly done by some coordinator processes).
    */
    static int rank_master;

    /// total number of MPI_COMM_WORLD processes spawned at program start
    static unsigned int num_processes;

    /// flag whether this process is the master process
    static bool is_master;
  }; // struct Process

  // define static member variables
  int Process::rank = MPI_PROC_NULL;
  int Process::rank_master = MPI_PROC_NULL;
  unsigned int Process::num_processes = 0;
  bool Process::is_master = false;

} // namespace FEAST

#endif // guard KERNEL_PROCESS_HPP
