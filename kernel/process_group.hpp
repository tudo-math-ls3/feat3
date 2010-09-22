#pragma once
#ifndef KERNEL_PROCESS_GROUP_HPP
#define KERNEL_PROCESS_GROUP_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/process.hpp>
#include <kernel/worker.hpp>

/**
* \brief group of processes sharing an MPI communicator
*/
class ProcessGroup
{

private:

protected:
  /* *****************
  * member variables *
  *******************/
  /**
  * \brief MPI group representing the processes of the group
  */
  MPI_Group _group;

  /**
  * \brief communicator shared by all processes of the group
  */
  MPI_Comm _comm;

  /**
  * \brief number of processes in this group
  */
  int _num_processes;

  /**
  * \brief rank of this process with respect to the ProcessGroup's communicator
  */
  int _rank;

  /**
  * \brief array of ranks the processes of this group have in the parent group
  *
  * <em>Dimension:</em> [#_num_processes]
  */
  // COMMENT_HILMAR, 15.9.2010:
  // Currently, this is simply a pointer to the part of the corresponding array
  // _group_ranks_world[my_group] in universe.hpp. Maybe it makes sense to copy the array and delete the one in
  // universe.hpp.
  int* _ranks_group_parent;

  /**
  * \brief pointer to the process group from which this process group has been spawned
  */
  ProcessGroup* _process_group_parent;

//    /**
//    * \brief process groups spawned from this process group
//    */
// COMMENT_HILMAR, 15.9.2010: not sure yet if this is needed
//    ProcessGroup** _process_group_child;

  /**
  * \brief id of the process group
  */
  const int _group_id;


public:
  /* *************
  * constructors *
  ***************/
  /**
  * \brief constructor for the case the MPI_Group and the MPI_Comm already exist
  *
  * This constructor is intended for creating a process group object containing all COMM_WORLD processes.
  */
  ProcessGroup(
    MPI_Comm comm,
    int num_processes)
    : _comm(comm),
      _num_processes(num_processes),
      _ranks_group_parent(nullptr),
      _process_group_parent(nullptr),
      _group_id(-1)
  {
    // get MPI_Group object for the given communicator
    int mpi_error_code = MPI_Comm_group(_comm, &_group);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_group");

    // and finally look up the local rank of this process w.r.t. the group's communicator
    mpi_error_code = MPI_Group_rank(_group, &_rank);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_rank");

    // since this is the world group of processes, the local and the global rank should be equal
    assert(Process::rank == _rank);
  }

  /**
  * \brief constructor for the case the MPI_Group and the corresponding communicator have to be created
  */
  ProcessGroup(
    int num_processes,
    int ranks_group_parent[],
    ProcessGroup* process_group_parent,
    const int group_id)
    : _num_processes(num_processes),
      _ranks_group_parent(ranks_group_parent),
      _process_group_parent(process_group_parent),
      _group_id(group_id)
  {
    int mpi_error_code = MPI_Group_incl(_process_group_parent->_group, _num_processes,
                                        _ranks_group_parent, &_group);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_incl");

    // Create the group communicator for, among others, collective operations.
    // It is essential that *all* processes in MPI_COMM_WORLD participate since MPI_COMM_WORLD
    // is a parent sub-universe, i.e., something to spawn from
    mpi_error_code = MPI_Comm_create(_process_group_parent->_comm, _group, &_comm);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_create");

    // and finally look up the local rank of this process w.r.t. the group's communicator
    mpi_error_code = MPI_Group_rank(_group, &_rank);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_rank");
  }

  /* ******************
  * getters & setters *
  ********************/
  /**
  * \brief getter for the MPI group
  *
  * Return a reference in order to avoid making a copy. Make this return reference constant so that the user
  * cannot change the object.
  */
  inline const MPI_Group& group() const
  {
    return _group;
  }

  /**
  * \brief getter for the MPI communicator
  *
  * Return a reference in order to avoid making a copy. Make this return reference constant so that the user
  * cannot change the object.
  */
  inline const MPI_Comm& comm() const
  {
    return _comm;
  }

  /**
  * \brief getter for the number of processes
  */
  inline int num_processes() const
  {
    return _num_processes;
  }

  /**
  * \brief getter for the pointer to the parent process group
  */
  inline ProcessGroup* process_group_parent() const
  {
    return _process_group_parent;
  }

  /**
  * \brief getter for the rank in the group communicator this process belongs to
  */
  inline int group_id() const
  {
    return _group_id;
  }

  /**
  * \brief getter for the rank in the group communicator this process belongs to
  */
  inline int rank() const
  {
    return _rank;
  }
}; // class ProcessGroup


/**
* \brief class describing a work group, i.e. a set of worker processes sharing the same MPI communicator and performing
*        the same task
*
* A WorkGroup object managing n processes typically consists of one Worker object (living on this process) and n-1
* RemoteProcess objects. WorkGroup objects are created by the load balancer. Example:
* The process group of a load balancer consists of 6 processes. One process (either that with local rank 0 or the
* dedicated one) reads the mesh and the solver configuration and decides that the coarse grid problem is to be
* treated by 2 processes (local ranks 0 and 1) and the fine grid problems by all 6 processes. Then 2 work
* groups are created: one consisting of the 2 processes with local ranks 0 and 1, and one consisting of all 6
* processes, each of the work groups having its own MPI communicator. In building the work groups, also the
* corresponding (Remote)Worker objects are created. Communication between different work groups or with the dedicated
* load balancer process is done via the enclosing ProcessGroup communicator.
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*
*/
class WorkGroup
  : public ProcessGroup
{

private:

  /* *****************
  * member variables *
  *******************/
  /**
  * \brief worker object living on this process
  */
  Worker* _worker;

  /**
  * \brief vector of remote workers in the work group
  *
  * Here, RemoteWorker objects are used (instead of Worker objects) since they exist on remote processes.
  */
  std::vector<RemoteWorker*> _remote_workers;


public:

  /* *************
  * constructors *
  ***************/
  /**
  * \brief constructor for the case the MPI_Group and the corresponding communicator have to be created
  */
  WorkGroup(
    const int num_processes,
    int ranks_group_parent[],
    ProcessGroup* process_group_parent,
    const int group_id)
    : ProcessGroup(num_processes, ranks_group_parent, process_group_parent, group_id)
  {
    _worker = new Worker(_comm, _rank, process_group_parent->comm(), _process_group_parent->rank());
  }
};

#endif // guard KERNEL_PROCESS_GROUP_HPP
