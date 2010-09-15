#pragma once
#ifndef KERNEL_PROCESS_GROUP_HPP
#define KERNEL_PROCESS_GROUP_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/base_header.hpp>


/**
* \brief group of processes sharing an MPI communicator
*/
class ProcessGroup
{
  /* *****************
   * private members *
   *******************/
  private:
    /* ******************
     * member variables *
     ********************/
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

    /**
    * \brief rank of this process within this process group
    */
    int _my_rank;

  /* ****************
   * public members *
   ******************/
  public:
    /* **************
     * constructors *
     ****************/
    /**
    * \brief constructor for the case the MPI_Group and the MPI_Comm already exist
    *
    * This constructor is intended for creating a process group object containing all COMM_WORLD processes.
    */
    ProcessGroup(
      const MPI_Comm comm,
      const int num_processes)
      : _comm(comm),
        _num_processes(num_processes),
        _ranks_group_parent(nullptr),
        _process_group_parent(nullptr),
        _group_id(-1)
    {
      // get MPI_Group object for the given communicator
      int mpi_error_code = MPI_Comm_group(_comm, &_group);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_group");

      // and finally look up the local rank of this process within the group
      mpi_error_code = MPI_Group_rank(_group, &_my_rank);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_rank");
    }

    /**
    * \brief constructor for the case the MPI_Group and the corresponding communicator have to be created
    */
    ProcessGroup(
      const int num_processes,
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

      // and finally look up the local rank of this process within the group
      mpi_error_code = MPI_Group_rank(_group, &_my_rank);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_rank");
    }

    /* *******************
     * getters & setters *
     *********************/
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
    inline int my_rank() const
    {
      return _my_rank;
    }
}; // class ProcessGroup


/**
* \brief Class describing a work group, i.e. a set of worker processes sharing the same MPI communicator.
*
* Only the load balancer creates objects of this class to maintain its work groups. Example:
* The process group of the load balancer consists of 6 processes. The load balancer reads the mesh and the solver
* configuration and decides that the coarse grid problem is to be treated by two processes (local ranks 0 and 1) and
* the fine grid problems by all six processes. Then it creates two work groups: one consisting of the two processes
* with local ranks 0 and 1, and one consisting of all six processes. The load balancer tells the waiting GroupProcess
* objects to create Worker objects. These worker objects then create their own communicator and store it. (The load
* balancer will only communicate to these processes via the ProcessGroup communicator, *not* via the WorkGroup
* communicator since it is simply not part of this communicator. So, the WorkGroup class doesn't define this
* work group communicator, only the worker processes do this!
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*
*/
class WorkGroup
  : public ProcessGroup
{
  private:
//    /**
//    * \brief array of workers in the work group
//    *
//    * Here, RemoteWorker objects are used (instead of Worker objects) since they exist on remote processes.
//    * (Note, that WorkGroup objects are only instantiated on load balancer processes.)
//    */
//    RemoteWorker* _workers;

  /* ****************
   * public members *
   ******************/
  public:
    /* **************
     * constructors *
     ****************/
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

    }
};

#endif // guard KERNEL_PROCESS_GROUP_HPP
