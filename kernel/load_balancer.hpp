#pragma once
#ifndef KERNEL_LOAD_BAL_HPP
#define KERNEL_LOAD_BAL_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>

/*
 * Note that only the load bal. processes are directly available to the user! All other processes are in an infinite
 * loop waiting for messages.
 * The user knows what his load balancers should do. He calls, e.g.,
 * if (load_bal.group_id == 0)
 * {
 *   load_bal.readMesh();
 *   ...
 * }
 * else
 * {
 *   do something else
 * }
 *
 * The load bal. with id 0 in the example above then
 *   - reads in the mesh
 *   - builds the necessary WorkGroup objects (e.g., one WorkGroup for the fine mesh problems and one for the
 *     coarse mesh problem)
 *   - tells each member of the work groups to create a Worker object
 *   - lets the workers of one work group build their own communicator
 *     (the worker with rank 0 in this communicator is automatically the coordinator who communicates with the master)
 *   - creates corresponding RemoteWorker objects (as members of the WorkGroup objects) representing the remote Worker
 *     objects
 *   - tells each Worker which (Remote)Worker objects he has to communicate with *in other workgroups* via the
 *     communicator the load balancer shares. E.g., the fine mesh workers have to send the restricted defect vector to
 *     a certain coarse mesh worker, while the coarse mesh worker has to send the coarse mesh correction to one or more
 *     fine mesh workers. Two such communicating workers in different work groups live either on the same process
 *     (internal communication = copy) or on different processes (external communication = MPI send/recv)!
 *     Example:
 *
 *     distribution of submeshes to workers A-G on different levels:
 *
 *     ---------------      ---------------      ---------------
 *     |             |      |      |      |      |      |      |
 *     |             |      |      |      |      |  D   |  G   |
 *     |             |      |      |      |      |      |      |
 *     |      A      |      |  B   |  C   |      ---------------
 *     |             |      |      |      |      |      |      |
 *     |             |      |      |      |      |  E   |  F   |
 *     |             |      |      |      |      |      |      |
 *     ---------------      ---------------      ---------------
 *       level 0               level 1              levels 2-L
 *
 *     * case a:
 *       process group rank:  0  1  2  3
 *             work group 2:  D  E  F  G         (four workers for the problems on level 2-L)
 *             work group 1:  B     C            (two workers for the problem on level 1)
 *             work group 0:  A                  (one worker for the coarse mesh problem on level 0)
 *
 *       Communications:
 *       A <--> B (internal, rank 0) A <--> C (external, ranks 0+2)
 *       B <--> D (internal, rank 0) B <--> E (external, ranks 0+1)
 *       C <--> F (internal, rank 2) C <--> G (external, ranks 2+3)
 *
 *     * case b:
 *       process group rank:  0  1  2  3  4
 *             work group 2:     D  E  F  G
 *             work group 1:     B     C
 *             work group 0:  A
 *
 *       Communications:
 *       A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+3)
 *       B <--> D (internal, rank 1) B <--> E (external, ranks 1+2)
 *       C <--> F (internal, rank 3) C <--> G (external, ranks 3+4)
 *
 *     * case c:
 *       process group rank:  0  1  2  3  4  5  6
 *             work group 2:           D  E  F  G
 *             work group 1:     B  C
 *             work group 0:  A
 *
 *       Communications:
 *       A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+2)
 *       B <--> D (external, ranks 1+3) B <--> E (external, ranks 1+4)
 *       C <--> F (external, ranks 2+5) C <--> G (external, ranks 2+6)
 *
 *   - sends corresponding parts of the mesh to the Worker objects
 */

/**
* \brief class defining a load balancer process
*
* @author Hilmar Wobker
* @author Dominik Goeddeke
*/
class LoadBalancer
  : public Process
{
  /* *****************
   * private members *
   *******************/
  private:

    /* ******************
     * member variables *
     ********************/
    /**
    * \brief process group the load balancer manages
    */
    ProcessGroup* _process_group;

    /**
    * \brief vector of work groups the load balancer manages
    */
    std::vector<WorkGroup> _work_groups;

    /**
    * \brief flag whether this process is a dedicated load balancer process
    *
    *
    */
    bool _dedicated_load_bal_process;

  /* ****************
   * public members *
   ******************/
  public:
    /* **************
     * constructors *
     ****************/
    /**
    * \brief constructor requiring six parameters
    */
    LoadBalancer(
      const int rank_world,
      const int rank_master,
      ProcessGroup* process_group,
      bool dedicated_load_bal_process)
      : Process(rank_world, rank_master),
        _process_group(process_group),
        _dedicated_load_bal_process(dedicated_load_bal_process)
    {
//      std::cout << "Loadbalancer = user entry point tut jetzt mal so als ob er was machen wuerde." << std::endl;
    }

    /* *******************
     * getters & setters *
     *********************/
    /**
    * \brief getter for the process group
    */
    inline ProcessGroup* process_group() const
    {
      return _process_group;
    }

    /**
    * \brief getter for the flag whether this a dedicated load balancer process
    */
    inline bool dedicated_load_bal_process() const
    {
      return _dedicated_load_bal_process;
    }

    /* ******************
     * member functions *
     ********************/
    /**
    * \brief dummy function in preparation of a function reading in a mesh file
    */
    void read_mesh()
    {
    }

    /**
    * \brief dummy function in preparation of a function for managing work groups
    *
    * This dummy function creates two work groups: one consisting of two workers responsible for the coarse grid
    * problem and one consisting of all the other workers responsible for the fine grid problem. Currently, everything
    * is hard-coded. Later, the user must be able to control the creation of work groups and even later the load
    * balancer has to apply clever strategies to create these work groups automatically so that the user doesn't have
    * to do anything.
    */
    void create_work_groups()
    {
      int num_processes = _process_group->num_processes();
      // set number of work groups manually to 2
      int num_work_groups = 2;
      // set number of workers manually to 2 and _num_processes - 3
      int num_workers_in_group[] = {2, num_processes - 3};

      // Partition the ranks of the process group communicator into groups, by simply enumerating the process
      // group ranks and assigning them consecutively to the requested number of processes. Note that the load
      // balancer is the last rank in the process group, i.e. _num_processes - 1.

      // 2-dim. array for storing the process group ranks building the work groups
      int** work_group_ranks;
      work_group_ranks = new int*[num_work_groups];
      // array for storing the work group id of each group process
      int work_group_id[num_processes - 1];
      // iterator for process group ranks, used to split them among the groups
      int iter_group_rank(-1);
      for(int igroup(0) ; igroup < num_work_groups ; ++igroup)
      {
        work_group_ranks[igroup] = new int[num_workers_in_group[igroup]];
        for(int j(0) ; j < num_workers_in_group[igroup] ; ++j)
        {
          // increase group rank
          ++iter_group_rank;
          // set group rank
          work_group_ranks[igroup][j] = iter_group_rank;
          // set group id for the current group rank
          work_group_id[iter_group_rank] = igroup;
        }
      }
      // final sanity check (rank assigned last (=sum of workers - 1) must be smaller than group rank of the load bal.
      assert(iter_group_rank < num_processes - 1);

      // - send work group id to each group process
      // - tell each group process to call
      //
//        mpi_error_code = MPI_Group_incl(_mpi_group, num_workers_in_group[igroup],
//                                        work_group_ranks[my_group], &process_group);
//        MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_incl");
//
      //   while the load balancer sets
//
//        // this is a special valid handle for an empty group that can be passed to operations like
//        // MPI_Comm_create that expect *all* processes in the parent communicator to participate.
//        process_group = MPI_GROUP_EMPTY;
//
      // - then tell each group process to call

//      mpi_error_code = MPI_Comm_create(MPI_COMM_WORLD, process_group, &group_comm);
//      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_create");

      // create work group for coarse grid problems with id 0
//      _work_groups[0] = new WorkGroup(0, 2);
      // create work group for fine grid problems with id 1
//      _work_groups[1] = new WorkGroup(1, );
    }
};

#endif // guard KERNEL_LOAD_BAL_HPP
