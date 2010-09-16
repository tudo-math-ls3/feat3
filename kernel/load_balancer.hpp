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
 *     coarse mesh problem) and creates corresponding MPI communicators ((the worker with rank 0 in this communicator
 *     is automatically the coordinator who communicates with the master)
 *   - tells each member of the work groups to create a Worker object (representing the Worker on this process itself)
*      and corresponding RemoteWorker objects (as members of the WorkGroup objects) representing the remote Worker
 *     objects
 *   - tells each Worker which (Remote)Worker objects he has to communicate with *in other workgroups* via the
 *     communicator they all share within the parent process group. E.g., the fine mesh workers have to send the
 *     restricted defect vector to a certain coarse mesh worker, while the coarse mesh worker has to send the coarse
 *     mesh correction to one or more fine mesh workers. Two such communicating workers in different work groups live
 *     either on the same process (internal communication = copy) or on different processes (external communication
 *     = MPI send/recv)!
 *
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
{
  /* *****************
   * private members *
   *******************/
  private:

    /* ******************
     * member variables *
     ********************/
    /**
    * \brief pointer to the process group the load balancer manages
    */
    ProcessGroup* _process_group;

    /**
    * \brief flag whether the load balancer's process group uses a dedicated load balancer process
    */
    bool _group_uses_dedicated_load_bal;

    /**
    * \brief flag whether this process is the dedicated load balancer process (if there is one)
    */
    bool _is_dedicated_load_bal;

    /**
    * \brief vector of work groups the load balancer manages
    */
    std::vector<WorkGroup*> _work_groups;

    /**
    * \brief number of work groups
    */
    int _num_work_groups;

    /**
    * \brief array of number of workers in each work group
    *
    * <em>Dimension:</em> [#_num_work_groups]
    */
    int* _num_workers_in_group;

    /**
    * \brief 2-dim. array for storing the process group ranks building the work groups
    *
    * <em>Dimension:</em> [#_num_work_groups][#_num_workers_in_group[group_id]]
    */
    int** _work_group_ranks;
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
      ProcessGroup* process_group,
      bool group_uses_dedicated_load_bal)
       : _process_group(process_group),
         _group_uses_dedicated_load_bal(group_uses_dedicated_load_bal)
    {
        // Inquire whether this process is a dedicated load balancer process. This is the case when the process group
        // uses a dedicated load bal. and when this process is the last in the process group.
        _is_dedicated_load_bal = _group_uses_dedicated_load_bal &&
                                 _process_group->rank() == _process_group->num_processes()-1;
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
    * \brief getter for the flag whether this process is a dedicated load balancer process
    */
    inline bool is_dedicated_load_bal() const
    {
      return _is_dedicated_load_bal;
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
      // shortcut to the number processes in the load balancer's process group
      int num_processes = _process_group->num_processes();
      // shortcut to the process group rank of this process
      int my_rank = _process_group->rank();

      // number of work groups, manually set to 2
      _num_work_groups = 2;
      // array of numbers of workers per work group
      _num_workers_in_group = new int[2];
      // set number of workers manually to 2 and remaining processes, resp.
      _num_workers_in_group[0] = 2;
      _num_workers_in_group[1] = num_processes - _num_workers_in_group[0];
      // if a dedicated load balancer process is used, decrease the number of workers in the second work group by 1
      if (_group_uses_dedicated_load_bal)
      {
        --_num_workers_in_group[1];
      }
      // assert that the number of processes in the second group is positive, i.e. that enough processes are available
      assert(_num_workers_in_group[1] > 0);

      // Partition the ranks of the process group communicator into groups, by simply enumerating the process
      // group ranks and assigning them consecutively to the requested number of processes. Note that the dedicated
      // load balancer, if required, is the last rank in the process group, i.e. _num_processes - 1.

      _work_group_ranks = new int*[_num_work_groups];
      // group id of this process
      int my_group(-1);
      // iterator for process group ranks, used to split them among the groups
      int iter_group_rank(-1);
      // now partition the ranks
      for(int igroup(0) ; igroup < _num_work_groups ; ++igroup)
      {
        _work_group_ranks[igroup] = new int[_num_workers_in_group[igroup]];
        for(int j(0) ; j < _num_workers_in_group[igroup] ; ++j)
        {
          // increase group rank
          ++iter_group_rank;
          // set group rank
          _work_group_ranks[igroup][j] = iter_group_rank;
          // inquire whether this process belongs to the current group
          if (my_rank == _work_group_ranks[igroup][j])
          {
            my_group = igroup;
          }
        }
      }
      // sanity check for the rank assigned last
      if (_group_uses_dedicated_load_bal)
      {
        // the dedicated load balancer has the last rank num_processes - 1
        assert(iter_group_rank == num_processes - 2);
      }
      else
      {
        assert(iter_group_rank == num_processes - 1);
      }

      // create WorkGroup objects including MPI groups and MPI communicators
      // (Exclude the dedicated load balancer process if there is one because it is only a member of the process group
      // but not of any work group we set up. So a special group is needed since otherwise the forking call below will
      // deadlock.)
      // It is not possible to set up all WorkGroups in one call, since the processes building the WorkGroups might
      // not be disjunct (see case a and case b in the example above). Hence, there are as many calls as there are
      // WorkGroups. All processes not belonging to the WorkGroup currently created call the MPI_Comm_create() function
      // with a dummy communicator.
      _work_groups.resize(_num_work_groups, nullptr);
      for(int igroup(0) ; igroup < _num_work_groups ; ++igroup)
      {
        if (my_group == igroup)
        {
          _work_groups[igroup] = new WorkGroup(_num_workers_in_group[my_group], _work_group_ranks[my_group],
                                               _process_group, my_group);
        }
        else
        {
          // *All* processes of the parent MPI group have to call the MPI_Comm_create() routine (otherwise the forking
          // will deadlock), so let all processes that are not part of the current work group call it with special
          // MPI_GROUP_EMPTY and dummy communicator.
          MPI_Comm dummy_comm;
          int mpi_error_code = MPI_Comm_create(_process_group->comm(), MPI_GROUP_EMPTY, &dummy_comm);
          MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_create");
        }
      }
    }
};

#endif // guard KERNEL_LOAD_BAL_HPP
