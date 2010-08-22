#pragma once
#ifndef KERNEL_LOAD_BAL_HPP
#define KERNEL_LOAD_BAL_HPP 1

#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <kernel/base_header.hpp>
#include <kernel/process.hpp>

using namespace std;

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

/* Class describing a work group, i.e. a set of worker processes sharing the same MPI communicator. Only the load
 * balancer creates objects of this class to maintain its work groups. Example:
 * The process group of the load balancer consists of 6 processes. The load balancer reads the mesh and the solver
 * configuration and decides that the coarse grid problem is to be treated by two processes (local ranks 0 and 1) and
 * the fine grid problems by all six processes. Then it creates two work groups: one consisting of the two processes
 * with local ranks 0 and 1, and one consisting of all six processes. The load balancer tells the waiting GroupProcess
 * objects to create Worker objects. These worker objects then create their own communicator and store it. (The load
 * balancer will only communicate to these processes via the ProcessGroup communicator, *not* via the WorkGroup
 * communicator since it is simply not part of this communicator. So, the WorkGroup class doesn't define this
 * work group communicator, only the worker processes do this!
 */
class WorkGroup
{
  // number of workers in the group
  int _num_workers;

  // array of workers in the work group. Here, RemoteWorker objects are used (instead of Worker objects) since they
  // exist on remote processes.
  RemoteWorker* _workers;

// BRAL: Instead of the array of RemoteWorker objects one could also simply store the ranks of the participating
// processes (with respect to the process group ranks). Not sure yet which is more appropriate.
// int* _ranks_local;

};


// class defining a load balancer process
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
    MPI_Comm _comm_local;
    int _rank_local;
    int _group_id;
    int* _ranks_group;

    // array of work groups the load balancer manages
    WorkGroup* _work_groups;

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
      int& rank_world,
      int& rank_master,
      MPI_Comm& comm_local,
      int& rank_local,
      int& group_id,
      int ranks_group[])
      : Process(rank_world, rank_master),
        _comm_local(comm_local),
        _rank_local(rank_local),
        _group_id(group_id),
        _ranks_group(ranks_group)
    {
    }

    int get_group_id()
    {
      return _group_id;
    }

    int get_rank_local()
    {
      return _rank_local;
    }

    // dummy routine
    void read_mesh()
    {
    }

    // dummy routine
    void create_work_groups()
    {
    }
};

#endif // guard KERNEL_LOAD_BAL_HPP
