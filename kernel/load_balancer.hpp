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
#include <kernel/base_mesh.hpp>

/**
* \brief class defining a load balancer
*
* Each initial process group is organised by a load balancer. It runs on all processes of the process group, which
* eases work flow organisation significantly. The user can choose to use a dedicated load balancer process which does
* not solve any linear algebra problem, but only performs organisation and scheduling tasks (reading the mesh, building
* graph structures, ...). The dedicated load balancer process is the last rank within its process group. If there is no
* such dedicated process, the coordinator of the process group (usually that with rank 0) performs these tasks.
*
* The user knows what each process group and its respective load balancer should do. He calls, e.g.,
* if(load_bal->group_id() == 0)
* {
*   load_bal->readMesh();
*   ...
* }
* else
* {
*   coffee_machine.start();
* }
*
* The load bal. with id 0 in the example above then
* 1) reads in the mesh (this is only done by the dedicated load balancer process or the coordinator, resp.)
* 2) builds the necessary WorkGroup objects (e.g., one WorkGroup for the fine mesh problems and one for the
*    coarse mesh problem) and creates corresponding MPI communicators. The worker with rank 0 in this communicator
*    is usually the coordinator which communicates with the master or the dedicated load balancer. The process
*    topologies of the work groups are then optimised by building corresponding graph structures. There are two cases:
*    a) There is a dedicated load balancer: The dedicated load balancer reads the mesh, creates work groups and and a
*       global graph structure for each work group. Then it distributes the relevant parts of the global graph to all
*       processes of the work groups, which then create their local graph structure and call MPI_Dist_graph_create(...)
*       to build the new MPI process topology in a distributed fashion.
*    b) There is no dedicated load balancer: Same as case a), but instead of the dedicated load balancer the coordinator
*       builds the global graph structure.
*    COMMENT_HILMAR: It might be more clever to also let the MPI implementation decide on which physical process the
*    dedicated load balancer should reside (instead of pinning it to the last rank in the process group). To improve
*    this is task of the ITMC.
*
* 3) tells each member of the work groups to create a Worker object (representing the Worker on this process itself)
*    and corresponding RemoteWorker objects (as members of the WorkGroup objects) representing the remote Worker
*    objects
* COMMENT_HILMAR: Probably, step 3) is skipped, i.e., there will be no extra Worker objects. Instead "the part of
* the WorkGroup living on this process" represents such a worker.
*
* 4) tells each work group which other work groups it has to communicate with via the communicator they all share
*    within the parent process group. E.g., the fine mesh work group has to send the restricted defect vector to the
*    coarse mesh work group, while the coarse mesh work group has to send the coarse mesh correction to the fine mesh
*    work group. Two such communicating work groups live either on the same process (internal communication = copy) or
*    on different processes (external communication = MPI send/recv). (See example below.)
*
* 5) sends corresponding parts of the mesh to the work groups
*
*     Example:
*
*     Distribution of submeshes to processes A-G on different levels (note that processes A-G are not necessarily
*     disjunct, i.e., several of them can refer to the same physical process, see cases a and b):
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
*     * case a, four physical processes:
*       process group rank:  0  1  2  3
*              WorkGroup 2:  D  E  F  G         (four WorkGroup processes for the problems on level 2-L)
*              WorkGroup 1:  B     C            (two WorkGroup processes for the problem on level 1)
*              WorkGroup 0:  A                  (one WorkGroup process for the coarse mesh problem on level 0)
*
*       Communication:
*       A <--> B (internal, rank 0) A <--> C (external, ranks 0+2)
*       B <--> D (internal, rank 0) B <--> E (external, ranks 0+1)
*       C <--> F (internal, rank 2) C <--> G (external, ranks 2+3)
*
*     * case b, five physical processes::
*       process group rank:  0  1  2  3  4
*              WorkGroup 2:     D  E  F  G
*              WorkGroup 1:     B     C
*              WorkGroup 0:  A
*
*       Communication:
*       A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+3)
*       B <--> D (internal, rank 1) B <--> E (external, ranks 1+2)
*       C <--> F (internal, rank 3) C <--> G (external, ranks 3+4)
*
*     * case c, seven physical processes:
*       process group rank:  0  1  2  3  4  5  6
*              WorkGroup 2:           D  E  F  G
*              WorkGroup 1:     B  C
*              WorkGroup 0:  A
*
*       Communication:
*       A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+2)
*       B <--> D (external, ranks 1+3) B <--> E (external, ranks 1+4)
*       C <--> F (external, ranks 2+5) C <--> G (external, ranks 2+6)
*
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
class LoadBalancer
{

private:

  /* *****************
  * member variables *
  *******************/
  /// pointer to the process group the load balancer manages
  ProcessGroup* _process_group;

  /// flag whether the load balancer's process group uses a dedicated load balancer process
  bool _group_has_dedicated_load_bal;

  /// flag whether this process is the dedicated load balancer process (if there is one)
  bool _is_dedicated_load_bal;

  /// vector of work groups the load balancer manages
  std::vector<WorkGroup*> _work_groups;

  /// number of work groups
  int _num_work_groups;

  /**
  * \brief array of number of workers in each work group
  *
  * Dimension: [#_num_work_groups]
  */
  int* _num_workers_in_group;

  /**
  * \brief 2-dim. array for storing the process group ranks building the work groups
  *
  * Dimension: [#_num_work_groups][#_num_workers_in_group[\a group_id]]
  */
  int** _work_group_ranks;

  /**
  * \brief base mesh the load balancer works with
  *
  * bla bla
  */
  BaseMesh* _base_mesh;

public:

  /* *************************
  * constructor & destructor *
  ***************************/
  /// constructor
  LoadBalancer(
    ProcessGroup* process_group,
    bool group_has_dedicated_load_bal)
    : _process_group(process_group),
      _group_has_dedicated_load_bal(group_has_dedicated_load_bal),
      _num_workers_in_group(nullptr),
      _work_group_ranks(nullptr),
      _base_mesh(nullptr)
  {
      // Inquire whether this process is a dedicated load balancer process. This is the case when the process group
      // uses a dedicated load bal. and when this process is the last in the process group.
      _is_dedicated_load_bal = _group_has_dedicated_load_bal &&
                               _process_group->rank() == _process_group->num_processes()-1;
  }

  /// destructor
  ~LoadBalancer()
  {
    if (_work_group_ranks != nullptr)
    {
      for(int igroup(0) ; igroup < _num_work_groups ; ++igroup)
      {
        delete [] _work_group_ranks[igroup];
      }
      delete [] _work_group_ranks;
      _work_group_ranks = nullptr;
    }
    for(unsigned int igroup(0) ; igroup < _work_groups.size() ; ++igroup)
    {
      delete _work_groups[igroup];
    }
    if (_num_workers_in_group != nullptr)
    {
      delete [] _num_workers_in_group;
    }
    if (_base_mesh != nullptr)
    {
      delete _base_mesh;
    }
  }

  /* ******************
  * getters & setters *
  ********************/
  /**
  * \brief getter for the process group
  *
  * \return ProcessGroup pointer #_process_group
  */
  inline ProcessGroup* process_group() const
  {
    return _process_group;
  }

  /* *****************
  * member functions *
  *******************/
  /**
  * \brief getter for the flag whether this process is a dedicated load balancer process
  *
  * \return boolean flag #_is_dedicated_load_bal
  */
  inline bool is_dedicated_load_bal() const
  {
    return _is_dedicated_load_bal;
  }

  /// dummy function in preparation of a function reading in a mesh file
  void read_mesh()
  {
    // the mesh is read by the dedicated load balancer if there is one, otherwise by the group coordinator
    if(_is_dedicated_load_bal || (!_group_has_dedicated_load_bal && _process_group->is_coordinator()))
    {
      _base_mesh = new BaseMesh();
      _base_mesh->read_mesh();
    }
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
    // create process topology (either the dedicated load balancer if there is one, otherwise the group coordinator)
    if(_is_dedicated_load_bal || (!_group_has_dedicated_load_bal && _process_group->is_coordinator()))
    {
      // ...
    }

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
    if(_group_has_dedicated_load_bal)
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
        if(my_rank == _work_group_ranks[igroup][j])
        {
          my_group = igroup;
        }
      }
    } // for(int igroup(0) ; igroup < _num_work_groups ; ++igroup)

    // sanity check for the rank assigned last
    if(_group_has_dedicated_load_bal)
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
      if(my_group == igroup)
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
    } // for(int igroup(0) ; igroup < _num_work_groups ; ++igroup)
  } // create_work_groups()
};

#endif // guard KERNEL_LOAD_BAL_HPP
