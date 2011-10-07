/* GENERAL_REMARK_BY_HILMAR:
 * In order to really start working on the load balancer, it is essential define some preliminary "statistics module"
 * providing some "pseudo" data for matrix patches / process patches. I think it is not really necessary that the
 * matrix patch / process patch stuff (issue00010) is already finished for that (although it would makes things
 * clearer and easier...) For a start, it should be sufficient to "fake" some matrix patches and statistics. At the
 * end the load balancer "only" needs a graph and node values...
 * Until now I only implemented two semi-hard-coded pseudo tests which show what kind of data the load balancer has to
 * prepare such that the manager can create work groups. These tests have of course nothing to do with real load
 * balancing...
 *
 * Also see COMMENT_HILMAR throughout this file.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */

#pragma once
#ifndef KERNEL_LOAD_BALANCER_HPP
#define KERNEL_LOAD_BALANCER_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>
#include <kernel/base_mesh/bm.hpp>

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

namespace FEAST
{
  /**
  * \brief class defining the load balancer
  *
  * \tparam space_dim_
  * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
  *
  * \tparam world_dim_
  * world dimension (determines the number of coordinates)
  *
  * \todo This class is only a dummy with some comments! No real load balancing functionalities implemented yet!
  *
  * \author Hilmar Wobker
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class LoadBalancer
  {

  private:

    /* *****************
    * member variables *
    *******************/

    /// pointer to the main process group
    ProcessGroup* _process_group;

    /// flag whether the main process group uses a dedicated load balancer process
    bool const _group_has_dedicated_load_bal;

    /// rank of the dedicated load balancer process (if there is one, otherwise MPI_PROC_NULL)
    int _rank_dedicated_load_balancer;

    /// underlying base mesh
    BaseMesh::BM<space_dim_, world_dim_>* _base_mesh;

    /// number of (ext.) work groups
    unsigned int _num_work_groups;

    ///array of numbers of processes per ext. work group
    unsigned int* _num_proc_in_ext_work_group;

    /**
    * \brief array indicating whether the ext. work groups contain an extra process for the coordinator
    *
    * This coordinator will not be a compute process in the corresponding work group.
    * Usually a boolean would do the trick, but we cannot send boolean arrays via MPI. (C++ bindings are deprecated,
    * hence we cannot use MPI::BOOL. MPI_LOGICAL cannot be used, since this is not equivalent to C++'s bool.)
    * So, we use unsigned char here and treat is as if it were boolean.
    */
    unsigned char* _group_contains_extra_coord;

    /// array defining the rank partitioning of ext. work groups
    int** _ext_work_group_ranks;

    /// array of Graphs representing the connectivity of work group processes
    Graph** _graphs;

    /// number of interlevel groups
    unsigned int _num_interlevel_groups;

    /**
    * \brief array of number of processes in each interlevel group
    *
    * Dimension: [#_num_interlevel_groups]
    */
    unsigned int* _num_proc_in_interlevel_group;

    /**
    * \brief array of positions of the root rank within the array _interlevel_group_ranks per interlevel group
    *
    * Dimension: [#_num_interlevel_groups]
    *
    * COMMENT_HILMAR: Maybe one could also store the rank of the root *always* in the first position of the array
    *   _interlevel_group_ranks and omit this extra array here... but I'm not sure if this is always possible.
    *   Furthermore, it would be more error-prone since one has to know that the root has to be at the first position.
    */
    unsigned int* _pos_of_root_in_interlevel_group;

    /**
    * \brief 2-dim. array of process group ranks building the interlevel groups
    *
    * Dimension: [#_num_interlevel_groups][#_num_proc_in_interlevel_group[\a group_id]]
    */
    int** _interlevel_group_ranks;

    /**
    * \brief array of flags whether graph has been created in this class
    *
    * Just used for the dummy test functions. Can be removed, when test dummy functions have been removed.
    */
    bool* _graph_has_been_created_here;

  public:

    /* ****************************
    * constructors and destructor *
    ******************************/
    /**
    * \brief CTOR for the case the load balancer lives on the coordinator process
    *
    * In this case, the base mesh can be simply passed as pointer
    */
    LoadBalancer(
      ProcessGroup* process_group,
      bool group_has_dedicated_load_bal,
      int rank_dedicated_load_balancer)
      : _process_group(process_group),
        _group_has_dedicated_load_bal(group_has_dedicated_load_bal),
        _rank_dedicated_load_balancer(rank_dedicated_load_balancer),
        _base_mesh(nullptr),
        _num_work_groups(0),
        _num_proc_in_ext_work_group(nullptr),
        _group_contains_extra_coord(nullptr),
        _ext_work_group_ranks(nullptr),
        _graphs(nullptr),
        _num_interlevel_groups(0),
        _num_proc_in_interlevel_group(nullptr),
        _pos_of_root_in_interlevel_group(nullptr),
        _interlevel_group_ranks(nullptr),
        _graph_has_been_created_here(nullptr)
    {
    }

    /// DTOR
    ~LoadBalancer()
    {
      if(_graphs != nullptr)
      {
        // this for loop is only necessary due to the dummy test functions define_work_groups_1/2()
        for(unsigned int i(0) ; i < _num_work_groups ; ++i)
        {
          if(_graphs[i] != nullptr && _graph_has_been_created_here[i])
          {
            delete _graphs[i];
          }
        }
        delete [] _graph_has_been_created_here;
        delete [] _graphs;
      }

      // Although created here, all arrays related to work groups and interlevel groups are not delete here, but
      // within the DTOR of class ManagerComp. The reason is that within the ManagerComp class these arrays have to
      // be deleted anyway on the non-coordinator processes. To avoid code duplication we delete the arrays on the
      // coordinator process (which point to the same arrays as the arrays in the LoadBalancer class) also in the DTOR
      // of ManagerComp class.
    }

    /* ******************
    * getters & setters *
    ********************/
    /// stores the base mesh
    inline void set_base_mesh(BaseMesh::BM<space_dim_, world_dim_>* base_mesh)
    {
      CONTEXT("LoadBalancer::set_base_mesh()");
      _base_mesh = base_mesh;
    }


    /**
    * \brief getter for the number of ext. work groups
    *
    * \return number of ext. work groups
    */
    inline unsigned int num_work_groups() const
    {
      CONTEXT("LoadBalancer::num_work_groups()");
      return _num_work_groups;
    }


    /**
    * \brief getter for the array of number of processes per ext. work group
    *
    * \return pointer to array of number of processes per ext. work group
    */
    inline unsigned int* num_proc_in_ext_work_group() const
    {
      CONTEXT("LoadBalancer::num_proc_in_ext_work_group()");
      return _num_proc_in_ext_work_group;
    }


    /**
    * \brief getter for the array indicating whether the ext. work groups contain an extra process for the coordinator
    *
    * \return pointer to array indicating whether the ext. work groups contain an extra process for the coordinator
    */
    inline unsigned char* group_contains_extra_coord() const
    {
      CONTEXT("LoadBalancer::group_contains_extra_coord()");
      return _group_contains_extra_coord;
    }


    /**
    * \brief getter for the 2D array of ext. work group ranks
    *
    * \return 2D array of ext. work group ranks
    */
    inline int** ext_work_group_ranks() const
    {
      CONTEXT("LoadBalancer::ext_work_group_ranks()");
      return _ext_work_group_ranks;
    }


    /**
    * \brief getter for the array of graph pointers
    *
    * \return array of graph pointers
    */
    inline Graph** graphs()
    {
      CONTEXT("LoadBalancer::graphs()");
      return _graphs;
    }


    /**
    * \brief getter for the number of interlevel groups
    *
    * \return number of interlevel groups
    */
    inline unsigned int num_interlevel_groups() const
    {
      CONTEXT("LoadBalancer::num_interlevel_groups()");
      return _num_interlevel_groups;
    }


    /**
    * \brief getter for the array of number of processes per interlevel group
    *
    * \return pointer to array of number of processes per interlevel group
    */
    inline unsigned int* num_proc_in_interlevel_group() const
    {
      CONTEXT("LoadBalancer::num_proc_in_interlevel_group()");
      return _num_proc_in_interlevel_group;
    }


    /**
    * \brief getter for the array of position of root rank per interlevel group
    *
    * \return pointer to array of position of root rank per interlevel group
    */
    inline unsigned int* pos_of_root_in_interlevel_group() const
    {
      CONTEXT("LoadBalancer::pos_of_root_in_interlevel_group()");
      return _pos_of_root_in_interlevel_group;
    }


    /**
    * \brief getter for the 2D array of interlevel group ranks
    *
    * \return 2D array of interlevel group ranks
    */
    inline int** interlevel_group_ranks() const
    {
      CONTEXT("LoadBalancer::interlevel_group_ranks()");
      return _interlevel_group_ranks;
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief dummy test function in preparation of a function for defining work groups
    *
    * This test function creates two work groups: one consisting of two workers responsible for the coarse grid
    * problem and one consisting of all the other workers responsible for the fine grid problem. Currently, everything
    * is hard-coded. Later, the user must be able to control the creation of work groups and even later the load
    * balancer has to apply clever strategies to create these work groups automatically so that the user doesn't have
    * to do anything.
    *
    * The following test code is semi hard-wired in the sense that we schedule one BMC to each processor.
    * Furthermore, two extra processes are required for testing the creation of a second work group.
    * Later, this has to be all done in some auto-magically way.
    *
    * In this dummy test example you can imagine we have something like the following situation:
    \verbatim
    -----------------        -----------------
    |               |        |12 |13 |14 |15 |
    |      1        |        |---+---+---+---|
    |               |        | 8 | 9 |10 |11 |
    |---------------|        |---+---+---+---|
    |               |        | 4 | 5 | 6 | 7 |
    |      0        |        |---+---+---+---|
    |               |        | 0 | 1 | 2 | 3 |
    -----------------        -----------------
    coarse grid (CG)          fine grid (FG)
        2 CGMPs                  16 PPs
    \endverbatim
    * The fine grid problem consists of 16 process patches (here PP = BMC) and is solved by n (in the figure n=16)
    * processes, i.e., the fine grid work group contains n=16 processes. The coarse grid problem consists of
    * 2 coarse grid matrix patches (CGMP) and is solved by 2 processes, i.e., the coarse grid work group contains
    * 2 processes.
    *
    * To optimise the communication between the coordinator of the main process group and the work groups, we add
    * this coordinator to a work group if it is not a compute process of this work group anyway. Hence, for each work
    * group, there are two possibilities:
    * The coordinator of the main process group is...
    *   1) ... not part of the work group:
    *     --> work group adds the coordinator as extra process
    *   2) ... is part of the work group:
    *     --> work group does not have to add the coordinator
    * The work group plus the coordinator process is called 'extended work group'. In the case the coordinator is part
    * of the work group anyway, the extended work group and the actual work group are equal.
    *
    * Depending on whether a dedicated load balancer is used or not, two different tests are performed (where n is the
    * number of base mesh cells):
    * 1) with dedicated load balancer process
    *    - work group for coarse grid: 2 processes: {0, 1}
    *    - work group for fine grid: n processes: {1, ..., n}
    *    - i.e. process 1 is in both work groups
    *    - coordinator process: 0
    *    - dedicated load balancer process: n+1
    * 2) without dedicated load balancer process
    *    - work group for coarse grid: 2 processes: {0, 1}
    *    - work group for fine grid: n processes: {2, ..., n+1}
    *    - i.e. the two work groups are disjunct
    *    - coordinator process: 0
    * This distinction is done only to be able to run both tests with n+2 processes in total. There is of course no
    * relation between overlapping of work groups on the one hand and use of a dedicated load balancer on the other
    * hand.
    *
    * We need two interlevel groups which connects each of the 2 processes in the coarse grid level work group to
    * 8 processes, resp., of the fine level work group. In the example above, the process dealing with CGMP 0 (1) builds
    * a interlevel group with the processes dealing with PPs 0-7 (8-15).
    *
    * COMMENT_HILMAR: We assume here that each process receives exactly one BMC and that the index of the cell in the
    *   graph structure equals the local rank within the work group. Later, there will be the matrix patch layer and
    *   the process patch layer, which both have their own connectivity structure. Here, we then actually need the
    *   connectivity graph of the process patch layer and a mapping between process patch index and work group ranks.
    *
    * \warning If this dummy function is removed / changed, also remove the deletion of _graphs[0] in the DTOR
    */
    void define_work_groups_2()
    {
      CONTEXT("LoadBalancer::define_work_groups_2()");

      // set number of ext. work groups manually to 2
      _num_work_groups = 2;
      // set number of interlevel groups manually to 2
      _num_interlevel_groups = 2;
      // allocate arrays
      _num_proc_in_ext_work_group = new unsigned int[_num_work_groups];
      _group_contains_extra_coord = new unsigned char[_num_work_groups];
      _ext_work_group_ranks = new int*[_num_work_groups];
      _graphs = new Graph*[_num_work_groups];
      _interlevel_group_ranks = new int*[_num_interlevel_groups];
      _num_proc_in_interlevel_group = new unsigned int[_num_interlevel_groups];
      _pos_of_root_in_interlevel_group = new unsigned int[_num_interlevel_groups];

      // shortcut to the number of processes in the manager's process group
      unsigned int num_processes = _process_group->num_processes();

      // shortcut to the number of cells in the base mesh
      unsigned int num_cells = unsigned int(_base_mesh->num_cells());

      // debug output
      Logger::log("num_processes: " + stringify(num_processes) + "\nnum_cells: " + stringify(num_cells) + "\n",
        Logger::master);
      // assert that the number of processes is n+2
      ASSERT(num_processes == num_cells + 2, "Number of processes " + stringify(num_processes)
             + " must be number of cells + 2, i.e., " + stringify(num_cells + 2) + ".");

      // set up the two test cases

      if(_group_has_dedicated_load_bal)
      {
        // test case 1 (n is the number of base mesh cells)
        // with dedicated load balancer process
        //  - work group for coarse grid: 2 processes: {0, 1}
        //  - work group for fine grid: n processes: {1, ..., n}
        //  - i.e. process 1 is in both work groups
        //  - coordinator process: 0
        //  - dedicated load balancer: n+1

        // the coordinator (rank 0) is at the same time a compute process of the first work group, so only the second
        // work group has to add an extra process
        _group_contains_extra_coord[0] = false;
        _group_contains_extra_coord[1] = true;

        // set number of processes per group
        _num_proc_in_ext_work_group[0] = 2;
        _num_proc_in_ext_work_group[1] = num_cells + 1;

        // partition the process group ranks into work groups

        // coarse grid work group
        _ext_work_group_ranks[0] = new int[_num_proc_in_ext_work_group[0]];
        _ext_work_group_ranks[0][0] = 0;
        _ext_work_group_ranks[0][1] = 1;

        // fine grid work group
        _ext_work_group_ranks[1] = new int[_num_proc_in_ext_work_group[1]];
        // set entries to {0, 1, ..., n}
        for(unsigned int i(0) ; i < _num_proc_in_ext_work_group[1] ; ++i)
        {
          _ext_work_group_ranks[1][i] = i;
        }

        // set number of processes per interlevel group: one for the coarse grid process and half of the fine grid
        // processes (if num_cells is odd, then the first interlevel group contains one process less than the second
        // interlevel group)
        _num_proc_in_interlevel_group[0] = 1 + num_cells/2;
        _num_proc_in_interlevel_group[1] = 1 + num_cells - num_cells/2;

        // set ranks for the interlevel work groups

        // first interlevel group connecting 0 and 1, ..., n/2
        _interlevel_group_ranks[0] = new int[_num_proc_in_interlevel_group[0]];
        // root of the interlevel group = coarse grid process
        _interlevel_group_ranks[0][0] = _ext_work_group_ranks[0][0];
        _pos_of_root_in_interlevel_group[0] = 0;
        // set entries to {1, 2, ..., n/2} (omit _ext_work_group_ranks[1][0] = 0 since this is only the extra
        // coordinator process and not a compute process)
        for(unsigned int i(1) ; i < _num_proc_in_interlevel_group[0] ; ++i)
        {
          _interlevel_group_ranks[0][i] = _ext_work_group_ranks[1][i];
        }

        // second interlevel group connecting 1 and n/2 + 1, ..., n
        _interlevel_group_ranks[1] = new int[_num_proc_in_interlevel_group[1]];
        // root of the interlevel group = coarse grid process = _ext_work_group_ranks[0][1] = 1
        _interlevel_group_ranks[1][0] = _ext_work_group_ranks[0][1];
        _pos_of_root_in_interlevel_group[1] = 0;
        // set entries to {n/2 + 1, ..., n}
        for(unsigned int i(1) ; i < _num_proc_in_interlevel_group[1] ; ++i)
        {
          _interlevel_group_ranks[1][i] = _ext_work_group_ranks[1][i + _num_proc_in_interlevel_group[0] - 1];
        }
        // COMMENT_HILMAR: Note that process with rank 1 is in both interlevel groups! This is not optimal (since
        //   when communication within the two groups is started simultaneously, they have to wait on each other), but
        //   it should work (i.e., MPI should be able to deal with this situation automatically without deadlocking).
      }
      else
      {
        // test case 2 (n is the number of base mesh cells)
        // without dedicated load balancer process
        //  - work group for coarse grid: 2 processes: {0, 1}
        //  - work group for fine grid: n processes: {2, ..., n+1}
        //  - i.e. the two work groups are disjunct
        //  - coordinator (and load balancer) process: 0

        // the coordinator (rank 0) is at the same time a compute process of the first work group, so only the second
        // work group has to add an extra process
        _group_contains_extra_coord[0] = false;
        _group_contains_extra_coord[1] = true;

        // set number of processes per group
        _num_proc_in_ext_work_group[0] = 2;
        _num_proc_in_ext_work_group[1] = num_cells + 1;

        // partition the process group ranks into work groups
        _ext_work_group_ranks[0] = new int[_num_proc_in_ext_work_group[0]];
        _ext_work_group_ranks[0][0] = 0;
        _ext_work_group_ranks[0][1] = 1;
        _ext_work_group_ranks[1] = new int[_num_proc_in_ext_work_group[1]];
        // set entries to {0, 2, ..., n+1}
        _ext_work_group_ranks[1][0] = 0;
        for(unsigned int i(1) ; i < _num_proc_in_ext_work_group[1] ; ++i)
        {
          _ext_work_group_ranks[1][i] = i+1;
        }

        // set number of processes per interlevel group: one for the coarse grid process and half of the fine grid
        // processes (if num_cells is odd, then the first interlevel group contains one process less than the second
        // interlevel group)
        _num_proc_in_interlevel_group[0] = 1 + num_cells/2;
        _num_proc_in_interlevel_group[1] = 1 + num_cells - num_cells/2;

        // set ranks for the interlevel work groups

        // first interlevel group connecting 0 and 2, ..., n/2 + 1
        _interlevel_group_ranks[0] = new int[_num_proc_in_interlevel_group[0]];
        // root of the interlevel group = coarse grid process
        _interlevel_group_ranks[0][0] = _ext_work_group_ranks[0][0];
        _pos_of_root_in_interlevel_group[0] = 0;
        // set entries to {2, ..., n/2 + 1} (omit _ext_work_group_ranks[1][0] = 0 since this is only the extra
        // coordinator process and not a compute process)
        for(unsigned int i(1) ; i < _num_proc_in_interlevel_group[0] ; ++i)
        {
          _interlevel_group_ranks[0][i] = _ext_work_group_ranks[1][i];
        }

        // second interlevel group connecting 1 and n/2 + 2, ..., n + 1
        _interlevel_group_ranks[1] = new int[_num_proc_in_interlevel_group[1]];
        // root of the interlevel group = coarse grid process
        _interlevel_group_ranks[1][0] = _ext_work_group_ranks[0][1];
        _pos_of_root_in_interlevel_group[1] = 0;
        // set entries to {n/2 + 2, ..., n + 1}
        for(unsigned int i(1) ; i < _num_proc_in_interlevel_group[1] ; ++i)
        {
          _interlevel_group_ranks[1][i] = _ext_work_group_ranks[1][i + _num_proc_in_interlevel_group[0] - 1];
        }
      }

      /* *********************************************************
      * create graph structures corresponding to the work groups *
      ***********************************************************/
      // build an artificial graph mimicing the distribution of the base mesh cells to two processors (e.g. for a mesh
      // of 16 BMCs: BMCs 0-7 on proc 1 and BMCs 8-15 on proc 2) which start an imagined coarse grid solver; this graph
      // will be used for the coarse grid work group
      Index* index = new Index[3];
      Index* neighbours = new Index[2];
      index[0] = 0;
      index[1] = 1;
      index[2] = 2;
      neighbours[0] = 1;
      neighbours[1] = 0;
      // Artificially create a graph object here. Usually, this comes from somewhere else.
      _graphs[0] = new Graph(2, index, neighbours);
      // arrays index and neighbours are copied within the Graph CTOR, hence they can be deallocated here
      delete [] index;
      delete [] neighbours;

      // get connectivity graph of the base mesh; this one will be used for the fine grid work group
      _graphs[1] = _base_mesh->graph();

      _graph_has_been_created_here = new bool[_num_work_groups];
      _graph_has_been_created_here[0] = true;
      _graph_has_been_created_here[1] = false;
    }


    /**
    * \brief second dummy test function in preparation of a function for defining work groups
    *
    * Same as above, but only one work group consisting of all processes of the compute process group and no
    * interlevel group are created.
    */
    void define_work_groups_1()
    {
      CONTEXT("LoadBalancer::define_work_groups_1()");

      // set number of ext. work groups manually to 1
      _num_work_groups = 1;

      // set number of interlevel groups manually to 0
      _num_interlevel_groups = 0;

      // allocate arrays
      _num_proc_in_ext_work_group = new unsigned int[_num_work_groups];
      _group_contains_extra_coord = new unsigned char[_num_work_groups];
      _ext_work_group_ranks = new int*[_num_work_groups];
      _graphs = new Graph*[_num_work_groups];

      // shortcut to the number of processes in the manager's process group
      unsigned int num_processes = _process_group->num_processes();

      // shortcut to the number of cells in the base mesh
      unsigned int num_cells = unsigned int(_base_mesh->num_cells());

      // debug output
      Logger::log("num_processes: " + stringify(num_processes) + "\nnum_cells: " + stringify(num_cells) + "\n",
        Logger::master);
      // if there is dedicated load balancer process, substract one process from the number of processes (to get the
      // actual number of compute processes)
      if(_group_has_dedicated_load_bal)
      {
        --num_processes;
      }
      // assert that the number of compute processes is n
      ASSERT(num_processes == num_cells, "Number of processes " + stringify(num_processes)
             + " must equal number of cells " + stringify(num_cells) + ".");

      // set up the test case

      // no extra coordinator process is needed
      _group_contains_extra_coord[0] = false;
      // set number of processes in the work group
      _num_proc_in_ext_work_group[0] = num_cells;
      // work group ranks
      _ext_work_group_ranks[0] = new int[_num_proc_in_ext_work_group[0]];
      // set entries to {0, ..., n-1}
      for(unsigned int i(0) ; i < _num_proc_in_ext_work_group[0] ; ++i)
      {
        _ext_work_group_ranks[0][i] = i;
      }

      // get connectivity graph of the base mesh; this one will be used for the fine grid work group
      _graphs[0] = _base_mesh->graph();

      _graph_has_been_created_here = new bool[_num_work_groups];
      _graph_has_been_created_here[0] = false;
    }
  };
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_LOAD_BALANCER_HPP
