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

namespace FEAST
{

  /**
  * \brief class defining a load balancer
  *
  * Each initial process group is organised by a load balancer. It runs on all processes of the process group, which
  * eases work flow organisation significantly. There is one coordinator process which is the only one knowing the
  * complete computational mesh. It is responsible for reading and distributing the mesh to the other processes, and for
  * organising partitioning and load balancing (collect and process matrix patch statistics, ...). This coordinator is
  * always the process with the largest rank within the process group.
  *
  * The user can choose if this coordinator process should also perform compute tasks (solving linear systems etc),
  * or if it should be a dedicated load balancing / coordinator process doing nothing else. This means, the coordinator
  * process and the dedicated load balancer process, if the latter exists, coincide.
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
  * The load bal. with id 0 in the example above then...
  * 1) ...reads in the mesh (this is only done by the dedicated load balancer process or the coordinator, resp.)
  * 2) ...defines which base mesh cells (BMCs) build a matrix patch (MP) and which MPs build a process patch (PP)

COMMENT_HILMAR: Currently, only perform the most simple case: BMC = MP = PP, i.e. one BMC per MPI job

  * 3) ...decides with how many processes (usually equals the number of virtual coarse grid matrix patches, CGMPs) to
  *    solve the global coarse grid problem and how to distribute the MPs to these processes/CGMPs.
  *    The standard case will be that the coarse grid problem is solved on one processor, i.e. one CGMP containing
  *    all MPs.
  * 4) ...builds the necessary WorkGroup objects (e.g., one WorkGroup for the fine mesh problems and one for the
  *    coarse mesh problem) and creates corresponding MPI communicators. The workers with the highest work group ranks
  *    are the coordinators of these work groups, which communicate with the master or the dedicated load balancer. The
  *    process topologies of the work groups are then optimised by building corresponding graph structures. There are
  *    two cases:
  *    a) There is a dedicated load balancer: The dedicated load balancer reads the mesh, creates work groups and a
  *       global graph structure for each work group. Then it distributes to each process of the work group the relevant
  *       parts of the global graph. Each work group process then creates its local graph structure and calls
  *       MPI_Dist_graph_create(...) to build the new MPI process topology in a distributed fashion.
  *    b) There is no dedicated load balancer: Same as case a), but instead of the dedicated load balancer the
  *       coordinator of the process group builds the global graph structure. In this case, it must be distinguished
  *       whether the coordinator is part of the work group or not.
  *    The load balancer does not create the work group objects directly, but "intermediate" groups, the so called
  *    ProcessSubgroup objects. These objects then contain the actual work groups. (See the description of class
  *    ProcessSubgroup.)
  *    COMMENT_HILMAR: It might be more clever to also let the MPI implementation decide on which physical process the
  *    dedicated load balancer should reside (instead of pinning it to the last rank in the process group). To improve
  *    this is task of the ITMC.
  * 5) ...tells each work group which other work groups it has to communicate with (via the communicator they all share
  *    within the parent process group). E.g., the fine mesh work group has to send the restricted defect vector to the
  *    coarse mesh work group, while the coarse mesh work group has to send the coarse mesh correction to the fine mesh
  *    work group. Two such communicating workers live either on the same process (internal communication = copy) or
  *    on different processes (external communication = MPI send/recv). (See example below.)
  *
  * 6) ...sends corresponding parts of the mesh to the work groups
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

    /// vector of work process subgroups the load balancer manages
    std::vector<ProcessSubgroup*> _subgroups;

    /// vector of graph structures representing the process topology within the work groups
    std::vector<Graph*> _graphs;

    /// number of work groups
    unsigned int _num_subgroups;

    /**
    * \brief array of number of workers in each work group
    *
    * Dimension: [#_num_subgroups]
    */
    unsigned int* _num_proc_in_subgroup;

    /**
    * \brief 2-dim. array for storing the process group ranks building the subgroups
    *
    * Dimension: [#_num_subgroups][#_num_proc_in_subgroup[\a group_id]]
    */
    int** _subgroup_ranks;

    /**
    * \brief boolean array indicating to which work groups this process belongs
    *
    * Dimension: [#_num_subgroups]
    */
    bool* _belongs_to_group;

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
        _num_proc_in_subgroup(nullptr),
        _subgroup_ranks(nullptr),
        _belongs_to_group(nullptr),
        _base_mesh(nullptr)
    {
    }

    /// destructor
    ~LoadBalancer()
    {
      if (_subgroup_ranks != nullptr)
      {
        for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)
        {
          delete [] _subgroup_ranks[igroup];
        }
        delete [] _subgroup_ranks;
        _subgroup_ranks = nullptr;
      }

      if (_belongs_to_group != nullptr)
      {
        delete [] _belongs_to_group;
        _belongs_to_group = nullptr;
      }

      for(unsigned int igroup(0) ; igroup < _subgroups.size() ; ++igroup)
      {
        delete _subgroups[igroup];
      }

      if (_num_proc_in_subgroup != nullptr)
      {
        delete [] _num_proc_in_subgroup;
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
    /// dummy function in preparation of a function reading in a mesh file
    void read_mesh(std::string const & mesh_file)
    {
      // the mesh is read by the process group coordinator
      if(_process_group->is_coordinator())
      {
        _base_mesh = new BaseMesh();
        // Hilmar's dummy routine
        //_base_mesh->dummy_read_mesh();
        // or Dominik's 'real' one
        try
        {
          Logger::log_master("Reading mesh file " + mesh_file + "...\n", Logger::SCREEN_FILE);
          _base_mesh->read_mesh(mesh_file);
        }
        catch (Exception& e)
        {
          ErrorHandler::exception_occured(e, ErrorHandler::CRITICAL);
        }
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
    *
    * To optimise the communication between the coordinator of the main process group and the work groups, we add
    * this coordinator to a work group if it is not a compute process of this work group anyway. Hence, for each work
    * group, there are three different possibilities:
    * 1) There is a dedicated load balancer process, which is automatically the coordinator of the main process group
    *    and belongs to no work group:
    *    --> work group adds the coordinator as extra process
    * 2) There is no dedicated load balancer process, and the coordinator of the main process group ...
    *   a) ... is not part of the work group:
    *     --> work group adds the coordinator as extra process
    *   b) ... is part of the work group:
    *     --> work group does not have to add the coordinator
    * Thus the 1-to-n or n-to-1 communication between coordinator and n work group processes can be performed via
    * MPI_Scatter() and MPI_Gather() (which always have to be called by all members of an MPI process group). This
    * is more efficient then using n calls of MPI_send() / MPI_recv() via the communicator of the main process group.
    */
    void create_subgroups()
    {
      // shortcut to the number of processes in the load balancer's process group
      unsigned int num_processes = _process_group->num_processes();

      /* ********************************************************************************
      * The coordinator is the only one knowing the base mesh, so only the coordinator  *
      * decides over the number of work groups and the process distribution to them.    *
      **********************************************************************************/

      // COMMENT_HILMAR:
      // The following test code is semi hard-wired in the sense that we schedule one BMC to each processor.
      // Furthermore, two extra processes are required for testing the creation of a second work group and the use
      // of a dedicated load balancing process.
      // Later, this has to be all done in some auto-magically way.

      // Two different tests can be performed (n is the number of base mesh cells):
      // 1) with dedicated load balancer process
      //    - work group for coarse grid: 2 processes: {0, 1}
      //    - work group for fine grid: n processes: {1, ..., n}
      //    - i.e. process 1 is in both work groups
      //    - dedicated load balancer and coordinator process: n+1
      // 2) without dedicated load balancer process
      //    - work group for coarse grid: 2 processes: {0, 1}
      //    - work group for fine grid: n processes: {2, ..., n+1}
      //    - i.e. the two work groups are disjunct
      //    - coordinator process: n+1
      // Both tests need n+2 processes in total. To choose the test, change the first entry of the boolean array
      // includes_dedicated_load_bal[] in universe_test.cpp.

      bool* group_contains_extra_coord;

      if(_process_group->is_coordinator())
      {
        unsigned int num_cells;
        num_cells = _base_mesh->num_cells();
        // debug output
        Logger::log_master("num_cells: " + StringUtils::stringify(num_cells) + "\n");
        // assert that the number of processes is n+2
        assert(num_processes == num_cells + 2);

        // set up the two test cases

        // number of work groups, manually set to 2
        _num_subgroups = 2;
        // array of numbers of processes per work group
        _num_proc_in_subgroup = new unsigned int[2];

        // Boolean array indicating whether the work groups contain an extra process for the coordinator (which will then
        // not be a compute process in this work group)
        group_contains_extra_coord = new bool[_num_subgroups];

        // allocate first dimension of the array for rank partitioning
        _subgroup_ranks = new int*[_num_subgroups];

        if(_group_has_dedicated_load_bal)
        {
          // test case 1 (n is the number of base mesh cells)
          // with dedicated load balancer process
          //  - work group for coarse grid: 2 processes: {0, 1}
          //  - work group for fine grid: n processes: {1, ..., n}
          //  - i.e. process 1 is in both work groups
          //  - dedicated load balancer and coordinator process: n+1

          // since there is a dedicated load balancer process, this has to be added to both work groups as extra
          // coordinator process
          group_contains_extra_coord[0] = true;
          group_contains_extra_coord[1] = true;

          // set number of processes per group
          _num_proc_in_subgroup[0] = 2 + 1;
          _num_proc_in_subgroup[1] = num_cells + 1;

          // partition the process group ranks into work groups
          // coarse grid work group
          _subgroup_ranks[0] = new int[_num_proc_in_subgroup[0]];
          _subgroup_ranks[0][0] = 0;
          _subgroup_ranks[0][1] = 1;
          _subgroup_ranks[0][2] = num_cells + 1;

          // fine grid work group
          _subgroup_ranks[1] = new int[_num_proc_in_subgroup[1]];
          // set entries to {1, ..., n+1}
          for(unsigned int i(0) ; i < _num_proc_in_subgroup[1] ; ++i)
          {
            _subgroup_ranks[1][i] = i+1;
          }
        }
        else
        {
          // test case 2 (n is the number of base mesh cells)
          // without dedicated load balancer process
          //  - work group for coarse grid: 2 processes: {0, 1}
          //  - work group for fine grid: n processes: {2, ..., n+1}
          //  - i.e. the two work groups are disjunct
          //  - coordinator process: n+1

          // the coordinator is at the same time a compute process of the second work group, so only the first work group
          // has to add an extra process
          group_contains_extra_coord[0] = true;
          group_contains_extra_coord[1] = false;

          // set number of processes per group
          _num_proc_in_subgroup[0] = 2 + 1;
          _num_proc_in_subgroup[1] = num_cells;

          // partition the process group ranks into work groups
          _subgroup_ranks[0] = new int[_num_proc_in_subgroup[0]];
          _subgroup_ranks[0][0] = 0;
          _subgroup_ranks[0][1] = 1;
          _subgroup_ranks[0][2] = num_cells + 1;
          _subgroup_ranks[1] = new int[_num_proc_in_subgroup[1]];
          // set entries to {2, ..., n+1}
          for(unsigned int i(0) ; i < _num_proc_in_subgroup[1] ; ++i)
          {
            _subgroup_ranks[1][i] = i+2;
          }
        }
      }  // if(_process_group->is_coordinator())
      /* *****************************
      * End of semi-hard-wired code. *
      *******************************/

      // now the coordinator broadcasts the relevant data to the other processes, that is:
      //   - _num_subgroups
      //   - _num_proc_in_subgroup
      //   - group_contains_extra_coord
      //   - _subgroup_ranks

      int mpi_error_code = MPI_Bcast(&_num_subgroups, 1, MPI_UNSIGNED, _process_group->rank_coord(),
                                     _process_group->comm());
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Bcast");

      if(!_process_group->is_coordinator())
      {
        _num_proc_in_subgroup = new unsigned int[_num_subgroups];
        group_contains_extra_coord = new bool[_num_subgroups];
      }

      mpi_error_code = MPI_Bcast(_num_proc_in_subgroup, _num_subgroups, MPI_UNSIGNED, _process_group->rank_coord(),
                                 _process_group->comm());
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Bcast");

      mpi_error_code = MPI_Bcast(group_contains_extra_coord, _num_subgroups, MPI_LOGICAL, _process_group->rank_coord(),
                                 _process_group->comm());
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Bcast");

      // let the non-coordinator processes allocate the array _subgroup_ranks for rank partitioning
      if(!_process_group->is_coordinator())
      {
        // debug output
//        std::cout << "Process " << _process_group->rank() << "received _num_subgroups=" << _num_subgroups << std::endl;

        _subgroup_ranks = new int*[_num_subgroups];
        for (unsigned int i(0) ; i<_num_subgroups ; ++i)
        {
          _subgroup_ranks[i] = new int[_num_proc_in_subgroup[i]];
          // debug output
//          std::cout << "Process " << _process_group->rank() << "received _num_proc_in_subgroup[" << i << "] = "
//                    << _num_proc_in_subgroup[i] << std::endl;
        }
      }

// COMMENT_HILMAR:
// The most elegant way to broadcasting the 2D array _subgroup_ranks is by defining a corresponding MPI datatype.
// The following code should do it, but it doesn't... :-)
// Until I find the bug, the workaround code below does the trick.

//      // create MPI datatype for sending the 2D array _subgroup_ranks
//      MPI_Aint base;
//      MPI_Address(_subgroup_ranks[0], &base);
//      MPI_Aint displacements[_num_subgroups];
//      for (unsigned int i(0) ; i<_num_subgroups ; ++i)
//      {
//        MPI_Address(_subgroup_ranks[i], &displacements[i]);
//        displacements[i] -= base;
//      }
//
//      MPI_Datatype iarray2D;
//      MPI_Type_hindexed(_num_subgroups, _num_proc_in_subgroup, displacements, MPI_UNSIGNED, &iarray2D);
//      MPI_Type_commit(&iarray2D);
//
//      // let the coordinator send the array to all non-coordinator processes
//      MPI_Bcast(_subgroup_ranks, 1, iarray2D, _process_group->rank_coord(), _process_group->comm());
//
//      MPI_Type_free(&iarray2D);

/* ****************
* workaround code *
******************/

      // let all processes create a 1D array which will hold a copy of the 2D array to be broadcast
      unsigned int total_length(0);
      int* _subgroup_ranks_1D;
      for (unsigned int i(0) ; i<_num_subgroups ; ++i)
      {
        total_length += _num_proc_in_subgroup[i];
      }
      _subgroup_ranks_1D = new int[total_length];

      // let the coordinator copy the 2D array into a 1D array
      if(_process_group->is_coordinator())
      {
        unsigned int pos(0);
        for (unsigned int i(0) ; i<_num_subgroups ; ++i)
        {
          for(unsigned int j(0) ; j < _num_proc_in_subgroup[i] ; ++j)
          {
            _subgroup_ranks_1D[pos] = _subgroup_ranks[i][j];
            ++pos;
          }
        }
      }

      // broadcast the 1D array
      mpi_error_code = MPI_Bcast(_subgroup_ranks_1D, total_length, MPI_INTEGER, _process_group->rank_coord(),
                                 _process_group->comm());
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Bcast");

      // let the non-coordinator processes copy the 1D array into the 2D array
      if(!_process_group->is_coordinator())
      {
        unsigned int pos(0);
        for (unsigned int i(0) ; i<_num_subgroups ; ++i)
        {
          for(unsigned int j(0) ; j < _num_proc_in_subgroup[i] ; ++j)
          {
            _subgroup_ranks[i][j] = _subgroup_ranks_1D[pos];
            ++pos;
          }
        }
      }
      // and delete the temporary 1D array again
      delete [] _subgroup_ranks_1D;

/* ***********************
* end of workaround code *
*************************/

//      // debug output
//      for (unsigned int i(0) ; i<_num_subgroups ; ++i)
//      {
//        std::string s = StringUtils::stringify(_subgroup_ranks[i][0]);
//        for(unsigned int j(1) ; j < _num_proc_in_subgroup[i] ; ++j)
//        {
//          s +=  " " + StringUtils::stringify(_subgroup_ranks[i][j]);
//        }
//        std::cout << "Process " << _process_group->rank() << "received _subgroup_ranks[" << i << "] = "
//                  << s << std::endl;
//      }

      /* *********************************
      * now begin creating the subgroups *
      ***********************************/

      // boolean array indicating to which work groups this process belongs
      _belongs_to_group = new bool[_num_subgroups];
      for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)
      {
        // intialise with false
        _belongs_to_group[igroup] = false;
        for(unsigned int j(0) ; j < _num_proc_in_subgroup[igroup] ; ++j)
        {
          if(_process_group->rank() == _subgroup_ranks[igroup][j])
          {
            _belongs_to_group[igroup] = true;
          }
        }
      } // for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)

      // create ProcessSubgroup objects including MPI groups and MPI communicators
      // It is not possible to set up all subgroups in one call, since the processes building the subgroups are
      // not necessarily disjunct. Hence, there are as many calls as there are subgroups. All processes not belonging
      // to the subgroup currently created call the MPI_Comm_create() function with a with dummy communicator and
      // dummy group.

      _subgroups.resize(_num_subgroups, nullptr);
      for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)
      {
        if(_belongs_to_group[igroup])
        {
          _subgroups[igroup] = new ProcessSubgroup(_num_proc_in_subgroup[igroup], _subgroup_ranks[igroup],
                                                   _process_group, igroup, group_contains_extra_coord[igroup]);
        }
        else
        {
          // *All* processes of the parent MPI group have to call the MPI_Comm_create() routine (otherwise the forking
          // will deadlock), so let all processes that are not part of the current work group call it with dummy
          // communicator and dummy group.  (The dummy group is necessary here since the other MPI_Group object is
          // hidden inside the ProcessSubgroup constructor above.)
          MPI_Comm dummy_comm;
          MPI_Group dummy_group;
          int rank_aux = _process_group->rank();
          int mpi_error_code = MPI_Group_incl(_process_group->group(), 1, &rank_aux, &dummy_group);
          MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_incl");
          mpi_error_code = MPI_Comm_create(_process_group->comm(), dummy_group, &dummy_comm);
          MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_create");

          // COMMENT_HILMAR: First, I used this simpler version:
          //   mpi_error_code = MPI_Comm_create(_process_group->comm(), MPI_GROUP_EMPTY, &dummy_comm);
          // It worked with OpenMPI 1.4.2 and MPICH2, but does not with OpenMPI 1.4.3. We are not quite sure yet, if that
          // is a bug in OpenMPI 1.4.3, or if this use of MPI_GROUP_EMPTY is incorrect.

          // Within the ProcessSubgroup constructor, another MPI communicator is created. So, do another dummy call.
          mpi_error_code = MPI_Comm_create(_process_group->comm(), dummy_group, &dummy_comm);
          MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_create");

          MPI_Comm_free(&dummy_comm);
          MPI_Group_free(&dummy_group);
        }
      } // for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)

      // delete aux. array
      delete [] group_contains_extra_coord;


      /* *********************************************************
      * create graph structures corresponding to the work groups *
      ***********************************************************/

      // let the coordinator create the process topology
      if(_process_group->is_coordinator())
      {
        _graphs.resize(_num_subgroups, nullptr);

        // build an artificial graph mimicing the distribution of the 16 base mesh cells to two processors
        // (e.g. BMCs 0-7 on proc 1 and BMCs 8-15 on proc 2) which start an imagined coarse grid solver; this graph will
        // be used for the coarse grid work group
        unsigned int* index = new unsigned int[3];
        unsigned int* neighbours = new unsigned int[2];
        index[0] = 0;
        index[1] = 1;
        index[2] = 2;
        neighbours[0] = 1;
        neighbours[1] = 0;
        _graphs[0] = new Graph(2, index, neighbours);
        _graphs[0]->print();

        // get connectivity graph of the base mesh; this one will be used for the fine grid work group
        _graphs[1] = _base_mesh->graph();
// COMMENT_HILMAR:
// We assume here that each process receives exactly one BMC and that the index of the cell in the graph structure
// equals the local rank within the work group.
      }

      /* ***************************************************************************
      * now let the coordinator send the relevant parts of the global graph to the *
      * corresponding work group members                                           *
      *****************************************************************************/

      for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)
      {
        if(_belongs_to_group[igroup])
        {
          unsigned int num_neighbours_local;
          unsigned int* neighbours_local;
          int root = _subgroups[igroup]->rank_coord();
          if(_subgroups[igroup]->is_coordinator())
          {
            /* **********************************
            * code for the sending root process *
            ************************************/

            unsigned int count = _graphs[igroup]->num_nodes();
            if(_subgroups[igroup]->contains_extra_coord())
            {
              // MPI_Scatter() counts the extra coordinator process as well (which only sends data). Arrays have to
              // be prepared with n+1 segments although only n processes receive data.
              ++count;
            }
            unsigned int num_neighbours[count];
            unsigned int* index = _graphs[igroup]->index();
            for(unsigned int i(0) ; i < _graphs[igroup]->num_nodes() ; ++i)
            {
              num_neighbours[i] = index[i+1] - index[i];
            }

            if(_subgroups[igroup]->contains_extra_coord())
            {
              // the extra coordinator process is the last in the work rank, so set the last entry of the array to zero
              num_neighbours[count-1] = 0;

              // send the number of neighbours to the non-root processes (use MPI_IN_PLACE to indicate that the root
              // does not receive/store any data)
              MPI_Scatter(num_neighbours, 1, MPI_UNSIGNED, MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, root,
                          _subgroups[igroup]->comm());
              // send the neighbours to the non-root processes (usually the neighbours array must also have n+1
              // segments, but since num_neighbours[count-1] == 0, the last segment is empty)
              MPI_Scatterv(_graphs[igroup]->neighbours(), reinterpret_cast<int*>(num_neighbours),
                           reinterpret_cast<int*>(index), MPI_UNSIGNED, MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, root,
                           _subgroups[igroup]->comm());
// COMMENT_HILMAR:
// Is it problematic to define num_neighbours as unsigned int* and then to cast it to int* when passed to
// MPI_Scatterv()? (This question also applies to other places in the code!)
// Should we better define it as int* right from the beginning?
            }
            else
            {
              // When there is no extra coordinator process, then the root is part of the compute work group and
              // also sends data to itself.

              // scatter the number of neighbours to the non-root processes and to the root process itself
              MPI_Scatter(num_neighbours, 1, MPI_UNSIGNED, &num_neighbours_local, 1, MPI_UNSIGNED, root,
                          _subgroups[igroup]->comm());
              neighbours_local = new unsigned int[num_neighbours_local];
              // scatter the neighbours to the non-root processes and to the root process itself
              MPI_Scatterv(_graphs[igroup]->neighbours(), reinterpret_cast<int*>(num_neighbours),
                           reinterpret_cast<int*>(index), MPI_UNSIGNED,neighbours_local, num_neighbours_local,
                           MPI_UNSIGNED, root,_subgroups[igroup]->comm());
            }
          }
          else
          {
            /* ******************************************
            * code for the receiving non-root processes *
            ********************************************/
            // receive the number of neighbours from the root process
            MPI_Scatter(nullptr, 0, MPI_DATATYPE_NULL, &num_neighbours_local, 1, MPI_UNSIGNED, root,
                        _subgroups[igroup]->comm());

            // receive the neighbours
            neighbours_local = new unsigned int[num_neighbours_local];
            MPI_Scatterv(nullptr, 0, nullptr, MPI_DATATYPE_NULL, neighbours_local, num_neighbours_local, MPI_UNSIGNED,
                         root, _subgroups[igroup]->comm());
          }

          if (!(_subgroups[igroup]->is_coordinator() && _subgroups[igroup]->contains_extra_coord()))
          {
            // now create distributed graph structure within the compute work groups (the array neighbours_local will be
            // deallocated in the destructor of the distributed graph object)
            _subgroups[igroup]->work_group()->set_graph_distributed(num_neighbours_local, neighbours_local);
          }
        } // if(_belongs_to_group[igroup])
      } // for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)

      // test local neighbourhood communication
      for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)
      {
        if (_belongs_to_group[igroup] &&
            !(_subgroups[igroup]->is_coordinator() && _subgroups[igroup]->contains_extra_coord()))
        {
          _subgroups[igroup]->work_group()->do_exchange();
        }
      }
    } // create_subgroups()
  };

} // namespace FEAST

#endif // guard KERNEL_LOAD_BAL_HPP
