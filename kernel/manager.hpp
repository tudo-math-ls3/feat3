#pragma once
#ifndef KERNEL_MANAGER_HPP
#define KERNEL_MANAGER_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/process.hpp>
#include <kernel/process_group.hpp>
#include <kernel/process_subgroup.hpp>
#include <kernel/base_mesh/file_parser.hpp>
#include <kernel/base_mesh/bm.hpp>

/// FEAST namespace
namespace FEAST
{
// COMMENT_HILMAR: entsprechend der Idee, den load balancer besser abzugrenzen, sollte diese ganze Klasse besser
// ProcessOrganiser, ProcessManager (oder so ähnlich) genannt werden. Die eigentliche load balancer Klasse wird dann
// (erstmal) nur auf dem coordinator der zugrundeliegenden Prozessgruppe instantiiert. Sie wird dann z.B. die
// Funktion define_work_groups() enthalten. (Waehrend die Funktion create_work_groups() hier in dieser Klasse bleiben
// sollte, da sie nur umsetzt, was der load balancer ihr vorgibt.)
  /**
  * \brief class defining a process group manager
  *
  * Each initial process group is organised by a manager. It runs on all processes of the process group, which
  * eases work flow organisation significantly. There is one coordinator process which is the only one knowing the
  * complete computational mesh. It is part of at least one work group, i.e., it performs compute tasks, but it is
  * also responsible for reading and distributing the mesh to the other processes, and for organising load balancing.
  * The coordinator has rank 0 within the process group.
  *
  * One process is responsible for performing load balancing (collect and process matrix patch statistics, compute
  * partitioning, define work groups, ...). This is either the coordinator process of the process group (being a compute
  * process at the same time), or a dedicated load balancer process (which does nothing else). If the user decides to
  * use such a decicated load balancer process, it is the process with the highest rank in the process group.
  *
  * The user knows what each process group and its respective manager should do. He calls, e.g.,
  \verbatim
  if(manager->group_id() == 0)
  {
    manager->readMesh();
    ...
  }
  else
  {
    manager->start_coffee_machine();
  }
  \endverbatim
  *
  * The manager with id 0 in the example above then...
  * -# ...reads in the mesh (this is only done by the coordinator)
  * -# ...lets its load balancer define which base mesh cells (BMCs) build a matrix patch (MP) and which MPs build a
  *    process patch (PP)
  *    COMMENT_HILMAR: Currently, only perform the most simple case: BMC = MP = PP, i.e. one BMC per MPI job
  * -# ...lets its load balancer define with how many processes (usually equals the number of virtual coarse grid
  *    matrix patches, CGMPs) to solve the global coarse grid problem and how to distribute the MPs to these
  *    processes/CGMPs. The standard case will be that the coarse grid problem is solved on one processor, i.e. one
  *    CGMP containing all MPs.
  * -# ...lets its load balancer define the necessary WorkGroup objects (e.g., one WorkGroup for the fine mesh problems
  *    and one for the coarse mesh problem) and creates corresponding MPI communicators. The workers with rank 0 within
  *    the work group (being itself a process group) are the coordinators of these work groups, which communicate with
  *    the master or the dedicated load balancer. The process topologies of the work groups are then optimised by
  *    building corresponding graph structures. There are two cases:
  *    -# There is a dedicated load balancer: The dedicated load balancer gets all necessary information from the
  *       coordinator (e.g., graph describing the neighourhood of the base mesh cells), defines MP and PP partitioning.
  *       Then it returns the necessary information to the coordinator of the process group which organises the actual
  *       setup of work groups and process distribution, e.g. by distributing to each worker process the relevant
  *       parts of the global graph. Each worker process then creates its local graph structure and calls
  *       MPI_Dist_graph_create(...) to build the new MPI process topology in a distributed fashion.
  *    -# There is no dedicated load balancer: Same as case a), but the tasks of the dedicated load balancer are
  *       performed by the coordinator of the process group.
  *    In both cases, one must distinguish whether the coordinator is part of the specific work group currently set up.
  *    The manager does not create the work group objects directly, but 'intermediate' groups, the so called
  *    ProcessSubgroup objects. These objects then contain the actual work groups and eventually additionally the
  *    coordinator process if it is not part of the work group anyway (see the description of class
  *    ProcessSubgroup for details).
  *    COMMENT_HILMAR: It might be more clever to also let the MPI implementation decide on which physical process the
  *    dedicated load balancer should reside (instead of pinning it to the last rank in the process group). To improve
  *    this is task of the ITMC.
  * -# ...tells each work group which other work groups it has to communicate with (via the communicator they all share
  *    within the parent process group). E.g., the fine mesh work group has to send the restricted defect vector to the
  *    coarse mesh work group, while the coarse mesh work group has to send the coarse mesh correction to the fine mesh
  *    work group. Two such communicating workers live either on the same process (internal communication = copy) or
  *    on different processes (external communication = MPI send/recv). (See example below.)
  * -# ...sends corresponding parts of the mesh to the work groups
  *
  * Example:
  *
  * Distribution of submeshes to processes A-G on different levels (note that processes A-G are not necessarily
  * disjunct, i.e., several of them can refer to the same physical process, see cases a and b):
  *
  \verbatim
  ---------------      ---------------      ---------------
  |             |      |      |      |      |      |      |
  |             |      |      |      |      |  D   |  G   |
  |             |      |      |      |      |      |      |
  |      A      |      |  B   |  C   |      ---------------
  |             |      |      |      |      |      |      |
  |             |      |      |      |      |  E   |  F   |
  |             |      |      |      |      |      |      |
  ---------------      ---------------      ---------------
    level 0               level 1              levels 2-L
  \endverbatim
  *
  * -# case a, four physical processes:
  \verbatim
  process group rank:  0  1  2  3
         WorkGroup 2:  D  E  F  G         (four WorkGroup processes for the problems on level 2-L)
         WorkGroup 1:  B     C            (two WorkGroup processes for the problem on level 1)
         WorkGroup 0:  A                  (one WorkGroup process for the coarse mesh problem on level 0)
  \endverbatim
  *
  *   Communication:
  *   A <--> B (internal, rank 0) A <--> C (external, ranks 0+2)
  *   B <--> D (internal, rank 0) B <--> E (external, ranks 0+1)
  *   C <--> F (internal, rank 2) C <--> G (external, ranks 2+3)
  *
  * -# case b, five physical processes::
  \verbatim
  process group rank:  0  1  2  3  4
         WorkGroup 2:     D  E  F  G
         WorkGroup 1:     B     C
         WorkGroup 0:  A
  \endverbatim
  *
  *   Communication:
  *   A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+3)
  *   B <--> D (internal, rank 1)    B <--> E (external, ranks 1+2)
  *   C <--> F (internal, rank 3)    C <--> G (external, ranks 3+4)
  *
  * -# case c, seven physical processes:
  \verbatim
  process group rank:  0  1  2  3  4  5  6
         WorkGroup 2:           D  E  F  G
         WorkGroup 1:     B  C
         WorkGroup 0:  A
  \endverbatim
  *
  *   Communication:
  *   A <--> B (external, ranks 0+1) A <--> C (external, ranks 0+2)
  *   B <--> D (external, ranks 1+3) B <--> E (external, ranks 1+4)
  *   C <--> F (external, ranks 2+5) C <--> G (external, ranks 2+6)
  *
  * \tparam space_dim_
  * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
  *
  * \tparam world_dim_
  * world dimension (determines the number of coordinates)
  *
  * \note Since a manager manages one (global) process group, it could also have been called 'ProcessGroupManager'.
  * But since 99 % of the code 'happens' within process groups (only application programmers dealing with multiphysics
  * or similar stuff have to organise different process groups), I omitted the leading 'ProcessGroup'.
  * Other names which came to my mind were 'ProgramFlowManager' (too long) and 'FlowManager' (misleading in a FE tool
  * dealing with CFD problems). (Hilmar)
  *
  * \author Hilmar Wobker
  * \author Dominik Goeddeke
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class Manager
  {

  private:

    /* *****************
    * member variables *
    *******************/
    /// process group the manager manages
    ProcessGroup* _process_group;

    /// flag whether the manager's process group uses a dedicated load balancer process
    bool const _group_has_dedicated_load_bal;

    /**
    * \brief rank of the load balancer process
    *
    * If the process group uses a dedicated load balancer process, then it has rank #_num_processes - 1. Otherwise,
    * it coincides with the coordinator process, i.e. it has rank 0.
    */
    int _rank_load_balancer;

    /// vector of process subgroups the manager manages
    std::vector<ProcessSubgroup*> _subgroups;

    /// vector of graph structures representing the process topology within the work groups
    std::vector<Graph*> _graphs;

    /// number of subgroups / work groups
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

    /// base mesh the manager works with
    BaseMesh::BM<space_dim_, world_dim_>* _base_mesh;

  public:

    /* *************************
    * constructor & destructor *
    ***************************/
    /// CTOR
    Manager(
      ProcessGroup* process_group,
      bool group_has_dedicated_load_bal)
      : _process_group(process_group),
        _group_has_dedicated_load_bal(group_has_dedicated_load_bal),
        _rank_load_balancer(MPI_PROC_NULL),
        _num_subgroups(0),
        _num_proc_in_subgroup(nullptr),
        _subgroup_ranks(nullptr),
        _belongs_to_group(nullptr),
        _base_mesh(nullptr)
    {
      CONTEXT("Manager::Manager()");
      // set rank of the load balancer process
      if(_group_has_dedicated_load_bal)
      {
        // if the group contains a dedicated load balancer, it is the process with the highest rank in the process group
        _rank_load_balancer = _process_group->num_processes()-1;
      }
      else
      {
        // otherwise the coordinator plays the role of the load balancer
        _rank_load_balancer = _process_group->rank_coord();
      }
    }

    /// DTOR
    ~Manager()
    {
      CONTEXT("Manager::~Manager()");
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

      while (!_subgroups.empty())
      {
        // call the destructor of the last element in the vector
        delete _subgroups.back();
        // delete the pointer from the vector
        _subgroups.pop_back();
      }

      if (_num_proc_in_subgroup != nullptr)
      {
        delete [] _num_proc_in_subgroup;
      }

      if (_base_mesh != nullptr)
      {
        delete _base_mesh;
      }

      // assume that the _graphs vector only holds copies of the graph object pointers and that the graph objects
      // themselves are destroyed where they were created
      _graphs.clear();

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
      CONTEXT("Manager::process_group()");
      return _process_group;
    }


    /**
    * \brief getter for the base mesh
    *
    * \return BaseMesh::BM<space_dim_, world_dim_> pointer #_base_mesh
    */
    inline BaseMesh::BM<space_dim_, world_dim_>* base_mesh() const
    {
      CONTEXT("Manager::base_mesh()");
      return _base_mesh;
    }


    /**
    * \brief getter for the flag whether the process group has a dedicated load balancer
    *
    * \return BaseMesh::BM<space_dim_, world_dim_> pointer #_base_mesh
    */
    inline bool group_has_dedicated_load_bal() const
    {
      CONTEXT("Manager::group_has_dedicated_load_bal()");
      return _group_has_dedicated_load_bal;
    }


    /**
    * \brief getter for rank of the load balancer process
    *
    * If the process group uses a dedicated load balancer process, then it has rank #_num_processes - 1. Otherwise,
    * it coincides with the coordinator process, i.e. it has rank 0.
    *
    * \return rank of the load balancer process
    */
    inline int rank_load_balancer() const
    {
      CONTEXT("Manager::rank_load_balancer()");
      return _rank_load_balancer;
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief checks whether this process is the load balancer of the process group
    *
    * \return true if this process is the load balancer of the process group, otherwise false
    */
    inline bool is_load_balancer() const
    {
      CONTEXT("ProcessGroup::is_load_balancer()");
      return _process_group->rank() == _rank_load_balancer;
    }


    /**
    * \brief checks whether this process is the dedicated load balancer of the process group
    *
    * \return true if this process is the dedicated load balancer of the process group, otherwise false
    */
    inline bool is_dedicated_load_balancer() const
    {
      CONTEXT("ProcessGroup::is_dedicated_load_balancer()");
      return group_has_dedicated_load_bal() && is_load_balancer();
    }


    /// read in a mesh file and set up base mesh
    void read_mesh(std::string const & mesh_file)
    {
      CONTEXT("Manager::read_mesh()");
      // the mesh is read by the process group coordinator
      if(_process_group->is_coordinator())
      {
        _base_mesh = new BaseMesh::BM<space_dim_, world_dim_>();
        BaseMesh::FileParser<space_dim_, world_dim_> parser;
        Logger::log_master("Reading mesh file " + mesh_file + "...\n", Logger::SCREEN_FILE);
        try
        {
          parser.parse(mesh_file, _base_mesh);
        }
        catch(Exception& e)
        {
          // abort the program
          ErrorHandler::exception_occured(e);
        }
        // set cell numbers (equal to indices since all cells are active)
        _base_mesh->set_cell_numbers();
        // create base mesh's graph structure
        _base_mesh->create_graph();
        // print base mesh
        std::string s = _base_mesh->print();
        Logger::log_master(s, Logger::SCREEN);
        Logger::log(s);
        // validate base mesh
        _base_mesh->validate(Logger::file);
      }
    }

    /**
    * \brief dummy function in preparation of a function for defining work groups
    *
    * As long as we cannot do this in an automatic way, the test driver calls a corresponding function, which sets
    * work groups manually.
    */
    void define_work_groups()
    {
      CONTEXT("Manager::define_work_groups()");
    }


    /**
    * \brief function that sets up subgroups and work groups basing on the provided information
    *
    * This function is called on all processes of the process group.
    *
    * \param[in] num_subgroups
    * number of subgroups
    *
    * \param[in] num_proc_in_subgroup
    * array of numbers of processes per work group
    *
    * \param[in] group_contains_extra_coord
    * Array indicating whether the subgroups contain an extra process for the coordinator (which will then
    * not be a compute process in the corresponding work group).
    * Usually a boolean would do the trick, but we cannot send boolean arrays via MPI. (C++ bindings are deprecated,
    * hence we cannot use MPI::BOOL. MPI_LOGICAL cannot be used, since this is not equivalent to C++'s bool.)
    * So, we use unsigned char here and treat is as if it were boolean.
    *
    * \param[in] subgroup_ranks
    * array for defining the rank partitioning
    *
    * \param[in] graphs
    * array of Graph pointers representing the connectivity of work group processes
    */
    void create_work_groups(
      unsigned int num_subgroups,
      unsigned int* num_proc_in_subgroup,
      unsigned char* group_contains_extra_coord,
      int** subgroup_ranks,
      Graph** graphs)
    {
      CONTEXT("Manager::create_work_groups()");
      // set data/pointers on the coordinator process (where the data is already available)
      if(_process_group->is_coordinator())
      {
        _num_subgroups = num_subgroups;
        _num_proc_in_subgroup = num_proc_in_subgroup;
        _subgroup_ranks = subgroup_ranks;
      }

      // now the coordinator broadcasts the relevant data to the other processes, that is:
      //   - _num_subgroups
      //   - _num_proc_in_subgroup
      //   - group_contains_extra_coord
      //   - _subgroup_ranks

      int mpi_error_code = MPI_Bcast(&_num_subgroups, 1, MPI_UNSIGNED, _process_group->rank_coord(),
                                     _process_group->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

      // allocate arrays on the non-coordinator processes
      if(!_process_group->is_coordinator())
      {
        ASSERT(_num_subgroups > 0, "Number of subgroups must be greater than 0.");
        _num_proc_in_subgroup = new unsigned int[_num_subgroups];
        group_contains_extra_coord = new unsigned char[_num_subgroups];
      }

      mpi_error_code = MPI_Bcast(_num_proc_in_subgroup, _num_subgroups, MPI_UNSIGNED, _process_group->rank_coord(),
                                 _process_group->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

      mpi_error_code = MPI_Bcast(group_contains_extra_coord, _num_subgroups, MPI_UNSIGNED_CHAR,
                                 _process_group->rank_coord(), _process_group->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

      // let the non-coordinator processes allocate the array _subgroup_ranks for rank partitioning
      if(!_process_group->is_coordinator())
      {
        _subgroup_ranks = new int*[_num_subgroups];
        for (unsigned int i(0) ; i<_num_subgroups ; ++i)
        {
          _subgroup_ranks[i] = new int[_num_proc_in_subgroup[i]];
        }
      }

// COMMENT_HILMAR:
// The most elegant way to broadcasting the 2D array _subgroup_ranks is by defining a corresponding MPI datatype.
// The following code should do it, but it doesn't... :-)
// Until I find the bug, the workaround code below does the trick.

//      // create MPI datatype for sending the 2D array _subgroup_ranks
//      MPI_Aint base;
//      MPI_Address(_subgroup_ranks[0], &base);
//      MPI_Aint* displacements = new MPI_Aint[_num_subgroups];
//      for (unsigned int i(0) ; i < _num_subgroups ; ++i)
//      {
//        MPI_Address(_subgroup_ranks[i], &displacements[i]);
//        displacements[i] -= base;
//      }
//
//      MPI_Datatype int_array_2d;
//      MPI_Type_create_hindexed(_num_subgroups, reinterpret_cast<int*>(_num_proc_in_subgroup), displacements,
//                               MPI_INTEGER, &int_array_2d);
//      MPI_Type_commit(&int_array_2d);
//
//      // let the coordinator send the array to all non-coordinator processes
//      MPI_Bcast(_subgroup_ranks, 1, int_array_2d, _process_group->rank_coord(), _process_group->comm());
//
//      MPI_Type_free(&int_array_2d);
//      delete [] displacements;

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
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

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
//        std::string s = stringify(_subgroup_ranks[i][0]);
//        for(unsigned int j(1) ; j < _num_proc_in_subgroup[i] ; ++j)
//        {
//          s +=  " " + stringify(_subgroup_ranks[i][j]);
//        }
//        s = "Process " + stringify(_process_group->rank()) + " received _subgroup_ranks[" + stringify(i) + "] = " + s;
//        _process_group->log_indiv_master(s);
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
      // to the subgroup currently created call the MPI_Comm_create() function with dummy communicator and dummy group.

      _subgroups.resize(_num_subgroups, nullptr);
      for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)
      {
        if(_belongs_to_group[igroup])
        {
          _subgroups[igroup] = new ProcessSubgroup(_num_proc_in_subgroup[igroup], _subgroup_ranks[igroup],
                                                   _process_group, igroup, (bool)group_contains_extra_coord[igroup]);
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
          validate_error_code_mpi(mpi_error_code, "MPI_Group_incl");
          mpi_error_code = MPI_Comm_create(_process_group->comm(), dummy_group, &dummy_comm);
          validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");
          MPI_Comm_free(&dummy_comm);

          // COMMENT_HILMAR: First, I used this simpler version:
          //   mpi_error_code = MPI_Comm_create(_process_group->comm(), MPI_GROUP_EMPTY, &dummy_comm);
          // It worked with OpenMPI 1.4.2 and MPICH2, but does not with OpenMPI 1.4.3. We are not quite sure yet, if
          // that is a bug in OpenMPI 1.4.3, or if this use of MPI_GROUP_EMPTY is incorrect.

          // Within the ProcessSubgroup constructor, another MPI communicator is created. So, do another dummy call.
          mpi_error_code = MPI_Comm_create(_process_group->comm(), dummy_group, &dummy_comm);
          validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");

          MPI_Comm_free(&dummy_comm);
          MPI_Group_free(&dummy_group);
        }
      } // for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)

      // delete aux. array (on all processes)
      delete [] group_contains_extra_coord;

      /* ***************************************************************************
      * now let the coordinator send the relevant parts of the global graph to the *
      * corresponding work group members                                           *
      *****************************************************************************/

      for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)
      {
        if(_belongs_to_group[igroup])
        {
          // number of neighbours of the graph node corresponding to this process (use unsigned int datatype here
          // instead of index_glob_t since MPI routines expect it)
          unsigned int num_neighbours_local;
          index_glob_t* neighbours_local(nullptr);
          int rank_coord = _subgroups[igroup]->rank_coord();
          if(_subgroups[igroup]->is_coordinator())
          {
            /* **********************************************
            * code for the sending coordinator process *
            ************************************************/

            // set the graph pointer for the current subgroup
            _graphs.push_back(graphs[igroup]);

            // since the MPI routines used below expect integer arrays, we have to copy two index_glob_t arrays
            // within the graph structures to corresponding int arrays
// COMMENT_HILMAR: Gibt es eine Moeglichkeit, das zu vermeiden? Ein reinterpret_cast<int*>(unsigned long) funzt nicht!
            unsigned int* num_neighbours_aux;
            unsigned int* index_aux;

            index_glob_t num_nodes;
            if(_subgroups[igroup]->contains_extra_coord())
            {
              // In case there is an extra coordinator process, we have to add one pseudo node to the graph and the
              // index array must be modified correspondingingly. This pseudo node corresponds to the coordinator
              // process itself which has to be included in the call of MPI_Scatterv(...). It has to appear at the
              // position in the arrays num_neighbours_aux[] and index_aux[] corresponding to the rank of the
              // coordinator. Although we know, that this is rank 0, we do not explicitly exploit this information here
              // since it might be that this is changed in future.
              num_nodes = _graphs[igroup]->num_nodes() + 1;
              index_aux = new unsigned int[num_nodes + 1];
              // copy the first part of the graphs's index array to the aux array, performing implicit cast from
              // index_glob_t to unsigned int
              for(index_glob_t i(0) ; i < (index_glob_t)rank_coord+1 ; ++i)
              {
                index_aux[i] = _graphs[igroup]->index()[i];
              }
              // insert the pseudo node
              index_aux[rank_coord+1] = index_aux[rank_coord];
              // copy the remaining part of the graphs's index array to the aux array, performing implicit cast from
              // index_glob_t to unsigned int
              for(index_glob_t i(rank_coord+1) ; i < num_nodes ; ++i)
              {
                index_aux[i+1] = _graphs[igroup]->index()[i];
              }
            }
            else
            {
              // in case there is no extra coordinator process, the number of neighbours equals the number of nodes
              // in the graph and the index array does not have to be modified
              num_nodes = _graphs[igroup]->num_nodes();
              index_aux = new unsigned int[num_nodes + 1];
              // copy the graphs's index array to the aux array, performing implicit cast from index_glob_t to
              // unsigned int
              for(index_glob_t i(0) ; i < num_nodes+1 ; ++i)
              {
                index_aux[i] = _graphs[igroup]->index()[i];
              }
            }

            // now determine the number of neighbours per node (eventually including the pseudo node for the extra
            // coordinator process)
            num_neighbours_aux = new unsigned int[num_nodes];
            for(index_glob_t i(0) ; i < num_nodes ; ++i)
            {
              num_neighbours_aux[i] = index_aux[i+1] - index_aux[i];
            }

            if(_subgroups[igroup]->contains_extra_coord())
            {
              // send the number of neighbours to the non-coordinator processes (use MPI_IN_PLACE to indicate that the
              // coordinator does not receive/store any data)
              MPI_Scatter(num_neighbours_aux, 1, MPI_UNSIGNED, MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                          rank_coord, _subgroups[igroup]->comm());
              // send the neighbours to the non-coordinator processes
              MPI_Scatterv(_graphs[igroup]->neighbours(), reinterpret_cast<int*>(num_neighbours_aux),
                           reinterpret_cast<int*>(index_aux), MPIType<index_glob_t>::value(), MPI_IN_PLACE, 0,
                           MPI_DATATYPE_NULL, rank_coord, _subgroups[igroup]->comm());
            }
            else
            {
              // When there is no extra coordinator process, then the coordinator is part of the compute work group and
              // also sends data to itself.

              // scatter the number of neighbours to the non-coordinator processes and to the coordinator process itself
              MPI_Scatter(reinterpret_cast<int*>(num_neighbours_aux), 1, MPI_UNSIGNED, &num_neighbours_local, 1,
                          MPI_INTEGER, rank_coord, _subgroups[igroup]->comm());
              neighbours_local = new index_glob_t[num_neighbours_local];
              // scatter the neighbours to the non-coordinator processes and to the coordinator process itself
              MPI_Scatterv(_graphs[igroup]->neighbours(), reinterpret_cast<int*>(num_neighbours_aux),
                           reinterpret_cast<int*>(index_aux), MPIType<index_glob_t>::value(), neighbours_local,
                           num_neighbours_local, MPIType<index_glob_t>::value(), rank_coord,
                           _subgroups[igroup]->comm());
            }
            // delete aux. arrays again
            delete [] num_neighbours_aux;
            delete [] index_aux;
          }
          else // !_subgroups[igroup]->is_coordinator()
          {
            /* *************************************************
            * code for the receiving non-coordinator processes *
            ***************************************************/
            // receive the number of neighbours from the coordinator process
            MPI_Scatter(nullptr, 0, MPI_DATATYPE_NULL, &num_neighbours_local, 1, MPI_INTEGER,
                        rank_coord, _subgroups[igroup]->comm());

            // receive the neighbours
            neighbours_local = new index_glob_t[num_neighbours_local];
            MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DATATYPE_NULL, neighbours_local,
                         num_neighbours_local, MPIType<index_glob_t>::value(), rank_coord,
                         _subgroups[igroup]->comm());
          }

          if (!(_subgroups[igroup]->is_coordinator() && _subgroups[igroup]->contains_extra_coord()))
          {
            // now create distributed graph structure within the compute work groups (the array neighbours_local is
            // copied in the constructor of the distributed graph object, hence it can be deallocated afterwards)
            _subgroups[igroup]->work_group()->set_graph_distributed(num_neighbours_local, neighbours_local);
            delete [] neighbours_local;
          }
        } // if(_belongs_to_group[igroup])
      } // for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)

// TODO: remove this test code here and call it from somewhere else
      // test local neighbourhood communication
      for(unsigned int igroup(0) ; igroup < _num_subgroups ; ++igroup)
      {
        if (_belongs_to_group[igroup] &&
            !(_subgroups[igroup]->is_coordinator() && _subgroups[igroup]->contains_extra_coord()))
        {
          _subgroups[igroup]->work_group()->do_exchange();
        }
      }
    } // create_work_groups()
  };
} // namespace FEAST

#endif // guard KERNEL_MANAGER_HPP
