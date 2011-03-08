#pragma once
#ifndef KERNEL_UNIVERSE_HPP
#define KERNEL_UNIVERSE_HPP 1

// includes, system
#include <iostream>
#include <stdlib.h>

// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/process.hpp>
#include <kernel/master.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/load_balancer.hpp>

/// FEAST namespace
namespace FEAST
{
  /**
  * \brief creates and manages the parallel universe
  *
  * At the beginning of the program the one and only object of this class is created via the static routine
  * Universe::create(...) which calls the constructor. A second call of the constructor is automatically prevented.
  * (The class is a so called 'singleton'.) The destruction of this object via the static routine Universe::destroy()
  * will also end the program. The universe object creates and manages all MPI processes. To this end, it sets up
  * MPI groups and communicators to partition the top level default universe MPI_COMM_WORLD into work units, roles,
  * and the whole lot for multiphysics applications. The normal user and kernel programmer should not have to interfere
  * with the universe for singlephysics applications.
  *
  * Idea: The soon-to-be-implemented global Feast_init() routine encapsulates the entire universe and makes physics-level
  *       subgroups available to the programmer.
  *
  * \tparam space_dim_
  * space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
  *
  * \tparam world_dim_
  * world dimension (determines the number of coordinates)
  *
  * \author Hilmar Wobker
  * \author Dominik Goeddeke
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class Universe
    : public InstantiationPolicy<Universe<space_dim_, world_dim_>, Singleton>
  {

  private:

    /* *****************
    * member variables *
    *******************/

    /// flag whether universe has already been created
    bool _universe_created;

    /// number p of processes in MPI_COMM_WORLD where p is specified via, e.g., 'mpirun -np p ...'
    unsigned int _num_processes;

    /// process group consisting of all processes
    ProcessGroup* _world_group;

    /// process group consisting of all processes but the master process (needed to cleanly finish the program)
    ProcessGroup* _world_group_without_master;

    /**
    * \brief process group responsible for one problem
    *
    * When the user wants to solve independent problems simultaneously (for example in some multiphysics application),
    * he can create process groups. The distribution of processes to such process groups cannot be changed in the course
    * of the program, i.e. the number of process groups and the number of processes per group is fixed. Hence, there
    * also will be no automatic load balancing across such process groups, only within each of the groups.
    */
    ProcessGroup* _process_group;

    /// number of process groups requested by the user via some top-level configuration file
    unsigned int _num_process_groups;

    /**
    * \brief array of number of processes in each process group (including eventual dedicated load balancer process)
    *
    * This array must be provided when more than one process group is used.\n
    * Dimension: [#_num_process_groups]
    */
    unsigned int* _num_processes_in_group;

    /**
    * \brief array of flags whether a dedicated load balancer process is needed in process group
    *
    * This array must be provided when more than one process group is used.\n
    * Dimension: [#_num_process_groups]
    */
    bool* _includes_dedicated_load_bal;

    /**
    * \brief 2-dim. array of MPI_COMM_WORLD ranks in top-level process group that each process unambigously belongs to
    *
    * Dimension: [#_num_process_groups][#_num_processes_in_group[group_id]]
    */
    int** _group_ranks_world;

    /**
    * \brief master process responsible for screen output and initial file IO
    *
    * Will be nullptr if not living on this process. The MPI_COMM_WORLD rank of the master process is the
    * last available one, i.e. _num_processes-1.
    */
    Master* _master;

    /**
    * \brief load balancer process
    *
    * Will be nullptr if not living on this process.
    */
    LoadBalancer<space_dim_, world_dim_>* _load_balancer;


    /**
    * \brief variable storing the base name of the log files
    *
    * If not set by the user, the file base name is set to 'feast'. The variable is passed to class Logger.
    */
    std::string _logfile_base_name;

    /* *************************
    * constructor & destructor *
    ***************************/
    /**
    * \brief CTOR
    *
    * This constructor is deliberately chosen to be private. Is called exactly once by each process in the whole life
    * time of the program.  It initialises the _universe_created variable and ensures that MPI is already initialised.
    */
    Universe()
      : _universe_created(false),
        _num_processes(0),
        _world_group(nullptr),
        _world_group_without_master(nullptr),
        _process_group(nullptr),
        _num_process_groups(0),
        _num_processes_in_group(nullptr),
        _includes_dedicated_load_bal(nullptr),
        _group_ranks_world(nullptr),
        _master(nullptr),
        _load_balancer(nullptr),
        _logfile_base_name("")
    {
      CONTEXT("Universe::Universe()");
      int mpi_is_initialised;
      MPI_Initialized(&mpi_is_initialised);
      if(!mpi_is_initialised)
      {
        throw InternalError("MPI is not initialised yet! Call init_mpi(argc,argv) first!");
      }
    }


    /**
    * \brief DTOR which automatically finalizes the MPI environment
    *
    * Just like the constructor, the destructor is deliberately chosen to be private.
    */
    ~Universe()
    {
      CONTEXT("Universe::~Universe()");
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief manages initial process distribution
    *
    * This function manages the initial process distribution into master and process groups.
    * Let there be
    * \li \c n processes with \c MPI_COMM_WORLD ranks <code>0, .., n-1</code>,
    * \li \c k process groups <code>0, .., k-1</code>, process group \c i comprising <code>0 < P_i < n-1</code>
    *     processes including the (eventually dedicated) load load balancer (<code>P_0 + ... + P_(k-1) = n-1)</code>
    *
    * Then the MPI_COMM_WORLD ranks are distributed as follows:
    * \li process group 0: ranks <code>0, ..., P_0-1</code> (where <code>P_0-1</code> is the coordinator
    *     (and eventually dedicated load balancer) process)
    * \li process group 1: ranks <code>P_0, ..., P_0+P_1-1</code>
    * \li process group 2: ranks <code>P_0 + P_1, ..., P_0 + P_1 + P_2 - 1</code>
    * \li ...
    * \li process group k-1: ranks <code>P_0 + ... + P_(k-2), ..., P_0 + ... + P_(k-1) - 1</code>
    * \li master: rank <code>P_0 + ... + P_(k-1) = n-1</code>
    */
    void _init()
    {
      CONTEXT("Universe::_init()");
      // get MPI_COMM_WORLD rank of this process
      int my_rank;
      int mpi_error_code = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      validate_error_code_mpi(mpi_error_code, "MPI_Comm_rank");

      // set MPI_COMM_WORLD rank of the master process
      int rank_master = _num_processes-1;

      // set static Process class variables for this MPI_COMM_WORLD process
      Process::rank = my_rank;
      Process::rank_master = rank_master;
      Process::num_processes = _num_processes;
      Process::is_master = (Process::rank == Process::rank_master);

      // open log files
      Logger::open_log_file(_logfile_base_name);

      // calculate the number of processes needed; and check whether there are enough MPI processes available:
      unsigned int num_processes_needed;

      // (1) the master needs one process
      num_processes_needed = 1;

      for(unsigned int igroup(0) ; igroup < _num_process_groups ; ++igroup)
      {
        // (2) add number of processes in each of the process groups (eventually including a dedicated load
        // balancer process)
        num_processes_needed += _num_processes_in_group[igroup];
      }

      // now check whether the number of available processes is sufficient
      try
      {
        if(_num_processes < num_processes_needed)
        {
          throw InternalError("Only " + stringify(_num_processes) + " processes available, but "
                              + stringify(num_processes_needed) + " processes needed!");
        }
        else if(_num_processes > num_processes_needed)
        {
          throw InternalError(stringify(_num_processes) + " processes available, and only "
                              + stringify(num_processes_needed) + " processes needed! FEAST does not "
                              + "support trailing orphaned processes.");
        }
        else
        {
          // all ok, let only one process comment on this
          if(Process::is_master)
          {
            std::string s = stringify(_num_processes) + " processes available and " + stringify(num_processes_needed)
                            + " needed.\n";
            Logger::log(s);
            std::cout << s;
          }
        }
      }
      catch (Exception& e)
      {
        // abort the program
        ErrorHandler::exception_occured(e);
      }

      // create ProcessGroup object representing the group of all COMM_WORLD processes
      _world_group = new ProcessGroup(MPI_COMM_WORLD, _num_processes);

      // create ProcessGroup object representing the group of all COMM_WORLD processes excluding the master process
      // Note that *all* processes of the parent MPI group have to call the MPI_Comm_create() routine (otherwise the
      // forking will deadlock), so let the master call it as well. Since the master does not belong to the group, the
      // new communicator is MPI_COMM_NULL on the master process.
      MPI_Group gr_without_master;
      MPI_Comm gr_comm;
      mpi_error_code = MPI_Group_excl(_world_group->group(), 1, &rank_master, &gr_without_master);
      validate_error_code_mpi(mpi_error_code, "MPI_Group_excl");
      mpi_error_code = MPI_Comm_create(_world_group->comm(), gr_without_master, &gr_comm);
      validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");
      if(!Process::is_master)
      {
        _world_group_without_master = new ProcessGroup(gr_comm, _num_processes-1);
      }

      // iterator for MPI_COMM_WORLD ranks, used to split them among the groups
      int iter_MPC_rank(-1);

      // partition global ranks into groups, including master and eventually dedicated load balancers, by simply
      // enumerating the global ranks and assigning them consecutively to the requested number of processes
      _group_ranks_world = new int*[_num_process_groups];

      // the group this process belongs to (all processes except the master will belong to a group)
      unsigned int my_group(0);

      // loop over all groups
      for(unsigned int igroup(0) ; igroup < _num_process_groups ; ++igroup)
      {
        // fill the rank array for this group
        _group_ranks_world[igroup] = new int[_num_processes_in_group[igroup]];
        for(unsigned int j(0) ; j < _num_processes_in_group[igroup] ; ++j)
        {
          ++iter_MPC_rank;
          _group_ranks_world[igroup][j] = iter_MPC_rank;
          // inquire whether this process belongs to the current group
          if(Process::rank == _group_ranks_world[igroup][j])
          {
            my_group = igroup;
          }
        }
      }
      // final sanity check (rank assigned last must be master rank minus 1)
      ASSERT(iter_MPC_rank == (int)_num_processes-2, "iter_MPC_rank " + stringify(iter_MPC_rank)
             + " must be equal to " + stringify((int)_num_processes-2) + ".");

      // Create ProcessGroup object. The constructor automatically calls the corresponding MPI routines for creating
      // MPI group and MPI communicator. Exclude the master because it is only a member of COMM_WORLD and not
      // of any group we set up.
      if(!Process::is_master)
      {
        _process_group = new ProcessGroup(_num_processes_in_group[my_group], _group_ranks_world[my_group],
                                          _world_group, my_group);
      }
      else
      {
        // *All* processes of the parent MPI group have to call the MPI_Comm_create() routine (otherwise the forking
        // will deadlock), so let the master call it with dummy communicator and dummy group. (The dummy group is
        // necessary here since the other MPI_Group object is hidden inside the ProcessGroup constructor above.)
        MPI_Comm dummy_comm;
        MPI_Group dummy_group;
        int mpi_error_code = MPI_Group_incl(_world_group->group(), 1, &Process::rank_master, &dummy_group);
        validate_error_code_mpi(mpi_error_code, "MPI_Group_incl");
        mpi_error_code = MPI_Comm_create(_world_group->comm(), dummy_group, &dummy_comm);
        validate_error_code_mpi(mpi_error_code, "MPI_Comm_create");
        // COMMENT_HILMAR: First, I used this simpler version:
        //   int mpi_error_code = MPI_Comm_create(_world_group->comm(), MPI_GROUP_EMPTY, &dummy_comm);
        // It worked with OpenMPI 1.4.2 and MPICH2, but does not with OpenMPI 1.4.3. We are not quite sure yet, if that
        // is a bug in OpenMPI 1.4.3, or if this use of MPI_GROUP_EMPTY is incorrect.
        MPI_Comm_free(&dummy_comm);
        MPI_Group_free(&dummy_group);
      }

      // all ok, now decide if this process is a regular one or the group's load-balancer or even the master

      // initialise the pointers as nullptr
      _load_balancer = nullptr;
      _master = nullptr;

      if(Process::is_master)
      {
        // for the last rank (_num_processes-1) create master process object (responsible for screen output)
        _master = new Master();

        // debug output
        std::string s = "Process " + stringify(Process::rank) + " is the MASTER OF THE UNIVERSE!\n";
        Logger::log(s);
        std::cout << s;

        // start the infinite service loop on the master, which waits for messages
        _master->service();
      }
      else
      {
        // create load balancer object in each process of the process group
        _load_balancer = new LoadBalancer<space_dim_, world_dim_>(_process_group,
                                                                  _includes_dedicated_load_bal[my_group]);
        // debug output
        Logger::log("Process " + stringify(Process::rank) + " belongs to process group " + stringify(my_group) + ".\n");
      }
    } // _init()


    /**
    * \brief cleanup counterpart for _init()/create()
    *
    * Deletes all dynamically allocated memory in _init() and create().
    */
    void _cleanup()
    {
      CONTEXT("Universe::_cleanup()");
      // open log files
      Logger::close_log_file();

      delete _world_group;
      for(unsigned int igroup(0) ; igroup < _num_process_groups ; ++igroup)
      {
        delete [] _group_ranks_world[igroup];
      }
      delete [] _group_ranks_world;
      if(Process::is_master)
      {
        delete _master;
      }
      else
      {
        delete _process_group;
        delete _load_balancer;
        delete _world_group_without_master;
      }
      delete [] _includes_dedicated_load_bal;
      delete [] _num_processes_in_group;
    }


  public:

    // set parent class as friend such that it can execute CTOR and DTOR of this class
    friend class InstantiationPolicy<Universe<space_dim_, world_dim_>, Singleton>;

    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the load balancer object
    *
    * \return pointer to LoadBalancer #_load_balancer
    */
    inline LoadBalancer<space_dim_, world_dim_>* load_balancer() const
    {
      CONTEXT("Universe::load_balancer()");
      ASSERT(_universe_created, "Universe must be created first by calling create(...).");
      return _load_balancer;
    }


    /**
    * \brief getter for the master process object
    *
    * \return pointer to Master #_master
    */
    inline Master* master() const
    {
      CONTEXT("Universe::master()");
      ASSERT(_universe_created, "Universe must be created first by calling create(...).");
      return _master;
    }


    /* ****************
    * member functions*
    ******************/
    /**
    * \brief function for creating a universe with exactly one process group
    *
    * When using this function, only the master process and one process group using no dedicated load balancer process
    * are created.
    *
    * \param[in] logfile_base_name
    * optional base name of log files (name pattern '<logfile_base_name><world rank of the process>.log'); if not
    * provided, it is set to 'feast'
    */
    void create(std::string const logfile_base_name = "feast")
    {
      CONTEXT("Universe::create()");
      if(!_universe_created)
      {
        ASSERT(logfile_base_name.size() > 0, "Logfile base name must not be empty.");
        _logfile_base_name = logfile_base_name;
        _universe_created = true;
        // get total number of processes
        int num_proc;
        int mpi_error_code = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        validate_error_code_mpi(mpi_error_code, "MPI_Comm_size");
        ASSERT(num_proc >= 1, "Number of processes must be at least 1.");
        _num_processes = (unsigned int) num_proc;

        // Set variables. Assume that one process group is needed using all available processes (minus one
        // master process) and no dedicated load balancer process.
        _num_process_groups = 1;
        _num_processes_in_group = new unsigned int[1];
        _num_processes_in_group[0] = _num_processes - 1;
        _includes_dedicated_load_bal = new bool[1];
        _includes_dedicated_load_bal[0] = false;

        // call the init routine
        _init();
      }
      else
      {
        throw InternalError("Universe has been created already!");
      }
    }


    /**
    * \brief function for creating a universe with several process groups
    *
    * When using this function, one can choose how many process groups are created. Additionally, one load balancer
    * per process group and one master process is created. The caller has to ensure that sufficient MPI processes are
    * available.
    *
    * \param[in] num_process_groups
    * number of process groups
    *
    * \param[in] num_processes_in_group
    * array of numbers of processes in each process group, dimension [\a num_process_groups]
    *
    * \param[in] includes_dedicated_load_bal
    * array of flags whether a dedicated load balancer process is needed in process group,
    * dimension [\a num_process_groups]
    *
    * \param[in] logfile_base_name
    * optional base name of log files (name pattern '<logfile_base_name><world rank of the process>.log'); if not
    * provided, it is set to 'feast'
    *
    * \todo When the number of parameters increases even more, than one should think about creating some struct which
    * contains all this data (something like "struct UniverseData"), so that this struct can be extended in future
    * without changing the interface of the routine. Then one could also merge the two versions of the create routine.
    */
    void create(
      unsigned int num_process_groups,
      unsigned int num_processes_in_group[],
      bool includes_dedicated_load_bal[],
      std::string const logfile_base_name = "feast")
    {
      CONTEXT("Universe::create()");
      if(!_universe_created)
      {
        ASSERT(logfile_base_name.size() > 0, "Logfile base name must not be empty.");
        _logfile_base_name = logfile_base_name;
        _universe_created = true;
        // get total number of processes
        int num_proc;
        int mpi_error_code = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        validate_error_code_mpi(mpi_error_code, "MPI_Comm_size");
        ASSERT(num_proc >= 1, "Number of processes must be at least 1.");
        _num_processes = (unsigned int) num_proc;
        // set variables
        _num_process_groups = num_process_groups;
        ASSERT(_num_process_groups >= 1, "Number of process groups must be at least 1.");
        _num_processes_in_group = num_processes_in_group;
        for(unsigned int igroup(0) ; igroup < _num_process_groups ; ++igroup)
        {
          ASSERT(_num_processes_in_group[igroup] >= 1, "Number of processes in process group " + stringify(igroup)
                 + " must be at least 1.");
        }
        _includes_dedicated_load_bal = includes_dedicated_load_bal;
        // call the init routine
        _init();
        // debug output
        Logger::log("Universe created on process " + stringify(Process::rank) + ".\n");
      }
      else
      {
        throw InternalError("Universe has been created already!");
      }
    }


    /**
    * \brief cleanup of the universe
    *
    * This function has to be called by all COMM_WORLD processes to cleanly end the program.
    */
    void destroy()
    {
      CONTEXT("Universe::destroy()");
      if(_universe_created)
      {
        if(!Process::is_master)
        {
          // wait for all non-master processes (the master is still in its infinite service loop)
          MPI_Barrier(_world_group_without_master->comm());

          // the coordinator of the process group _world_group_without_master tells the master to stop its service loop
          if(_world_group_without_master->is_coordinator())
          {
            // the coordinator inits a new message with corresponding ID
            Comm::init(ServiceIDs::MASTER_FINISH_SERVICE);
            // send message
            Comm::send();
            // now the master ends its service loops and automatically calls Universe::destroy()
          }
        }
        // clean up dynamically allocated memory
        _cleanup();
      }
      // debug output
      Logger::log("Universe destroyed on process " + stringify(Process::rank) + ".\n");
    }
  };
} // namespace FEAST

#endif // guard KERNEL_UNIVERSE_HPP
