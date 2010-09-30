#pragma once
#ifndef KERNEL_UNIVERSE_HPP
#define KERNEL_UNIVERSE_HPP 1

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cassert>

// includes, Feast
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/base_header.hpp>
#include <kernel/process.hpp>
#include <kernel/master.hpp>
#include <kernel/load_balancer.hpp>

// TODO: All classes in this file/compilation unit should be namespace Feast.
using namespace Feast;

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
* \author Hilmar Wobker
* \author Dominik Goeddeke
*/
class Universe
{

private:

  /* *****************
  * member variables *
  *******************/
  /**
  * \brief pointer to the one and only instance of the universe
  *
  * Will be set within the routine Universe::create().
  */
  static Universe* _universe;

  /// total number of processes in MPI_COMM_WORLD
  int _num_processes;

  /// process group consisting of all processes
  ProcessGroup* _world_group;

  /// process group responsible for one problem
  ProcessGroup* _process_group;

  /// number of process groups requested by the user via some top-level configuration file
  const int _num_process_groups;

  /**
  * \brief array of number of processes in each process group (including eventual dedicated load balancer process)
  *
  * This array must be provided when more than one process group is used.
  * <em>Dimension:</em> [#_num_process_groups]
  */
  int const * _num_processes_in_group;

  /**
  * \brief array of flags whether a dedicated load balancer process is needed in process group
  *
  * This array must be provided when more than one process group is used.
  * <em>Dimension:</em> [#_num_process_groups]
  */
  bool const * _includes_dedicated_load_bal;

  /// number of processes actually needed by the user, based on some top-level configuration file
  int _num_processes_needed;

  /**
  * \brief 2-dim. array of MPI_COMM_WORLD ranks in top-level process group that each process unambigously belongs to
  *
  * <em>Dimension:</em> [#_num_process_groups][#_num_processes_in_group[group_id]]
  */
  int ** _group_ranks_world;

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
  * Will be nullptr if not living on this process
  */
  LoadBalancer* _load_balancer;


  /* *************************
  * constructor & destructor *
  ***************************/
  /**
  * \brief constructor for creating the one and only instance of the Universe class
  *
  * The constructor is deliberately chosen to be private: That way one can prevent users to create
  * more than one instance of the Universe class ('singleton').
  *
  * \param[in] num_processes
  * total number of MPI processes
  *
  * \param[in] num_process_groups
  * number of process groups
  *
  * \param[in] num_processes_in_group
  * array of numbers of processes in work groups
  * (\em dimension [Universe::Universe#num_process_groups])
  *
  * \param[in] includes_dedicated_load_bal
  * array of flags whether dedicated load balancer required in work groups
  * (\em dimension [Universe::Universe#num_process_groups])
  */
  Universe(
    const int num_processes,
    const int num_process_groups,
    const int num_processes_in_group[],
    const bool includes_dedicated_load_bal[])
    : _num_processes(num_processes),
      _num_process_groups(num_process_groups),
      _num_processes_in_group(num_processes_in_group),
      _includes_dedicated_load_bal(includes_dedicated_load_bal)
  {
    _init();
//    // debug output
//    std::cout << "Process " << Process::rank << " created its part of the Universe!" << std::endl;
  }


  /**
  * \brief destructor which automatically finalizes the MPI environment
  *
  * Just as the constructor, the destructor is deliberately chosen to be private.
  */
  ~Universe()
  {
    // clean up dynamically allocated memory
    _cleanup();
//    // debug output
//    std::cout << "Process " << Process::rank << " now destroys its part of the Universe!" << std::endl;
  }

  /// copy constructor, intentionally undefined to prevent object instantiation
  Universe(const Universe &);

  /// assignment constructor, intentionally undefined to prevent object instantiation
  Universe & operator=(const Universe &);


  /* *****************
  * member functions *
  *******************/
  /**
  * \brief initialises MPI, returns the total number of processes
  *
  * \param[in] argc
  * argument count passed to the main() method
  *
  * \param[in] argv
  * arguments passed to the main() method
  *
  * \param[out] num_processes
  * number p of available MPI processes as specified via 'mpirun -np p ...'
  */
  static void _init_mpi(
    int& argc,
    char* argv[],
    int& num_processes)
  {
    // init MPI
    int mpi_error_code = MPI_Init(&argc, &argv);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Init");

    // get total number of processes
    mpi_error_code = MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_size");
  }


  /**
  * \brief manages initial process distribution
  *
  * This function manages the initial process distribution into master, load balancers and process groups.
  * Let there be
  * \li \c n processes with \c MPI_COMM_WORLD ranks <code>0, .., n-1</code>,
  * \li \c k process groups <code>0, .., k-1</code>, process group \c i comprising <code>0 < P_i < n-1</code>
  *     processes and one load balancer.
  *
  * Then the MPI_COMM_WORLD ranks are distributed as follows:
  * \li process group 0: ranks <code>0, ..., P_0-1</code> plus rank <code>P_0</code> for load balancer 0
  * \li process group 1: ranks <code>P_0+1, ..., P_0+P_1+1-1</code> plus rank <code>P_0 + P_1 + 1</code> for
  *     load balancer 1
  * \li process group 2: ranks <code>P_0 + P_1 + 2, ..., P_0 + P_1 + P_2 + 2 - 1</code> plus rank
  *     <code>P_0 + P_1 + P_2 + 2</code> for load balancer 2
  * \li ...
  * \li process group k-1: ranks <code>P_0 + ... + P_(k-2) + k-1, .., P_0 + ... + P_(k-1) + k-1 - 1</code> plus rank
  *     rank <code>P_0 + ... + P_(k-1) + k-1</code> for load balancer k-1
  * \li master: rank <code>P_0 + ... + P_(k-1) + k = n_needed-1</code> where <code>n_needed <= n</code>
  */
  void _init()
  {
    // get MPI_COMM_WORLD rank of this process
    int my_rank;
    int mpi_error_code = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_rank");

    // set MPI_COMM_WORLD rank of the master process
    int rank_master = _num_processes-1;

    // set static Process class variables for this MPI_COMM_WORLD process
    Process::rank = my_rank;
    Process::rank_master = rank_master;
    Process::is_master = (Process::rank == Process::rank_master);

    // calculate the number of processes needed; and check whether there are enough MPI processes available:

    // (1) the master needs one process
    _num_processes_needed = 1;

    for(int igroup(0) ; igroup < _num_process_groups ; ++igroup)
    {
      // (2) add number of processes in each of the process groups (eventually including a dedicated load
      // balancer process)
      _num_processes_needed += _num_processes_in_group[igroup];
    }

    // now check whether the number of available processes is sufficient
    if(_num_processes < _num_processes_needed)
    {
      MPIUtils::abort("Error! Only " + StringUtils::stringify(_num_processes) + " processes available, but "
        + StringUtils::stringify(_num_processes_needed) + " processes needed!");
    }
    else if(_num_processes > _num_processes_needed)
    {
      MPIUtils::abort("Error!  " + StringUtils::stringify(_num_processes) + " processes available, and only "
        + StringUtils::stringify(_num_processes_needed) + " processes needed!. Since Feast does not "
        + "support trailing orphaned processes, this is an error.");
    }
    else
    {
      // all ok, let only one process comment on this
      if(Process::is_master)
      {
        std::cout << _num_processes << " processes available and " << _num_processes_needed
                  << " needed." << std::endl;
      }
    }

    // create ProcessGroup object representing the group of all COMM_WORLD processes
    _world_group = new ProcessGroup(MPI_COMM_WORLD, _num_processes);

    // iterator for MPI_COMM_WORLD ranks, used to split them among the groups
    int iter_MPC_rank(-1);

    // partition global ranks into groups, including master and eventually dedicated load balancers, by simply
    // enumerating the global ranks and assigning them consecutively to the requested number of processes
    _group_ranks_world = new int*[_num_process_groups];

    // the group this process belongs to (all processes except the master will belong to a group)
    int my_group(-1);

    // loop over all groups
    for(int igroup(0) ; igroup < _num_process_groups ; ++igroup)
    {
      // fill the rank array for this group
      _group_ranks_world[igroup] = new int[_num_processes_in_group[igroup]];
      for(int j(0) ; j < _num_processes_in_group[igroup] ; ++j)
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
    assert(iter_MPC_rank == _num_processes-2);

    // create ProcessGroup object. The constructor automatically calls the corresponding MPI routines for creating
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
      // will deadlock), so let the master call it with special MPI_GROUP_EMPTY and dummy communicator.
      MPI_Comm dummy_comm;
      int mpi_error_code = MPI_Comm_create(_world_group->comm(), MPI_GROUP_EMPTY, &dummy_comm);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_create");
    }

    // all ok, now decide if this process is a regular one or the group's load-balancer or even the master

    // initialise the pointers as nullptr
    _load_balancer = nullptr;
    _master = nullptr;

    if(Process::is_master)
    {
      // for the last rank (_num_processes-1) create master process object (responsible for screen output)
      _master = new Master();

      std::cout << "Process " << Process::rank << " is the MASTER OF THE UNIVERSE!" << std::endl;

      // start some dummy wait routine on the master
//      _master->wait();

      // start the infinite service loop on the master, which waits for messages
      _master->service();
    }
    else
    {
      // create load balancer object in each process of the process group
      _load_balancer = new LoadBalancer(_process_group, _includes_dedicated_load_bal[my_group]);
    }
  } // _init()


  /**
  * \brief cleanup counterpart for _init().
  *
  * Deletes all dynamically allocated memory in _init().
  */
  void _cleanup()
  {
    delete _world_group;
    for(int igroup(0) ; igroup < _num_process_groups ; ++igroup)
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
    }
  }


public:

  /* ****************
  * member functions*
  ******************/
  /**
  * \brief function for creating a universe with exactly one process group
  *
  * When using this function, only the master process and one process group using an dedicated load balancer process
  * are created.
  *
  * \param[in] argc
  * argument count passed to the main() method
  *
  * \param[in] argv
  * arguments passed to the main() method
  *
  * \return Universe pointer #_universe
  */
  static Universe* create(
    int argc,
    char* argv[])
  {
    if(_universe == nullptr)
    {
      int num_processes;
      _init_mpi(argc, argv, num_processes);

      // Assume that one process group is needed using all available processes (minus one master process).
      // The dedicated load balancer process is included in this number.
      int num_process_groups = 1;
      int num_processes_in_group[1] = {num_processes-1};
      bool includes_dedicated_load_bal[1] = {true};

      // create the one and only instance of the Universe class
      _universe = new Universe(num_processes, num_process_groups, num_processes_in_group,
                               includes_dedicated_load_bal);
      return _universe;
    }
    else
    {
      MPIUtils::abort("Universe can be created only once!");
      return nullptr;
    }
  }


  /**
  * \brief function for creating a universe with several process groups
  *
  * When using this function, one can choose how many process groups are created. Additionally, one load balancer
  * per process group and one master process is created. The caller has to ensure that sufficient MPI processes are
  * available.
  *
  * \param[in] argc
  * argument count passed to the main() method
  *
  * \param[in] argv
  * arguments passed to the main() method
  *
  * \param[in] num_process_groups
  * number of process groups
  *
  * \param[in] num_processes_in_group
  * array of numbers of processes in each process group (\em dimension [Universe::create#num_process_groups])
  *
  * \param[in] includes_dedicated_load_bal
  * array of flags whether a dedicated load balancer process is needed in process group
  * (\em dimension [Universe::create#num_process_groups])
  *
  * \return Universe pointer #_universe
  */
  static Universe* create(
    int argc,
    char* argv[],
    const int num_process_groups,
    const int num_processes_in_group[],
    const bool includes_dedicated_load_bal[])
  {
    if(_universe == nullptr)
    {
      int num_processes;
      // init MPI and get number of available processes and the rank of this processor
      _init_mpi(argc, argv, num_processes);
      // create the one and only instance of the Universe class
      _universe = new Universe(num_processes, num_process_groups, num_processes_in_group,
                               includes_dedicated_load_bal);
      return _universe;
    }
    else
    {
      MPIUtils::abort("Universe can be created only once!");
      return nullptr;
    }
  }


  /// function for destroying the universe and finalizing MPI at the end of the program
  static void destroy()
  {
    if(_universe != nullptr)
    {
      // delete the one and only instance of Universe
      delete _universe;
      _universe = nullptr;
      // shut down MPI
      int mpi_is_initialised;
      MPI_Initialized(&mpi_is_initialised);
      if(mpi_is_initialised)
      {
        MPI_Finalize();
      }
    }
    else
    {
      MPIUtils::abort("There is no universe to destroy!");
    }
  }


  /* ******************
  * getters & setters *
  ********************/
  /**
  * \brief getter for the load balancer object
  *
  * \return pointer to LoadBalancer #_load_balancer
  */
  inline LoadBalancer* load_balancer() const
  {
    return _load_balancer;
  }

  /**
  * \brief getter for the master process object
  *
  * \return pointer to Master #_master
  */
  inline Master* master() const
  {
    return _master;
  }
};
// initialise static member
Universe* Universe::_universe = nullptr;

#endif // guard KERNEL_UNIVERSE_HPP
