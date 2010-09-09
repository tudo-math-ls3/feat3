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
 * @author Hilmar Wobker
 * @author Dominik Goeddeke
 */
class Universe
{
  /* *****************
   * private members *
   *******************/
  private:

    /* ******************
     * member variables *
     ********************/

    /**
     * \brief pointer to the one and only instance of the universe
     *
     * Will be set within the routine Universe::create().
     */
    static Universe* _universe;

    /**
     * \brief total number of processes in MPI_COMM_WORLD
     */
    int _num_processes;

    /**
     * \brief MPI_COMM_WORLD rank of this process
     */
    int _my_rank;

    /**
     * \brief number of process groups requested by the user via some top-level configuration file
     */
    int _num_process_groups;

    /**
     * \brief array of number of processes in each process groups (excl. load balancer process)
     *
     * This array must be provided when more than one process group is used.
     * <em>Dimension:</em> [#_num_process_groups]
     */
    int* _num_processes_in_group;

    /**
     * \brief number of processes actually needed by the user, based on some top-level configuration file
     */
    int _num_processes_needed;

    /**
     * \brief array of ranks in top-level process group that each process unambigously belongs to, per top-level group
     *
     * <em>Dimension:</em> [#_num_process_groups][#_num_processes_in_group[group_id]+1]
     */
    int** _group_ranks;

    /**
     * \brief master process responsible for screen output and initial file IO
     *
     * Will be nullptr if not living on this process.
     */
    Master* _master;

    /**
     * \brief load balancer process
     *
     * Will be nullptr if not living on this process which implicitly includes top-level group membership
     */
    LoadBalancer* _load_balancer;

    /**
     * \brief group process
     *
     * Will be nullptr if not living on this process which implicitly includes top-level group membership
     */
    GroupProcess* _group_process;

    /**
     * \brief handle to group associated with MPI_COMM_WORLD communicator.
     *
     * Needs to be stored explicitly.
     */
    MPI_Group _mpi_world_group;


    /* *************
     * constructor *
     ***************/
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
     * array of numbers of processes in each work group (dimension [#num_process_groups])
     */
    Universe(
      int num_processes,
      int my_rank,
      int num_process_groups,
      int num_processes_in_group[])
      : _num_processes(num_processes),
        _my_rank(my_rank),
        _num_process_groups(num_process_groups),
        _num_processes_in_group(num_processes_in_group)
    {
      _init();
      std::cout << "Process " << _my_rank << " created its part of the Universe!" << std::endl;
    }


    /* ************
     * destructor *
     **************/
    /**
     * \brief destructor which automatically finalizes the MPI environment
     *
     * Just as the constructor, the destructor is deliberately chosen to be private.
     */
    ~Universe()
    {
      // clean up dynamically allocated memory
      _cleanup();
      std::cout << "Process " << _my_rank << " now destroys its part of the Universe!" << std::endl;
    }


    /* ******************
     * member functions *
     ********************/
    /**
     * \brief initialises MPI, returns the total number of processes and the rank of this processor
     *
     * \param[in] argc
     * argument count passed to the main() method
     *
     * \param[in] argv
     * arguments passed to the main() method
     *
     * \param[out] num_processes
     * number p of available MPI processes as specified via 'mpirun -np p ...'
     *
     * \param[out] my_rank
     * MPI_COMM_WORLD rank of this MPI process
     */
    static void _init_mpi(
      int& argc,
      char* argv[],
      int& num_processes,
      int& my_rank)
    {
      // init MPI
      int mpi_error_code = MPI_Init(&argc, &argv);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Init");

      // get total number of processes
      mpi_error_code = MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_size");

      // rank of this process
      mpi_error_code = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_rank");
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

      // top-level (physics-level) process group that this process belongs to
      MPI_Group process_group;
      // communicator corresponding to this process'es (TODO GRAMMAR) group, aka process group
      MPI_Comm group_comm;

      // calculate the number of processes needed; and check whether there are enough MPI processes available:
      // (1) the master needs one process
      _num_processes_needed = 1;

      // (2) add number of processes in each of the process groups
      for (int igroup(0) ; igroup < _num_process_groups ; ++igroup)
      {
        _num_processes_needed += _num_processes_in_group[igroup];
      }
      // (3) add one load balancer process per process group
      _num_processes_needed += _num_process_groups;

      // now check whether the number of available processes is sufficient
      if (_num_processes < _num_processes_needed)
      {
        MPIUtils::abort("Error! Only " + StringUtils::stringify(_num_processes) + " processes available, but "
          + StringUtils::stringify(_num_processes_needed) + " processes needed!");
      }
      else if (_num_processes > _num_processes_needed)
      {
        MPIUtils::abort("Error!  " + StringUtils::stringify(_num_processes) + " processes available, and only "
          + StringUtils::stringify(_num_processes_needed) + " processes needed!. Since Feast does not "
          + "support trailing orphaned processes, this is an error.");
      }
      else
      {
        // all ok, let only one process comment on this
        if (_my_rank == 0)
        {
          std::cout << _num_processes << " processes available and " << _num_processes_needed
                    << " needed." << std::endl;
        }
      }

      // iterator for MPI_COMM_WORLD ranks, used to split the global ranks among the groups
      int iter_MPC_rank(0);

      // extract world group handle and store it for later use,
      // in particular for forking off new groups from the global group/comm
      int mpi_error_code = MPI_Comm_group(MPI_COMM_WORLD, &_mpi_world_group);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_group");

      // partition global ranks into groups, including master and loadbalancer(s), by simply enumerating
      // the global ranks and assigning them consecutively to the requested number of processes
      _group_ranks = new int*[_num_process_groups];
      int my_group;
      // loop over all groups
      for(int igroup(0) ; igroup < _num_process_groups ; ++igroup)
      {
        _group_ranks[igroup] = new int[_num_processes_in_group[igroup]+1];
        // fill the rank array for this group
        for(int j(0) ; j < _num_processes_in_group[igroup]+1 ; ++j)
        {
          _group_ranks[igroup][j] = iter_MPC_rank;
          // std::cout << _my_rank << ": igroup=" << igroup << ", j=" << j << ", ranks[j]=" << _group_ranks[j] << std::endl;
          ++iter_MPC_rank;
          // inquire whether this process belongs to the current group
          if (_my_rank == _group_ranks[igroup][j])
          {
            my_group = igroup;
            // std::cout << _my_rank << " belongs to group " << igroup << std::endl;
          }
        }
      }
      // final sanity check
      assert(iter_MPC_rank == _num_processes-1);

      // set MPI_COMM_WORLD rank of the master process
      // TODO include in doxygen: Master is the very last available process
      int rank_master = _num_processes-1;

      // create MPI groups (exclude the master because it is only a member of COMM_WORLD and not
      // of any group we set up, so a special group is needed since otherwise the forking call below
      // will deadlock)
      if (_my_rank != rank_master)
      {
        mpi_error_code = MPI_Group_incl(_mpi_world_group, _num_processes_in_group[my_group]+1, _group_ranks[my_group],
                                        &process_group);
        MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_incl");
      }
      else
      {
        // this is a special valid handle for an empty group that can be passed to operations like
        // MPI_Comm_create that expect *all* processes in the parent communicator to participate.
        process_group = MPI_GROUP_EMPTY;
      }

      // Create the group communicator for, among others, collective operations.
      // It is essential that *all* processes in MPI_COMM_WORLD participate since MPI_COMM_WORLD
      // is a parent sub-universe, i.e., something to spawn from
      mpi_error_code = MPI_Comm_create(MPI_COMM_WORLD, process_group, &group_comm);
      MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Comm_create");

      // rank of this process within the process group it will be scheduled to
      int rank_local;

      if (_my_rank != rank_master)
      {
        // and finally look up the local rank of this process within the group
        mpi_error_code = MPI_Group_rank(process_group, &rank_local);
        MPIUtils::validate_mpi_error_code(mpi_error_code, "MPI_Group_rank");
      }

      // all ok, now decide if this process is a regular one or the group's load-balancer or even the master

      // initialise the pointers as nullptr
      _group_process = nullptr;
      _load_balancer = nullptr;
      _master = nullptr;

      if (_my_rank == rank_master)
      {
        // for the last rank (_num_processes-1) create master process object (responsible for screen output)
        _master = new Master(_my_rank);

        // start some infinite loop in the Master object, which waits for messages
        _master->wait();
      }
      else
      {
        int rank_load_bal = _group_ranks[my_group][_num_processes_in_group[my_group]];
        if (_my_rank == rank_load_bal)
        {
          // create load balancer object for the last rank in the current group
          _load_balancer = new LoadBalancer(_my_rank, rank_master, group_comm, rank_local, my_group,
                                            _group_ranks[my_group]);
        }
        else
        {
          // create group process object for all the other ranks in the current group and
          // start infinite listener loop
          _group_process = new GroupProcess(_my_rank, rank_master, group_comm, rank_local, my_group, rank_load_bal);
          _group_process->wait();
        }
      }
    } // _init()

    /**
     * \brief cleanup counterpart for _init().
     *
     * Deletes all dynamically allocated memory in _init().
     */
    void _cleanup()
    {
      for(int igroup(0) ; igroup < _num_process_groups ; ++igroup)
      {
        delete [] _group_ranks[igroup];
      }
      delete [] _group_ranks;
    }


  /* ****************
   * public members *
   ******************/
  public:

    /* *****************
     * member functions*
     *******************/

    /**
     * \brief function for creating a universe with exactly one process group
     *
     * When using this function, only the master process and one process group with load balancer are created.
     *
     *
     * \param[in] argc
     * argument count passed to the main() method
     *
     * \param[in] argv
     * arguments passed to the main() method
     */
    static Universe* create(
      int argc,
      char* argv[])
    {
      if (_universe == nullptr)
      {
        int num_processes;
        int my_rank;
        // init MPI and get number of available processes and the rank of this processor
        _init_mpi(argc, argv, num_processes, my_rank);

        // assume that one process group is needed using all available processes (minus one master and
        // one load balancer)
        int num_process_groups = 1;
        int num_processes_in_group[1] = {num_processes-2};

        // create the one and only instance of the Universe class
        _universe = new Universe(num_processes, my_rank, num_process_groups, num_processes_in_group);
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
     * array of numbers of processes in each work group (dimension [#num_process_groups])
     */
    static Universe* create(
      int argc,
      char* argv[],
      int num_process_groups,
      int num_processes_in_group[])
    {
      if (_universe == nullptr)
      {
        int num_processes;
        int my_rank;
        // init MPI and get number of available processes and the rank of this processor
        _init_mpi(argc, argv, num_processes, my_rank);
        // create the one and only instance of the Universe class
        _universe = new Universe(num_processes, my_rank, num_process_groups, num_processes_in_group);
        return _universe;
      }
      else
      {
        MPIUtils::abort("Universe can be created only once!");
        return nullptr;
      }
    }


    /**
     * \brief function for destroying the universe and finalizing MPI at the end of the program
     */
    static void destroy()
    {
      if (_universe != nullptr)
      {
        // delete the one and only instance of Universe
        delete _universe;
        _universe = nullptr;
        // shut down MPI
        int mpi_is_initialised;
        MPI_Initialized(&mpi_is_initialised);
        if (mpi_is_initialised)
        {
          MPI_Finalize();
        }
      }
      else
      {
        MPIUtils::abort("There is no universe to destroy!");
      }
    }


    /* ***********
     * accessors *
     *************/
    /**
     * \brief accessor for the load balancer object
     */
    inline LoadBalancer* get_load_balancer() const
    {
      return _load_balancer;
    }

    /**
     * \brief accessor for the group process object
     */
    inline GroupProcess* get_group_process() const
    {
      return _group_process;
    }

    /**
     * \brief accessor for the master process object
     */
    inline Master* get_master() const
    {
      return _master;
    }

};
// initialise static member
Universe* Universe::_universe = nullptr;

#endif // guard KERNEL_UNIVERSE_HPP
