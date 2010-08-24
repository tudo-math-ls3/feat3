#pragma once
#ifndef KERNEL_MYCLASS_HPP
#define KERNEL_MYCLASS_HPP 1

#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cassert>
#include <kernel/util/string_utils.hpp>
#include <kernel/base_header.hpp>
#include <kernel/process.hpp>
#include <kernel/load_balancer.hpp>

using namespace Feast;

// In order to increase performance one can ignore the status object in some MPI functions
// (see Section 3.2.6 in the MPI 2.2 standard).
// In debug mode, we use a real status object, otherwise the special MPI status MPI_STATUS_IGNORE.
// Requires that the MPI_Status object *always* has the name 'status'.
#ifdef NDEBUG
#  define FEAST_MPI_STATUS MPI_STATUS_IGNORE
#else
#  define FEAST_MPI_STATUS &status
#endif

/**
 * \brief creates and manages the parallel universe
 *
 * At the beginning of the program exactly one object of this class has to be created. The destruction of this object
 * will also end the program. The universe object creates and manages all MPI processes.
 * @Hilmar: Waere es nicht vielleicht schlau, diese Klasse gleich aus Prinzip statisch zu machen, so dass sie nur einmal
 *          existiert und beim Programmstart automatisch schon da ist? Dann kann man von ueberall per
 *          Universe::bral() Funktionen aufrufen, der Benutzer braucht diese Klasse nicht zu instantiieren, und der
 *          automatisch ausgefuehrte Konstruktor baut dann Prozessgruppen usw. basierend auf Konfigurationsdateien.
 *          So koennte man immer noch dem Benutzer volle Kontrolle ueber "mehr als eine Physik" geben, denn dafuer muss
 *          man das Universum nicht explizit erschaffen (Ich liebe diese MPI-Sprache nebenbei...)
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
     * \brief indicates whether the universe has already been created
     *
     * This flag will be set to \c true as soon as a constructor is called. Each following constructor call aborts the
     * program.
     */
    static bool _universe_created;

    /**
     * \brief total number of processes in MPI_COMM_WORLD
     */
    int _num_processes;

    /**
     * \brief number of processes actually needed by the user
     */
    int _num_processes_needed;

    /**
     * \brief number of process groups
     */
    int _num_process_groups;

    /**
     * \brief array of number of processes in process groups
     *
     * This array must be provided when more than one process group is used.
     * <em>Dimension:</em> [#_num_process_groups]
     */
    int* _num_processes_in_group;

    /**
     * \brief master process responsible for screen output
     *
     * Will be null if not living on this process.
     */
    Master* _master;

    /**
     * \brief load balancer process
     *
     * Will be null if not living on this process.
     */
    LoadBalancer* _load_balancer;

    /**
     * \brief group process
     *
     * Will be null if not living on this process.
     */
    GroupProcess* _group_process;

    /* ******************
     * member functions *
     ********************/

    /**
     * \brief inits MPI and inquires the total number of processes
     *
     * \param[in] argc argument count passed to the main() method
     * \param[in] argv arguments passed to the main() method
     */
    void _init_mpi(
      int& argc,
      char* argv[])
    {
      // init MPI
      int mpi_error_code = MPI_Init(&argc, &argv);
      _validate_mpi_error_code(mpi_error_code, "MPI_Init");
      // get total number of processes
      mpi_error_code = MPI_Comm_size(MPI_COMM_WORLD, &_num_processes);
      _validate_mpi_error_code(mpi_error_code, "MPI_Comm_size");
    }

    /**
     * \brief evaluates the return value of MPI routines
     * @Hilmar: Ich wuerde das vielleicht nicht ins Universum packen, sondern irgendwo hin, wo man bequem von ueberall
     *          aus auf zugreifen kann. Oder kennt jedes Objekt in FEAST2 (bzw in diesem Prototypen) das Universum?
     *
     * \param[in] error_code
     * MPI error code
     *
     * \param[in] mpi_function_name
     * name of the calling routine
     */
    void _validate_mpi_error_code(
      int error_code,
      std::string mpi_function_name)
    {
      if (error_code != MPI_SUCCESS)
      {
        std::cerr << mpi_function_name << " failed with error code " << error_code <<". Aborting." << std::endl;
        int mpi_is_initialised;
        MPI_Initialized(&mpi_is_initialised);
        if (mpi_is_initialised)
        {
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        exit(1);
      }
    }

    /**
     * \brief aborts the program
     *
     * This function aborts the program and especially the MPI universe.
     *
     * \param[in] msg
     * message explaining the reason for the abortion
     */
    void _abort(std::string msg)
    {
      std::cerr << msg << " Aborting program..." << std::endl;
      int mpi_is_initialised;
      MPI_Initialized(&mpi_is_initialised);
      if (mpi_is_initialised)
      {
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      exit(1);
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

#ifndef NDEBUG
      MPI_Status status;
#endif
      // world group
      MPI_Group MPI_GROUP_WORLD;
      MPI_Group process_group;
      MPI_Comm comm_group;

      // rank of this process
      int my_rank;
      int mpi_error_code = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      _validate_mpi_error_code(mpi_error_code, "MPI_Comm_rank");

      // set MPI_COMM_WORLD rank of the master process
      int rank_master = _num_processes-1;

      // calculate the number of processes needed and check whether there are enough MPI processes available

      // the master needs one process
      _num_processes_needed = 1;

      // add number of processes in each of the process groups
      for (int i(0) ; i < _num_process_groups ; ++i)
      {
        _num_processes_needed += _num_processes_in_group[i];
      }
      // add one load balancer process per process group
      _num_processes_needed += _num_process_groups;

      // now check whether the number of available processes is sufficient
      if (_num_processes < _num_processes_needed)
      {
        _abort("Error! Only " + StringUtils::stringify(_num_processes) + " processes available, but "
               + StringUtils::stringify(_num_processes_needed) + " processes needed!");
      }
      else if (_num_processes > _num_processes_needed)
      {
        if (my_rank == 0)
        {
          std::cout << "Warning! "<< _num_processes << " processes available, but only " << _num_processes_needed
               << " processes needed!" << std::endl;
        }
      }
      else
      {
        if (my_rank == 0)
        {
          std::cout << _num_processes << " processes available and " << _num_processes_needed
               << " needed." << std::endl;
        }
      }

      // get world group
      MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

      // iterator for MPI_COMM_WORLD ranks
      int iter_MPC_rank(0);
      // rank of this process within the process group it will be scheduled to
      int rank_local;

      // initialise the pointers as nullptr
      _group_process = nullptr;
      _load_balancer = nullptr;
      _master = nullptr;

      // create MPI_COMM_WORLD rank arrays describing the process distribution to process groups and load balancers
      for(int igroup(0) ; igroup < _num_process_groups ; ++igroup)
      {
        // int array for the MPI_COMM_WORLD ranks of the process group (including the load balancer process)
        int* ranks = new int[_num_processes_in_group[igroup]+1];
        // flag whether this process belongs to the current group
        bool belongs_to_group = false;
        // fill the rank array
        for(int j(0) ; j < _num_processes_in_group[igroup]+1 ; ++j)
        {
          ranks[j] = iter_MPC_rank;
          std::cout << my_rank << ": igroup=" << igroup << ", j=" << j << ", ranks[j]=" << ranks[j] << std::endl;
          ++iter_MPC_rank;
          // inquire whether this process belongs to the current group
          if (my_rank == ranks[j])
          {
            std::cout << my_rank << " belongs to group " << igroup << std::endl;
            belongs_to_group = true;
          }
        }
        if (belongs_to_group)
        {
          // each process whose rank is in ranks now calls the corresponding MPI routines to build a sub group
          MPI_Group_incl(MPI_GROUP_WORLD, _num_processes_in_group[igroup]+1, ranks, &process_group);
          MPI_Group_rank(process_group, &rank_local);
          MPI_Comm_create(MPI_COMM_WORLD, process_group, &comm_group);

          int rank_load_bal = ranks[_num_processes_in_group[igroup]];
          if (my_rank == rank_load_bal)
          {
            // create load balancer object for the last rank in the current group
            _load_balancer = new LoadBalancer(my_rank, rank_master, comm_group, rank_local, igroup, ranks);
            // @hilmar: Debug-Ausgabe in diesem Konstruktor erscheint nie auf dem Bildschirm
          }
          else
          {
            // create group process object for all the other ranks in the current group
            _group_process = new GroupProcess(my_rank, rank_master, rank_load_bal, comm_group, rank_local);
            _group_process->wait();
            // @hilmar: Debug-Ausgabe in dieser Endlosschleife erscheint nie auf dem Bildschirm
          }
        }
      }
      assert(iter_MPC_rank == _num_processes_needed-1);

      // for the last rank (_num_processes-1) create master process object (responsible for screen output)
      if (my_rank == rank_master)
      {
        _master = new Master(my_rank);

        // start some infinite loop in the Master object, which waits for messages
        _master->wait();
        // @hilmar: Nur diese Bildschirmausgabe erscheint
      }
    }

    /**
     * \brief cleanup counterpart for _init().
     *
     * Deletes all dynamically allocated memory in _init().
     */
    void _cleanup()
    {
      // TODO: implement
    }

  /* ****************
   * public members *
   ******************/
  public:

    /* **************
     * constructors *
     ****************/

    /**
     * \brief simple constructor for creating a universe with exactly one work group of processes
     *
     * When using this constructor, only the master process and one work group without load balancer is created.
     * The constructor may only be called once.
     *
     * \param[in] argc
     * argument count passed to the main() method
     *
     * \param[in] argv
     * arguments passed to the main() method
     */
     // @Hilmar: Hab irgendwo gelesen, dass moeglichst alle elementaren member so initialisiert werden sollen, hab aber
     //          vergessen warum. Code entsprechend geaendert
    Universe(
      int& argc,
      char* argv[])
      : _num_process_groups(1)
    {
      if (_universe_created)
      {
        _abort("Only *one* universe can be created!");
      }
      _universe_created = true;
      std::cout << "Universe created!" << std::endl;
      _init_mpi(argc, argv);
      _num_processes_in_group = new int[_num_process_groups];
      _num_processes_in_group[0] = _num_processes - 1;
      _init();
    }

    /**
     * \brief constructor for creating a universe with several work groups
     *
     * When using this constructor, one can choose how many process groups are created. Additionally, one load balancer
     * per process group and one master process is created. The caller has to ensure that sufficient MPI processes are
     * available. The constructor may only be called once.
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
    Universe(
      int& argc,
      char* argv[],
      int& num_process_groups,
      int num_processes_in_group[])
      : _num_process_groups(num_process_groups),
        _num_processes_in_group(num_processes_in_group)
    {
      if (_universe_created)
      {
        _abort("Only *one* universe can be created!");
      }
      _universe_created = true;
      std::cout << "Universe created!" << std::endl;
      _init_mpi(argc, argv);

      _init();
    }

    /* ************
     * destructor *
     **************/
    /**
     * \brief destructor which automatically finalizes the MPI environment
     */
    ~Universe()
    {
      // clean up dynamically allocated memory
      _cleanup();

      // shut down MPI
      int mpi_is_initialised;
      MPI_Initialized(&mpi_is_initialised);
      if (mpi_is_initialised)
      {
        MPI_Finalize();
      }
      std::cout << "Universe destroyed!" << std::endl;
    }


    /* *****************
     * member functions*
     *******************/
    /**
     * \brief getter function for the load balancer object
     */
    LoadBalancer* get_load_balancer()
    {
      return _load_balancer;
    }

    /**
     * \brief getter function for the group process object
     */
    GroupProcess* get_group_process()
    {
      return _group_process;
    }

    /**
     * \brief getter function for the master process object
     */
    Master* get_master()
    {
      return _master;
    }

};
// initialise static member
bool Universe::_universe_created = false;

#endif // guard KERNEL_UNIVERSE_HPP
