// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/util/assertion.hpp>
#include <test_system/test_system.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>

using namespace TestSystem;
using namespace FEAST;

/**
* \brief testing creation of work groups based on a 2D mesh
*
* \test
* This test creates a universe, reads a 2D mesh, builds a base mesh and creates work groups.
*
* \tparam Tag_
* description missing
*
* \tparam DT_
* description missing
*
* \tparam space_dim_
* space dimension (must be <= world_dim_; it is < world_dim_, e.g., when doing FE on 2D surfaces in a 3D world)
*
* \tparam world_dim_
* world dimension (determines the number of coordinates)
*
* \author Hilmar Wobker
*/
template<
  typename Tag_,
  typename DT_,
  unsigned char space_dim_,
  unsigned char world_dim_>
class LoggerTest
  : public TaggedTest<Tag_, DT_>
{

public:

  /**
  * \brief CTOR
  *
  * \param[in] mesh_file
  * name of the mesh file to be read
  */
  LoggerTest()
    : TaggedTest<Tag_, DT_>("logger_test")
  {
  }

  /**
  * \brief sets the number of processes to use to 10 (2 + 3 + 4 for 3 process groups + 1 master)
  *
  * For the hard-wired example creating three process groups we need 10 processes in total (2, 3 and 4 processes, resp.,
  * for the process groups and 1 for the master).
  */
  unsigned long mpi_proc_count() const
  {
    return 10;
  }


  /// main routine
  void run() const
  {
    CONTEXT("LoggerTest::run()");
    // init MPI
    init_mpi();

    // set number of process groups
    unsigned int num_process_groups(3);

    // array of numbers of processes in process groups (must be provided when num_process_groups > 1)
    unsigned int* num_processes_in_group = new unsigned int[num_process_groups];

    // set number of processes for the first process groups
    num_processes_in_group[0] = 2;
    num_processes_in_group[1] = 3;
    num_processes_in_group[2] = 4;

    // array of flags whether a dedicated load balancer process is needed in process groups
    // (must be provided when num_process_groups > 1)
    bool* includes_dedicated_load_bal = new bool[num_process_groups];
    includes_dedicated_load_bal[0] = false;
    includes_dedicated_load_bal[1] = false;
    includes_dedicated_load_bal[2] = false;

    // set shortcut to the one and only instance of Universe (since this is the first call of
    // Universe<space_dim_, world_dim_>::instance(), it also calls the constructor of the Universe singleton class)
    Universe<space_dim_, world_dim_>* universe = Universe<space_dim_, world_dim_>::instance();

    // create universe, let the outer test system catch eventual exceptions
    universe->create(num_process_groups, num_processes_in_group, includes_dedicated_load_bal, "logger_test");

    // Get process objects. Note that on each process only one of the following two exists (the other one is the
    // null pointer).
    LoadBalancer<space_dim_, world_dim_>* load_balancer = universe->load_balancer();
    Master* master = universe->master();

    if(load_balancer != nullptr)
    {
      ProcessGroup* process_group = load_balancer->process_group();
      unsigned int group_id = process_group->group_id();

      int rank_process_group = process_group->rank();

      std::string prefix(std::string("Group " + stringify(group_id) + ", proc " + stringify(rank_process_group)));

      // let the coordinators perform some extra log output
      if(process_group->is_coordinator())
      {
        PrettyPrinter pp(25 + 10*(group_id+1), '#', prefix + " ");
        pp.add_line_sep();
        pp.add_line_centered("Logger test");
        pp.add_line_sep();
        pp.add_line("bla blub");
        pp.add_line("long long long long long long long long long long");
        pp.add_line_sep();
        // print it like this...
        Logger::log(pp.block());
        // ... or like this...
        pp.print(Logger::file);
        // send it to the master for screen and file output
        Logger::log_master(pp.block());

        // test log_master_array()
        std::vector<std::string> messages(3);
        messages[0] = prefix + ": Testing this...\n";
        messages[1] = prefix + ": ...vector logging...\n";
        messages[2] = prefix + ": ...feature!\n";
        Logger::log_master_array(messages, Logger::SCREEN_FILE);
        Logger::log(messages);
      }

      // let all processes write something to their log file
      std::string s("Process with world rank " + stringify(Process::rank) + " has local rank "
                    + stringify(rank_process_group) + " in group " + stringify(group_id) + ".\n");
      Logger::log(s);

      // perform actions depending on the group id
      if(group_id == 0)
      {
        process_group->log_indiv_master(prefix + ": NI!\n", Logger::SCREEN_FILE);
      }
      else if(group_id == 1)
      {
        process_group->log_indiv_master(prefix + ": Ecky Ecky Ecky F'tang F'tang\n", Logger::SCREEN);
      }
      else if(group_id == 2)
      {
        process_group->log_indiv_master(prefix + ": A SHRUBBERY!\n", Logger::FILE);
      }
      else
      {
        throw InternalError("group_id " + stringify(group_id) + " is not valid!");
      }

// TODO: add some TEST_CHECK(...)

      // Everything done, call universe destruction routine.
      universe->destroy();
    }
    else if(master != nullptr)
    {
      // This branch is entered when the infinite service loop of the master has been finished.
      // This, however, usually happens only at program end.
      // Everything done, call universe destruction routine.
      universe->destroy();
    }
    else
    {
      // This branch must not be entered. Throw InternalError which is caught by outer test system.
      throw InternalError("Process with rank " + stringify(Process::rank)
                          + " has no particular role, this should not happen.");
    }
  } // run()
}; // LoggerTest

// create test instance, using space and world dimension 1
LoggerTest<Nil, Nil, 1, 1> logger_test;
