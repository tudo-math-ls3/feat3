/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR in this file.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/foundation/logger.hpp>
#include <test_system/test_system.hpp>
#ifdef PARALLEL
#include <kernel/util/mpi_utils.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/foundation/process.hpp>
#include <kernel/foundation/universe.hpp>
#endif // PARALLEL

using namespace FEAST;
using namespace FEAST::TestSystem;

#ifdef PARALLEL
/* ********************************************************************************************* */
/*   -  P A R A L L E L   T E S T   D R I V E R  *  P A R A L L E L   T E S T   D R I V E R   -  */
/* ********************************************************************************************* */
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
    : TaggedTest<Tag_, DT_>("logger_test (parallel)")
  {
  }

  /**
  * \brief sets the number of processes needed to 10 (2 + 3 + 4 for 3 process groups + 1 master)
  *
  * For the hard-wired example creating three process groups we need 10 processes in total (2, 3 and 4 processes, resp.,
  * for the process groups and 1 for the master).
  *
  * \return number of processes needed
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
    Manager<space_dim_, world_dim_>* manager = universe->manager();
    Master* master = universe->master();

    if(manager != nullptr)
    {
      ProcessGroup* process_group = manager->process_group_main();
      unsigned int group_id = process_group->group_id();

      int rank_process_group = process_group->rank();

      String prefix(String("Group " + stringify(group_id) + ", proc " + stringify(rank_process_group)));

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
        // send it to the master for screen and file output
        Logger::log(pp.block(), Logger::master);
      }

      // let all processes write something to their log file
      String s("Process with world rank " + stringify(Process::rank) + " has local rank "
                    + stringify(rank_process_group) + " in group " + stringify(group_id) + ".\n");
      Logger::log(s);

      // perform actions depending on the group id
      /*if(group_id == 0)
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
      }*/

      // Everything done, call universe destruction routine.
      universe->destroy();

      // COMMENT_HILMAR:
      // It would be nice to have some TEST_CHECK(...) here... But it's hard to check whether some screen output
      // appeared as expected. Concerning logfile output... maybe one could read the content of the file and compare
      // it to what is expected to be in the file...
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

#else
/* ************************************************************************************* */
/*   -  S E R I A L   T E S T   D R I V E R  *  S E R I A L   T E S T   D R I V E R   -  */
/* ************************************************************************************* */
class LoggerTest
  : public TaggedTest<Nil, Nil>
{
public:
  /// CTOR
  LoggerTest()
    : TaggedTest<Nil, Nil>("logger_test (serial)")
  {
  }

  /// main routine
  void run() const
  {
    CONTEXT("LoggerTest::run()");

    // open two log files
    Logger::open("./logger_test_0", 0);
    Logger::open("./logger_test_1", 1);

    // print something
    Logger::log("This should appear on the screen and nowhere else\n", Logger::screen);
    Logger::log("This should appear on the screen and in log file 0\n", Logger::local);
    Logger::log("This should appear in log file 1\n", Logger::local_file_1);
    Logger::log("This should appear in log file 0 and 1\n", Logger::local_file_0|Logger::local_file_1);
    Logger::log("This should not appear anywhere\n", Logger::none);

    // close log files
    Logger::close_all();
  }
}; // LoggerTest

LoggerTest logger_test;

#endif // PARALLEL
