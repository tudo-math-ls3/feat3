/* GENERAL_REMARK_BY_HILMAR:
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
// includes, Feast
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <kernel/util/assertion.hpp>
#include <test_system/test_system.hpp>
#include <kernel/comm.hpp>
#include <kernel/graph.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
* \brief testing creation of work groups based on a 1D mesh
*
* \test
* This test creates a universe, reads a 1D mesh (currently hard coded), builds a base mesh and creates work groups.
*
* \author Hilmar Wobker
*/
template <typename Tag_, typename DT_, unsigned char world_dim_, unsigned char space_dim_>
class WorkGroupTest1D
  : public TaggedTest<Tag_, DT_>
{

public:

  /// CTOR
  WorkGroupTest1D()
    : TaggedTest<Tag_, DT_>("work_group_test_1d")
  {
  }

  /**
  * \brief sets the number of processes needed to 4 (3 BMCs + 1 master)
  *
  * For this hard-wired example using the fixed 1D mesh defined in the file parser with 3 base mesh cells we need 4
  * processes in total (3 for the one and only process group and 1 for the master process).
  *
  * \return number of processes needed
  */
  unsigned long mpi_proc_count() const
  {
    return 4;
  }


  /// main routine
  void run() const
  {
    CONTEXT("WorkGroupTest1D::run()");
    // init MPI
    init_mpi();

    // set shortcut to the one and only instance of Universe (since this is the first call of
    // Universe<space_dim_, world_dim_>::instance(), it also calls the constructor of the Universe singleton class)
    Universe<space_dim_, world_dim_>* universe = Universe<space_dim_, world_dim_>::instance();

    // create universe consisting of one process group without dedicated load balancer,
    // let the outer test system catch eventual exceptions
    universe->create("work_group_test_1d");

    // Get process objects. Note that on each process only one of the following two exists (the other one is the
    // null pointer).
    Manager<space_dim_, world_dim_>* manager = universe->manager();
    Master* master = universe->master();

    if(manager != nullptr)
    {
      // get pointer to the process group this process belongs to
      ProcessGroup* process_group = manager->process_group_main();

      // debug output
      int rank_process_group = process_group->rank();
      std::string s("Process " + stringify(Process::rank) + " is manager process with local rank "
                    + stringify(rank_process_group) + ".\n");
      Logger::log(s);

      // let the manager "read" the mesh, which currently means: create a hard-wired mesh consisting of 3 edges
      manager->read_mesh("dummy");

      // let the load balancer define the work groups and interlevel groups (currently hard-coded)
      manager->define_work_groups_1();

      // let the manager create the work groups and interlevel groups
      manager->create_work_groups();

      // perform a hard-coded test checking whether work group and interlevel group communication works
      unsigned int test_result = manager->test_communication();

      // Everything done, call universe destruction routine.
      universe->destroy();

      // check whether the test was succesful (Don't try to do this check *before* destruction of the universe! If
      // the test failed (i.e., test_result > 0), then the program deadlocks.)
      TEST_CHECK(test_result == 0);

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
}; // WorkGroupTest1D

// create test instance, using space and world dimension 1
WorkGroupTest1D<Nil, Nil, 1, 1> work_group_test_1d;
#endif // PARALLEL
