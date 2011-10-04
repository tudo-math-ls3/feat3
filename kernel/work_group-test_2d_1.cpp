/* GENERAL_REMARK_BY_HILMAR:
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
// includes, Feast
#include <kernel/base_header.hpp>
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
class WorkGroupTest2D1
  : public TaggedTest<Tag_, DT_>
{

private:

  /// name of the mesh file to be read
  std::string _mesh_file;

public:

  /**
  * \brief CTOR
  *
  * \param[in] mesh_file
  * name of the mesh file to be read
  */
  WorkGroupTest2D1(std::string const mesh_file)
    : TaggedTest<Tag_, DT_>("work_group_test_2d_1"),
      _mesh_file(mesh_file)
  {
  }


  /**
  * \brief sets the number of processes needed to 17 (16 BMCs + 1 master)
  *
  * For the semi-hard-wired example using a mesh with n base mesh cells we need n+1 processes in total
  * (n for the one and only process group and 1 for the master process).
  *
  * \return number of processes needed
  *
  * \todo Instead of returning a hard-coded 17 here, get the number of BMCs from the name of the mesh file or somehow
  * else.
  */
  unsigned long mpi_proc_count() const
  {
    return 17;
  }


  /// main routine
  void run() const
  {
    CONTEXT("WorkGroupTest2D1::run()");
    // init MPI
    init_mpi();

    // set shortcut to the one and only instance of Universe (since this is the first call of
    // Universe<space_dim_, world_dim_>::instance(), it also calls the constructor of the Universe singleton class)
    Universe<space_dim_, world_dim_>* universe = Universe<space_dim_, world_dim_>::instance();

    // create universe consisting of one process group without dedicated load balancer,
    // let the outer test system catch eventual exceptions
    universe->create("work_group_test_2d_1");

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

      // let the manager read the mesh file and create a base mesh (this is only done by the coordinator process)
      manager->read_mesh(_mesh_file);

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
}; // WorkGroupTest2D1

// create test instance, using space and world dimension 2
WorkGroupTest2D1<Nil, Nil, 2, 2> work_group_test_2d_1(stringify(FEAST_SRC_DIR) + "/data/meshes/test_16bmc.feast");
