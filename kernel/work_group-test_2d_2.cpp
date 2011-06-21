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
class WorkGroupTest2D2
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
  WorkGroupTest2D2(std::string const mesh_file)
    : TaggedTest<Tag_, DT_>("work_group_test_2d_2"),
      _mesh_file(mesh_file)
  {
  }


  /**
  * \brief sets the number of processes to use to 21
  *
  * For the semi-hard-wired example using a mesh with n base mesh cells we need n+5 processes in total. Two
  * 'main process groups' (for two completely independent problems) are created. Each process group has its own
  * manager. Only the first process group (consisting of n+2 processes does some useful things, the second one
  * (consisting of 2 processes) is idle.
  *
  * \return number of processes needed
  *
  * \todo Instead of returning a hard-coded 21 here, get the number of BMCs from the name of the mesh file or somehow
  * else. (Also store it in some variable, it is needed in run().)
  */
  unsigned long mpi_proc_count() const
  {
    return 21;
  }


  /// main routine
  void run() const
  {
    CONTEXT("WorkGroupTest2D2::run()");
    // init MPI
    init_mpi();

    // set number of process groups
    unsigned int num_process_groups(2);

    // TODO: change this hard-coded 16 to the number determined in mpi_proc_count().
    unsigned int num_cells = 16;
    ASSERT(num_cells > 0, "Number of cells must not be zero.");

    // array of numbers of processes in process groups (must be provided when num_process_groups > 1)
    unsigned int* num_processes_in_group = new unsigned int[num_process_groups];

    // set number of processes for the first process group to num_cells + 2
    num_processes_in_group[0] = num_cells + 2;
    // set number of processes for the second process group to 2
    num_processes_in_group[1] = 2;

    // array of flags whether a dedicated load balancer process is needed in process groups
    // (must be provided when num_process_groups > 1)
    // Set the first entry either to false or to true to test two different configurations. (You don't have to change
    // the number of processes for that.) (See description of the example in routine define_work_groups().)
    bool* includes_dedicated_load_bal = new bool[num_process_groups];
    includes_dedicated_load_bal[0] = false;
    includes_dedicated_load_bal[1] = false;

    // set shortcut to the one and only instance of Universe (since this is the first call of
    // Universe<space_dim_, world_dim_>::instance(), it also calls the constructor of the Universe singleton class)
    Universe<space_dim_, world_dim_>* universe = Universe<space_dim_, world_dim_>::instance();

    // create universe, let the outer test system catch eventual exceptions
    universe->create(num_process_groups, num_processes_in_group, includes_dedicated_load_bal, "work_group_test_2d_2");

    // Get process objects. Note that on each process only one of the following two exists (the other one is the
    // null pointer).
    Manager<space_dim_, world_dim_>* manager = universe->manager();
    Master* master = universe->master();

    if(manager != nullptr)
    {
      // get pointer to the process group this process belongs to
      ProcessGroup* process_group = manager->process_group_main();
      unsigned int group_id = process_group->group_id();

      // debug output
      int rank_process_group = process_group->rank();
      std::string s("Process " + stringify(Process::rank) + " is the manager with local rank "
                    + stringify(rank_process_group) + " in group " + stringify(group_id) + ".\n");
      Logger::log(s);

      // perform actions depending on the group id
      if(group_id == 0)
      {
        // let the manager read the mesh file and create a base mesh (this is only done by the coordinator process)
        manager->read_mesh(_mesh_file);

        // let the load balancer define the work groups and interlevel groups (currently hard-coded)
        manager->define_work_groups_2();

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
      else
      {
        // the second process group does something else, programmed by the application outside the kernel...
        // Everything done, call universe destruction routine.
        universe->destroy();
     }
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
}; // WorkGroupTest2D2

// create test instance, using space and world dimension 2
WorkGroupTest2D2<Nil, Nil, 2, 2> work_group_test_2d_2(stringify(FEAST_SRC_DIR) + "/data/meshes/test_16bmc.feast");
