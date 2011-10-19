/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR in this file.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
// includes, Feast
#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#ifdef PARALLEL
#include <kernel/util/mpi_utils.hpp>
#include <kernel/foundation/comm.hpp>
#include <kernel/foundation/process.hpp>
#include <kernel/foundation/universe.hpp>

#include <kernel/foundation/cell_subdivision.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
* \brief testing subdivision of an edge in 1D
*
* \test
* This test creates a universe, reads a 1D mesh (currently hard coded), builds a base mesh and subdivides an edge.
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
class Cell1DEdgeTestSubdivision
  : public TaggedTest<Tag_, DT_>
{

public:

  /// CTOR
  Cell1DEdgeTestSubdivision()
    : TaggedTest<Tag_, DT_>("cell_1d_edge-test_subdivision")
  {
  }


  /**
  * \brief sets the number of processes needed to 2 (1 manager + 1 master)
  *
  * \return number of processes needed
  */
  unsigned long mpi_proc_count() const
  {
    return 2;
  }


  /// main routine
  void run() const
  {
    // init MPI
    init_mpi();

    // set shortcut to the one and only instance of Universe (since this is the first call of
    // Universe<space_dim_, world_dim_>::instance(), it also calls the constructor of the Universe singleton class)
    Universe<space_dim_, world_dim_>* universe = Universe<space_dim_, world_dim_>::instance();

    // create universe consisting of one process group without dedicated load balancer,
    // let the outer test system catch eventual exceptions
    universe->create("cell_1d_edge_test_subdivision");

    // Get process objects. Note that on each process only one of the following two exists (the other one is the
    // null pointer).
    Manager<space_dim_, world_dim_>* manager = universe->manager();
    Master* master = universe->master();

    if(manager != nullptr)
    {
      // get pointer to the process group this process belongs to
      ProcessGroup* process_group = manager->process_group_main();

      // let the manager "read" the mesh, which currently means: create a hard-wired mesh consisting of 3 edges
      manager->read_mesh("dummy");

      // get pointer to the base mesh
      BaseMesh::BM<space_dim_, world_dim_>* bm = manager->base_mesh();

      // let the coordinator subdivide a cell
      if(process_group->is_coordinator())
      {
       // subdivide cell 1
        Logger::log("******************\n");
        Logger::log("Subdividing cell 1\n");
        Logger::log("******************\n");
        BaseMesh::SubdivisionData<space_dim_, space_dim_, world_dim_>* subdiv_data
          = new BaseMesh::SubdivisionData<space_dim_, space_dim_, world_dim_>(BaseMesh::CONFORM_SAME_TYPE);
        bm->cell(1)->subdivide(subdiv_data);

        // add created cells and subcells to the corresponding base mesh vectors
        bm->add_created_items(bm->cell(1)->subdiv_data());
        // set cell numbers (now they differ from indices)
        bm->set_cell_numbers();
        // print base mesh
        Logger::log(bm->print(), Logger::local_file_0);
       // validate base mesh
        Logger::log(bm->validate(), Logger::local_file_0);

        // COMMENT_HILMAR: neighbourhood update

        Logger::log("!!! Neighbourhood update after subdivision not implemented yet!!!\n");
        Logger::log("!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!\n");

// COMMENT_HILMAR: add some TEST_CHECK(...)

      }

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

}; // Cell1DEdgeTestSubdivision

// create test instance, using space and world dimension 1
Cell1DEdgeTestSubdivision<Nil, Nil, 1, 1> cell_1d_edge_test_subdivision;
#endif // PARALLEL
