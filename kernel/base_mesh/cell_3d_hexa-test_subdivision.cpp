// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <test_system/test_system.hpp>
#include <kernel/comm.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>
#include <kernel/base_mesh/cell_subdivision.hpp>

/**
* \brief testing subdivision of a hexa in 3D
*
* \test
* This test creates a universe, reads a 3D mesh (currently hard coded), builds a base mesh and subdivides a hexa.
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
class Cell3DHexaTestSubdivision
  : public TaggedTest<Tag_, DT_>
{

public:

  /// CTOR
  Cell3DHexaTestSubdivision()
    : TaggedTest<Tag_, DT_>("cell_3d_hexa_test_subdivision")
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
    universe->create("cell_3d_hexa-test_subdivision");

    // Get process objects. Note that on each process only one of the following two exists (the other one is the
    // null pointer).
    Manager<space_dim_, world_dim_>* manager = universe->manager();
    Master* master = universe->master();

    if(manager != nullptr)
    {
      // get pointer to the process group this process belongs to
      ProcessGroup* process_group = manager->process_group_main();

      // let the manager "read" the mesh, which currently means: create a hard-wired mesh consisting of 4 hexas
      manager->read_mesh("dummy");

      // get pointer to the base mesh
      BaseMesh::BM<space_dim_, world_dim_>* bm = manager->base_mesh();

      // let the coordinator subdivide a cell
      if(process_group->is_coordinator())
      {
       // subdivide cell 0
        Logger::log("******************\n");
        Logger::log("Subdividing cell 0\n");
        Logger::log("******************\n");
        BaseMesh::SubdivisionData<space_dim_, space_dim_, world_dim_>* subdiv_data
          = new BaseMesh::SubdivisionData<space_dim_, space_dim_, world_dim_>(BaseMesh::NONCONFORM_SAME_TYPE);
        bm->cell(0)->subdivide(subdiv_data);

        // add created cells and subcells to the corresponding base mesh vectors
        bm->add_created_items(bm->cell(0)->subdiv_data());
        // set cell numbers (now they differ from indices)
        bm->set_cell_numbers();
        // print base mesh
        bm->print(Logger::file);
       // validate base mesh
        bm->validate(Logger::file);

        // TODO: neighbourhood update

        Logger::log("!!! Neighbourhood update after subdivision not implemented yet!!!\n");
        Logger::log("!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!\n");

// add some TEST_CHECK(...)

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

}; // Cell3DHexaTestSubdivision

// create test instance, using space and world dimension 3
Cell3DHexaTestSubdivision<Nil, Nil, 3, 3> cell_3d_hexa_test_subdivision;
