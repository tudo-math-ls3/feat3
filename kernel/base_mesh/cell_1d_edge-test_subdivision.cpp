// includes, Feast
#include <kernel/base_header.hpp>
#include <kernel/util/string_utils.hpp>
#include <kernel/util/mpi_utils.hpp>
#include <test_system/test_system.hpp>
#include <kernel/comm.hpp>
#include <kernel/process.hpp>
#include <kernel/universe.hpp>
#include <kernel/base_mesh/cell_subdivision.hpp>

using namespace TestSystem;
using namespace FEAST;

/**
* \brief testing subdivision of an edge in 1D
*
* \test
* This test creates a universe, reads a 1D mesh (currently hard coded), builds a base mesh and subdivides an edge.
*
* \author Hilmar Wobker
*/
template <typename Tag_, typename DT_, unsigned char world_dim_, unsigned char space_dim_>
class Cell1DEdgeTestSubdivision
  : public TaggedTest<Tag_, DT_>
{

public:

  /// CTOR
  Cell1DEdgeTestSubdivision()
    : TaggedTest<Tag_, DT_>("cell_1d_edge-test_subdivision")
  {
  }


  /// sets the number of processes to use to 2 (1 load balancer + 1 master)
  unsigned long mpi_proc_count() const
  {
    return 2;
  }


  /// main routine
  void run() const
  {
    // init MPI
    MPIUtils::init_MPI();

    // set shortcut to the one and only instance of Universe (since this is the first call of
    // Universe<space_dim_, world_dim_>::instance(), it also calls the constructor of the Universe singleton class)
    Universe<space_dim_, world_dim_>* universe = Universe<space_dim_, world_dim_>::instance();

    // create universe, let the outer test system catch eventual exceptions
    universe->create("cell_1d_edge_test_subdivision");

    // Get process objects. Note that on each process only one of the following two exists (the other one is the
    // null pointer).
    LoadBalancer<space_dim_, world_dim_>* load_balancer = universe->load_balancer();
    Master* master = universe->master();

    if(load_balancer != nullptr)
    {
      ProcessGroup* process_group = load_balancer->process_group();

      // let the load balancer "read" the mesh, which currently means: create a hard-wired mesh consisting of 3 edges
      load_balancer->read_mesh("dummy");

      // get pointer to the base mesh
      BaseMesh::BM<space_dim_, world_dim_>* bm = load_balancer->base_mesh();

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

}; // Cell1DEdgeTestSubdivision

// create test instance, using space and world dimension 1
Cell1DEdgeTestSubdivision<Nil, Nil, 1, 1> cell_1d_edge_test_subdivision;
