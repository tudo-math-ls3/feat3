/* GENERAL_REMARK_BY_HILMAR:
 * See COMMENT_HILMAR in this file.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
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
* \brief testing subdivision of a quad in 2D
*
* \test
* This test creates a universe, reads a 2D mesh, builds a base mesh and subdivides a quad.
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
class Cell2DQuadTestSubdivision
  : public TaggedTest<Tag_, DT_>
{

private:

  /// name of the mesh file to be read
  std::string _mesh_file;

public:

  /// CTOR
  Cell2DQuadTestSubdivision(std::string const mesh_file)
    : TaggedTest<Tag_, DT_>("cell_2d_quad_test_subdivision"),
      _mesh_file(mesh_file)
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
    universe->create("cell_2d_quad-test_subdivision");

    // Get process objects. Note that on each process only one of the following two exists (the other one is the
    // null pointer).
    Manager<space_dim_, world_dim_>* manager = universe->manager();
    Master* master = universe->master();

    if(manager != nullptr)
    {
      // get pointer to the process group this process belongs to
      ProcessGroup* process_group = manager->process_group_main();

      // let the manager read the mesh
      manager->read_mesh(_mesh_file);

      // get pointer to the base mesh
      BaseMesh::BM<space_dim_, world_dim_>* bm = manager->base_mesh();

      // let the coordinator subdivide a cell
      if(process_group->is_coordinator())
      {
       // subdivide cell 0
        Logger::log("******************\n");
        Logger::log("Subdividing cell 2\n");
        Logger::log("******************\n");
        BaseMesh::SubdivisionData<space_dim_, space_dim_, world_dim_>* subdiv_data
          = new BaseMesh::SubdivisionData<space_dim_, space_dim_, world_dim_>(BaseMesh::NONCONFORM_SAME_TYPE);
        bm->cell(2)->subdivide(subdiv_data);

        // add created cells and subcells to the corresponding base mesh vectors
        bm->add_created_items(bm->cell(2)->subdiv_data());
        // set cell numbers (now they differ from indices)
        bm->set_cell_numbers();
        // print base mesh
        bm->print(Logger::file);
       // validate base mesh
        bm->validate(Logger::file);

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

}; // Cell2DQuadTestSubdivision

// create test instance, using space and world dimension 2
Cell2DQuadTestSubdivision<Nil, Nil, 2, 2>
  cell_2d_quad_test_subdivision(stringify(FEAST_SRC_DIR) + "/data/meshes/test_16bmc.feast");
