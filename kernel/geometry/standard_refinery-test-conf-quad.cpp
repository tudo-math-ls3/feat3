#include <test_system/test_system.hpp>
#include <kernel/archs.hpp>
#include <kernel/geometry/test_aux/tetris_quad.hpp>
#include <kernel/geometry/test_aux/standard_quad.hpp>
#include <kernel/geometry/standard_refinery.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;
using namespace FEAST::Geometry::Conformal;

typedef ConformalMeshPolicy<Shape::Quadrilateral> RootMeshPolicy;
typedef ConformalSubMeshPolicy<Shape::Quadrilateral> SubMeshPolicy;

typedef ConformalMesh<RootMeshPolicy> RootMesh;
typedef ConformalSubMesh<SubMeshPolicy> SubMesh;
typedef CellSubSet<Shape::Quadrilateral> SubSet;

typedef StandardRefinery<RootMesh> RootMeshRefinery;
typedef StandardRefinery<SubMesh> SubMeshRefinery;
typedef StandardRefinery<SubSet> SubSetRefinery;


/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Peter Zajac
 */
class StandardRefineryTestConfQuad
  : public TestSystem::TaggedTest<Archs::None, Nil>
{
public:
  StandardRefineryTestConfQuad() :
    TestSystem::TaggedTest<Archs::None, Nil>("standard_refinery-test-conf-quad")
  {
  }

  virtual void run() const
  {
    quad_std_test();
    quad_tetris_test();
  }

  void quad_tetris_test() const
  {
    // create a 2D tetris mesh
    RootMesh* quad_mesh_coarse = TestAux::create_tetris_mesh_2d();

    // create an edge submesh
    SubMesh* edge_submesh_coarse = TestAux::create_tetris_edge_submesh_2d();

    // create a quad submesh
    SubMesh* quad_submesh_coarse = TestAux::create_tetris_quad_submesh_2d();

    // create a cell sub-set
    SubSet* cell_subset_coarse = TestAux::create_tetris_quad_cellsubset_2d();

    // create refineries
    RootMeshRefinery* quad_mesh_refinery = new RootMeshRefinery(*quad_mesh_coarse);
    SubMeshRefinery* edge_submesh_refinery = new SubMeshRefinery(*edge_submesh_coarse);
    SubMeshRefinery* quad_submesh_refinery = new SubMeshRefinery(*quad_submesh_coarse);
    SubSetRefinery* cell_subset_refinery = new SubSetRefinery(*cell_subset_coarse);

    // refine the meshes
    RootMesh* quad_mesh_fine = quad_mesh_refinery->refine();
    SubMesh* edge_submesh_fine = edge_submesh_refinery->refine(*quad_mesh_coarse);
    SubMesh* quad_submesh_fine = quad_submesh_refinery->refine(*quad_mesh_coarse);
    SubSet* cell_subset_fine = cell_subset_refinery->refine(*quad_mesh_coarse);

    // validate refined meshes
    try
    {
      TestAux::validate_refined_tetris_mesh_2d(*quad_mesh_fine);
      TestAux::validate_refined_tetris_edge_submesh_2d(*edge_submesh_fine);
      TestAux::validate_refined_tetris_quad_submesh_2d(*quad_submesh_fine);
      TestAux::validate_refined_tetris_quad_cellsubset_2d(*cell_subset_fine);
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    // clean up
    delete cell_subset_fine;
    delete quad_submesh_fine;
    delete edge_submesh_fine;
    delete quad_mesh_fine;
    delete cell_subset_refinery;
    delete quad_submesh_refinery;
    delete edge_submesh_refinery;
    delete quad_mesh_refinery;
    delete cell_subset_coarse;
    delete quad_submesh_coarse;
    delete edge_submesh_coarse;
    delete quad_mesh_coarse;
  }

  void quad_std_test() const
  {
    RootMesh* quad_mesh_coarse;
    RootMeshRefinery* quad_mesh_refinery;
    RootMesh* quad_mesh_fine;
    try
    {
      for(Index i(0); i < 4; ++i) //loop over all possible orientations (max. 8)
      {
        // create a 2D quad element mesh
        quad_mesh_coarse = TestAux::create_quad_mesh_2d(i);

        // create refineries
        quad_mesh_refinery = new RootMeshRefinery(*quad_mesh_coarse);

        // refine the meshes
        quad_mesh_fine = quad_mesh_refinery->refine();

        // validate refined meshes
        TestAux::validate_refined_quad_mesh_2d(*quad_mesh_fine,i);

        // clean up
        delete quad_mesh_fine;
        delete quad_mesh_refinery;
        delete quad_mesh_coarse;
      }
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  }

} standard_refinery_test_conf_quad;
