#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/standard_tria.hpp>
#include <kernel/geometry/standard_refinery.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;

typedef ConformalMeshPolicy<Shape::Triangle> RootMeshPolicy;
typedef ConformalSubMeshPolicy<Shape::Triangle> SubMeshPolicy;

typedef ConformalMesh<RootMeshPolicy> RootMesh;
typedef ConformalSubMesh<SubMeshPolicy> SubMesh;

typedef StandardRefinery<RootMesh> RootMeshRefinery;
typedef StandardRefinery<SubMesh> SubMeshRefinery;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Constantin Christof
 */

class StandardRefineryTestConfTriangle
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  StandardRefineryTestConfTriangle() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("standard_refinery-test-conf-triangle")
  {
  }

  virtual void run() const
  {
    tria_std_test();
    tria_patch_test();
  }

  void tria_std_test() const
  {
    RootMesh* triangle_mesh_coarse;
    RootMeshRefinery* triangle_mesh_refinery;
    RootMesh* triangle_mesh_fine;

    try
    {
      for(Index i(0); i < 4; ++i) //loop over all possible orientations
      {
        // create a 2D triangle element mesh
        triangle_mesh_coarse = TestAux::create_tria_mesh_2d(i);

        // create refineries
        triangle_mesh_refinery = new RootMeshRefinery(*triangle_mesh_coarse);

        // refine the meshes
        triangle_mesh_fine = triangle_mesh_refinery->refine();

        // validate refined meshes
        TestAux::validate_refined_tria_mesh_2d(*triangle_mesh_fine,i);

        // clean up
        delete triangle_mesh_fine;
        delete triangle_mesh_refinery;
        delete triangle_mesh_coarse;
      }
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  }

  void tria_patch_test() const
  {
    RootMesh* triangle_mesh_coarse;
    RootMeshRefinery* triangle_mesh_refinery;
    RootMesh* triangle_mesh_fine;

    // create an edge submesh
    SubMesh* edge_submesh_coarse = TestAux::create_patch_edge_submesh_2d();

    // create a tria submesh
    SubMesh* tria_submesh_coarse = TestAux::create_patch_tria_submesh_2d();

    try
    {

        // create a 2D triangle element mesh
        triangle_mesh_coarse = TestAux::create_patch_tria_mesh_2d();

        // create refineries
        triangle_mesh_refinery = new RootMeshRefinery(*triangle_mesh_coarse);
        SubMeshRefinery* edge_submesh_refinery = new SubMeshRefinery(*edge_submesh_coarse);
        SubMeshRefinery* tria_submesh_refinery = new SubMeshRefinery(*tria_submesh_coarse);

        // refine the meshes
        triangle_mesh_fine = triangle_mesh_refinery->refine();
        SubMesh* edge_submesh_fine = edge_submesh_refinery->refine(*triangle_mesh_coarse);
        SubMesh* tria_submesh_fine = tria_submesh_refinery->refine(*triangle_mesh_coarse);

        // validate refined meshes
        TestAux::validate_refined_patch_tria_mesh_2d(*triangle_mesh_fine);
        TestAux::validate_refined_patch_edge_submesh_2d(*edge_submesh_fine);
        TestAux::validate_refined_patch_tria_submesh_2d(*tria_submesh_fine);

        // clean up
        delete triangle_mesh_fine;
        delete tria_submesh_fine;
        delete edge_submesh_fine;
        delete triangle_mesh_refinery;
        delete tria_submesh_refinery;
        delete edge_submesh_refinery;
        delete triangle_mesh_coarse;
        delete tria_submesh_coarse;
        delete edge_submesh_coarse;

    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  }
} standard_refinery_test_conf_triangle;
