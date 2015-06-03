#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/standard_tria.hpp>
#include <kernel/geometry/conformal_mesh.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;

typedef ConformalMesh<Shape::Triangle> RootMesh;
typedef MeshPart<RootMesh> SubMesh;

typedef StandardRefinery<RootMesh> RootMeshRefinery;
typedef StandardRefinery<SubMesh, RootMesh> SubMeshRefinery;

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
    // loop over all possible orientations
    for(int i(0); i < 4; ++i)
    {
      // create a 2D triangle element mesh
      RootMesh* triangle_mesh_coarse = TestAux::create_tria_mesh_2d(i);

      try
      {
        // create refineries
        RootMeshRefinery triangle_mesh_refinery(*triangle_mesh_coarse);

        // refine the meshes
        RootMesh triangle_mesh_fine(triangle_mesh_refinery);

        // validate refined meshes
        TestAux::validate_refined_tria_mesh_2d(triangle_mesh_fine, i);
      }
      catch(const String& msg)
      {
        TEST_CHECK_MSG(false, msg);
      }

      // clean up
      delete triangle_mesh_coarse;
    }
  }

  void tria_patch_test() const
  {
    // create a 2D triangle element mesh
    RootMesh* triangle_mesh_coarse = TestAux::create_patch_tria_mesh_2d();

    // create an edge submesh
    SubMesh* edge_submesh_coarse = TestAux::create_patch_edge_submesh_2d();

    // create a tria submesh
    SubMesh* tria_submesh_coarse = TestAux::create_patch_tria_submesh_2d();

    try
    {
      // create refineries
      RootMeshRefinery triangle_mesh_refinery(*triangle_mesh_coarse);
      SubMeshRefinery edge_submesh_refinery(*edge_submesh_coarse, *triangle_mesh_coarse);
      SubMeshRefinery tria_submesh_refinery(*tria_submesh_coarse, *triangle_mesh_coarse);

      // refine the meshes
      RootMesh triangle_mesh_fine(triangle_mesh_refinery);
      SubMesh edge_submesh_fine(edge_submesh_refinery);
      SubMesh tria_submesh_fine(tria_submesh_refinery);

      // validate refined meshes
      TestAux::validate_refined_patch_tria_mesh_2d(triangle_mesh_fine);
      TestAux::validate_refined_patch_edge_submesh_2d(edge_submesh_fine);
      TestAux::validate_refined_patch_tria_submesh_2d(tria_submesh_fine);
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    // clean up
    delete triangle_mesh_coarse;
    delete tria_submesh_coarse;
    delete edge_submesh_coarse;
  }
} standard_refinery_test_conf_triangle;
