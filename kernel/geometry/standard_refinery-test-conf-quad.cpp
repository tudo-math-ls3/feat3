#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/tetris_quad.hpp>
#include <kernel/geometry/standard_refinery.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;
using namespace FEAST::Geometry::Conformal;

typedef ConformalMeshPolicy<Shape::Quadrilateral> RootMeshPolicy;
typedef ConformalSubMeshPolicy<Shape::Quadrilateral> SubMeshPolicy;

typedef ConformalMesh<RootMeshPolicy> RootMesh;
typedef ConformalSubMesh<SubMeshPolicy>  SubMesh;

typedef StandardRefinery<RootMesh> RootMeshRefinery;
typedef StandardRefinery<SubMesh> SubMeshRefinery;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Peter Zajac
 */
class StandardRefineryTestConfQuad
  : public TestSystem::TaggedTest<Nil, Nil>
{
public:
  StandardRefineryTestConfQuad() :
    TestSystem::TaggedTest<Nil, Nil>("standard_refinery-test-conf-quad")
  {
  }

  virtual void run() const
  {
    // create a 2D tetris mesh
    RootMesh* quad_mesh_coarse = TestAux::create_tetris_mesh_2d();

    // create an edge submesh
    SubMesh* edge_submesh_coarse = TestAux::create_tetris_edge_submesh_2d();

    // create a quad submesh
    SubMesh* quad_submesh_coarse = TestAux::create_tetris_quad_submesh_2d();

    // create refineries
    RootMeshRefinery* quad_mesh_refinery = new RootMeshRefinery(*quad_mesh_coarse);
    SubMeshRefinery* edge_submesh_refinery = new SubMeshRefinery(*edge_submesh_coarse);
    SubMeshRefinery* quad_submesh_refinery = new SubMeshRefinery(*quad_submesh_coarse);

    // refine the meshes
    RootMesh* quad_mesh_fine = quad_mesh_refinery->refine();
    SubMesh* edge_submesh_fine = edge_submesh_refinery->refine(*quad_mesh_coarse);
    SubMesh* quad_submesh_fine = quad_submesh_refinery->refine(*quad_mesh_coarse);

    // validate refined meshes
    try
    {
      TestAux::validate_refined_tetris_mesh_2d(*quad_mesh_fine);
      TestAux::validate_refined_tetris_edge_submesh_2d(*edge_submesh_fine);
      TestAux::validate_refined_tetris_quad_submesh_2d(*quad_submesh_fine);
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    // clean up
    delete quad_submesh_fine;
    delete edge_submesh_fine;
    delete quad_mesh_fine;
    delete quad_submesh_refinery;
    delete edge_submesh_refinery;
    delete quad_mesh_refinery;
    delete quad_submesh_coarse;
    delete edge_submesh_coarse;
    delete quad_mesh_coarse;
  }
} standard_refinery_test_conf_quad;
