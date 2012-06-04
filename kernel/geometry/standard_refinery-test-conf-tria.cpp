#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/standard_tria.hpp>
#include <kernel/geometry/standard_refinery.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;
using namespace FEAST::Geometry::Conformal;

typedef ConformalMeshPolicy<Shape::Triangle> RootMeshPolicy;
typedef ConformalMesh<RootMeshPolicy> RootMesh;
typedef StandardRefinery<RootMesh> RootMeshRefinery;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Constantin Christof
 */

class StandardRefineryTestConfTriangle
  : public TestSystem::TaggedTest<Nil, Nil>
{
public:
  StandardRefineryTestConfTriangle() :
    TestSystem::TaggedTest<Nil, Nil>("standard_refinery-test-conf-triangle")
  {
  }

  virtual void run() const
  {
    triangle_refinementtest();
    triangle_quad_refinementtest();
  }

  void triangle_refinementtest() const
  {
    RootMesh* triangle_mesh_coarse;
    RootMeshRefinery* triangle_mesh_refinery;
    RootMesh* triangle_mesh_fine;

    try
    {
      for(Index i(0); i < 4; ++i) //loop over all possible orientations
      {
        // create a 2D triangle element mesh
        triangle_mesh_coarse = TestAux::create_trianglerefinement_mesh_2d(i);

        // create refineries
        triangle_mesh_refinery = new RootMeshRefinery(*triangle_mesh_coarse);

        // refine the meshes
        triangle_mesh_fine = triangle_mesh_refinery->refine();

        // validate refined meshes
        TestAux::validate_refined_trianglerefinement_mesh_2d(*triangle_mesh_fine,i);

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

  void triangle_quad_refinementtest() const
  {
    RootMesh* triangle_mesh_coarse;
    RootMeshRefinery* triangle_mesh_refinery;
    RootMesh* triangle_mesh_fine;

    try
    {

        // create a 2D triangle element mesh
        triangle_mesh_coarse = TestAux::create_triangle_refinement_mesh_2d();

        // create refineries
        triangle_mesh_refinery = new RootMeshRefinery(*triangle_mesh_coarse);

        // refine the meshes
        triangle_mesh_fine = triangle_mesh_refinery->refine();

        // validate refined meshes
        TestAux::validate_refined_triangle_refinement_mesh_2d(*triangle_mesh_fine);

        // clean up
        delete triangle_mesh_fine;
        delete triangle_mesh_refinery;
        delete triangle_mesh_coarse;

    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  }

} standard_refinery_test_conf_triangle;
