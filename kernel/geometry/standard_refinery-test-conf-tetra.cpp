#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/standard_tetra.hpp>
#include <kernel/geometry/standard_refinery.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;
using namespace FEAST::Geometry::Conformal;

typedef ConformalMeshPolicy<Shape::Tetrahedron> TetRootMeshPolicy;
typedef ConformalMesh<TetRootMeshPolicy> TetRootMesh;
typedef StandardRefinery<TetRootMesh> TetRootMeshRefinery;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Constantin Christof
 */

class StandardRefineryTestConfTetrahedron
  : public TestSystem::TaggedTest<Nil, Nil>
{
public:
  StandardRefineryTestConfTetrahedron() :
    TestSystem::TaggedTest<Nil, Nil>("standard_refinery-test-conf-tetrahedron")
  {
  }

  virtual void run() const
  {
    tetrahedron_refinementtest();
    tetrahedron_refinementtest2();
  }

  void tetrahedron_refinementtest() const
  {

    TetRootMesh* tetrahedron_mesh_coarse;
    TetRootMeshRefinery* tetrahedron_mesh_refinery;
    TetRootMesh* tetrahedron_mesh_fine;

    try
    {
      for(Index i(0); i < 3; ++i) //loop over all possible orientations
      {
        // create a 3D tetrahedron element mesh
        tetrahedron_mesh_coarse = TestAux::create_tetrahedronrefinement_mesh_3d(i);

        // create refineries
        tetrahedron_mesh_refinery = new TetRootMeshRefinery(*tetrahedron_mesh_coarse);

        // refine the meshes
        tetrahedron_mesh_fine = tetrahedron_mesh_refinery->refine();

        // validate refined meshes
        TestAux::validate_refined_tetrahedronrefinement_mesh_3d(*tetrahedron_mesh_fine,i);

        // clean up
        delete tetrahedron_mesh_fine;
        delete tetrahedron_mesh_refinery;
        delete tetrahedron_mesh_coarse;
      }
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  }

  void tetrahedron_refinementtest2() const
  {

    TetRootMesh* tetrahedron_mesh_coarse;
    TetRootMeshRefinery* tetrahedron_mesh_refinery;
    TetRootMesh* tetrahedron_mesh_fine;

    try
    {
      // create a 3D tetrahedron element mesh
      tetrahedron_mesh_coarse = TestAux::create_bigtetrahedronrefinement_mesh_3d();

      // create refineries
      tetrahedron_mesh_refinery = new TetRootMeshRefinery(*tetrahedron_mesh_coarse);

      // refine the meshes
      tetrahedron_mesh_fine = tetrahedron_mesh_refinery->refine();

      // validate refined meshes
      TestAux::validate_refined_bigtetrahedronrefinement_mesh_3d(*tetrahedron_mesh_fine);

      // clean up
      delete tetrahedron_mesh_fine;
      delete tetrahedron_mesh_refinery;
      delete tetrahedron_mesh_coarse;
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  }

} standard_refinery_test_conf_tetrahedron;
