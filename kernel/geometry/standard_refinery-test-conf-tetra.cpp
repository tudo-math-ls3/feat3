#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/standard_tetra.hpp>
#include <kernel/geometry/standard_refinery.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;
using namespace FEAST::Geometry::Conformal;

typedef ConformalMeshPolicy<Shape::Tetrahedron> TetRootMeshPolicy;
typedef ConformalSubMeshPolicy<Shape::Tetrahedron> SubMeshPolicy;

typedef ConformalMesh<TetRootMeshPolicy> TetRootMesh;
typedef ConformalSubMesh<SubMeshPolicy> SubMesh;

typedef StandardRefinery<TetRootMesh> TetRootMeshRefinery;
typedef StandardRefinery<SubMesh> SubMeshRefinery;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Constantin Christof
 */

class StandardRefineryTestConfTetrahedron
  : public TestSystem::TaggedTest<Archs::None, Nil>
{
public:
  StandardRefineryTestConfTetrahedron() :
    TestSystem::TaggedTest<Archs::None, Nil>("standard_refinery-test-conf-tetrahedron")
  {
  }

  virtual void run() const
  {
    tetra_std_test();
    tetra_block_test();
  }

  void tetra_std_test() const
  {

    TetRootMesh* tetrahedron_mesh_coarse;
    TetRootMeshRefinery* tetrahedron_mesh_refinery;
    TetRootMesh* tetrahedron_mesh_fine;

    try
    {
      for(Index i(0); i < 3; ++i) //loop over all possible orientations
      {
        // create a 3D tetrahedron element mesh
        tetrahedron_mesh_coarse = TestAux::create_tetra_mesh_3d(i);

        // create refineries
        tetrahedron_mesh_refinery = new TetRootMeshRefinery(*tetrahedron_mesh_coarse);

        // refine the meshes
        tetrahedron_mesh_fine = tetrahedron_mesh_refinery->refine();

        // validate refined meshes
        TestAux::validate_refined_tetra_mesh_3d(*tetrahedron_mesh_fine,i);

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

  void tetra_block_test() const
  {

    TetRootMesh* tetrahedron_mesh_coarse;
    TetRootMeshRefinery* tetrahedron_mesh_refinery;
    TetRootMesh* tetrahedron_mesh_fine;

    try
    {
      // create a 3D tetrahedron element mesh
      tetrahedron_mesh_coarse = TestAux::create_block_tetra_mesh_3d();

      // create tria submesh
      SubMesh* tria_submesh_coarse = TestAux::create_block_tria_submesh_3d();

      // create refineries
      tetrahedron_mesh_refinery = new TetRootMeshRefinery(*tetrahedron_mesh_coarse);
      SubMeshRefinery* tria_submesh_refinery = new SubMeshRefinery(*tria_submesh_coarse);

      // refine the meshes
      tetrahedron_mesh_fine = tetrahedron_mesh_refinery->refine();
      SubMesh* tria_submesh_fine = tria_submesh_refinery->refine(*tetrahedron_mesh_coarse);

      // validate refined meshes
      TestAux::validate_refined_block_tetra_mesh_3d(*tetrahedron_mesh_fine);
      TestAux::validate_refined_block_tria_submesh_3d(*tria_submesh_fine);

      // clean up
      delete tetrahedron_mesh_fine;
      delete tria_submesh_fine;
      delete tetrahedron_mesh_refinery;
      delete tria_submesh_refinery;
      delete tetrahedron_mesh_coarse;
      delete tria_submesh_coarse;
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  }

} standard_refinery_test_conf_tetrahedron;
