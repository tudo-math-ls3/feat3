#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/standard_hexa.hpp>
#include <kernel/geometry/test_aux/tetris_hexa.hpp>
#include <kernel/geometry/standard_refinery.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;
using namespace FEAST::Geometry::Conformal;

typedef ConformalMeshPolicy<Shape::Hexahedron> HexaMeshPolicy;
typedef ConformalSubMeshPolicy<Shape::Hexahedron> SubMeshPolicy;

typedef ConformalMesh<HexaMeshPolicy> HexaMesh;
typedef ConformalSubMesh<SubMeshPolicy> SubMesh;

typedef StandardRefinery<HexaMesh> HexaMeshRefinery;
typedef StandardRefinery<SubMesh> SubMeshRefinery;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Constantin Christof
 */
class StandardRefineryTestConfHexa
  : public TestSystem::TaggedTest<Nil, Nil>
{
public:
  StandardRefineryTestConfHexa() :
    TestSystem::TaggedTest<Nil, Nil>("standard_refinery-test-conf-hexa")
  {
  }

  virtual void run() const
  {
    hexa_refinementtest();
    hexa_tetristest();
  }

  void hexa_refinementtest() const
  {
    HexaMesh* hexa_mesh_coarse;
    HexaMeshRefinery* hexa_mesh_refinery;
    HexaMesh* hexa_mesh_fine;
    try
    {
      for(Index i(0); i < 3; ++i) //loop over all possible orientations
      {
        // create a 3D hex-element mesh
        hexa_mesh_coarse = TestAux::create_hexarefinement_mesh_3d(i);

        // create refineries
        hexa_mesh_refinery = new HexaMeshRefinery(*hexa_mesh_coarse);

        // refine the meshes
        hexa_mesh_fine = hexa_mesh_refinery->refine();

        // validate refined meshes
        TestAux::validate_refined_hexarefinement_mesh_3d(*hexa_mesh_fine,i);

        // clean up
        delete hexa_mesh_fine;
        delete hexa_mesh_refinery;
        delete hexa_mesh_coarse;
      }
    }
    catch(const String& msg)
        {
          TEST_CHECK_MSG(false, msg);
        }
  }

  void hexa_tetristest() const
  {
    // create a 3D tetris mesh
    HexaMesh* hexa_mesh_coarse = TestAux::create_tetris_mesh_3d();

    // create a quad submesh
    SubMesh* quad_submesh_coarse = TestAux::create_tetris_quad_submesh_3d();

    // create refineries
    HexaMeshRefinery* hexa_mesh_refinery = new HexaMeshRefinery(*hexa_mesh_coarse);
    SubMeshRefinery* quad_submesh_refinery = new SubMeshRefinery(*quad_submesh_coarse);

    // refine the meshes
    HexaMesh* hexa_mesh_fine = hexa_mesh_refinery->refine();
    SubMesh* quad_submesh_fine = quad_submesh_refinery->refine(*hexa_mesh_coarse);

    // validate refined meshes
    try
    {
      TestAux::validate_refined_tetris_mesh_3d(*hexa_mesh_fine);
      TestAux::validate_refined_tetris_quad_submesh_3d(*quad_submesh_fine);
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    // clean up
    delete hexa_mesh_fine;
    delete quad_submesh_fine;
    delete hexa_mesh_refinery;
    delete quad_submesh_refinery;
    delete hexa_mesh_coarse;
    delete quad_submesh_coarse;
  }

} standard_refinery_test_conf_hexa;
