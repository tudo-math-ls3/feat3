#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/standard_hexa.hpp>
#include <kernel/geometry/test_aux/tetris_hexa.hpp>
#include <kernel/geometry/conformal_mesh.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

typedef ConformalMesh<Shape::Hexahedron> HexaMesh;
typedef MeshPart<HexaMesh> SubMesh;

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
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  StandardRefineryTestConfHexa() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("standard_refinery-test-conf-hexa")
  {
  }

  virtual ~StandardRefineryTestConfHexa()
  {
  }

  virtual void run() const override
  {
    hexa_std_test();
    hexa_tetris_test();
  }

  void hexa_std_test() const
  {
    // loop over all possible orientations
    for(int i(0); i < 3; ++i)
    {
      // create a 3D hex-element mesh
      HexaMesh* hexa_mesh_coarse = TestAux::create_hexa_mesh_3d(i);

      try
      {
        // create refineries
        HexaMeshRefinery hexa_mesh_refinery(*hexa_mesh_coarse);

        // refine the meshes
        HexaMesh hexa_mesh_fine(hexa_mesh_refinery);

        // validate refined meshes
        TestAux::validate_refined_hexa_mesh_3d(hexa_mesh_fine, i);
      }
      catch(const String& msg)
      {
        TEST_CHECK_MSG(false, msg);
      }

      // clean up
      delete hexa_mesh_coarse;
    }
  }

  void hexa_tetris_test() const
  {
    // create a 3D tetris mesh
    HexaMesh* hexa_mesh_coarse = TestAux::create_tetris_mesh_3d();

    // create a quad submesh
    SubMesh* quad_submesh_coarse = TestAux::create_tetris_quad_submesh_3d();

    try
    {
      // create refineries
      HexaMeshRefinery hexa_mesh_refinery(*hexa_mesh_coarse);
      SubMeshRefinery quad_submesh_refinery(*quad_submesh_coarse, *hexa_mesh_coarse);

      // refine the meshes
      HexaMesh hexa_mesh_fine(hexa_mesh_refinery);
      SubMesh quad_submesh_fine(quad_submesh_refinery);

      // validate refined meshes
      TestAux::validate_refined_tetris_mesh_3d(hexa_mesh_fine);
      TestAux::validate_refined_tetris_quad_submesh_3d(quad_submesh_fine);
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    // clean up
    delete hexa_mesh_coarse;
    delete quad_submesh_coarse;
  }

} standard_refinery_test_conf_hexa;
