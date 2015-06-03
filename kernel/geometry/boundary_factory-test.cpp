#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/tetris_quad.hpp>
#include <kernel/geometry/test_aux/standard_quad.hpp>
#include <kernel/geometry/boundary_factory.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;

typedef ConformalMesh<Shape::Quadrilateral> RootMesh;
typedef MeshPart<RootMesh> SubSet;


typedef StandardRefinery<RootMesh> RootMeshRefinery;
typedef BoundaryFactory<RootMesh> RootBndFactory;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Peter Zajac
 */
class BoundaryFactoryTest
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  BoundaryFactoryTest() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("BoundaryFactoryTest")
  {
  }

  virtual void run() const
  {
    quad_tetris_test();
  }

  void quad_tetris_test() const
  {
    // create a 2D tetris mesh
    RootMesh* quad_mesh_coarse = TestAux::create_tetris_mesh_2d();

    {
      // create refinery
      RootMeshRefinery quad_mesh_refinery(*quad_mesh_coarse);

      // refine the mesh
      RootMesh quad_mesh_fine(quad_mesh_refinery);

      // create boundary factories
      RootBndFactory bnd_factory_coarse(*quad_mesh_coarse);
      RootBndFactory bnd_factory_fine(quad_mesh_fine);

      // create boundary cellsets
      SubSet bnd_cells_coarse(bnd_factory_coarse);
      SubSet bnd_cells_fine(bnd_factory_fine);

      // validate refined meshes
      try
      {
        TestAux::validate_tetris_quad_boundary_cellsubset_2d(bnd_cells_coarse);
        TestAux::validate_refined_tetris_quad_boundary_cellsubset_2d(bnd_cells_fine);
      }
      catch(const String& msg)
      {
        TEST_CHECK_MSG(false, msg);
      }
    }

    // clean up
    delete quad_mesh_coarse;
  }

} boundary_factory_test;
