// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/tetris_quad.hpp>
#include <kernel/geometry/test_aux/standard_quad.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

typedef ConformalMesh<Shape::Quadrilateral> RootMesh;
typedef MeshPart<RootMesh> SubMesh;

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
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  StandardRefineryTestConfQuad() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("standard_refinery-test-conf-quad")
  {
  }

  virtual ~StandardRefineryTestConfQuad()
  {
  }

  virtual void run() const override
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
    SubMesh* cell_subset_coarse = TestAux::create_tetris_quad_cellsubset_2d();

    // create refineries
    RootMeshRefinery quad_mesh_refinery(*quad_mesh_coarse);
    SubMeshRefinery edge_submesh_refinery(*edge_submesh_coarse, *quad_mesh_coarse);
    SubMeshRefinery quad_submesh_refinery(*quad_submesh_coarse, *quad_mesh_coarse);
    SubMeshRefinery cell_subset_refinery(*cell_subset_coarse, *quad_mesh_coarse);

    // refine the meshes
    RootMesh quad_mesh_fine(quad_mesh_refinery);
    SubMesh edge_submesh_fine(edge_submesh_refinery);
    SubMesh quad_submesh_fine(quad_submesh_refinery);
    SubMesh cell_subset_fine(cell_subset_refinery);

    // validate refined meshes
    try
    {
      TestAux::validate_refined_tetris_mesh_2d(quad_mesh_fine);
      TestAux::validate_refined_tetris_edge_submesh_2d(edge_submesh_fine);
      TestAux::validate_refined_tetris_quad_submesh_2d(quad_submesh_fine);
      TestAux::validate_refined_tetris_quad_cellsubset_2d(cell_subset_fine);
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    // clean up
    delete cell_subset_coarse;
    delete quad_submesh_coarse;
    delete edge_submesh_coarse;
    delete quad_mesh_coarse;
  }

  void quad_std_test() const
  {
    // loop over all possible orientations (max. 8)
    for(int i(0); i < 4; ++i)
    {
      // create a 2D quad element mesh
      RootMesh* quad_mesh_coarse = TestAux::create_quad_mesh_2d(i);

      try
      {
        // create refinery
        RootMeshRefinery quad_mesh_refinery(*quad_mesh_coarse);

        // refine the mesh
        RootMesh quad_mesh_fine(quad_mesh_refinery);

        // validate refined mesh
        TestAux::validate_refined_quad_mesh_2d(quad_mesh_fine,i);
      }
      catch(const String& msg)
      {
        TEST_CHECK_MSG(false, msg);
      }

      // clean up
      delete quad_mesh_coarse;
    }
  }

} standard_refinery_test_conf_quad;
