// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/standard_tetra.hpp>
#include <kernel/geometry/conformal_mesh.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

typedef ConformalMesh<Shape::Tetrahedron> RootMesh;
typedef MeshPart<RootMesh> SubMesh;

typedef StandardRefinery<RootMesh> RootMeshRefinery;
typedef StandardRefinery<SubMesh> SubMeshRefinery;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Constantin Christof
 */

class StandardRefineryTestConfTetrahedron
  : public UnitTest
{
public:
  StandardRefineryTestConfTetrahedron() :
    UnitTest("standard_refinery-test-conf-tetrahedron")
  {
  }

  virtual ~StandardRefineryTestConfTetrahedron()
  {
  }

  virtual void run() const override
  {
    tetra_std_test();
    big_tetra_mesh_test();
  }

  void tetra_std_test() const
  {
    // loop over all possible orientations
    for(int i(0); i < 3; ++i)
    {
      // create a 3D tetrahedron element mesh
      RootMesh* tetrahedron_mesh_coarse = TestAux::create_tetra_mesh_3d(i);

      try
      {
        // create refineries
        RootMeshRefinery tetrahedron_mesh_refinery(*tetrahedron_mesh_coarse);

        // refine the meshes
        RootMesh tetrahedron_mesh_fine(tetrahedron_mesh_refinery);

        // validate refined meshes
        TestAux::validate_refined_tetra_mesh_3d(tetrahedron_mesh_fine, i);
      }
      catch(const String& msg)
      {
        TEST_CHECK_MSG(false, msg);
      }

      // clean up
      delete tetrahedron_mesh_coarse;
    }
  }

  void big_tetra_mesh_test() const
  {
    // create a 3D tetrahedron element mesh
    RootMesh* tetrahedron_mesh_coarse = TestAux::create_big_tetra_mesh_3d();

    // create tria submesh
    SubMesh* tria_submesh_coarse = TestAux::create_tria_submesh_3d();

    try
    {
      // create refineries
      RootMeshRefinery tetrahedron_mesh_refinery(*tetrahedron_mesh_coarse);
      SubMeshRefinery tria_submesh_refinery(*tria_submesh_coarse, *tetrahedron_mesh_coarse);

      // refine the meshes
      RootMesh tetrahedron_mesh_fine(tetrahedron_mesh_refinery);
      SubMesh tria_submesh_fine(tria_submesh_refinery);

      // validate refined meshes
      TestAux::validate_refined_big_tetra_mesh_3d(tetrahedron_mesh_fine);
      TestAux::validate_refined_tria_submesh_3d(tria_submesh_fine);
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    // clean up
    delete tetrahedron_mesh_coarse;
    delete tria_submesh_coarse;
  }

} standard_refinery_test_conf_tetrahedron;
