// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>
#include <kernel/geometry/structured_mesh.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

typedef StructuredMesh<2> StructMesh;

typedef StandardRefinery<StructMesh> StructMeshRefinery;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 *
 * \author Peter Zajac
 */
class StandardRefineryTestStructQuad
  : public UnitTest
{
public:
  StandardRefineryTestStructQuad() :
    UnitTest("standard_refinery-test-struct-quad")
  {
  }

  virtual ~StandardRefineryTestStructQuad()
  {
  }

  virtual void run() const override
  {
    // allocate a 2x3 structured mesh
    Index num_slices[] =
    {
      2, 3
    };
    StructMesh* mesh_coarse = new StructMesh(num_slices);

    // validate entity counts; a 2x3 mesh has 12 vertices, 17 edges and 6 quads
    TEST_CHECK_MSG(mesh_coarse->get_num_entities(0) == 12, "invalid coarse mesh vertex count");
    TEST_CHECK_MSG(mesh_coarse->get_num_entities(1) == 17, "invalid coarse mesh edge count");
    TEST_CHECK_MSG(mesh_coarse->get_num_entities(2) == 6, "invalid coarse mesh quad count");

    // set vertex array
    static Real vtx[] =
    {
      0.0, 0.0, 1.0, 0.0, 2.0, 0.0,
      0.0, 1.0, 1.0, 1.0, 2.0, 1.0,
      0.0, 2.0, 1.0, 2.0, 2.0, 2.0,
      0.0, 3.0, 1.0, 3.0, 2.0, 3.0
    };
    TestAux::copy_vtx(mesh_coarse->get_vertex_set(), vtx);

    // create a refinery
    StructMeshRefinery* refinery = new StructMeshRefinery(*mesh_coarse);

    // refine the mesh
    StructMesh* mesh_fine = new StructMesh(*refinery);

    // validate slice counts
    TEST_CHECK_MSG(mesh_fine->get_num_slices(0) == 4, "invalid fine mesh X-slice count");
    TEST_CHECK_MSG(mesh_fine->get_num_slices(1) == 6, "invalid fine mesh Y-slice count");

    // validate entity counts
    TEST_CHECK_MSG(mesh_fine->get_num_entities(0) == 35, "invalid fine mesh vertex count");
    TEST_CHECK_MSG(mesh_fine->get_num_entities(1) == 58, "invalid fine mesh edge count");
    TEST_CHECK_MSG(mesh_fine->get_num_entities(2) == 24, "invalid fine mesh quad count");

    // validate vertex coordinates
    static Real wtx[] =
    {
      0.0, 0.0, 0.5, 0.0, 1.0, 0.0, 1.5, 0.0, 2.0, 0.0,
      0.0, 0.5, 0.5, 0.5, 1.0, 0.5, 1.5, 0.5, 2.0, 0.5,
      0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.5, 1.0, 2.0, 1.0,
      0.0, 1.5, 0.5, 1.5, 1.0, 1.5, 1.5, 1.5, 2.0, 1.5,
      0.0, 2.0, 0.5, 2.0, 1.0, 2.0, 1.5, 2.0, 2.0, 2.0,
      0.0, 2.5, 0.5, 2.5, 1.0, 2.5, 1.5, 2.5, 2.0, 2.5,
      0.0, 3.0, 0.5, 3.0, 1.0, 3.0, 1.5, 3.0, 2.0, 3.0
    };
    TEST_CHECK_MSG(TestAux::comp_vtx(mesh_fine->get_vertex_set(), wtx), "vertex refinement failure");

    // clean up
    delete mesh_fine;
    delete refinery;
    delete mesh_coarse;
  }
} standard_refinery_test_struct_quad;
