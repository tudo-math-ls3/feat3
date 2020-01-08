// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>
#include <kernel/geometry/structured_mesh.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

typedef StructuredMesh<3> StructMesh3d;
typedef StandardRefinery<StructMesh3d> StructMeshRefinery3d;

/**
 * \brief Test class for the StandardRefinery class template.
 *
 * \test Tests the StandardRefinery class templates.
 */
class StandardRefineryTestStructHexa
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  StandardRefineryTestStructHexa() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("standard_refinery-test-struct-hexa")
  {
  }

  virtual ~StandardRefineryTestStructHexa()
  {
  }

  virtual void run() const override
  {

    // allocate a 3x1x2 structured hexa-mesh
    Index num_slices_3d[] =
    {
      3, 1, 2
    };
    StructMesh3d mesh_coarse(num_slices_3d);

    // validate entity counts; a 3x1x2 mesh has 24 vertices, 46 edges, 29 quads and 6 hexas
    TEST_CHECK_MSG(mesh_coarse.get_num_entities(0) == 24, "invalid coarse mesh vertex count");
    TEST_CHECK_MSG(mesh_coarse.get_num_entities(1) == 46, "invalid coarse mesh edge count");
    TEST_CHECK_MSG(mesh_coarse.get_num_entities(2) == 29, "invalid coarse mesh quad count");
    TEST_CHECK_MSG(mesh_coarse.get_num_entities(3) == 6, "invalid coarse mesh hexa count");

    // set vertex array
    static Real vtx[] =
    {
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 4.0, 1.0, 0.0,
      0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 4.0, 0.0, 1.0,
      0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 4.0, 1.0, 1.0,
      0.0, 0.0, 3.0, 1.0, 0.0, 3.0, 2.0, 0.0, 3.0, 4.0, 0.0, 3.0,
      0.0, 1.0, 3.0, 1.0, 1.0, 3.0, 2.0, 1.0, 3.0, 4.0, 1.0, 3.0,

    };
    TestAux::copy_vtx(mesh_coarse.get_vertex_set(), vtx);

    // create a refinery
    StructMeshRefinery3d refinery(mesh_coarse);

    // refine the mesh
    StructMesh3d mesh_fine(refinery);

    // validate slice counts
    TEST_CHECK_MSG(mesh_fine.get_num_slices(0) == 6, "invalid fine mesh X-slice count");
    TEST_CHECK_MSG(mesh_fine.get_num_slices(1) == 2, "invalid fine mesh Y-slice count");
    TEST_CHECK_MSG(mesh_fine.get_num_slices(2) == 4, "invalid fine mesh Z-slice count");

    // validate entity counts
    TEST_CHECK_MSG(mesh_fine.get_num_entities(0) == 105, "invalid fine mesh vertex count");
    TEST_CHECK_MSG(mesh_fine.get_num_entities(1) == 244, "invalid fine mesh edge count");
    TEST_CHECK_MSG(mesh_fine.get_num_entities(2) == 188, "invalid fine mesh quad count");
    TEST_CHECK_MSG(mesh_fine.get_num_entities(3) == 48, "invalid fine mesh hexa count");

    // validate vertex coordinates
    static Real wtx[] =
    {
      0.0, 0.0, 0.0, 0.5, 0.0, 0.0,
      1.0, 0.0, 0.0, 1.5, 0.0, 0.0,
      2.0, 0.0, 0.0, 3.0, 0.0, 0.0,
      4.0, 0.0, 0.0,
      0.0, 0.5, 0.0, 0.5, 0.5, 0.0,
      1.0, 0.5, 0.0, 1.5, 0.5, 0.0,
      2.0, 0.5, 0.0, 3.0, 0.5, 0.0,
      4.0, 0.5, 0.0,
      0.0, 1.0, 0.0, 0.5, 1.0, 0.0,
      1.0, 1.0, 0.0, 1.5, 1.0, 0.0,
      2.0, 1.0, 0.0, 3.0, 1.0, 0.0,
      4.0, 1.0, 0.0,

      0.0, 0.0, 0.5, 0.5, 0.0, 0.5,
      1.0, 0.0, 0.5, 1.5, 0.0, 0.5,
      2.0, 0.0, 0.5, 3.0, 0.0, 0.5,
      4.0, 0.0, 0.5,
      0.0, 0.5, 0.5, 0.5, 0.5, 0.5,
      1.0, 0.5, 0.5, 1.5, 0.5, 0.5,
      2.0, 0.5, 0.5, 3.0, 0.5, 0.5,
      4.0, 0.5, 0.5,
      0.0, 1.0, 0.5, 0.5, 1.0, 0.5,
      1.0, 1.0, 0.5, 1.5, 1.0, 0.5,
      2.0, 1.0, 0.5, 3.0, 1.0, 0.5,
      4.0, 1.0, 0.5,

      0.0, 0.0, 1.0, 0.5, 0.0, 1.0,
      1.0, 0.0, 1.0, 1.5, 0.0, 1.0,
      2.0, 0.0, 1.0, 3.0, 0.0, 1.0,
      4.0, 0.0, 1.0,
      0.0, 0.5, 1.0, 0.5, 0.5, 1.0,
      1.0, 0.5, 1.0, 1.5, 0.5, 1.0,
      2.0, 0.5, 1.0, 3.0, 0.5, 1.0,
      4.0, 0.5, 1.0,
      0.0, 1.0, 1.0, 0.5, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.5, 1.0, 1.0,
      2.0, 1.0, 1.0, 3.0, 1.0, 1.0,
      4.0, 1.0, 1.0,

      0.0, 0.0, 2.0, 0.5, 0.0, 2.0,
      1.0, 0.0, 2.0, 1.5, 0.0, 2.0,
      2.0, 0.0, 2.0, 3.0, 0.0, 2.0,
      4.0, 0.0, 2.0,
      0.0, 0.5, 2.0, 0.5, 0.5, 2.0,
      1.0, 0.5, 2.0, 1.5, 0.5, 2.0,
      2.0, 0.5, 2.0, 3.0, 0.5, 2.0,
      4.0, 0.5, 2.0,
      0.0, 1.0, 2.0, 0.5, 1.0, 2.0,
      1.0, 1.0, 2.0, 1.5, 1.0, 2.0,
      2.0, 1.0, 2.0, 3.0, 1.0, 2.0,
      4.0, 1.0, 2.0,

      0.0, 0.0, 3.0, 0.5, 0.0, 3.0,
      1.0, 0.0, 3.0, 1.5, 0.0, 3.0,
      2.0, 0.0, 3.0, 3.0, 0.0, 3.0,
      4.0, 0.0, 3.0,
      0.0, 0.5, 3.0, 0.5, 0.5, 3.0,
      1.0, 0.5, 3.0, 1.5, 0.5, 3.0,
      2.0, 0.5, 3.0, 3.0, 0.5, 3.0,
      4.0, 0.5, 3.0,
      0.0, 1.0, 3.0, 0.5, 1.0, 3.0,
      1.0, 1.0, 3.0, 1.5, 1.0, 3.0,
      2.0, 1.0, 3.0, 3.0, 1.0, 3.0,
      4.0, 1.0, 3.0
    };
    TEST_CHECK_MSG(TestAux::comp_vtx(mesh_fine.get_vertex_set(), wtx), "vertex refinement failure");
  }
} standard_refinery_test_struct_hexa;
