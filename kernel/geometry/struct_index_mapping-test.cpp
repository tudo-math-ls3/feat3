// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/validate_structured_meshes.hpp>
#include <kernel/geometry/structured_mesh.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

/**
 * \brief Test class for the struct_index_mapping formulas
 *
 * \test Tests the struct_index_mapping formulas.
 *
 * \author Constantin Christof
 */

class IndexMappingTest
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  IndexMappingTest() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("index_mapping-test")
  {
  }

  // test
  virtual void run() const override
  {
    // 2d mesh
    Index num_slices2d[2] =
    {
      6, 2
    };

    // 3d mesh
    Index num_slices3d[3] =
    {
      4, 2, 3
    };

    // create meshes
    StructuredMesh<2>* mesh2d = new StructuredMesh<2>(num_slices2d);
    StructuredMesh<3>* mesh3d = new StructuredMesh<3>(num_slices3d);

    // validate
    try
    {
      TestAux::validate_structured_mesh_2d(mesh2d);
      TestAux::validate_structured_mesh_3d(mesh3d);
    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    delete mesh2d;
    delete mesh3d;

 }

} index_mapping_test;
