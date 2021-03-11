// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/parsed_hit_test_factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>         // for RefinedUnitCubeFactor
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart
#include <kernel/util/string.hpp>                          // for String

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;
#ifdef FEAT_HAVE_FPARSER
/**
 * \brief Test class for the ParsedHitTestFactoryTest class.
 *
 * \test Tests the ParsedHitTestFactoryTest class.
 *
 * \author Gesa Pottbrock
 */

class ParsedHitTestFactoryTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  ParsedHitTestFactoryTest() :
    TaggedTest<Archs::None, Archs::None>("parsed_hit_test_factory-test")
  {
  }

  virtual ~ParsedHitTestFactoryTest()
  {
  }

  void test_0() const
  {
    typedef Shape::Quadrilateral ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::MeshPart<MeshType> BoundaryType;
    typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactoryType;
    typedef Geometry::StandardRefinery<MeshType> RefineryType;

    // some parameters for this special test
    Index level(6);

    // create a meshfactory
    MeshFactoryType mesh_factory(level);
    // create a mesh and a refined mesh
    MeshType mesh (mesh_factory);
    RefineryType refinery(mesh);
    MeshType mesh_fine(refinery);

    // create two ParsedHitTestFactories
    ParsedHitTestFactory< MeshType,int(2)> hit_test( mesh, "0.3-sqrt((x-0.5)^2+(y-0.5)^2)");
    BoundaryType sphere(hit_test);
    ParsedHitTestFactory< MeshType,int(2)> hit_test_fine( mesh_fine, "0.3-sqrt((x-0.5)^2+(y-0.5)^2)");
    BoundaryType sphere_fine(hit_test_fine);

    // number of nodes, edges, quads
    Index nv(sphere.get_num_entities(0)), ne(sphere.get_num_entities(1)), nq(sphere.get_num_entities(2));
    // check if the number of nodes fits
    TEST_CHECK_EQUAL(nv + ne + nq, sphere_fine.get_num_entities(0));
    // helper vector to compare the nodes caught by the hit test of the mesh and refined mesh
    std::vector<Index> nodes(nv + ne + nq);

    // compute the nodes that should be caught by the hit test of the refined mesh
    for (Index i(0); i < sphere.get_num_entities(0); ++i)
    {
      nodes[i] = sphere.get_target_set<0>()[i];
    }
    for (Index i(0); i < sphere.get_num_entities(1); ++i)
    {
      nodes[i + nv] = sphere.get_target_set<1>()[i] + mesh.get_num_entities(0);
    }
    for (Index i(0); i < sphere.get_num_entities(2); ++i)
    {
      nodes[i + nv + ne] = sphere.get_target_set<2>()[i] + mesh.get_num_entities(0) + mesh.get_num_entities(1);
    }
    // compare the node numbers
    for (Index i(0); i < sphere_fine.get_num_entities(0); ++i)
    {
      TEST_CHECK_EQUAL(nodes[i], sphere_fine.get_target_set<0>()[i]);
    }
  } // test_0

  virtual void run() const override
  {
    // run test #0
    test_0();
  }
} hit_test_factory_test;

#endif // defined(FEAT_HAVE_FPARSER)
