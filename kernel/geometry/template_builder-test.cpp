// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/raw_refinement_templates.hpp>
#include <kernel/geometry/refinement_types.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/geometry/template_builder.hpp>
#include <kernel/geometry/template_sets.hpp>
#include <kernel/geometry/templates/schneiders_data.hpp>
#include <kernel/geometry/templates/sun_zhao_ma_data.hpp>
#include <kernel/geometry/templates/sun_zhao_ma_expansion_data.hpp>
#include <kernel/geometry/templates/two_refinement_data.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/string.hpp>
#include <test_system/test_system.hpp>

#include <kernel/geometry/adaptive_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;
using namespace FEAT::Geometry::Intern;

class TemplateBuilderTest : public UnitTest
{
public:
  TemplateBuilderTest() : UnitTest("TemplateBuilderTest")
  {
  }

  ~TemplateBuilderTest() override = default;

  using RawEdge = RawEntity<Shape::Hypercube<1>>;
  using RawFace = RawEntity<Shape::Hypercube<2>>;
  using RawCell = RawEntity<Shape::Hypercube<3>>;

  using EdgeType = StandardRefinementType<Shape::Hypercube<1>>;
  using FaceType = StandardRefinementType<Shape::Hypercube<2>>;
  using CellType = StandardRefinementType<Shape::Hypercube<3>>;

  void run() const override
  {
    test_utils();
    test_search_space();
    test_template_construction();

    test_orientation_mapping<2, 2>();
    test_orientation_mapping<2, 1>();
    test_orientation_mapping<2, 0>();
    test_orientation_mapping<1, 1>();
    test_orientation_mapping<1, 0>();

    // Check raw dimensions are correct
    TEST_CHECK_EQUAL(RawEdge::num_coords, 1);
    TEST_CHECK_EQUAL(RawEdge::num_vertices, 2);

    TEST_CHECK_EQUAL(RawFace::num_coords, 2);
    TEST_CHECK_EQUAL(RawFace::num_vertices, 4);

    TEST_CHECK_EQUAL(RawCell::num_coords, 3);
    TEST_CHECK_EQUAL(RawCell::num_vertices, 8);
  }

  void test_utils() const
  {
    // Check is_internal_vertex
    TEST_CHECK(is_internal_vertex(Tiny::Vector<Real, 3>{0.25, 0.5, 0.75}));
    TEST_CHECK(!is_internal_vertex(Tiny::Vector<Real, 3>{0.0, 0.5, 0.75}));

    // Test is_same_raw_entity
    {
      RawEntity<Shape::Hypercube<2>> a(
        Tiny::Vector<Real, 2>{0.0, 0.0},
        Tiny::Vector<Real, 2>{1.0, 0.0},
        Tiny::Vector<Real, 2>{0.0, 1.0},
        Tiny::Vector<Real, 2>{1.0, 1.0});

      RawEntity<Shape::Hypercube<2>> b(
        Tiny::Vector<Real, 2>{0.5, 0.5},
        Tiny::Vector<Real, 2>{1.0, 0.0},
        Tiny::Vector<Real, 2>{0.0, 1.0},
        Tiny::Vector<Real, 2>{1.0, 1.0});

      TEST_CHECK(is_same_raw_entity(a, a));
      TEST_CHECK(!is_same_raw_entity(a, b));

      RawEntity<Shape::Hypercube<1>, 3> edge_a(
        Tiny::Vector<Real, 3>{0.5, 0.5, 0.0},
        Tiny::Vector<Real, 3>{0.0, 0.5, 0.0});
      RawEntity<Shape::Hypercube<1>, 3> edge_b(
        Tiny::Vector<Real, 3>{0.0, 0.5, 0.0},
        Tiny::Vector<Real, 3>{0.0, 0.5, 0.0});

      TEST_CHECK(!is_same_raw_entity(edge_a, edge_b));
      TEST_CHECK(!is_same_raw_entity(edge_b, edge_a));
    }

    // Test bitset rotations
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_2d(std::bitset<4>(0b0001)), std::bitset<4>(0b0010));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_2d(std::bitset<4>(0b0010)), std::bitset<4>(0b1000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_2d(std::bitset<4>(0b1000)), std::bitset<4>(0b0100));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_2d(std::bitset<4>(0b0100)), std::bitset<4>(0b0001));

    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_xaxis(std::bitset<8>(0b00000001)), std::bitset<8>(0b00000100));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_xaxis(std::bitset<8>(0b00000100)), std::bitset<8>(0b01000000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_xaxis(std::bitset<8>(0b01000000)), std::bitset<8>(0b00010000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_xaxis(std::bitset<8>(0b00010000)), std::bitset<8>(0b00000001));

    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_xaxis(std::bitset<8>(0b00000010)), std::bitset<8>(0b00001000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_xaxis(std::bitset<8>(0b00001000)), std::bitset<8>(0b10000000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_xaxis(std::bitset<8>(0b10000000)), std::bitset<8>(0b00100000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_xaxis(std::bitset<8>(0b00100000)), std::bitset<8>(0b00000010));

    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_yaxis(std::bitset<8>(0b00000001)), std::bitset<8>(0b00010000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_yaxis(std::bitset<8>(0b00010000)), std::bitset<8>(0b00100000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_yaxis(std::bitset<8>(0b00100000)), std::bitset<8>(0b00000010));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_yaxis(std::bitset<8>(0b00000010)), std::bitset<8>(0b00000001));

    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_yaxis(std::bitset<8>(0b00000100)), std::bitset<8>(0b01000000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_yaxis(std::bitset<8>(0b01000000)), std::bitset<8>(0b10000000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_yaxis(std::bitset<8>(0b10000000)), std::bitset<8>(0b00001000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_yaxis(std::bitset<8>(0b00001000)), std::bitset<8>(0b00000100));

    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_zaxis(std::bitset<8>(0b00000001)), std::bitset<8>(0b00000010));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_zaxis(std::bitset<8>(0b00000010)), std::bitset<8>(0b00001000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_zaxis(std::bitset<8>(0b00001000)), std::bitset<8>(0b00000100));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_zaxis(std::bitset<8>(0b00000100)), std::bitset<8>(0b00000001));

    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_zaxis(std::bitset<8>(0b00010000)), std::bitset<8>(0b00100000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_zaxis(std::bitset<8>(0b00100000)), std::bitset<8>(0b10000000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_zaxis(std::bitset<8>(0b10000000)), std::bitset<8>(0b01000000));
    TEST_CHECK_EQUAL(Geometry::Intern::rotate_arraylike_zaxis(std::bitset<8>(0b01000000)), std::bitset<8>(0b00010000));

    // Test 2d template rotations
    {
      using V = typename RawTemplate<Shape::Hypercube<2>>::VertexType;

      RawTemplate<Shape::Hypercube<2>> tmplt;
      tmplt.add_entity(V{0.0, 0.0}, V{0.5, 0.0}, V{0.0, 0.5}, V{0.5, 0.5});

      auto rot_tmplt = rotate_template_2d(tmplt);
      auto rot_entity = rot_tmplt.entities[0];

      TEST_CHECK(is_same_raw_entity(
        rot_entity,
        RawEntity<Shape::Hypercube<2>>(V{0.5, 0.0}, V{1.0, 0.0}, V{0.5, 0.5}, V{1.0, 0.5})));
    }

    // Test 3d template rotations
    {
      using V = typename RawTemplate<Shape::Hypercube<3>>::VertexType;

      RawTemplate<Shape::Hypercube<3>> tmplt;
      tmplt.add_entity(
        V{0.0, 0.0, 0.0},
        V{0.5, 0.0, 0.0},
        V{0.0, 0.5, 0.0},
        V{0.5, 0.5, 0.0},
        V{0.0, 0.0, 0.5},
        V{0.5, 0.0, 0.5},
        V{0.0, 0.5, 0.5},
        V{0.5, 0.5, 0.5});

      auto rot_tmplt_x = rotate_template_xaxis(tmplt);
      auto rot_tmplt_y = rotate_template_yaxis(tmplt);
      auto rot_tmplt_z = rotate_template_zaxis(tmplt);

      auto rot_entity_x = rot_tmplt_x.entities[0];
      auto rot_entity_y = rot_tmplt_y.entities[0];
      auto rot_entity_z = rot_tmplt_z.entities[0];

      TEST_CHECK(is_same_raw_entity(
        rot_entity_x,
        RawEntity<Shape::Hypercube<3>>(
          V{0.0, 0.5, 0.0},
          V{0.5, 0.5, 0.0},
          V{0.0, 1.0, 0.0},
          V{0.5, 1.0, 0.0},
          V{0.0, 0.5, 0.5},
          V{0.5, 0.5, 0.5},
          V{0.0, 1.0, 0.5},
          V{0.5, 1.0, 0.5})));

      TEST_CHECK(is_same_raw_entity(
        rot_entity_y,
        RawEntity<Shape::Hypercube<3>>(
          V{0.0, 0.0, 0.5},
          V{0.5, 0.0, 0.5},
          V{0.0, 0.5, 0.5},
          V{0.5, 0.5, 0.5},
          V{0.0, 0.0, 1.0},
          V{0.5, 0.0, 1.0},
          V{0.0, 0.5, 1.0},
          V{0.5, 0.5, 1.0})));

      TEST_CHECK(is_same_raw_entity(
        rot_entity_z,
        RawEntity<Shape::Hypercube<3>>(
          V{0.5, 0.0, 0.0},
          V{1.0, 0.0, 0.0},
          V{0.5, 0.5, 0.0},
          V{1.0, 0.5, 0.0},
          V{0.5, 0.0, 0.5},
          V{1.0, 0.0, 0.5},
          V{0.5, 0.5, 0.5},
          V{1.0, 0.5, 0.5})));
    }

    // Test RawTemplate
    {
      using V = typename RawTemplate<Shape::Hexahedron>::VertexType;

      RawTemplate<Shape::Hexahedron> tmplt;

      RawEntity<Shape::Hexahedron> expected(
        V{0.0, 0.0, 0.0},
        V{0.5, 0.0, 0.0},
        V{0.0, 0.5, 0.0},
        V{0.5, 0.5, 0.0},
        V{0.0, 0.0, 0.5},
        V{0.5, 0.0, 0.5},
        V{0.0, 0.5, 0.5},
        V{0.5, 0.5, 0.5});

      tmplt
        .add_entity(
          V{0.0, 0.0, 0.0},
          V{0.5, 0.0, 0.0},
          V{0.0, 0.5, 0.0},
          V{0.5, 0.5, 0.0},
          V{0.0, 0.0, 0.5},
          V{0.5, 0.0, 0.5},
          V{0.0, 0.5, 0.5},
          V{0.5, 0.5, 0.5})
        .axis_aligned(V{0.0, 0.0, 0.0}, V{0.5, 0.5, 0.5})
        .grid({2, 2, 2}, V{0.5, 0.5, 0.5});

      // .add_entity entity should be the same as expected
      TEST_CHECK(is_same_raw_entity(tmplt.entities[0], expected));

      // .axis_aligned entity should be the same as expected
      TEST_CHECK(is_same_raw_entity(tmplt.entities[1], expected));

      // first .grid entity should be the same as expected
      TEST_CHECK(is_same_raw_entity(tmplt.entities[2], expected));

      RawTemplate<Shape::Hexahedron> unit_cube_tmplt;
      unit_cube_tmplt.axis_aligned(V{0.0, 0.0, 0.0}, V{1.0, 1.0, 1.0});

      unit_cube_tmplt.recurse(0, tmplt);

      // Recursively inserted tmplt into the unit cube, unit_cube_tmplt should now be identical to tmplt
      TEST_CHECK_EQUAL(tmplt.entities.size(), unit_cube_tmplt.entities.size());
      for(std::size_t i(0); i < tmplt.entities.size(); ++i)
      {
        TEST_CHECK(is_same_raw_entity(tmplt.entities[i], unit_cube_tmplt.entities[i]));
      }

      RawTemplate<Shape::Hexahedron> quarter;
      quarter.axis_aligned(V{0.0, 0.0, 0.0}, V{0.5, 0.5, 0.5});

      RawTemplate<Shape::Hexahedron> quarter_copy = quarter;

      quarter.recurse(0, quarter_copy);

      TEST_CHECK(is_same_raw_entity(
        quarter.entities[0],
        RawEntity<Shape::Hexahedron>(
          V{0.0, 0.0, 0.0},
          V{0.25, 0.0, 0.0},
          V{0.0, 0.25, 0.0},
          V{0.25, 0.25, 0.0},
          V{0.0, 0.0, 0.25},
          V{0.25, 0.0, 0.25},
          V{0.0, 0.25, 0.25},
          V{0.25, 0.25, 0.25})));
    }

    // Test is_internal functions
    {
      // Check topology vertex
      TEST_CHECK(!is_internal_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 0.0}));
      // Check vertex on edge
      TEST_CHECK(!is_internal_vertex(Tiny::Vector<Real, 3>{0.5, 0.0, 0.0}));
      // Check vertex on face
      TEST_CHECK(!is_internal_vertex(Tiny::Vector<Real, 3>{0.5, 0.5, 0.0}));
      // Check internal vertex
      TEST_CHECK(is_internal_vertex(Tiny::Vector<Real, 3>{0.5, 0.5, 0.5}));

      using Entity = RawEntity<Shape::Hypercube<1>, 3>;
      using V = typename Entity::VertexType;

      // Outer edge shouldn't be internal
      TEST_CHECK(!is_internal_entity(Entity(V{0.0, 0.0, 0.0}, V{1.0, 0.0, 0.0})));
      // Edge with all internal vertices should be internal
      TEST_CHECK(is_internal_entity(Entity(V{0.25, 0.25, 0.25}, V{0.75, 0.75, 0.75})));
      // Edge with all internal vertices, but crossing the interior, should be internal
      TEST_CHECK(is_internal_entity(Entity(V{0.0, 0.0, 0.0}, V{1.0, 1.0, 1.0})));
    }

    // Test vertex projections
    {
      Tiny::Vector<Real, 1> left{0.0};
      Tiny::Vector<Real, 1> center{0.5};
      Tiny::Vector<Real, 1> right{1.0};

      using Embedder = VertexEmbedder<Shape::Quadrilateral, 1>;

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{0.0, 0.0}, Embedder::embed(0, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{0.5, 0.0}, Embedder::embed(0, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{1.0, 0.0}, Embedder::embed(0, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{0.0, 1.0}, Embedder::embed(1, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{0.5, 1.0}, Embedder::embed(1, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{1.0, 1.0}, Embedder::embed(1, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{0.0, 0.0}, Embedder::embed(2, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{0.0, 0.5}, Embedder::embed(2, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{0.0, 1.0}, Embedder::embed(2, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{1.0, 0.0}, Embedder::embed(3, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{1.0, 0.5}, Embedder::embed(3, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 2>{1.0, 1.0}, Embedder::embed(3, right)));
    }

    {
      Tiny::Vector<Real, 1> left{0.0};
      Tiny::Vector<Real, 1> center{0.5};
      Tiny::Vector<Real, 1> right{1.0};

      using Embedder = VertexEmbedder<Shape::Hexahedron, 1>;

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 0.0}, Embedder::embed(0, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.5, 0.0, 0.0}, Embedder::embed(0, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 0.0}, Embedder::embed(0, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 0.0}, Embedder::embed(1, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.5, 1.0, 0.0}, Embedder::embed(1, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 0.0}, Embedder::embed(1, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 1.0}, Embedder::embed(2, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.5, 0.0, 1.0}, Embedder::embed(2, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 1.0}, Embedder::embed(2, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 1.0}, Embedder::embed(3, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.5, 1.0, 1.0}, Embedder::embed(3, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 1.0}, Embedder::embed(3, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 0.0}, Embedder::embed(4, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.5, 0.0}, Embedder::embed(4, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 0.0}, Embedder::embed(4, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 0.0}, Embedder::embed(5, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.5, 0.0}, Embedder::embed(5, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 0.0}, Embedder::embed(5, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 1.0}, Embedder::embed(6, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.5, 1.0}, Embedder::embed(6, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 1.0}, Embedder::embed(6, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 1.0}, Embedder::embed(7, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.5, 1.0}, Embedder::embed(7, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 1.0}, Embedder::embed(7, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 0.0}, Embedder::embed(8, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 0.5}, Embedder::embed(8, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 1.0}, Embedder::embed(8, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 0.0}, Embedder::embed(9, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 0.5}, Embedder::embed(9, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 1.0}, Embedder::embed(9, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 0.0}, Embedder::embed(10, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 0.5}, Embedder::embed(10, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 1.0}, Embedder::embed(10, right)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 0.0}, Embedder::embed(11, left)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 0.5}, Embedder::embed(11, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 1.0}, Embedder::embed(11, right)));
    }

    {
      Tiny::Vector<Real, 2> bl{0.0, 0.0};
      Tiny::Vector<Real, 2> br{1.0, 0.0};
      Tiny::Vector<Real, 2> center{0.5, 0.5};
      Tiny::Vector<Real, 2> tl{0.0, 1.0};
      Tiny::Vector<Real, 2> tr{1.0, 1.0};

      using Embedder = VertexEmbedder<Shape::Hexahedron, 2>;

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 0.0}, Embedder::embed(0, bl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 0.0}, Embedder::embed(0, br)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.5, 0.5, 0.0}, Embedder::embed(0, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 0.0}, Embedder::embed(0, tl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 0.0}, Embedder::embed(0, tr)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 1.0}, Embedder::embed(1, bl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 1.0}, Embedder::embed(1, br)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.5, 0.5, 1.0}, Embedder::embed(1, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 1.0}, Embedder::embed(1, tl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 1.0}, Embedder::embed(1, tr)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 0.0}, Embedder::embed(2, bl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 0.0}, Embedder::embed(2, br)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.5, 0.0, 0.5}, Embedder::embed(2, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 1.0}, Embedder::embed(2, tl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 1.0}, Embedder::embed(2, tr)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 0.0}, Embedder::embed(3, bl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 0.0}, Embedder::embed(3, br)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.5, 1.0, 0.5}, Embedder::embed(3, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 1.0}, Embedder::embed(3, tl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 1.0}, Embedder::embed(3, tr)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 0.0}, Embedder::embed(4, bl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 0.0}, Embedder::embed(4, br)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.5, 0.5}, Embedder::embed(4, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 0.0, 1.0}, Embedder::embed(4, tl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{0.0, 1.0, 1.0}, Embedder::embed(4, tr)));

      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 0.0}, Embedder::embed(5, bl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 0.0}, Embedder::embed(5, br)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.5, 0.5}, Embedder::embed(5, center)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 0.0, 1.0}, Embedder::embed(5, tl)));
      TEST_CHECK(is_same_vertex(Tiny::Vector<Real, 3>{1.0, 1.0, 1.0}, Embedder::embed(5, tr)));
    }
  }

  void test_search_space() const
  {
    // Test topology
    {
      TemplateSearchSpace<Shape::Quadrilateral> space;

      space.add_topology();

      const std::array<Tiny::Vector<Real, 2>, 4> expected_vertices{
        Tiny::Vector<Real, 2>{0.0, 0.0},
        Tiny::Vector<Real, 2>{1.0, 0.0},
        Tiny::Vector<Real, 2>{0.0, 1.0},
        Tiny::Vector<Real, 2>{1.0, 1.0},
      };

      const auto& vertices = space.template entries<0>();
      const auto& edges = space.template entries<1>();

      TEST_CHECK_EQUAL(vertices.size(), 4);
      TEST_CHECK_EQUAL(edges.size(), 4);

      for(int i(0); i < 4; i++)
      {
        TEST_CHECK(is_same_vertex(vertices[i].raw_entity.coords[0], expected_vertices[i]));
        TEST_CHECK_EQUAL(vertices[i].reference.index, Index(i));
        TEST_CHECK_EQUAL(vertices[i].reference.source, EntitySource::ParentTopology);
      }

      for(int i(0); i < 4; i++)
      {
        using FaceMapping = Geometry::Intern::FaceIndexMapping<Shape::Quadrilateral, 1, 0>;
        TEST_CHECK(is_same_vertex(edges[i].raw_entity.coords[0], expected_vertices[FaceMapping::map(i, 0)]));
        TEST_CHECK(is_same_vertex(edges[i].raw_entity.coords[1], expected_vertices[FaceMapping::map(i, 1)]));

        TEST_CHECK_EQUAL(edges[i].reference.index, Index(i));
        TEST_CHECK_EQUAL(edges[i].reference.source, EntitySource::ParentTopology);
      }
    }

    {
      TemplateSearchSpace<Shape::Hexahedron> space;

      space.add_topology();

      const std::array<Tiny::Vector<Real, 3>, 8> expected_vertices{
        Tiny::Vector<Real, 3>{0.0, 0.0, 0.0},
        Tiny::Vector<Real, 3>{1.0, 0.0, 0.0},
        Tiny::Vector<Real, 3>{0.0, 1.0, 0.0},
        Tiny::Vector<Real, 3>{1.0, 1.0, 0.0},
        Tiny::Vector<Real, 3>{0.0, 0.0, 1.0},
        Tiny::Vector<Real, 3>{1.0, 0.0, 1.0},
        Tiny::Vector<Real, 3>{0.0, 1.0, 1.0},
        Tiny::Vector<Real, 3>{1.0, 1.0, 1.0},
      };

      const auto& vertices = space.template entries<0>();
      const auto& edges = space.template entries<1>();
      const auto& faces = space.template entries<2>();

      TEST_CHECK_EQUAL(vertices.size(), 8);
      TEST_CHECK_EQUAL(edges.size(), 12);
      TEST_CHECK_EQUAL(faces.size(), 6);

      for(int i(0); i < 8; i++)
      {
        TEST_CHECK(is_same_vertex(vertices[i].raw_entity.coords[0], expected_vertices[i]));
        TEST_CHECK_EQUAL(vertices[i].reference.index, Index(i));
        TEST_CHECK_EQUAL(vertices[i].reference.source, EntitySource::ParentTopology);
      }

      for(int i(0); i < 12; i++)
      {
        using FaceMapping = Geometry::Intern::FaceIndexMapping<Shape::Hexahedron, 1, 0>;
        TEST_CHECK(is_same_vertex(edges[i].raw_entity.coords[0], expected_vertices[FaceMapping::map(i, 0)]));
        TEST_CHECK(is_same_vertex(edges[i].raw_entity.coords[1], expected_vertices[FaceMapping::map(i, 1)]));

        TEST_CHECK_EQUAL(edges[i].reference.index, Index(i));
        TEST_CHECK_EQUAL(edges[i].reference.source, EntitySource::ParentTopology);
      }

      for(int i(0); i < 6; i++)
      {
        using FaceMapping = Geometry::Intern::FaceIndexMapping<Shape::Hexahedron, 2, 0>;

        TEST_CHECK(is_same_vertex(faces[i].raw_entity.coords[0], expected_vertices[FaceMapping::map(i, 0)]));
        TEST_CHECK(is_same_vertex(faces[i].raw_entity.coords[1], expected_vertices[FaceMapping::map(i, 1)]));
        TEST_CHECK(is_same_vertex(faces[i].raw_entity.coords[2], expected_vertices[FaceMapping::map(i, 2)]));
        TEST_CHECK(is_same_vertex(faces[i].raw_entity.coords[3], expected_vertices[FaceMapping::map(i, 3)]));

        TEST_CHECK_EQUAL(faces[i].reference.index, Index(i));
        TEST_CHECK_EQUAL(faces[i].reference.source, EntitySource::ParentTopology);
      }
    }

    // Test boundary elements
    {
      TemplateSearchSpace<Shape::Hexahedron> space;

      RawTemplate<Shape::Quadrilateral>& full_refinement =
        TwoRefinementData::raw_faces()[StandardRefinementType<Shape::Quadrilateral>(0b1111)];

      space.add_boundary_entity(0, full_refinement);

      const auto& vertices = space.template entries<0>();
      const auto& edges = space.template entries<1>();
      const auto& faces = space.template entries<2>();

      TEST_CHECK_EQUAL(vertices.size(), 1);
      TEST_CHECK_EQUAL(edges.size(), 4);
      TEST_CHECK_EQUAL(faces.size(), 4);

      TEST_CHECK(is_same_vertex(vertices[0].raw_entity.coords[0], Tiny::Vector<Real, 3>{0.5, 0.5, 0.0}));
      TEST_CHECK_EQUAL(vertices[0].reference.index, 0);
      TEST_CHECK_EQUAL(vertices[0].reference.entity, 0);
      TEST_CHECK_EQUAL(vertices[0].reference.source, EntitySource::BoundaryFace);

      std::vector<bool> seen_edges(4, false);
      {
        EntityReference entry = space.search(RawEntity<Shape::Hypercube<1>, 3>(
          Tiny::Vector<Real, 3>{0.5, 0.5, 0.0},
          Tiny::Vector<Real, 3>{0.0, 0.5, 0.0}));

        TEST_CHECK_EQUAL(entry.source, EntitySource::BoundaryFace);
        TEST_CHECK_EQUAL(entry.entity, 0);

        TEST_CHECK_EQUAL(seen_edges[entry.index], false);
        seen_edges[entry.index] = true;
      }

      {
        EntityReference entry = space.search(RawEntity<Shape::Hypercube<1>, 3>(
          Tiny::Vector<Real, 3>{0.5, 0.5, 0.0},
          Tiny::Vector<Real, 3>{1.0, 0.5, 0.0}));

        TEST_CHECK_EQUAL(entry.source, EntitySource::BoundaryFace);
        TEST_CHECK_EQUAL(entry.entity, 0);
        TEST_CHECK_EQUAL(seen_edges[entry.index], false);
        seen_edges[entry.index] = true;
      }

      {
        EntityReference entry = space.search(RawEntity<Shape::Hypercube<1>, 3>(
          Tiny::Vector<Real, 3>{0.5, 0.5, 0.0},
          Tiny::Vector<Real, 3>{0.5, 0.0, 0.0}));

        TEST_CHECK_EQUAL(entry.source, EntitySource::BoundaryFace);
        TEST_CHECK_EQUAL(entry.entity, 0);
        TEST_CHECK_EQUAL(seen_edges[entry.index], false);
        seen_edges[entry.index] = true;
      }

      {
        EntityReference entry = space.search(RawEntity<Shape::Hypercube<1>, 3>(
          Tiny::Vector<Real, 3>{0.5, 0.5, 0.0},
          Tiny::Vector<Real, 3>{0.5, 1.0, 0.0}));

        TEST_CHECK_EQUAL(entry.source, EntitySource::BoundaryFace);
        TEST_CHECK_EQUAL(entry.entity, 0);
        TEST_CHECK_EQUAL(seen_edges[entry.index], false);
        seen_edges[entry.index] = true;
      }

      std::vector<bool> seen_faces(4, false);
      {
        EntityReference entry = space.search(RawEntity<Shape::Hypercube<2>, 3>(
          Tiny::Vector<Real, 3>{0.0, 0.0, 0.0},
          Tiny::Vector<Real, 3>{0.5, 0.0, 0.0},
          Tiny::Vector<Real, 3>{0.0, 0.5, 0.0},
          Tiny::Vector<Real, 3>{0.5, 0.5, 0.0}));

        TEST_CHECK_EQUAL(entry.source, EntitySource::BoundaryFace);
        TEST_CHECK_EQUAL(entry.entity, 0);

        TEST_CHECK_EQUAL(seen_faces[entry.index], false);
        seen_edges[entry.index] = true;
      }

      {
        EntityReference entry = space.search(RawEntity<Shape::Hypercube<2>, 3>(
          Tiny::Vector<Real, 3>{0.5, 0.0, 0.0},
          Tiny::Vector<Real, 3>{1.0, 0.0, 0.0},
          Tiny::Vector<Real, 3>{0.5, 0.5, 0.0},
          Tiny::Vector<Real, 3>{1.0, 0.5, 0.0}));

        TEST_CHECK_EQUAL(entry.source, EntitySource::BoundaryFace);
        TEST_CHECK_EQUAL(entry.entity, 0);

        TEST_CHECK_EQUAL(seen_faces[entry.index], false);
        seen_edges[entry.index] = true;
      }

      {
        EntityReference entry = space.search(RawEntity<Shape::Hypercube<2>, 3>(
          Tiny::Vector<Real, 3>{0.5, 0.5, 0.0},
          Tiny::Vector<Real, 3>{1.0, 0.5, 0.0},
          Tiny::Vector<Real, 3>{0.5, 1.0, 0.0},
          Tiny::Vector<Real, 3>{1.0, 1.0, 0.0}));

        TEST_CHECK_EQUAL(entry.source, EntitySource::BoundaryFace);
        TEST_CHECK_EQUAL(entry.entity, 0);

        TEST_CHECK_EQUAL(seen_faces[entry.index], false);
        seen_edges[entry.index] = true;
      }

      {
        EntityReference entry = space.search(RawEntity<Shape::Hypercube<2>, 3>(
          Tiny::Vector<Real, 3>{0.0, 0.5, 0.0},
          Tiny::Vector<Real, 3>{0.5, 0.5, 0.0},
          Tiny::Vector<Real, 3>{0.0, 1.0, 0.0},
          Tiny::Vector<Real, 3>{0.5, 1.0, 0.0}));

        TEST_CHECK_EQUAL(entry.source, EntitySource::BoundaryFace);
        TEST_CHECK_EQUAL(entry.entity, 0);

        TEST_CHECK_EQUAL(seen_faces[entry.index], false);
        seen_edges[entry.index] = true;
      }
    }
  }

  void test_template_construction() const
  {
    TemplateBuilder<SchneidersData> builder;

    // Check edge templates
    TEST_CHECK_EQUAL((builder.template get_template<1>(EdgeType(0)).get_vertex_coefficients().size()), 0);
    TEST_CHECK_EQUAL((builder.template get_template<1>(EdgeType(1)).get_vertex_coefficients().size()), 1);
    TEST_CHECK_EQUAL((builder.template get_template<1>(EdgeType(2)).get_vertex_coefficients().size()), 1);
    TEST_CHECK_EQUAL((builder.template get_template<1>(EdgeType(3)).get_vertex_coefficients().size()), 2);

    const Real tol = Math::pow(Math::eps<Real>(), Real(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS((builder.template get_template<1>(EdgeType(1)).get_vertex_coefficients()[0][0]), 2.0 / 3.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((builder.template get_template<1>(EdgeType(1)).get_vertex_coefficients()[0][1]), 1.0 / 3.0, tol);

    TEST_CHECK_EQUAL_WITHIN_EPS((builder.template get_template<1>(EdgeType(2)).get_vertex_coefficients()[0][0]), 1.0 / 3.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((builder.template get_template<1>(EdgeType(2)).get_vertex_coefficients()[0][1]), 2.0 / 3.0, tol);

    TEST_CHECK_EQUAL_WITHIN_EPS((builder.template get_template<1>(EdgeType(3)).get_vertex_coefficients()[0][0]), 2.0 / 3.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((builder.template get_template<1>(EdgeType(3)).get_vertex_coefficients()[0][1]), 1.0 / 3.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((builder.template get_template<1>(EdgeType(3)).get_vertex_coefficients()[1][0]), 1.0 / 3.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((builder.template get_template<1>(EdgeType(3)).get_vertex_coefficients()[1][1]), 2.0 / 3.0, tol);

    TEST_CHECK_EQUAL(
      (builder.template get_template<1>(EdgeType(3)).get_topologies<1>()[0].get_references<0>()[0].source),
      EntitySource::ParentTopology);
    TEST_CHECK_EQUAL(
      (builder.template get_template<1>(EdgeType(3)).get_topologies<1>()[0].get_references<0>()[1].source),
      EntitySource::Sibling);

    // Check orientation mappings
    //// Check vertex mapping
    {
      TEST_CHECK_EQUAL((builder.template correct_for_orientation<1, 0>(EdgeType(3), 1, 0).first), 1);
      TEST_CHECK_EQUAL((builder.template correct_for_orientation<1, 0>(EdgeType(3), 1, 1).first), 0);
    }

    //// Check edge mapping
    {
      TEST_CHECK_EQUAL((builder.template correct_for_orientation<1, 1>(EdgeType(1), 1, 0).first), 1);
      TEST_CHECK_EQUAL((builder.template correct_for_orientation<1, 1>(EdgeType(1), 1, 1).first), 0);

      TEST_CHECK_EQUAL((builder.template correct_for_orientation<1, 1>(EdgeType(3), 1, 0).first), 2);
      TEST_CHECK_EQUAL((builder.template correct_for_orientation<1, 1>(EdgeType(3), 1, 1).first), 1);
      TEST_CHECK_EQUAL((builder.template correct_for_orientation<1, 1>(EdgeType(3), 1, 2).first), 0);
    }

    // Check face templates
    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(0)).get_vertex_coefficients().size(), 0);
    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(1)).get_vertex_coefficients().size(), 1);
    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(3)).get_vertex_coefficients().size(), 4);

    // Check type 1 face template
    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(1)).get_vertex_coefficients().size(), 1);
    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(1)).get_topologies<1>().size(), 3);
    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(1)).get_topologies<2>().size(), 3);

    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(0)).get_topologies<1>().size(), 0);
    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(1)).get_topologies<1>().size(), 3);
    TEST_CHECK_EQUAL(builder.template get_template<2>(FaceType(3)).get_topologies<1>().size(), 10);

    // Check cell templates
    TEST_CHECK_EQUAL(builder.template get_template<3>(CellType(0)).get_vertex_coefficients().size(), 0);
    TEST_CHECK_EQUAL(builder.template get_template<3>(CellType(1)).get_vertex_coefficients().size(), 1);
    TEST_CHECK_EQUAL(builder.template get_template<3>(CellType(2)).get_vertex_coefficients().size(), 1);
    TEST_CHECK_EQUAL(builder.template get_template<3>(CellType(3)).get_vertex_coefficients().size(), 4);
  }

  template<int dim_, int codim_>
  void test_orientation_mapping() const
  {
    TemplateBuilder<SchneidersData> builder;

    const int num_orientations = orientations<Shape::Hypercube<dim_>>();

    for(const auto& entry : SchneidersData::template raw_templates<dim_>())
    {
      auto type = entry.first;
      auto tmplt = builder.template get_template<dim_>(type);

      for(int orientation(0); orientation < num_orientations; ++orientation)
      {
        const Index size = tmplt.template num_entities<codim_>();

        // The mapping should be a bijection
        std::vector<int> seen(size, 0);

        // Check no element is mapped to twice
        for(Index i(0); i < size; ++i)
        {
          Index mapped = builder.template correct_for_orientation<dim_, codim_>(type, orientation, i).first;
          TEST_CHECK_MSG(
            seen[mapped] == 0,
            "test_orientation_mapping<" + stringify(dim_) + ", " + stringify(codim_) + ">: Element " +
              stringify(mapped) + " is mapped to twice for type " + stringify(type) + " and orientation " +
              stringify(orientation));
          seen[mapped] += 1;
        }

        // Check all elements are mapped to
        for(Index i(0); i < size; ++i)
        {
          TEST_CHECK_EQUAL(seen[i], 1);
        }
      }
    }
  }
};

static const TemplateBuilderTest template_builder_test;

template<typename TemplateSet_, typename Shape_>
class TemplateSetTest : public UnitTest
{
public:
  using TemplateSet = TemplateSet_;

  using CMeshType = ConformalMesh<Shape_>;

  using AMeshType = AdaptiveMesh<TemplateSet, Shape_>;

  using VertexType = typename AMeshType::VertexType;

  TemplateSetTest() : UnitTest("TemplateSetTest")
  {
  }

  ~TemplateSetTest() override = default;

  void run() const override
  {
    test_all_templates(1);
    test_all_templates(2);
    test_shared_elements();
  }

  void test_all_templates(Index level) const
  {
    Geometry::UnitCubeFactory<CMeshType> factory;
    CMeshType foundation(factory);

    Index max_markings = (1UL << (Index)Shape::FaceTraits<Shape_, 0>::count);
    for(std::uint64_t marking = 1; marking < max_markings; marking++)
    {
      AMeshType a_mesh(foundation);

      Geometry::SubdivisionLevels levels(foundation.get_num_vertices());
      for(Index i = 0; i < foundation.get_num_vertices(); ++i)
      {
        levels[i] = (marking & (1ULL << (std::uint64_t)(i))) > 0 ? level : 0;
      }
      a_mesh.adapt(levels);
      // a_mesh.fill_neighbors();

      CMeshType exported = a_mesh.to_conformal_mesh(Geometry::Layer{a_mesh.num_layers() - 1});

      TEST_CHECK(exported.get_num_vertices() > 0);
    }
  }

  void test_shared_elements() const
  {
    static constexpr int dim = Shape_::dimension;
    using SharedShape = typename Shape::FaceTraits<Shape_, dim - 1>::ShapeType;
    static constexpr int num_shared_verts = Shape::FaceTraits<SharedShape, 0>::count;

    Geometry::RefinedUnitCubeFactory<CMeshType> factory(1);
    CMeshType foundation(factory);

    // We are choosing face 1 of cell 0 as our shared face.
    // This is a shared element for both 2D and 3D meshes.
    Index facet = foundation.template get_index_set<dim, dim - 1>()(0, 1);
    IndexTuple<num_shared_verts> shared_facet = foundation.template get_index_set<dim - 1, 0>()[facet];

    Index num_markings = (1UL << (Index)num_shared_verts);

    for(Index marking(1); marking < num_markings; marking++)
    {
      Geometry::SubdivisionLevels sdls(foundation.get_num_vertices(), 0);

      for(int i(0); i < num_shared_verts; i++)
      {
        sdls[shared_facet[i]] = (marking & (1ULL << (Index)i)) > 0 ? 2 : 0;
      }

      AMeshType a_mesh(foundation);
      a_mesh.adapt(sdls, AMeshType::ImportBehaviour::All);

      CMeshType exported = a_mesh.to_conformal_mesh(Layer{2});

      TEST_CHECK(exported.get_num_vertices() > 0);
    }
  }
};

static const TemplateSetTest<StandardTemplateSet<TwoRefinementData>, Shape::Hypercube<2>>
  template_set_test_two_refinement_2d;
static const TemplateSetTest<StandardTemplateSet<TwoRefinementData>, Shape::Hypercube<3>>
  template_set_test_tworefinement_3d;

static const TemplateSetTest<StandardTemplateSet<SchneidersData>, Shape::Hypercube<2>> template_set_test_schneiders_2d;
static const TemplateSetTest<StandardTemplateSet<SchneidersData>, Shape::Hypercube<3>> template_set_test_schneiders_3d;

static const TemplateSetTest<IsolatedPointTemplateSet<SunZhaoMaData>, Shape::Hypercube<2>>
  template_set_test_sun_zhao_ma_2d;
static const TemplateSetTest<IsolatedPointTemplateSet<SunZhaoMaData>, Shape::Hypercube<3>>
  template_set_test_sun_zhao_ma_3d;

static const TemplateSetTest<IsolatedPointTemplateSet<SunZhaoMaExpansionData>, Shape::Hypercube<2>>
  template_set_test_sun_zhao_ma_expansion_2d;
static const TemplateSetTest<IsolatedPointTemplateSet<SunZhaoMaExpansionData>, Shape::Hypercube<3>>
  template_set_test_sun_zhao_ma_expansion_3d;
