// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/shape_convert_factory.hpp>
#include <kernel/geometry/test_aux/standard_quad.hpp>
#include <kernel/geometry/test_aux/standard_hexa.hpp>
#include <kernel/geometry/test_aux/standard_tetra.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

typedef ConformalMesh<Shape::Triangle> TriaMesh;
typedef ConformalMesh<Shape::Quadrilateral> QuadMesh;
typedef ConformalMesh<Shape::Tetrahedron> TetraMesh;
typedef ConformalMesh<Shape::Hexahedron> HexaMesh;

/**
 * \brief Test class for the ShapeConvertFactory class template
 *
 * \test Tests the ShapeConvertFactory class template.
 *
 * \author Peter Zajac
 */
class ShapeConvertTest
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  ShapeConvertTest() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("ShapeConvertTest")
  {
  }

  virtual ~ShapeConvertTest()
  {
  }

  virtual void run() const override
  {
    // test quad->tria conversion
    unit_quad_to_tria();
    // test hexa->tetra conversion
    unit_hexa_to_tetra();
    // test tria->quad conversion is used in test tetra->hexa

    // test tetra->hexa conversion
    unit_tetra_to_hexa();
  }

  void unit_quad_to_tria() const
  {
    // create a unit-quad mesh
    UnitCubeFactory<QuadMesh> factory;
    QuadMesh quad_mesh(factory);

    // convert quad to tri
    ShapeConvertFactory<TriaMesh> convertor(quad_mesh);
    TriaMesh tria_mesh(convertor);

    // check counts
    TEST_CHECK_EQUAL(tria_mesh.get_num_entities(0), Index(5)); // 5 vertices
    TEST_CHECK_EQUAL(tria_mesh.get_num_entities(1), Index(8)); // 8 edges
    TEST_CHECK_EQUAL(tria_mesh.get_num_entities(2), Index(4)); // 4 triangles

    // check vertices
    static const Real vtx[] =
    {
      0.0, 0.0,
      1.0, 0.0,
      0.0, 1.0,
      1.0, 1.0,
      0.5, 0.5
    };

    // set up vertices-at-edge array
    static const Index v_e[] =
    {
      0, 1,
      2, 3,
      0, 2,
      1, 3,
      0, 4,
      1, 4,
      2, 4,
      3, 4
    };

    // set up vertices-at-edge array
    static const Index v_t[] =
    {
      4, 0, 1,
      4, 3, 2,
      4, 2, 0,
      4, 1, 3
    };

    // set up edges-at-tri array
    static const Index e_t[] =
    {
      0, 5, 4,
      1, 6, 7,
      2, 4, 6,
      3, 7, 5
    };

    // validate
    TEST_CHECK_MSG(TestAux::comp_vtx(tria_mesh.get_vertex_set(), vtx), "Vertex conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tria_mesh.get_index_set<1,0>(), v_e), "Vertex-At-Edge conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tria_mesh.get_index_set<2,0>(), v_t), "Vertex-At-Tri conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tria_mesh.get_index_set<2,1>(), e_t), "Edge-At-Tri conversion failure");
  }

  void unit_hexa_to_tetra() const
  {
    // create a Hexa mesh
    HexaMesh* hexa_mesh = TestAux::create_hexa_mesh_3d(0);

    // convert hexa to tetra
    ShapeConvertFactory<TetraMesh> convertor(*hexa_mesh);
    TetraMesh tetra_mesh(convertor);

    delete hexa_mesh;

    // check counts
    TEST_CHECK_EQUAL(tetra_mesh.get_num_entities(0), Index(8+6+1)); // 15 vertices
    TEST_CHECK_EQUAL(tetra_mesh.get_num_entities(1), Index(12+4*6+8+6)); // 50 edges
    TEST_CHECK_EQUAL(tetra_mesh.get_num_entities(2), Index(4*6+4*6+12)); // 60 triangles
    TEST_CHECK_EQUAL(tetra_mesh.get_num_entities(3), Index(4*6)); // 24 tetras


    // check vertices
    static const Real vtx[3*15] =
    {
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      1.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
      1.0, 0.0, 1.0,
      0.0, 1.0, 1.0,
      1.0, 1.0, 1.0,
      0.5, 0.5, 0.0,
      0.5, 0.5, 1.0,
      0.5, 0.0, 0.5,
      0.5, 1.0, 0.5,
      0.0, 0.5, 0.5,
      1.0, 0.5, 0.5,
      0.5, 0.5, 0.5
    };

    // set up vertices-at-edge array
    static const Index v_e[] =
    {
      0, 1, // 0
      2, 3,
      4, 5,
      6, 7,
      0, 2,
      1, 3,// 5
      4, 6,
      5, 7,
      0, 4,
      1, 5,
      2, 6,// 10
      3, 7,
      0, 8,
      1, 8,
      2, 8,
      3, 8,// 15
      4, 9,
      5, 9,
      6, 9,
      7, 9,
      0, 10,// 20
      1, 10,
      4, 10,
      5, 10,
      2, 11,
      3, 11,// 25
      6, 11,
      7, 11,
      0, 12,
      2, 12,
      4, 12,// 30
      6, 12,
      1, 13,
      3, 13,
      5, 13,
      7, 13,// 35
      0, 14,
      1, 14,
      2, 14,
      3, 14,
      4, 14,// 40
      5, 14,
      6, 14,
      7, 14,
      8, 14,
      9, 14,// 45
      10, 14,
      11, 14,
      12, 14,
      13, 14
    };

    // set up vertices-at-triangle array
    static const Index v_t[] =
    {
      8, 0, 1,//  0
      8, 3, 2,
      8, 2, 0,
      8, 1, 3,
      9, 4, 5,
      9, 7, 6,// 5
      9, 6, 4,
      9, 5, 7,
      10, 0, 1,
      10, 5, 4,
      10, 4, 0,// 10
      10, 1, 5,
      11, 2, 3,
      11, 7, 6,
      11, 6, 2,
      11, 3, 7,// 15
      12, 0, 2,
      12, 6, 4,
      12, 4, 0,
      12, 2, 6,
      13, 1, 3,// 20
      13, 7, 5,
      13, 5, 1,
      13, 3, 7,
      // coarse edge triangles
      14, 0, 1,
      14, 3, 2,// 25
      14, 4, 5,
      14, 7, 6,
      14, 2, 0,
      14, 1, 3,
      14, 6, 4,// 30
      14, 5, 7,
      14, 4, 0,
      14, 1, 5,
      14, 6, 2,
      14, 3, 7,// 35
      // face triangles
      14, 8, 0,
      14, 8, 1,
      14, 8, 2,
      14, 8, 3,
      14, 9, 4,// 40
      14, 9, 5,
      14, 9, 6,
      14, 9, 7,
      14, 10, 0,
      14, 10, 1,// 45
      14, 10, 4,
      14, 10, 5,
      14, 11, 2,
      14, 11, 3,
      14, 11, 6,// 50
      14, 11, 7,
      14, 12, 0,
      14, 12, 2,
      14, 12, 4,
      14, 12, 6,// 55
      14, 13, 1,
      14, 13, 3,
      14, 13, 5,
      14, 13, 7
    };



    // set up vertex-at-tetra array
    static const Index v_tetra[] =
    {
    //face 0
      0, 1, 8, 14,
      3, 2, 8, 14,
      2, 0, 8, 14,
      1, 3, 8, 14,
    //face 1
      5, 4, 9, 14,
      6, 7, 9, 14,
      4, 6, 9, 14,
      7, 5, 9, 14,
    //face 2
      1, 0, 10, 14,
      4, 5, 10, 14,
      0, 4, 10, 14,
      5, 1, 10, 14,
    //face 3
      2, 3, 11, 14,
      7, 6, 11, 14,
      6, 2, 11, 14,
      3, 7, 11, 14,
    //face 4
      0, 2, 12, 14,
      6, 4, 12, 14,
      4, 0, 12, 14,
      2, 6, 12, 14,
    //face 5
      3, 1, 13, 14,
      5, 7, 13, 14,
      1, 5, 13, 14,
      7, 3, 13, 14
    };

    // set up edge-at-triangle array
    static const Index e_t[] =
    {
      0, 13, 12,
      1, 14, 15,
      4, 12, 14,
      5, 15, 13,
      2, 17, 16,
      3, 18, 19,
      6, 16, 18,
      7, 19, 17,
      0, 21, 20,
      2, 22, 23,
      8, 20, 22,
      9, 23, 21,
      1, 25, 24,
      3, 26, 27,
      10, 24, 26,
      11, 27, 25,
      4, 29, 28,
      6, 30, 31,
      8, 28, 30,
      10, 31, 29,
      5, 33, 32,
      7, 34, 35,
      9, 32, 34,
      11, 35, 33,
      0, 37, 36,
      1, 38, 39,
      2, 41, 40,
      3, 42, 43,
      4, 36, 38,
      5, 39, 37,
      6, 40, 42,
      7, 43, 41,
      8, 36, 40,
      9, 41, 37,
      10, 38, 42,
      11, 43, 39,

      12, 36, 44,
      13, 37, 44,
      14, 38, 44,
      15, 39, 44,

      16, 40, 45,
      17, 41, 45,
      18, 42, 45,
      19, 43, 45,

      20, 36, 46,
      21, 37, 46,
      22, 40, 46,
      23, 41, 46,

      24, 38, 47,
      25, 39, 47,
      26, 42, 47,
      27, 43, 47,

      28, 36, 48,
      29, 38, 48,
      30, 40, 48,
      31, 42, 48,

      32, 37, 49,
      33, 39, 49,
      34, 41, 49,
      35, 43, 49
    };

    // set up triangle-at-tetra array
    static const Index t_tetra[] =
    {
      37, 36, 24, 0,
      38, 39, 25, 1,
      36, 38, 28, 2,
      39, 37, 29, 3,

      40, 41, 26, 4,
      43, 42, 27, 5,
      42, 40, 30, 6,
      41, 43, 31, 7,

      44, 45, 24, 8,
      47, 46, 26, 9,
      46, 44, 32, 10,
      45, 47, 33, 11,

      49, 48, 25, 12,
      50, 51, 27, 13,
      48, 50, 34, 14,
      51, 49, 35, 15,

      53, 52, 28, 16,
      54, 55, 30, 17,
      52, 54, 32, 18,
      55, 53, 34, 19,

      56, 57, 29, 20,
      59, 58, 31, 21,
      58, 56, 33, 22,
      57, 59, 35, 23
    };

    // set up edge-at-tetra array
    static const Index e_tetra[] =
    {
      0, 12, 36, 13, 37, 44,
      1, 15, 39, 14, 38, 44,
      4, 14, 38, 12, 36, 44,
      5, 13, 37, 15, 39, 44,

      2, 17, 41, 16, 40, 45,
      3, 18, 42, 19, 43, 45,
      6, 16, 40, 18, 42, 45,
      7, 19, 43, 17, 41, 45,

      0, 21, 37, 20, 36, 46,
      2, 22, 40, 23, 41, 46,
      8, 20, 36, 22, 40, 46,
      9, 23, 41, 21, 37, 46,

      1, 24, 38, 25, 39, 47,
      3, 27, 43, 26, 42, 47,
      10, 26, 42, 24, 38, 47,
      11, 25, 39, 27, 43, 47,

      4, 28, 36, 29, 38, 48,
      6, 31, 42, 30, 40, 48,
      8, 30, 40, 28, 36, 48,
      10, 29, 38, 31, 42, 48,

      5, 33, 39, 32, 37, 49,
      7, 34, 41, 35, 43, 49,
      9, 32, 37, 34, 41, 49,
      11, 35, 43, 33, 39, 49
    };

    // validate
    TEST_CHECK_MSG(TestAux::comp_vtx(tetra_mesh.get_vertex_set(), vtx), "Vertex conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tetra_mesh.get_index_set<1,0>(), v_e), "Vertex-At-Edge conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tetra_mesh.get_index_set<2,0>(), v_t), "Vertex-At-Tri conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tetra_mesh.get_index_set<3,0>(), v_tetra), "Vertex-At-Tetra conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tetra_mesh.get_index_set<2,1>(), e_t), "Edge-At-Tri conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tetra_mesh.get_index_set<3,1>(), e_tetra), "Edge-At-Tetra conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(tetra_mesh.get_index_set<3,2>(), t_tetra), "Tri-At-Tetra conversion failure");
  }

  void unit_tetra_to_hexa() const
  {
    // create a Tetra mesh
    TetraMesh* tetra_mesh = TestAux::create_tetra_mesh_3d(0);

    // convert tetra to hexa
    ShapeConvertFactory<HexaMesh> convertor(*tetra_mesh);
    HexaMesh hexa_mesh(convertor);

    delete tetra_mesh;

    // check counts
    TEST_CHECK_EQUAL(hexa_mesh.get_num_entities(0), Index(4+6+4+1)); // 15 vertices
    TEST_CHECK_EQUAL(hexa_mesh.get_num_entities(1), Index(2*6+4*3+4)); // 28 edges
    TEST_CHECK_EQUAL(hexa_mesh.get_num_entities(2), Index(4*3+6)); // 18 triangles
    TEST_CHECK_EQUAL(hexa_mesh.get_num_entities(3), Index(4)); // 4 tetras


    // check vertices
    static const Real vtx[3*15] =
    {
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
      0.5, 0.0, 0.0,
      0.0, 0.5, 0.0,
      0.0, 0.0, 0.5,
      0.5, 0.5, 0.0,
      0.5, 0.0, 0.5,
      0.0, 0.5, 0.5,
      1.0/3.0, 1.0/3.0, 1.0/3.0,
      0.0, 1.0/3.0, 1.0/3.0,
      1.0/3.0, 0.0, 1.0/3.0,
      1.0/3.0, 1.0/3.0, 0.0,
      1.0/4.0, 1.0/4.0, 1.0/4.0
    };

    // set up vertices-at-edge array
    static const Index v_e[] =
    {
      // edges
      0, 4, // 0
      4, 1,
      0, 5,
      5, 2,
      0, 6,
      6, 3, // 5
      1, 7,
      7, 2,
      1, 8,
      8, 3,
      2, 9, // 10
      9, 3,
      // faces
      9, 10,
      8, 10,
      7, 10,
      9, 11, // 15
      6, 11,
      5, 11,
      8, 12,
      6, 12,
      4, 12, // 20
      7, 13,
      5, 13,
      4, 13,
      //
      10, 14,
      11, 14, // 25
      12, 14,
      13, 14
    };

    // set up vertices-at-quad array
    static const Index v_q[] =
    {
      10, 8, 7, 1, // 0
      10, 7, 9, 2,
      10, 9, 8, 3,
      11, 6, 5, 0,
      11, 5, 9, 2,
      11, 9, 6, 3, // 5
      12, 6, 4, 0,
      12, 4, 8, 1,
      12, 8, 6, 3,
      13, 5, 4, 0,
      13, 4, 7, 1, // 10
      13, 7, 5, 2,
      //
      14, 12, 13, 4,
      14, 11, 13, 5,
      14, 11, 12, 6,
      14, 10, 13, 7, // 15
      14, 10, 12, 8,
      14, 10, 11, 9
    };



    // set up vertex-at-hexa array
    static const Index v_hexa[] =
    {
      14, 11, 12, 6, 13, 5, 4, 0,
      14, 10, 12, 8, 13, 7, 4, 1,
      14, 10, 11, 9, 13, 7, 5, 2,
      14, 10, 11, 9, 12, 8, 6, 3
    };

    // set up edge-at-quad array
    static const Index e_q[] =
    {
      13, 6, 14, 8,
      14, 10, 12, 7,
      12, 9, 13, 11,
      16, 2, 17, 4,
      17, 10, 15, 3,
      15, 5, 16, 11,
      19, 0, 20, 4,
      20, 8, 18, 1,
      18, 5, 19, 9,
      22, 0, 23, 2,
      23, 6, 21, 1,
      21, 3, 22, 7,
      26, 23, 27, 20,
      25, 22, 27, 17,
      25, 19, 26, 16,
      24, 21, 27, 14,
      24, 18, 26, 13,
      24, 15, 25, 12
    };

    // set up quad-at-hexa array
    static const Index q_hexa[] =
    {
      14, 9, 13, 6, 12, 3,
      16, 10, 15, 7, 12, 0,
      17, 11, 15, 4, 13, 1,
      17, 8, 16, 5, 14, 2
    };

    // set up edge-at-hexa array
    static const Index e_hexa[] =
    {
      25, 19, 22, 0, 26, 16, 23, 2, 27, 17, 20, 4,
      24, 18, 21, 1, 26, 13, 23, 6, 27, 14, 20, 8,
      24, 15, 21, 3, 25, 12, 22, 7, 27, 14, 17, 10,
      24, 15, 18, 5, 25, 12, 19, 9, 26, 13, 16, 11
    };

    // validate
    TEST_CHECK_MSG(TestAux::comp_vtx(hexa_mesh.get_vertex_set(), vtx), "Vertex conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(hexa_mesh.get_index_set<1,0>(), v_e), "Vertex-At-Edge conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(hexa_mesh.get_index_set<2,0>(), v_q), "Vertex-At-Quad conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(hexa_mesh.get_index_set<3,0>(), v_hexa), "Vertex-At-Hexa conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(hexa_mesh.get_index_set<2,1>(), e_q), "Edge-At-Quad conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(hexa_mesh.get_index_set<3,1>(), e_hexa), "Edge-At-Hexa conversion failure");
    TEST_CHECK_MSG(TestAux::comp_idx(hexa_mesh.get_index_set<3,2>(), q_hexa), "Quad-At-Hexa conversion failure");
  }

} shape_convert_test;
