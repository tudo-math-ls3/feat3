// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "kernel/geometry/export_vtk.hpp"
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part_ops.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/shape.hpp>
#include <test_system/test_system.hpp>

#include <bitset>
#include <memory>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

const static class MeshPartOpsTest : public UnitTest
{
  using DataType = Real;
  using MeshType = ConformalMesh<Shape::Hexahedron, 3, DataType>;
  using MeshPartType = MeshPart<MeshType>;

  using MeshPartOps = MeshPartOperations<MeshType>;

  using AttributeSetType = AttributeSet<DataType>;

public:
  MeshPartOpsTest() : UnitTest("mesh_part_ops-test")
  {
  }

  MeshPartOpsTest(const MeshPartOpsTest&) = default;
  MeshPartOpsTest(MeshPartOpsTest&&) = delete;
  MeshPartOpsTest& operator=(const MeshPartOpsTest&) = delete;
  MeshPartOpsTest& operator=(MeshPartOpsTest&&) = delete;

  ~MeshPartOpsTest() override = default;

  void run() const override
  {
    MeshType mesh = create_test_mesh();

    test_face_intersection(mesh);
    test_cell_intersection(mesh);
    test_union(mesh);
    test_disjoint_union(mesh);
    test_difference(mesh);
  }

  void test_face_intersection(const MeshType& mesh) const
  {
    const MeshPartType left = make_mesh_part(mesh, std::bitset<3>(0b100), DataType(1));
    const MeshPartType middle = make_mesh_part(mesh, std::bitset<3>(0b010), DataType(3));

    const MeshPartType face_part =
      MeshPartOps::meshpart_intersection(left, middle, [&](DataType /*l*/, DataType r) { return r; });

    // Target sets of part
    const auto& vertex_target_set = face_part.get_target_set<0>();
    const auto& edge_target_set = face_part.get_target_set<1>();
    const auto& face_target_set = face_part.get_target_set<2>();

    // Face intersection should contain exactly one face
    TEST_CHECK_EQUAL(face_part.get_num_entities(0), 4);
    TEST_CHECK_EQUAL(face_part.get_num_entities(1), 4);
    TEST_CHECK_EQUAL(face_part.get_num_entities(2), 1);
    TEST_CHECK_EQUAL(face_part.get_num_entities(3), 0);

    // Face in intersection is right face of cell 0
    const Index face_idx = mesh.get_index_set<3, 2>()(0, 5);

    // Face of part and assumed face of intersection match
    TEST_CHECK_EQUAL(face_idx, face_target_set[0]);

    // Index sets of universe
    const auto& v_at_f_u = mesh.get_index_set<2, 0>();
    const auto& v_at_e_u = mesh.get_index_set<1, 0>();
    const auto& e_at_f_u = mesh.get_index_set<2, 1>();

    // Index sets of part
    const auto& v_at_f_p = face_part.get_index_set<2, 0>();
    const auto& v_at_e_p = face_part.get_index_set<1, 0>();
    const auto& e_at_f_p = face_part.get_index_set<2, 1>();

    // Sizes of index sets should match mesh part size
    TEST_CHECK(v_at_f_p.get_num_entities() == 1);
    TEST_CHECK(e_at_f_p.get_num_entities() == 1);
    TEST_CHECK(v_at_e_p.get_num_entities() == 4);

    // Mesh part index sets should be consistent with underlying mesh
    for(int i(0); i < 4; i++)
    {
      TEST_CHECK(v_at_f_u(face_idx, i) == vertex_target_set[v_at_f_p(0, i)]);
      TEST_CHECK(e_at_f_u(face_idx, i) == edge_target_set[e_at_f_p(0, i)]);
    }
    for(Index i(0); i < 4; i++)
    {
      Index part_edge = e_at_f_p(0, int(i));
      Index universe_edge = edge_target_set[part_edge];

      TEST_CHECK_EQUAL(v_at_e_u(universe_edge, 0), vertex_target_set[v_at_e_p(part_edge, 0)]);
      TEST_CHECK_EQUAL(v_at_e_u(universe_edge, 1), vertex_target_set[v_at_e_p(part_edge, 1)]);
    }
  }

  void test_cell_intersection(const MeshType& mesh) const
  {
    const MeshPartType universe = make_mesh_part(mesh, std::bitset<3>(0b111), DataType(2));
    const MeshPartType left = make_mesh_part(mesh, std::bitset<3>(0b100), DataType(1));

    const MeshPartType intersection_part =
      MeshPartOps::meshpart_intersection(left, universe, [&](DataType l, DataType r) { return l + r; });

    // Intersection result should be one cell
    TEST_CHECK_EQUAL(intersection_part.get_num_entities(0), 8);
    TEST_CHECK_EQUAL(intersection_part.get_num_entities(1), 12);
    TEST_CHECK_EQUAL(intersection_part.get_num_entities(2), 6);
    TEST_CHECK_EQUAL(intersection_part.get_num_entities(3), 1);

    // Intersection result should contain vertices 0, 1, 4, 5, 8, 9, 12, 13
    {
      const auto& vertex_target_set = intersection_part.get_target_set<0>();
      const Index* begin = vertex_target_set.get_indices();
      const Index* end = begin + vertex_target_set.get_num_entities();
      for(Index i : {0, 1, 4, 5, 8, 9, 12, 13})
      {
        TEST_CHECK(std::find(begin, end, i) != end);
      }
    }

    // Intersection result should contain cell 0
    TEST_CHECK_EQUAL(intersection_part.get_target_set<3>()[0], 0);

    // Intersection should have added attributes
    const AttributeSetType* attribute_set = intersection_part.find_attribute("attribute");
    for(Index i(0); i < intersection_part.get_num_entities(0); i++)
    {
      TEST_CHECK_EQUAL(attribute_set->operator()(i, 0), DataType(3));
    }
  }

  void test_union(const MeshType& mesh) const
  {
    const MeshPartType left = make_mesh_part(mesh, std::bitset<3>(0b100), DataType(1));
    const MeshPartType middle = make_mesh_part(mesh, std::bitset<3>(0b010), DataType(3));

    const MeshPartType union_part =
      MeshPartOps::meshpart_union(left, middle, [&](DataType l, DataType r) { return l + r; });

    // Union should be two adjacent cells
    TEST_CHECK_EQUAL(union_part.get_num_entities(0), 12);
    TEST_CHECK_EQUAL(union_part.get_num_entities(1), 20);
    TEST_CHECK_EQUAL(union_part.get_num_entities(2), 11);
    TEST_CHECK_EQUAL(union_part.get_num_entities(3), 2);

    // Intersection should have attribute
    const AttributeSetType* attribute_set = union_part.find_attribute("attribute");
    const TargetSet& vertex_target_set = union_part.get_target_set<0>();
    for(Index i(0); i < union_part.get_num_entities(0); i++)
    {
      Index u = vertex_target_set[i];
      if(u == 0 || u == 4 || u == 8 || u == 12)
      {
        // Left vertices should have kept their value
        TEST_CHECK_EQUAL(attribute_set->operator()(i, 0), DataType(1));
      }
      else if(u == 2 || u == 6 || u == 10 || u == 14)
      {
        // Right vertices should have kept their value
        TEST_CHECK_EQUAL(attribute_set->operator()(i, 0), DataType(3));
      }
      else
      {
        // All other vertices should have been added
        TEST_CHECK_EQUAL(attribute_set->operator()(i, 0), DataType(4));
      }
    }
  }

  void test_disjoint_union(const MeshType& mesh) const
  {
    const MeshPartType left = make_mesh_part(mesh, std::bitset<3>(0b100), DataType(1));
    const MeshPartType right = make_mesh_part(mesh, std::bitset<3>(0b001), DataType(3));

    const MeshPartType union_part = MeshPartOps::meshpart_union(
      left,
      right,
      [&](DataType l, DataType r)
      {
        // No common elements; merge function should never be called
        TEST_CHECK(false);
        return l + r;
      });

    // Union should be two disjoint cells
    TEST_CHECK_EQUAL(union_part.get_num_entities(0), 16);
    TEST_CHECK_EQUAL(union_part.get_num_entities(1), 24);
    TEST_CHECK_EQUAL(union_part.get_num_entities(2), 12);
    TEST_CHECK_EQUAL(union_part.get_num_entities(3), 2);
  }


  void test_difference(const MeshType& mesh) const
  {
    const MeshPartType left_middle = make_mesh_part(mesh, std::bitset<3>(0b110), DataType(1), false);
    const MeshPartType middle = make_mesh_part(mesh, std::bitset<3>(0b010), DataType(3), false);

    const MeshPartType difference_part =
      MeshPartOps::meshpart_difference(left_middle, middle);

    TEST_CHECK_EQUAL(difference_part.get_num_entities(0), 4);
    TEST_CHECK_EQUAL(difference_part.get_num_entities(1), 8);
    TEST_CHECK_EQUAL(difference_part.get_num_entities(2), 5);
    TEST_CHECK_EQUAL(difference_part.get_num_entities(3), 1);

    // All attribute values should be from left mesh part
    const AttributeSetType* attribute_set = difference_part.find_attribute("attribute");
    for(Index i(0); i < difference_part.get_num_entities(0); i++)
    {
      TEST_CHECK_EQUAL(attribute_set->operator()(i, 0), DataType(1));
    }
  }

private:
  static MeshPartType make_mesh_part(const MeshType& mesh, std::bitset<3> cell_mask, DataType attribute, bool create_topology=true)
  {
    const std::array<Index, 4> size{
      0, // deduced
      0, // deduced
      0, // deduced
      cell_mask.count(),
    };

    MeshPartType result = MeshPartType(size.data());
    auto& target_set = result.get_target_set<3>();

    Index next = 0;
    for(Index i(0); i < 3; i++)
    {
      // Invert mask to match left-to-right order of cells in test mesh
      if(cell_mask[2 - i])
      {
        target_set[next++] = i;
      }
    }

    result.deduct_target_sets_from_top<3>(mesh.get_index_set_holder());

    if(create_topology)
    {
      result.deduct_topology(mesh.get_index_set_holder());
    }

    auto attribute_set = std::make_unique<AttributeSetType>(result.get_num_entities(0));

    for(Index i(0); i < result.get_num_entities(0); i++)
    {
      attribute_set->operator()(i, 0) = attribute;
    }
    result.add_attribute(std::move(attribute_set), "attribute");

    return result;
  }
  /**
   * \brief Create test mesh
   *
   * The test mesh is a simple row of 3 cells,
   * which is enough for testing boolean operations on mesh parts.
   */
  static MeshType create_test_mesh()
  {
    using VertexType = typename MeshType::VertexType;

    const std::array<Index, 4> mesh_size{16, 0, 0, 3};
    MeshType mesh(mesh_size.data());

    auto& vertices = mesh.get_vertex_set();
    vertices[0] = VertexType{0.0, 0.0, 0.0};
    vertices[1] = VertexType{1.0, 0.0, 0.0};
    vertices[2] = VertexType{2.0, 0.0, 0.0};
    vertices[3] = VertexType{3.0, 0.0, 0.0};

    vertices[4] = VertexType{0.0, 1.0, 0.0};
    vertices[5] = VertexType{1.0, 1.0, 0.0};
    vertices[6] = VertexType{2.0, 1.0, 0.0};
    vertices[7] = VertexType{3.0, 1.0, 0.0};

    vertices[8] = VertexType{0.0, 0.0, 1.0};
    vertices[9] = VertexType{1.0, 0.0, 1.0};
    vertices[10] = VertexType{2.0, 0.0, 1.0};
    vertices[11] = VertexType{3.0, 0.0, 1.0};

    vertices[12] = VertexType{0.0, 1.0, 1.0};
    vertices[13] = VertexType{1.0, 1.0, 1.0};
    vertices[14] = VertexType{2.0, 1.0, 1.0};
    vertices[15] = VertexType{3.0, 1.0, 1.0};

    auto& cells = mesh.get_index_set<3, 0>();

    for(Index i = 0; i < 3; i++)
    {
      cells[i][0] = i;
      cells[i][1] = i + 1;
      cells[i][2] = i + 4;
      cells[i][3] = i + 4 + 1;

      cells[i][4] = i + 8;
      cells[i][5] = i + 1 + 8;
      cells[i][6] = i + 4 + 8;
      cells[i][7] = i + 4 + 1 + 8;
    }

    mesh.deduct_topology_from_top();

    ExportVTK<MeshType> exporter(mesh);
    exporter.write("mesh_part_ops-test");

    return mesh;
  }
} mesh_part_ops_test;
