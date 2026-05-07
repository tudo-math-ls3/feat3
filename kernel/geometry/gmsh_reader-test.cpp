// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/conformal_mesh.hpp>   // for ConformalMesh
#include <kernel/geometry/gmsh_reader.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/shape.hpp>
#include <sstream>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

template<typename T, class... StreamArgs>
inline std::basic_ostream<StreamArgs...>& operator<=(std::basic_ostream<StreamArgs...>& out, T const& data)
{
  out.write(reinterpret_cast<char const*>(&data), sizeof(T));
  return out;
}

class MshReaderTest : public UnitTest
{
public:
  MshReaderTest() : UnitTest("MshReaderTest")
  {
  }

  ~MshReaderTest() override = default;

  void run() const override
  {
    test_hex(hex(true));
    test_hex(hex(false));
    test_tet(tet(true));
    test_tet(tet(false));
    test_quad(quad(true));
    test_quad(quad(false));
    test_triangle(triangle(true));
    test_triangle(triangle(false));
  }

  void test_triangle(const std::stringstream& stream) const
  {
    using MeshType = ConformalMesh<Shape::Triangle, 3>;
    using VertexType = typename MeshType::VertexType;

    GmshFileReader reader(stream);
    reader.read_root_markup();

    TEST_CHECK_EQUAL(reader.get_meshtype_string(), "conformal:simplex:2:3");

    MeshAtlas<MeshType> atlas;
    auto node = reader.parse(atlas);

    const Real tol = TestSystem::tol<Real>();

    // Vertices should match nodes of file
    const auto& vertices = node->get_mesh()->get_vertex_set();
    TEST_CHECK_EQUAL(node->get_mesh()->get_num_vertices(), 3);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[0] - VertexType{0.0, 0.0, 0.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[1] - VertexType{1.0, 0.0, 0.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[2] - VertexType{0.0, 1.0, 0.0}).norm_euclid(), 0, tol);

    // NOTE: The indices below are correct, the .msh files are one-based, but FEAT is zero-based
    // Indices should be flipped to match FEAT zig-zag convention
    const auto& v_at_c = node->get_mesh()->get_index_set<2, 0>();
    TEST_CHECK_EQUAL(node->get_mesh()->get_num_elements(), 1);
    TEST_CHECK_EQUAL(v_at_c(0, 0), 0);
    TEST_CHECK_EQUAL(v_at_c(0, 1), 1);
    TEST_CHECK_EQUAL(v_at_c(0, 2), 2);
  }

  void test_quad(const std::stringstream& stream) const
  {
    using MeshType = ConformalMesh<Shape::Quadrilateral, 3>;
    using VertexType = typename MeshType::VertexType;

    GmshFileReader reader(stream);
    reader.read_root_markup();

    TEST_CHECK_EQUAL(reader.get_meshtype_string(), "conformal:hypercube:2:3");

    MeshAtlas<MeshType> atlas;
    auto node = reader.parse(atlas);

    const Real tol = TestSystem::tol<Real>();

    // Vertices should match nodes of file
    const auto& vertices = node->get_mesh()->get_vertex_set();
    TEST_CHECK_EQUAL(node->get_mesh()->get_num_vertices(), 4);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[0] - VertexType{-1.0, -1.0, 0.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[1] - VertexType{1.0, -1.0, 0.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[2] - VertexType{1.0, 1.0, 0.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[3] - VertexType{-1.0, 1.0, 0.0}).norm_euclid(), 0, tol);

    // NOTE: The indices below are correct, the .msh files are one-based, but FEAT is zero-based
    // Indices should be flipped to match FEAT zig-zag convention
    const auto& v_at_c = node->get_mesh()->get_index_set<2, 0>();
    TEST_CHECK_EQUAL(node->get_mesh()->get_num_elements(), 1);
    TEST_CHECK_EQUAL(v_at_c(0, 0), 0);
    TEST_CHECK_EQUAL(v_at_c(0, 1), 1);
    TEST_CHECK_EQUAL(v_at_c(0, 2), 3);
    TEST_CHECK_EQUAL(v_at_c(0, 3), 2);
  }

  void test_hex(const std::stringstream& stream) const
  {
    using MeshType = ConformalMesh<Shape::Hexahedron>;
    using VertexType = typename MeshType::VertexType;

    GmshFileReader reader(stream);
    reader.read_root_markup();

    TEST_CHECK_EQUAL(reader.get_meshtype_string(), "conformal:hypercube:3:3");

    MeshAtlas<MeshType> atlas;
    auto node = reader.parse(atlas);

    const Real tol = TestSystem::tol<Real>();

    // Vertices should match nodes of file
    const auto& vertices = node->get_mesh()->get_vertex_set();
    TEST_CHECK_EQUAL(node->get_mesh()->get_num_vertices(), 8);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[0] - VertexType{-1.0, -1.0, -1.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[1] - VertexType{1.0, -1.0, -1.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[2] - VertexType{1.0, 1.0, -1.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[3] - VertexType{-1.0, 1.0, -1.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[4] - VertexType{-1.0, -1.0, 1.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[5] - VertexType{1.0, -1.0, 1.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[6] - VertexType{1.0, 1.0, 1.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[7] - VertexType{-1.0, 1.0, 1.0}).norm_euclid(), 0, tol);

    // NOTE: The indices below are correct, the .msh files are one-based, but FEAT is zero-based
    // Indices should be flipped to match FEAT zig-zag convention
    const auto& v_at_c = node->get_mesh()->get_index_set<3, 0>();
    TEST_CHECK_EQUAL(node->get_mesh()->get_num_elements(), 1);
    TEST_CHECK_EQUAL(v_at_c(0, 0), 0);
    TEST_CHECK_EQUAL(v_at_c(0, 1), 1);
    TEST_CHECK_EQUAL(v_at_c(0, 2), 3);
    TEST_CHECK_EQUAL(v_at_c(0, 3), 2);
    TEST_CHECK_EQUAL(v_at_c(0, 4), 4);
    TEST_CHECK_EQUAL(v_at_c(0, 5), 5);
    TEST_CHECK_EQUAL(v_at_c(0, 6), 7);
    TEST_CHECK_EQUAL(v_at_c(0, 7), 6);
  }

  void test_tet(const std::stringstream& stream) const
  {
    using MeshType = ConformalMesh<Shape::Tetrahedron>;
    using VertexType = typename MeshType::VertexType;

    GmshFileReader reader(stream);
    reader.read_root_markup();

    TEST_CHECK_EQUAL(reader.get_meshtype_string(), "conformal:simplex:3:3");

    MeshAtlas<MeshType> atlas;
    auto node = reader.parse(atlas);

    const Real tol = TestSystem::tol<Real>();

    // Vertices should match nodes of file
    const auto& vertices = node->get_mesh()->get_vertex_set();
    TEST_CHECK_EQUAL(node->get_mesh()->get_num_vertices(), 4);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[0] - VertexType{0.0, 0.0, 0.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[1] - VertexType{1.0, 0.0, 0.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[2] - VertexType{0.0, 1.0, 0.0}).norm_euclid(), 0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS((vertices[3] - VertexType{0.0, 0.0, 1.0}).norm_euclid(), 0, tol);

    // NOTE: The indices below are correct, the .msh files are one-based, but FEAT is zero-based
    const auto& v_at_c = node->get_mesh()->get_index_set<3, 0>();
    TEST_CHECK_EQUAL(node->get_mesh()->get_num_elements(), 1);
    TEST_CHECK_EQUAL(v_at_c(0, 0), 0);
    TEST_CHECK_EQUAL(v_at_c(0, 1), 1);
    TEST_CHECK_EQUAL(v_at_c(0, 2), 2);
    TEST_CHECK_EQUAL(v_at_c(0, 3), 3);
  }

private:
  struct Node
  {
    int id;
    double x, y, z;
  };

  template<std::size_t len>
  struct Element
  {
    int id;
    int tag1;
    int tag2;
    std::array<int, len> nodes;
  };

  template<std::size_t num_nodes, std::size_t num_elements, std::size_t nodes_per_element>
  static std::stringstream make_msh(
    int element_type,
    const std::array<Node, num_nodes>& nodes,
    const std::array<Element<nodes_per_element>, num_elements>& elements,
    bool binary)
  {
    std::stringstream mts(binary ? std::ios::in | std::ios::out | std::ios::binary : std::ios::in | std::ios::out);

    // MeshFormat
    mts << "$MeshFormat\n";
    mts << "2.2 " << (binary ? 1 : 0) << " 8\n";

    if(binary)
    {
      mts <= 1;
      mts << "\n";
    }

    mts << "$EndMeshFormat\n";

    // Nodes
    mts << "$Nodes\n";
    mts << num_nodes << "\n";

    if(binary)
    {
      for(const auto& n : nodes)
      {
        mts <= n.id;
        mts <= n.x;
        mts <= n.y;
        mts <= n.z;
      }
    }
    else
    {
      for(Index i(0); i < num_nodes; i++)
      {
        const auto& n = nodes[i];
        mts << n.id << " " << n.x << " " << n.y << " " << n.z;

        if(i < num_nodes - 1)
        {
          mts << "\n";
        }
      }
    }

    mts << "\n$EndNodes\n";

    // Elements
    mts << "$Elements\n";
    mts << num_elements << "\n";

    if(binary)
    {
      // Chunk header
      mts <= element_type;
      mts <= (int)num_elements;
      mts <= 2; // num_tags

      for(const auto& element : elements)
      {
        mts <= element.id;
        mts <= element.tag1;
        mts <= element.tag2;

        for(int i : element.nodes)
        {
          mts <= i;
        }
      }
    }
    else
    {
      for(Index i(0); i < num_elements; i++)
      {
        const auto& element = elements[i];
        mts << element.id << " " << element_type << " 2 " << element.tag1 << " " << element.tag2 << " ";

        for(int n : element.nodes)
        {
          mts << n << " ";
        }

        if(i < num_elements - 1)
        {
          mts << "\n";
        }
      }
    }

    mts << "\n$EndElements\n";

    return mts;
  }

  static std::stringstream tet(bool binary)
  {
    constexpr std::array<Node, 4> nodes{{{1, 0, 0, 0}, {2, 1, 0, 0}, {3, 0, 1, 0}, {4, 0, 0, 1}}};

    const std::array<Element<4>, 1> elements{{{1, 0, 1, {1, 2, 3, 4}}}};

    return make_msh(4, nodes, elements, binary);
  }

  static std::stringstream hex(bool binary)
  {
    constexpr std::array<Node, 8> nodes{{
      {1, -1, -1, -1},
      {2, 1, -1, -1},
      {3, 1, 1, -1},
      {4, -1, 1, -1},
      {5, -1, -1, 1},
      {6, 1, -1, 1},
      {7, 1, 1, 1},
      {8, -1, 1, 1},
    }};

    const std::array<Element<8>, 1> elements{{{1, 0, 1, {1, 2, 3, 4, 5, 6, 7, 8}}}};

    return make_msh(5, nodes, elements, binary);
  }

  static std::stringstream quad(bool binary)
  {
    constexpr std::array<Node, 4> nodes{{
      {1, -1, -1, 0},
      {2, 1, -1, 0},
      {3, 1, 1, 0},
      {4, -1, 1, 0},
    }};

    const std::array<Element<4>, 1> elements{{{1, 0, 1, {1, 2, 3, 4}}}};

    return make_msh(3, nodes, elements, binary);
  }

  static std::stringstream triangle(bool binary)
  {
    constexpr std::array<Node, 3> nodes{{
      {1, 0, 0, 0},
      {2, 1, 0, 0},
      {3, 0, 1, 0},
    }};

    const std::array<Element<3>, 1> elements{{{1, 0, 1, {1, 2, 3}}}};

    return make_msh(2, nodes, elements, binary);
  }
}; // class MshReaderTest

static const MshReaderTest msh_reader_test;
