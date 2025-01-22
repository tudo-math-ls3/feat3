// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "kernel/geometry/intern/refinement_field.hpp"
#include "kernel/geometry/template_builder.hpp"
#include <kernel/geometry/adaptive_mesh.hpp>
#include <kernel/geometry/adaptive_mesh_layer.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/intern/adaptive_mesh_algorithms.hpp>
#include <kernel/geometry/intern/adaptive_mesh_storage.hpp>
#include <kernel/geometry/refinement_types.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/geometry/template_sets.hpp>
#include <kernel/geometry/templates/schneiders_data.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/type_traits.hpp>
#include <test_system/test_system.hpp>

#include <memory>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;

/**
 * \brief Unit tests for the AdaptiveMeshStorage class
 *
 * \author Markus Muegge
 */
template<typename DT_>
class AdaptiveMeshStorageTest : public UnitTest
{
public:
  using ShapeType = Shape::Quadrilateral;
  using TemplateSet = SchneidersTemplates;
  using Vertex = Tiny::Vector<DT_, 2>;

  using Storage = Geometry::Intern::AdaptiveMeshStorage<ShapeType, TemplateSet, Vertex>;

  using VertKey = typename Storage::template ElementRefByDim<0>;
  using EdgeKey = typename Storage::template ElementRefByDim<1>;
  using FaceKey = typename Storage::template ElementRefByDim<2>;
  using CellKey = typename Storage::template ElementRefByDim<3>;

  using AdaptiveVertex = typename Storage::AdaptiveVertexType;
  using AdaptiveEdge = typename Storage::template AdaptiveElementByDim<1>;
  using AdaptiveFace = typename Storage::template AdaptiveElementByDim<2>;

  using EdgeType = StandardRefinementType<Shape::Hypercube<1>>;
  using FaceType = StandardRefinementType<Shape::Quadrilateral>;

  AdaptiveMeshStorageTest() : UnitTest("AdaptiveMeshStorageTest", Type::Traits<DT_>::name())
  {
  }

  ~AdaptiveMeshStorageTest() override = default;

  /**
   * \brief Insert a type 1 face into the storage.
   *
   * \note Topologies are left empty. Only children are added
   */
  static FaceKey insert_test_tree(Storage& storage)
  {

    AdaptiveFace root_elem(FaceType(1), Layer{0}, {});

    root_elem.children.template by_dim<0>()[0] = storage.insert(Vertex{1.0 / 3.0, 1.0 / 3.0}, Layer{1});

    root_elem.children.template by_dim<1>()[0] = storage.insert(AdaptiveEdge(EdgeType(0), Layer{1}, {}));
    root_elem.children.template by_dim<1>()[1] = storage.insert(AdaptiveEdge(EdgeType(0), Layer{1}, {}));
    root_elem.children.template by_dim<1>()[2] = storage.insert(AdaptiveEdge(EdgeType(0), Layer{1}, {}));

    root_elem.children.template by_dim<2>()[0] = storage.insert(AdaptiveFace(FaceType(0), Layer{1}, {}));
    root_elem.children.template by_dim<2>()[1] = storage.insert(AdaptiveFace(FaceType(0), Layer{1}, {}));
    root_elem.children.template by_dim<2>()[2] = storage.insert(AdaptiveFace(FaceType(0), Layer{1}, {}));

    return storage.insert(root_elem);
  }

  void run() const override
  {
    test_insertion();
    test_num_entities();
    test_type_retrieval();
    test_erasing();
    test_iteration();
    test_indexing();
  }

  void test_insertion() const
  {
    Storage storage;

    const Vertex test_vert = Vertex{1.0, 2.0};
    AdaptiveVertex vertex(test_vert, Layer{0});
    VertKey vert_key = storage.insert(vertex);

    // Ensure we can retrieve the vertex
    TEST_CHECK_EQUAL(storage[vert_key].layer.idx, 0);

    AdaptiveEdge edge_a(EdgeType(0), Layer{0}, Geometry::Intern::ElementTopology<Shape::Hypercube<1>>{});
    EdgeKey key_a = storage.insert(edge_a);

    AdaptiveEdge edge_b(EdgeType(1), Layer{1}, Geometry::Intern::ElementTopology<Shape::Hypercube<1>>{});
    EdgeKey key_b = storage.insert(edge_b);

    // Ensure we can retirieve entities
    TEST_CHECK_EQUAL(storage[key_a].layer.idx, 0);
    TEST_CHECK_EQUAL(storage[key_b].layer.idx, 1);

    // Ensure layers are reported correctly
    TEST_CHECK_EQUAL(storage.num_layers(), 2);
  }

  void test_num_entities() const
  {
    Storage storage;

    insert_test_tree(storage);

    // Totals should equal elements added by insert_test_tree
    TEST_CHECK_EQUAL(storage.template num_total_entities<0>(), 1);
    TEST_CHECK_EQUAL(storage.template num_total_entities<1>(), 3);
    TEST_CHECK_EQUAL(storage.template num_total_entities<2>(), 4);

    // There should only be the root face on layer 0
    TEST_CHECK_EQUAL(storage.num_entities(Layer{0}, 0), 0);
    TEST_CHECK_EQUAL(storage.num_entities(Layer{0}, 1), 0);
    TEST_CHECK_EQUAL(storage.num_entities(Layer{0}, 2), 1);
    TEST_CHECK_EQUAL(storage.num_entities(Layer{0}, 3), 0);

    // There should only be the root face on layer 0
    TEST_CHECK_EQUAL(storage.num_entities(Layer{1}, 0), 1);
    TEST_CHECK_EQUAL(storage.num_entities(Layer{1}, 1), 3);
    TEST_CHECK_EQUAL(storage.num_entities(Layer{1}, 2), 3);
    TEST_CHECK_EQUAL(storage.num_entities(Layer{1}, 3), 0);
  }

  void test_type_retrieval() const
  {
    Storage storage;

    AdaptiveEdge edge_a(EdgeType(0), Layer{0}, Geometry::Intern::ElementTopology<Shape::Hypercube<1>>{});
    EdgeKey key_a = storage.insert(edge_a);

    AdaptiveEdge edge_b(EdgeType(0b10), Layer{0}, Geometry::Intern::ElementTopology<Shape::Hypercube<1>>{});
    EdgeKey key_b = storage.insert(edge_b);

    TEST_CHECK_EQUAL(storage.type(key_a), EdgeType(0));
    TEST_CHECK_EQUAL(storage.type(key_b), EdgeType(0b10));
  }

  void test_erasing() const
  {
    Storage storage;

    FaceKey root = insert_test_tree(storage);

    std::array<Index, 4> erased = storage.erase(root);

    // Ensure we erased the correct amount of elements
    TEST_CHECK_EQUAL(erased[0], 1);
    TEST_CHECK_EQUAL(erased[1], 3);
    TEST_CHECK_EQUAL(erased[2], 4);
    TEST_CHECK_EQUAL(erased[3], 0);

    // We should have erased all elements in the storage
    TEST_CHECK_EQUAL(storage.template num_total_entities<0>(), 0);
    TEST_CHECK_EQUAL(storage.template num_total_entities<1>(), 0);
    TEST_CHECK_EQUAL(storage.template num_total_entities<2>(), 0);
  }

  void test_iteration() const
  {
    Storage storage;

    struct Visitor
    {
      std::vector<String> visited;

      void operator()(VertKey /*key*/)
      {
        visited.emplace_back("Vertex");
      }

      void operator()(EdgeKey /*key*/)
      {
        visited.emplace_back("Edge");
      }

      void operator()(FaceKey /*key*/)
      {
        visited.emplace_back("Face");
      }

      void operator()(CellKey /*key*/)
      {
        visited.emplace_back("Cell");
      }
    };

    FaceKey root = insert_test_tree(storage);
    insert_test_tree(storage);

    {
      Visitor visitor{};
      storage.walk_subtree(root, visitor, Layer{1}, Storage::IterationOrder::PreOrder);

      std::array<String, 8> expected_pre_order = {"Face", "Vertex", "Edge", "Edge", "Edge", "Face", "Face", "Face"};

      TEST_CHECK_EQUAL(visitor.visited.size(), expected_pre_order.size());
      for(Index i = 0; i < expected_pre_order.size(); i++)
      {
        TEST_CHECK_EQUAL(visitor.visited[i], expected_pre_order[i]);
      }
    }

    {
      Visitor visitor{};
      storage.walk_tree(visitor, Layer{1}, Storage::IterationOrder::PreOrder);

      std::array<String, 16> expected_pre_order = {
        "Face",
        "Vertex",
        "Edge",
        "Edge",
        "Edge",
        "Face",
        "Face",
        "Face",
        "Face",
        "Vertex",
        "Edge",
        "Edge",
        "Edge",
        "Face",
        "Face",
        "Face",
      };

      TEST_CHECK_EQUAL(visitor.visited.size(), expected_pre_order.size());
      for(Index i = 0; i < expected_pre_order.size(); i++)
      {
        TEST_CHECK_EQUAL(visitor.visited[i], expected_pre_order[i]);
      }
    }

    {
      Visitor visitor{};
      storage.walk_subtree(root, visitor, Layer{1}, Storage::IterationOrder::PostOrder);

      std::array<String, 8> expected_post_order = {"Vertex", "Edge", "Edge", "Edge", "Face", "Face", "Face", "Face"};

      TEST_CHECK_EQUAL(visitor.visited.size(), expected_post_order.size());
      for(Index i = 0; i < expected_post_order.size(); i++)
      {
        TEST_CHECK_EQUAL(visitor.visited[i], expected_post_order[i]);
      }
    }

    {
      Visitor visitor{};
      storage.walk_tree(visitor, Layer{1}, Storage::IterationOrder::PostOrder);

      std::array<String, 16> expected_post_order = {
        "Vertex",
        "Edge",
        "Edge",
        "Edge",
        "Face",
        "Face",
        "Face",
        "Face",
        "Vertex",
        "Edge",
        "Edge",
        "Edge",
        "Face",
        "Face",
        "Face",
        "Face",
      };

      TEST_CHECK_EQUAL(visitor.visited.size(), expected_post_order.size());
      for(Index i = 0; i < expected_post_order.size(); i++)
      {
        TEST_CHECK_EQUAL(visitor.visited[i], expected_post_order[i]);
      }
    }
  }

  void test_indexing() const
  {
    Storage storage;

    // Test vertex indexing

    // NOTE: We use the first vertex coordinate to identify vertices
    std::array<Vertex, 6> vertices = {
      Vertex{0.0, 0.0},
      Vertex{1.0, 0.0},
      Vertex{2.0, 0.0},
      Vertex{3.0, 0.0},
      Vertex{4.0, 0.0},
      Vertex{5.0, 0.0},
    };

    std::array<VertKey, 6> vert_keys;
    vert_keys[0] = storage.insert(vertices[0], Layer{0});
    vert_keys[1] = storage.insert(vertices[1], Layer{0});
    vert_keys[2] = storage.insert(vertices[2], Layer{0});

    vert_keys[3] = storage.insert(vertices[3], Layer{1});
    vert_keys[4] = storage.insert(vertices[4], Layer{1});
    vert_keys[5] = storage.insert(vertices[5], Layer{1});

    storage.reindex();

    // NOTE: In the following we will use lists of seen coords to ensure all
    // vertices that were added were assigned a index and are available via the
    // index-based accessors, without imposing a order on the returned
    // vertices. At the time of writing this the order is identical, but there
    // is no requirement forcing this.

    // Test layer 0 indices
    {
      // Ensure keys and indices are consistent
      for(Index i = 0; i < 3; i++)
      {
        AdaptiveVertex key_vertex = storage[vert_keys[i]];
        AdaptiveVertex idx_vertex = storage.get_vertex_by_index(Layer{0}, storage.get_index(vert_keys[i]));

        TEST_CHECK_EQUAL(key_vertex.vertex[0], idx_vertex.vertex[0]);
      }

      // Ensure all elements are reachable via indices
      std::vector<Real> seen_coords;

      for(Index i = 0; i < storage.num_entities(Layer{0}, 0); i++)
      {
        typename Storage::AdaptiveVertexType& vertex = storage.get_vertex_by_index(Layer{0}, i);
        seen_coords.push_back(vertex.vertex[0]);
      }

      // We should have seen the first three vertices
      // All vertices are permanent, so we should have seen all vertices
      TEST_CHECK_EQUAL(seen_coords.size(), 3);

      std::sort(seen_coords.begin(), seen_coords.end());
      TEST_CHECK_EQUAL(seen_coords[0], 0.0);
      TEST_CHECK_EQUAL(seen_coords[1], 1.0);
      TEST_CHECK_EQUAL(seen_coords[2], 2.0);
    }

    // Test layer 1 indices
    {
      // Ensure keys and indices are consistent
      for(Index i = 0; i < 6; i++)
      {
        AdaptiveVertex key_vertex = storage[vert_keys[i]];
        AdaptiveVertex idx_vertex = storage.get_vertex_by_index(Layer{1}, storage.get_index(vert_keys[i]));

        TEST_CHECK_EQUAL(key_vertex.vertex[0], idx_vertex.vertex[0]);
      }

      // Ensure all elements are reachable via indices
      std::vector<Real> seen_coords;

      for(Index i = 0; i < storage.num_entities(Layer{1}, 0); i++)
      {
        AdaptiveVertex& vertex = storage.get_vertex_by_index(Layer{1}, i);
        seen_coords.push_back(vertex.vertex[0]);
      }

      // All vertices are permanent, so we should have seen all vertices
      TEST_CHECK_EQUAL(seen_coords.size(), 6);

      std::sort(seen_coords.begin(), seen_coords.end());
      TEST_CHECK_EQUAL(seen_coords[0], 0.0);
      TEST_CHECK_EQUAL(seen_coords[1], 1.0);
      TEST_CHECK_EQUAL(seen_coords[2], 2.0);
      TEST_CHECK_EQUAL(seen_coords[3], 3.0);
      TEST_CHECK_EQUAL(seen_coords[4], 4.0);
      TEST_CHECK_EQUAL(seen_coords[5], 5.0);
    }

    // Test element indexing

    // NOTE: We use the index to identify edges
    std::array<AdaptiveEdge, 6> edges = {
      AdaptiveEdge(EdgeType(0), Layer{0}, {}),
      AdaptiveEdge(EdgeType(1), Layer{0}, {}),
      AdaptiveEdge(EdgeType(1), Layer{0}, {}),

      AdaptiveEdge(EdgeType(2), Layer{1}, {}),
      AdaptiveEdge(EdgeType(2), Layer{1}, {}),
      AdaptiveEdge(EdgeType(2), Layer{1}, {}),
    };

    for(Index i(0); i < 6; ++i)
    {
      edges[i].index = i;
    }

    std::array<EdgeKey, 6> edge_keys;
    edge_keys[0] = storage.insert(edges[0]);
    edge_keys[1] = storage.insert(edges[1]);
    edge_keys[2] = storage.insert(edges[2]);
    edge_keys[3] = storage.insert(edges[3]);
    edge_keys[4] = storage.insert(edges[4]);
    edge_keys[5] = storage.insert(edges[5]);

    storage.reindex();

    // Test layer 0 edges
    {
      // Ensure keys and indices are consistent
      for(Index i = 0; i < 3; i++)
      {
        AdaptiveEdge key_edge = storage[edge_keys[i]];
        AdaptiveEdge idx_edge = storage.template get_by_index<1>(Layer{0}, storage.get_index(edge_keys[i]));

        TEST_CHECK_EQUAL(key_edge.type, idx_edge.type);
      }

      // Ensure all elements are reachable via indices
      std::vector<EdgeType> seen_edges;

      for(Index i = 0; i < storage.num_entities(Layer{0}, 1); i++)
      {
        AdaptiveEdge& edge = storage.template get_by_index<1>(Layer{0}, i);
        seen_edges.push_back(edge.type);
      }

      // We should have seen all Layer 0 edges
      TEST_CHECK_EQUAL(seen_edges.size(), 3);

      TEST_CHECK_EQUAL(seen_edges[0], EdgeType(0));
      TEST_CHECK_EQUAL(seen_edges[1], EdgeType(1));
      TEST_CHECK_EQUAL(seen_edges[2], EdgeType(1));
    }

    // Test layer 1 edges
    {
      // Ensure keys and indices are consistent
      for(Index i : {0, 3, 4, 5})
      {
        AdaptiveEdge key_edge = storage[edge_keys[i]];
        AdaptiveEdge idx_edge = storage.template get_by_index<1>(Layer{1}, storage.get_index(edge_keys[i]));

        TEST_CHECK_EQUAL(key_edge.type, idx_edge.type);
      }

      // Ensure all elements are reachable via indices
      std::vector<EdgeType> seen_edges;

      for(Index i = 0; i < storage.num_entities(Layer{1}, 1); i++)
      {
        AdaptiveEdge& edge = storage.template get_by_index<1>(Layer{1}, i);
        seen_edges.push_back(edge.type);
      }

      // We should have seen Layer 1 edges, plus the permanent edge of layer 0
      TEST_CHECK_EQUAL(seen_edges.size(), 4);

      TEST_CHECK_EQUAL(seen_edges[0], EdgeType(0));
      TEST_CHECK_EQUAL(seen_edges[1], EdgeType(2));
      TEST_CHECK_EQUAL(seen_edges[2], EdgeType(2));
      TEST_CHECK_EQUAL(seen_edges[3], EdgeType(2));
    }
  }
};

static const AdaptiveMeshStorageTest<float> adaptive_mesh_storage_test_float;
static const AdaptiveMeshStorageTest<double> adaptive_mesh_storage_test_double;

/**
 * \brief Unit tests for adaptive mesh refinement
 *
 * \author Markus Muegge
 */
template<typename DT_>
class AdaptiveRefinementTest : public UnitTest
{
public:
  using TemplateSet2D = SchneidersTemplates;
  using TemplateSet3D = SchneidersTemplates;
  using Mesh2D = ConformalMesh<Shape::Quadrilateral, 2, DT_>;
  using Mesh3D = ConformalMesh<Shape::Hexahedron, 3, DT_>;
  using AdaptiveMesh2D = AdaptiveMesh<TemplateSet2D, Shape::Quadrilateral, 2, DT_>;
  using AdaptiveMesh3D = AdaptiveMesh<TemplateSet3D, Shape::Hexahedron, 3, DT_>;
  using AdaptiveMeshLayer2D = AdaptiveMeshLayer<AdaptiveMesh2D>;
  using AdaptiveMeshLayer3D = AdaptiveMeshLayer<AdaptiveMesh3D>;

  AdaptiveRefinementTest() : UnitTest("AdaptiveRefinementTest", Type::Traits<DT_>::name())
  {
  }

  ~AdaptiveRefinementTest() override = default;

  template<typename MeshType_, int face_dim>
  void check_congruency(const MeshType_& mesh, Index cell) const
  {
    using ShapeType = typename MeshType_::ShapeType;
    constexpr int dim = ShapeType::dimension;

    static_assert(face_dim < dim);

    constexpr int num_faces = Shape::FaceTraits<ShapeType, face_dim>::count;

    for(int i = 0; i < num_faces; i++)
    {
      int c = Geometry::Intern::congruency<MeshType_, dim, face_dim>(mesh, cell, i);
      TEST_CHECK_NOT_EQUAL(c, -1);
    }
  }

  template<typename MeshType_, int face_dim>
  void check_well_formed(const MeshType_& mesh, Index cell) const
  {
    using ShapeType = typename MeshType_::ShapeType;
    constexpr int dim = ShapeType::dimension;
    static_assert(face_dim < dim);

    using FaceShape = typename Shape::FaceTraits<ShapeType, face_dim>::ShapeType;
    constexpr int face_vertices = Shape::FaceTraits<FaceShape, 0>::count;

    constexpr int num_faces = Shape::FaceTraits<ShapeType, face_dim>::count;
    auto& face_at_cell = mesh.template get_index_set<dim, face_dim>();
    auto& v_at_face = mesh.template get_index_set<face_dim, 0>();

    for(int face_index = 0; face_index < num_faces; face_index++)
    {
      Index face = face_at_cell[cell][face_index];
      for(int i = 0; i < face_vertices; i++)
      {
        for(int j = 0; j < i; j++)
        {
          TEST_CHECK_NOT_EQUAL(v_at_face[face][i], v_at_face[face][j]);
        }
      }
    }
  }

  void check_2d_mesh_consistency(const Mesh2D& mesh) const
  {
    for(Index cell = 0; cell < mesh.get_num_elements(); cell++)
    {
      check_well_formed<Mesh2D, 1>(mesh, cell);
      check_congruency<Mesh2D, 1>(mesh, cell);
    }
  }

  void check_3d_mesh_consistency(const Mesh3D& mesh) const
  {
    for(Index cell = 0; cell < mesh.get_num_elements(); cell++)
    {
      check_well_formed<Mesh3D, 1>(mesh, cell);
      check_congruency<Mesh3D, 1>(mesh, cell);
      check_well_formed<Mesh3D, 2>(mesh, cell);
      check_congruency<Mesh3D, 2>(mesh, cell);
    }
  }

  virtual void run() const override
  {
    test_mesh_creation();
    test_simple();
    test_neighbors();
    test_stats();
    test_export_and_refine();
    test_import();
    test_level_transfer();
    test_mesh_part_projection();
    test_coarse_fine_mapping_type_1();
    test_coarse_fine_mapping_type_15();
    test_adaptive_mesh_2d();
    test_adaptive_mesh_3d();
    test_orientation_handling();
  }

  void test_mesh_creation() const
  {
    Geometry::UnitCubeFactory<Mesh3D> mesh_factory;
    Mesh3D foundation_mesh(mesh_factory);

    Mesh3D refined_mesh =
      AdaptiveMesh3D::create_refined_mesh(foundation_mesh, [](Index) { return 1; });

    TEST_CHECK_EQUAL(refined_mesh.get_num_elements(), 27);
  }

  void test_simple() const
  {
    Geometry::UnitCubeFactory<Mesh3D> mesh_factory;
    Mesh3D base_mesh(mesh_factory);

    Geometry::SubdivisionLevels levels(base_mesh.get_num_vertices(), 0);
    levels[3] = 2;

    AdaptiveMesh3D a_mesh(base_mesh);
    a_mesh.adapt(levels, AdaptiveMesh3D::ImportBehaviour::All);
  }

  void test_neighbors() const
  {
    Geometry::RefinedUnitCubeFactory<Mesh3D> mesh_factory(1);
    Mesh3D base_mesh(mesh_factory);

    Geometry::SubdivisionLevels levels(base_mesh.get_num_vertices());
    auto a_mesh = std::make_shared<AdaptiveMesh3D>(base_mesh);
    a_mesh->adapt(levels, AdaptiveMesh3D::ImportBehaviour::All);

    AdaptiveMeshLayer3D layer(a_mesh, Geometry::Layer{0});

    layer.fill_neighbours();

    auto neighbors = layer.get_neighbors();

    TEST_CHECK_EQUAL(neighbors(0, 1), 4);
    TEST_CHECK_EQUAL(neighbors(0, 3), 2);
    TEST_CHECK_EQUAL(neighbors(0, 5), 1);

    TEST_CHECK_EQUAL(neighbors(7, 0), 3);
    TEST_CHECK_EQUAL(neighbors(7, 2), 5);
    TEST_CHECK_EQUAL(neighbors(7, 4), 6);
  }

  void test_stats() const
  {
    Geometry::UnitCubeFactory<Mesh3D> mesh_factory;
    Mesh3D base_mesh(mesh_factory);

    Geometry::SubdivisionLevels levels(base_mesh.get_num_vertices(), 0);

    AdaptiveMesh3D a_mesh(base_mesh);
    AdaptionStats stats = a_mesh.adapt(levels, AdaptiveMesh3D::ImportBehaviour::All);

    TEST_CHECK_EQUAL(stats.added.at(0), 8);
    TEST_CHECK_EQUAL(stats.added.at(1), 12);
    TEST_CHECK_EQUAL(stats.added.at(2), 6);
    TEST_CHECK_EQUAL(stats.added.at(3), 1);

    TEST_CHECK_EQUAL(stats.kept.at(0), 0);
    TEST_CHECK_EQUAL(stats.kept.at(1), 0);
    TEST_CHECK_EQUAL(stats.kept.at(2), 0);
    TEST_CHECK_EQUAL(stats.kept.at(3), 0);

    TEST_CHECK_EQUAL(stats.removed.at(0), 0);
    TEST_CHECK_EQUAL(stats.removed.at(1), 0);
    TEST_CHECK_EQUAL(stats.removed.at(2), 0);
    TEST_CHECK_EQUAL(stats.removed.at(3), 0);

    for(Index i = 0; i < 8; i++)
    {
      levels[i] = 1;
    }
    stats = a_mesh.adapt(levels, AdaptiveMesh3D::ImportBehaviour::All);

    TEST_CHECK_EQUAL(stats.added.at(0), 56);
    TEST_CHECK_EQUAL(stats.added.at(1), 144 + 12);
    TEST_CHECK_EQUAL(stats.added.at(2), 108 + 6);
    TEST_CHECK_EQUAL(stats.added.at(3), 27 + 1);

    TEST_CHECK_EQUAL(stats.kept.at(0), 8);
    TEST_CHECK_EQUAL(stats.kept.at(1), 0);
    TEST_CHECK_EQUAL(stats.kept.at(2), 0);
    TEST_CHECK_EQUAL(stats.kept.at(3), 0);

    TEST_CHECK_EQUAL(stats.removed.at(0), 0);
    TEST_CHECK_EQUAL(stats.removed.at(1), 12);
    TEST_CHECK_EQUAL(stats.removed.at(2), 6);
    TEST_CHECK_EQUAL(stats.removed.at(3), 1);

    stats = a_mesh.adapt(levels, AdaptiveMesh3D::ImportBehaviour::All);

    TEST_CHECK_EQUAL(stats.added.at(0), 0);
    TEST_CHECK_EQUAL(stats.added.at(1), 0);
    TEST_CHECK_EQUAL(stats.added.at(2), 0);
    TEST_CHECK_EQUAL(stats.added.at(3), 0);

    TEST_CHECK_EQUAL(stats.kept.at(0), 56 + 8);
    TEST_CHECK_EQUAL(stats.kept.at(1), 144 + 12);
    TEST_CHECK_EQUAL(stats.kept.at(2), 108 + 6);
    TEST_CHECK_EQUAL(stats.kept.at(3), 27 + 1);

    TEST_CHECK_EQUAL(stats.removed.at(0), 0);
    TEST_CHECK_EQUAL(stats.removed.at(1), 0);
    TEST_CHECK_EQUAL(stats.removed.at(2), 0);
    TEST_CHECK_EQUAL(stats.removed.at(3), 0);

    for(Index i = 0; i < 8; i++)
    {
      levels[i] = 0;
    }
    stats = a_mesh.adapt(levels, AdaptiveMesh3D::ImportBehaviour::All);

    TEST_CHECK_EQUAL(stats.added.at(0), 0);
    TEST_CHECK_EQUAL(stats.added.at(1), 12);
    TEST_CHECK_EQUAL(stats.added.at(2), 6);
    TEST_CHECK_EQUAL(stats.added.at(3), 1);

    TEST_CHECK_EQUAL(stats.kept.at(0), 8);
    TEST_CHECK_EQUAL(stats.kept.at(1), 0);
    TEST_CHECK_EQUAL(stats.kept.at(2), 0);
    TEST_CHECK_EQUAL(stats.kept.at(3), 0);

    TEST_CHECK_EQUAL(stats.removed.at(0), 56);
    TEST_CHECK_EQUAL(stats.removed.at(1), 144 + 12);
    TEST_CHECK_EQUAL(stats.removed.at(2), 108 + 6);
    TEST_CHECK_EQUAL(stats.removed.at(3), 27 + 1);

    levels[0] = 1;
    stats = a_mesh.adapt(levels, AdaptiveMesh3D::ImportBehaviour::All);

    TEST_CHECK_EQUAL(stats.kept.at(0), 8);
    TEST_CHECK_EQUAL(stats.kept.at(1), 9);
    TEST_CHECK_EQUAL(stats.kept.at(2), 3);
    TEST_CHECK_EQUAL(stats.kept.at(3), 0);
  }

  void test_export_and_refine() const
  {
    Geometry::RefinedUnitCubeFactory<Mesh2D> mesh_factory(2);
    Mesh2D base_mesh(mesh_factory);

    Geometry::SubdivisionLevels levels(base_mesh.get_num_vertices(), 0);
    levels[0] = 2;
    levels[1] = 2;
    levels[3] = 2;

    AdaptiveMesh2D a_mesh(base_mesh);
    a_mesh.adapt(levels, AdaptiveMesh2D::ImportBehaviour::All);

    Index num_faces = a_mesh.get_num_entities(Layer{2}, 2);

    Mesh2D refined_mesh = a_mesh.to_conformal_mesh(Layer{2});
    for(int i = 0; i < 2; ++i)
    {
      Geometry::StandardRefinery<Mesh2D> factory(refined_mesh);
      refined_mesh = factory.make();
    }

    TEST_CHECK_EQUAL(refined_mesh.get_num_entities(2), 4 * 4 * num_faces);
  }

  void test_import() const
  {
    Geometry::RefinedUnitCubeFactory<Mesh2D> mesh_factory(2);
    Mesh2D base_mesh(mesh_factory);

    Geometry::SubdivisionLevels levels(base_mesh.get_num_vertices(), 0);
    auto a_mesh = std::make_shared<AdaptiveMesh2D>(base_mesh);
    a_mesh->adapt(levels, AdaptiveMesh2D::ImportBehaviour::All);

    AdaptiveMeshLayer2D layer(a_mesh, Layer{0});

    TEST_CHECK_EQUAL(layer.get_num_entities(0), base_mesh.get_num_entities(0));
    TEST_CHECK_EQUAL(layer.get_num_entities(1), base_mesh.get_num_entities(1));
    TEST_CHECK_EQUAL(layer.get_num_entities(2), base_mesh.get_num_entities(2));
  }

  void test_level_transfer() const
  {
    UnitCubeFactory<Mesh2D> mesh_factory;
    Mesh2D mesh(mesh_factory);

    AdaptiveMesh2D adaptive_mesh(mesh);

    SubdivisionLevels levels(mesh.get_num_vertices());
    for(Index i(0); i < mesh.get_num_vertices(); ++i)
    {
      levels[i] = 1;
    }
    adaptive_mesh.adapt(levels);

    {
      SubdivisionLevels level_deltas(16);
      // Further refine the 4 inner vertices
      level_deltas[12] = 1;
      level_deltas[13] = 1;
      level_deltas[14] = 1;
      level_deltas[15] = 1;

      // Levels should transfer to regular vertices
      auto transferred = adaptive_mesh.transfer_sdls(level_deltas);
      TEST_CHECK_EQUAL(transferred[0], 1);
      TEST_CHECK_EQUAL(transferred[1], 1);
      TEST_CHECK_EQUAL(transferred[2], 1);
      TEST_CHECK_EQUAL(transferred[3], 1);

      // Mark an edge vertex with a delta of 2
      level_deltas[4] = 2;

      transferred = adaptive_mesh.transfer_sdls(level_deltas);
      TEST_CHECK_EQUAL(transferred[0], 2);
      TEST_CHECK_EQUAL(transferred[1], 1);
      TEST_CHECK_EQUAL(transferred[2], 1);
      TEST_CHECK_EQUAL(transferred[3], 1);
    }

    for(Index i(0); i < mesh.get_num_vertices(); ++i)
    {
      levels[i] = 2;
    }
    adaptive_mesh.adapt(levels);

    {
      SubdivisionLevels level_deltas(adaptive_mesh.get_num_entities(Layer{2}, 0));

      for(Index i(0); i < adaptive_mesh.get_num_entities(Layer{2}, 0); i++)
      {
        level_deltas[i] = 2;
      }

      auto transferred = adaptive_mesh.transfer_sdls(level_deltas);
      for(Index i(0); i < adaptive_mesh.get_num_entities(Layer{0}, 0); i++)
      {
        TEST_CHECK_EQUAL(transferred[i], 2);
      }
    }
  }

  void test_mesh_part_projection() const
  {
    {
      RefinedUnitCubeFactory<Mesh2D> mesh_factory(2);
      Mesh2D mesh(mesh_factory);

      AdaptiveMesh2D adaptive_mesh(mesh);
      SubdivisionLevels levels(mesh.get_num_vertices(), 1);
      adaptive_mesh.adapt(levels);

      BoundaryFactory<Mesh2D> boundary_factory(mesh);
      MeshPart<Mesh2D> reg_boundary(boundary_factory);

      MeshPart<AdaptiveMeshLayer2D> adaptive_boundary =
        adaptive_mesh.template project_meshpart<AdaptiveMeshLayer2D>(Layer{1}, reg_boundary);

      const Index reg_vertices = reg_boundary.template get_target_set<0>().get_num_entities();
      const Index reg_edges = reg_boundary.template get_target_set<1>().get_num_entities();

      const Index adaptive_vertices = adaptive_boundary.template get_target_set<0>().get_num_entities();
      const Index adaptive_edges = adaptive_boundary.template get_target_set<1>().get_num_entities();

      // Each regular edge is split into three.
      // There should be 3 times as many adaptive edges as regular edges
      TEST_CHECK_EQUAL(3 * reg_edges, adaptive_edges);

      // The adaptive boundary contains the same vertices as the regular boundary
      // plus 2 extra vertices for each refined edge
      TEST_CHECK_EQUAL(reg_vertices + 2 * reg_edges, adaptive_vertices);
    }

    {
      RefinedUnitCubeFactory<Mesh3D> mesh_factory(2);
      Mesh3D mesh(mesh_factory);

      AdaptiveMesh3D adaptive_mesh(mesh);
      SubdivisionLevels levels(mesh.get_num_vertices(), 1);
      adaptive_mesh.adapt(levels);

      BoundaryFactory<Mesh3D> boundary_factory(mesh);
      MeshPart<Mesh3D> reg_boundary(boundary_factory);

      MeshPart<AdaptiveMeshLayer3D> adaptive_boundary =
        adaptive_mesh.template project_meshpart<AdaptiveMeshLayer3D>(Layer{1}, reg_boundary);

      const Index reg_vertices = reg_boundary.template get_target_set<0>().get_num_entities();
      const Index reg_edges = reg_boundary.template get_target_set<1>().get_num_entities();
      const Index reg_faces = reg_boundary.template get_target_set<2>().get_num_entities();

      const Index adaptive_vertices = adaptive_boundary.template get_target_set<0>().get_num_entities();
      const Index adaptive_edges = adaptive_boundary.template get_target_set<1>().get_num_entities();
      const Index adaptive_faces = adaptive_boundary.template get_target_set<2>().get_num_entities();
      // Each regular face is split into nine parts
      // There should be 9 times as many adaptive edges as regular edges
      TEST_CHECK_EQUAL(9 * reg_faces, adaptive_faces);

      // Each regular edge is split into three.
      // Each split of a regular face introduces 12 new edges.
      TEST_CHECK_EQUAL(3 * reg_edges + 12 * reg_faces, adaptive_edges);

      // The adaptive boundary contains the same vertices as the regular boundary
      // plus 2 extra vertices for each refined edge
      // plus 4 extra vertices for each refined face
      TEST_CHECK_EQUAL(reg_vertices + 2 * reg_edges + 4 * reg_faces, adaptive_vertices);
    }
  }

  void test_coarse_fine_mapping_type_1() const
  {
    RefinedUnitCubeFactory<Mesh2D> mesh_factory(0);
    Mesh2D mesh(mesh_factory);

    auto adaptive_mesh = std::make_shared<AdaptiveMesh2D>(mesh);
    SubdivisionLevels levels(mesh.get_num_vertices(), 0);
    levels[0] = 3;

    adaptive_mesh->adapt(levels);

    AdaptiveMeshLayer layer_0 = AdaptiveMeshLayer(adaptive_mesh, Layer{0});
    AdaptiveMeshLayer layer_1 = AdaptiveMeshLayer(adaptive_mesh, Layer{1});
    AdaptiveMeshLayer layer_2 = AdaptiveMeshLayer(adaptive_mesh, Layer{2});
    AdaptiveMeshLayer layer_3 = AdaptiveMeshLayer(adaptive_mesh, Layer{3});

    AdaptiveChildMapping cf_mapping_0_1(layer_0, layer_1);
    AdaptiveChildMapping cf_mapping_1_2(layer_1, layer_2);
    AdaptiveChildMapping cf_mapping_2_3(layer_2, layer_3);

    TEST_CHECK_EQUAL(cf_mapping_0_1.get_num_nodes_domain(), 1);
    TEST_CHECK_EQUAL(cf_mapping_0_1.get_num_nodes_image(), 3);

    TEST_CHECK_EQUAL(cf_mapping_1_2.get_num_nodes_domain(), 3);
    TEST_CHECK_EQUAL(cf_mapping_1_2.get_num_nodes_image(), 5);

    TEST_CHECK_EQUAL(cf_mapping_2_3.get_num_nodes_domain(), 5);
    TEST_CHECK_EQUAL(cf_mapping_2_3.get_num_nodes_image(), 7);

    auto check_bijection = [&](AdaptiveChildMapping<AdaptiveMeshLayer2D>& mapping)
    {
      std::vector<Index> image_nodes;
      for(Index c = 0; c < mapping.get_num_nodes_domain(); c++)
      {
        auto it = mapping.image_begin(c);
        while(it != mapping.image_end(c))
        {
          image_nodes.push_back(*it);
          ++it;
        }
      }

      TEST_CHECK_EQUAL(image_nodes.size(), mapping.get_num_nodes_image());

      std::set<Index> seen;

      for(auto idx : image_nodes)
      {
        seen.insert(idx);
      }

      TEST_CHECK_EQUAL(seen.size(), mapping.get_num_nodes_image());
    };

    check_bijection(cf_mapping_0_1);
    check_bijection(cf_mapping_1_2);
    check_bijection(cf_mapping_2_3);
  }

  void test_coarse_fine_mapping_type_15() const
  {
    RefinedUnitCubeFactory<Mesh2D> mesh_factory(0);
    Mesh2D mesh(mesh_factory);

    auto adaptive_mesh = std::make_shared<AdaptiveMesh2D>(mesh);
    SubdivisionLevels levels(mesh.get_num_vertices());

    levels[0] = 1;
    levels[1] = 1;
    levels[2] = 1;
    levels[3] = 1;

    adaptive_mesh->adapt(levels);

    AdaptiveMeshLayer layer_0 = AdaptiveMeshLayer(adaptive_mesh, Layer{0});
    AdaptiveMeshLayer layer_1 = AdaptiveMeshLayer(adaptive_mesh, Layer{1});

    AdaptiveChildMapping cf_mapping(layer_0, layer_1);

    TEST_CHECK_EQUAL(cf_mapping.get_num_nodes_domain(), 1);
    TEST_CHECK_EQUAL(cf_mapping.get_num_nodes_image(), 9);

    std::vector<Index> image_nodes;

    auto it = cf_mapping.image_begin(0);
    while(it != cf_mapping.image_end(0))
    {
      Index value = *it;
      image_nodes.push_back(value);
      ++it;
    }

    for(Index i = 0; i < 9; i++)
    {
      TEST_CHECK(std::find(image_nodes.begin(), image_nodes.end(), i) != image_nodes.end());
    }
  }

  void test_adaptive_mesh_3d() const
  {
    RefinedUnitCubeFactory<Mesh3D> mesh_factory(0);
    Mesh3D mesh(mesh_factory);

    AdaptiveMesh3D adaptive_mesh(mesh);
    SubdivisionLevels levels(mesh.get_num_vertices());

    levels[0] = 2;
    levels[1] = 1;
    levels[2] = 0;
    levels[3] = 0;
    levels[4] = 0;
    levels[5] = 0;
    levels[6] = 0;
    levels[7] = 0;

    adaptive_mesh.adapt(levels);

    TEST_CHECK_EQUAL(adaptive_mesh.template num_total_entities<3>(), 16);
  }

  void test_adaptive_mesh_2d() const
  {
    RefinedUnitCubeFactory<Mesh2D> mesh_factory(0);
    Mesh2D mesh(mesh_factory);

    AdaptiveMesh2D adaptive_mesh(mesh);
    SubdivisionLevels levels(mesh.get_num_vertices());

    levels[0] = 2;
    levels[1] = 1;
    levels[2] = 0;
    levels[3] = 0;

    adaptive_mesh.adapt(levels);

    // Total entity counts should include all objects, including those of base mesh
    TEST_CHECK_EQUAL(adaptive_mesh.template num_total_entities<0>(), 15);
    TEST_CHECK_EQUAL(adaptive_mesh.template num_total_entities<1>(), 28);
    TEST_CHECK_EQUAL(adaptive_mesh.template num_total_entities<2>(), 11);

    // Level entity counts should only include elements of that level
    TEST_CHECK_EQUAL(adaptive_mesh.get_num_entities(Layer{2}, 0), 15);
    TEST_CHECK_EQUAL(adaptive_mesh.get_num_entities(Layer{2}, 1), 23);
    TEST_CHECK_EQUAL(adaptive_mesh.get_num_entities(Layer{2}, 2), 9);

    // Writing to mesh should work
    Mesh2D refined_mesh = adaptive_mesh.to_conformal_mesh(Layer{2});

    check_2d_mesh_consistency(refined_mesh);

    levels[0] = 1;
    levels[1] = 1;
    levels[2] = 1;
    levels[3] = 0;
    adaptive_mesh.adapt(levels);

    TEST_CHECK_EQUAL(adaptive_mesh.template num_total_entities<0>(), 14);
    TEST_CHECK_EQUAL(adaptive_mesh.template num_total_entities<1>(), 25);
    TEST_CHECK_EQUAL(adaptive_mesh.template num_total_entities<2>(), 9);

    TEST_CHECK_EQUAL(adaptive_mesh.get_num_entities(Layer{1}, 0), 14);
    TEST_CHECK_EQUAL(adaptive_mesh.get_num_entities(Layer{1}, 1), 21);
    TEST_CHECK_EQUAL(adaptive_mesh.get_num_entities(Layer{1}, 2), 8);

    refined_mesh = adaptive_mesh.to_conformal_mesh(Layer{1});

    check_2d_mesh_consistency(refined_mesh);
  }

  static void test_orientation_handling()
  {
    RefinedUnitCubeFactory<Mesh3D> mesh_factory(1);
    Mesh3D mesh(mesh_factory);
    auto& v_at_c = mesh.template get_index_set<3, 0>();

    {
      // Swap orientation of some edges
      auto& v_at_e = mesh.template get_index_set<1, 0>();
      auto& e_at_c = mesh.template get_index_set<3, 1>();
      std::swap(v_at_e[e_at_c[0][0]][0], v_at_e[e_at_c[0][0]][1]);
      std::swap(v_at_e[e_at_c[0][3]][0], v_at_e[e_at_c[0][3]][1]);
      std::swap(v_at_e[e_at_c[0][4]][0], v_at_e[e_at_c[0][4]][1]);
      std::swap(v_at_e[e_at_c[0][7]][0], v_at_e[e_at_c[0][7]][1]);
      std::swap(v_at_e[e_at_c[0][9]][0], v_at_e[e_at_c[0][9]][1]);
      std::swap(v_at_e[e_at_c[0][10]][0], v_at_e[e_at_c[0][10]][1]);
    }

    {
      // Swap orientation of some faces
      auto& v_at_f = mesh.template get_index_set<2, 0>();
      auto& e_at_f = mesh.template get_index_set<2, 1>();
      auto& f_at_c = mesh.template get_index_set<3, 2>();
      Index a, b, c, d, face;

      // Produce an orientation 2 face
      face = f_at_c[0][0];
      a = v_at_f[face][0];
      b = v_at_f[face][1];
      c = v_at_f[face][2];
      d = v_at_f[face][3];
      v_at_f[face][0] = b;
      v_at_f[face][1] = d;
      v_at_f[face][2] = a;
      v_at_f[face][3] = c;
      a = e_at_f[face][0];
      b = e_at_f[face][1];
      c = e_at_f[face][2];
      d = e_at_f[face][3];
      e_at_f[face][0] = d;
      e_at_f[face][1] = c;
      e_at_f[face][2] = a;
      e_at_f[face][3] = b;

      // Produce an orientation 5 face
      face = f_at_c[0][3];
      a = v_at_f[face][0];
      b = v_at_f[face][1];
      c = v_at_f[face][2];
      d = v_at_f[face][3];
      v_at_f[face][0] = b;
      v_at_f[face][1] = a;
      v_at_f[face][2] = d;
      v_at_f[face][3] = c;
      a = e_at_f[face][0];
      b = e_at_f[face][1];
      c = e_at_f[face][2];
      d = e_at_f[face][3];
      e_at_f[face][0] = a;
      e_at_f[face][1] = b;
      e_at_f[face][2] = d;
      e_at_f[face][3] = c;
    }

    // Adapt for all possible 3D template types, ensure no assertions are hit
    for(std::uint64_t type = 0; type < 256; type++)
    {
      SubdivisionLevels levels(mesh.get_num_vertices());

      for(int i = 0; i < 8; i++)
      {
        levels[v_at_c[0][i]] = (type & (1ULL << (std::uint64_t)i)) > 0 ? 1 : 0;
      }

      AdaptiveMesh3D adaptive_mesh(mesh);
      adaptive_mesh.adapt(levels);
    }
  }
};

static const AdaptiveRefinementTest<float> adaptive_refinement_test_float;
static const AdaptiveRefinementTest<double> adaptive_refinement_test_double;

/**
 * \brief Unit tests for the mock mesh interface
 *
 * \author Markus Muegge
 */
template<typename DT_>
class AdaptiveMeshLayerTest : public UnitTest
{
public:
  using ShapeType = Shape::Hexahedron;
  using TemplateSet = SchneidersTemplates;

  using FoundationMeshType = ConformalMesh<ShapeType, 3, DT_>;
  using AdaptiveMeshType = AdaptiveMesh<TemplateSet, ShapeType, 3, DT_>;
  using AdaptiveMeshLayerType = AdaptiveMeshLayer<AdaptiveMeshType>;

  AdaptiveMeshLayerTest() : UnitTest("AdaptiveMeshLayerTest", Type::Traits<DT_>::name())
  {
  }

  ~AdaptiveMeshLayerTest() override = default;

  void run() const override
  {
    Geometry::UnitCubeFactory<FoundationMeshType> factory;
    FoundationMeshType foundation_mesh(factory);
    auto adaptive_mesh = std::make_shared<AdaptiveMeshType>(foundation_mesh);

    Geometry::SubdivisionLevels sdls(8, 2);
    adaptive_mesh->adapt(sdls, AdaptiveMeshType::ImportBehaviour::All);

    AdaptiveMeshLayerType layer_one(adaptive_mesh, Geometry::Layer{1});
    AdaptiveMeshLayerType layer_two(adaptive_mesh, Geometry::Layer{2});

    // Test get_num_entities
    TEST_CHECK_EQUAL(layer_one.get_num_entities(0), 4 * 4 * 4);
    TEST_CHECK_EQUAL(layer_two.get_num_entities(0), 10 * 10 * 10);

    TEST_CHECK_EQUAL(layer_one.get_num_entities(3), 27);
    TEST_CHECK_EQUAL(layer_two.get_num_entities(3), 27 * 27);

    // Test vertex set
    auto& vertex_set_one = layer_one.get_vertex_set();
    auto& vertex_set_two = layer_two.get_vertex_set();

    TEST_CHECK_EQUAL(vertex_set_one.num_coords, 3);
    TEST_CHECK_EQUAL(vertex_set_two.num_coords, 3);

    TEST_CHECK_EQUAL(vertex_set_one.get_num_vertices(), 4 * 4 * 4);
    TEST_CHECK_EQUAL(vertex_set_two.get_num_vertices(), 10 * 10 * 10);

    // Since vertices are permanent, vertex sets should return same vertex at index 0
    auto& vertex_one = vertex_set_one[0];
    auto& vertex_two = vertex_set_two[0];

    TEST_CHECK_EQUAL(vertex_one[0], vertex_two[0]);
    TEST_CHECK_EQUAL(vertex_one[1], vertex_two[1]);
    TEST_CHECK_EQUAL(vertex_one[2], vertex_two[2]);

    // Test index set holder
    auto& index_set_holder = layer_one.get_index_set_holder();

    // Test index set
    const auto& v_at_c = index_set_holder.template get_index_set<3, 0>();
    const auto& e_at_c = index_set_holder.template get_index_set<3, 1>();
    const auto& f_at_c = index_set_holder.template get_index_set<3, 2>();

    TEST_CHECK_EQUAL(v_at_c.num_indices, 8);
    TEST_CHECK_EQUAL(v_at_c.get_num_entities(), 27);
    TEST_CHECK_EQUAL(v_at_c.get_index_bound(), 64);

    TEST_CHECK_EQUAL(e_at_c.num_indices, 12);
    TEST_CHECK_EQUAL(e_at_c.get_num_entities(), 27);
    TEST_CHECK_EQUAL(e_at_c.get_index_bound(), 144);

    TEST_CHECK_EQUAL(f_at_c.num_indices, 6);
    TEST_CHECK_EQUAL(f_at_c.get_num_entities(), 27);
    TEST_CHECK_EQUAL(f_at_c.get_index_bound(), 108);

    // Indices should be identical across the adaptive mesh, index sets, and index tuples
    for(Index i(0); i < v_at_c.get_num_entities(); i++)
    {
      auto vertices = v_at_c[i];
      auto edges = e_at_c[i];
      auto faces = f_at_c[i];

      for(int j(0); j < v_at_c.num_indices; j++)
      {
        const Index a = adaptive_mesh->template get_face_index<3, 0>(Geometry::Layer{1}, i, j);
        const Index b = v_at_c(i, j);
        const Index c = vertices[j];
        TEST_CHECK_EQUAL(a, b);
        TEST_CHECK_EQUAL(b, c);
      }
      for(int j(0); j < e_at_c.num_indices; j++)
      {
        const Index a = adaptive_mesh->template get_face_index<3, 1>(Geometry::Layer{1}, i, j);
        const Index b = e_at_c(i, j);
        const Index c = edges[j];
        TEST_CHECK_EQUAL(a, b);
        TEST_CHECK_EQUAL(b, c);
      }
      for(int j(0); j < f_at_c.num_indices; j++)
      {
        const Index a = adaptive_mesh->template get_face_index<3, 2>(Geometry::Layer{1}, i, j);
        const Index b = f_at_c(i, j);
        const Index c = faces[j];
        TEST_CHECK_EQUAL(a, b);
        TEST_CHECK_EQUAL(b, c);
      }
    }
  }
};

static const AdaptiveMeshLayerTest<float> adaptive_mesh_layer_test_float;
static const AdaptiveMeshLayerTest<double> adaptive_mesh_layer_test_double;

/// Unit-test class for testing the StandardTemplateSet and StandardRefinementType
class StandardTemplateSetTest : public UnitTest
{
public:
  StandardTemplateSetTest() : UnitTest("StandardTemplateSetTest")
  {
  }

  ~StandardTemplateSetTest() override = default;

  void run() const override
  {
    test_refinement_type_calculation();
    test_subdivisionlevel_spreading();
  }

  void test_refinement_type_calculation() const
  {
    // No bits set for type 0
    Geometry::Intern::RefinementFieldTuple<std::uint64_t, 4> tuple;

    tuple = {0, 0, 0, 0};
    TEST_CHECK_EQUAL(StandardRefinementType<Shape::Quadrilateral>(tuple).to_number(), 0b0000);
    // Order of bits is inverted, due to lsb being on the right
    tuple = {0, 0, 0, 1};
    TEST_CHECK_EQUAL(StandardRefinementType<Shape::Quadrilateral>(tuple).to_number(), 0b1000);
    // Anything above 0 should be set
    tuple = {2, 3, 0, 1};
    TEST_CHECK_EQUAL(StandardRefinementType<Shape::Quadrilateral>(tuple).to_number(), 0b1011);
    // All bits can be set
    tuple = {1, 1, 1, 1};
    TEST_CHECK_EQUAL(StandardRefinementType<Shape::Quadrilateral>(tuple).to_number(), 0b1111);
  }

  void test_subdivisionlevel_spreading() const
  {
    Geometry::Intern::RefinementFieldTuple<std::uint64_t, 8> ref_field_tuple{2, 2, 2, 2, 2, 2, 0, 0};

    // Topology vertices should be decremented
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::ParentTopology, 0, 0, 0),
        ref_field_tuple),
      1);

    // Topology vertices should not be decremented below 0
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::ParentTopology, 7, 0, 0),
        ref_field_tuple),
      0);

    // Siblings should receive minimum of decremented subdivision levels
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::Sibling, 0, 0, 0),
        ref_field_tuple),
      0);

    // Boundary edge vertices should receive minimum of decremented adjacent vertices
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::BoundaryEdge, 0, 0, 0),
        ref_field_tuple),
      1);
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::BoundaryEdge, 0, 0, 3),
        ref_field_tuple),
      0);

    // Boundary face vertices should receive minimum of decremented adjacent vertices
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::BoundaryFace, 0, 0, 0),
        ref_field_tuple),
      1);
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::BoundaryFace, 0, 0, 1),
        ref_field_tuple),
      0);
  }
};

static const StandardTemplateSetTest standard_template_set_test;

/// Unit-test class for testing the IsolatedPointTemplateSet and IsolatedPointRefinementType
class IsolatedPointTemplateSetTest : public UnitTest
{
public:
  IsolatedPointTemplateSetTest() : UnitTest("StandardTemplateSetTest")
  {
  }

  ~IsolatedPointTemplateSetTest() override = default;

  void run() const override
  {
    test_refinement_type_calculation();
    test_subdivisionlevel_spreading();
  }

  void test_refinement_type_calculation() const
  {
    using Marking = Geometry::IsolatedPointVertexMarking;
    // No bits set for type 0
    Geometry::Intern::RefinementFieldTuple<Marking, 4> tuple;

    tuple = {Marking{0, false}, Marking{0, false}, Marking{0, false}, Marking{0, false}};
    TEST_CHECK_EQUAL(IsolatedPointRefinementType<Shape::Quadrilateral>(tuple).to_number(), 0b0000);
    // Order of bits is inverted, due to lsb being on the right
    tuple = {Marking{0, false}, Marking{0, false}, Marking{0, false}, Marking{1, false}};
    TEST_CHECK_EQUAL(IsolatedPointRefinementType<Shape::Quadrilateral>(tuple).to_number(), 0b1000);
    // Anything above 0 should be set
    tuple = {Marking{2, false}, Marking{3, false}, Marking{0, false}, Marking{1, false}};
    TEST_CHECK_EQUAL(IsolatedPointRefinementType<Shape::Quadrilateral>(tuple).to_number(), 0b1011);
    // All bits can be set
    tuple = {Marking{1, false}, Marking{1, false}, Marking{1, false}, Marking{1, false}};
    TEST_CHECK_EQUAL(IsolatedPointRefinementType<Shape::Quadrilateral>(tuple).to_number(), 0b1111);
  }

  void test_subdivisionlevel_spreading() const
  {
    using Marking = Geometry::IsolatedPointVertexMarking;

    Geometry::Intern::RefinementFieldTuple<Marking, 8> ref_field_tuple{
      Marking{2, false},
      Marking{2, false},
      Marking{2, false},
      Marking{2, false},
      Marking{2, false},
      Marking{2, false},
      Marking{0, false},
      Marking{0, false}};

    // Topology vertices should be decremented
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::ParentTopology, 0, 0, 0),
        ref_field_tuple),
      Marking(1, false));

    // Topology vertices should not be decremented below 0
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::ParentTopology, 7, 0, 0),
        ref_field_tuple),
      Marking(0, false));

    // Siblings should receive minimum of decremented subdivision levels
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::Sibling, 0, 0, 0),
        ref_field_tuple),
      Marking(0, false));

    // Boundary edge vertices should receive minimum of decremented adjacent vertices
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::BoundaryEdge, 0, 0, 0),
        ref_field_tuple),
      Marking(1, false));
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::BoundaryEdge, 0, 0, 3),
        ref_field_tuple),
      Marking(0, false));

    // Boundary face vertices should receive minimum of decremented adjacent vertices
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::BoundaryFace, 0, 0, 0),
        ref_field_tuple),
      Marking(1, false));
    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::BoundaryFace, 0, 0, 1),
        ref_field_tuple),
      Marking(0, false));

    Geometry::Intern::RefinementFieldTuple<Marking, 8> isolated_ref_field_tuple{
      Marking{2, true},
      Marking{0, false},
      Marking{0, false},
      Marking{0, false},
      Marking{0, false},
      Marking{0, false},
      Marking{0, false},
      Marking{0, false}};

    TEST_CHECK_EQUAL(
      Geometry::Intern::spread_refinement_field<Shape::Hexahedron>(
        EntityReference(EntitySource::ParentTopology, 0, 0, 0),
        isolated_ref_field_tuple),
      Marking(1, true));
  }
};

static const IsolatedPointTemplateSetTest isolated_point_template_set_test;
