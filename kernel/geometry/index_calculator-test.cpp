#include <test_system/test_system.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/test_aux/index_calculator_meshes.hpp>
#include <kernel/geometry/test_aux/tetris_hexa.hpp>
#include <kernel/geometry/test_aux/standard_tetra.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;
using namespace FEAT::Geometry::TestAux;

/**
 * \brief Test class for the IndexCalculator class template and the IndexTree class template.
 *
 * \test Tests the IndexCalculator class template an the IndexTree class template.
 *
 * \author Constantin Christof
 */

class IndexCalculatorTest
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  IndexCalculatorTest() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("index_calculator-test")
  {
  }

  virtual ~IndexCalculatorTest()
  {
  }

  // run the tests
  virtual void run() const override
  {
    hexa_index_tree_test();
    tetra_index_tree_test();
    hexa_index_calculator_test();
    tetra_index_calculator_test();
  }

  void hexa_index_tree_test() const
  {

    try
    {
      // typedef index set type
      typedef IndexSet<2> IndexSetTypeEV;
      typedef IndexSet<4> IndexSetTypeQV;

      // create mesh and index-trees
      HexaMesh* mesh = create_tetris_mesh_3d();
      HEdgeIndexTree* edge_tree;
      QuadIndexTree* quad_tree;

      // fetch the quad-vertex- and the edge-vertex-index set
      const IndexSetTypeEV& index_set_e_v = mesh->get_index_set<1,0>();
      const IndexSetTypeQV& index_set_q_v = mesh->get_index_set<2,0>();

      // create index-tree objects
      Index vertex_number = mesh->get_num_entities(0);
      edge_tree = new HEdgeIndexTree(vertex_number);
      quad_tree = new QuadIndexTree(vertex_number);

      // parsing
      edge_tree->parse(index_set_e_v);
      quad_tree->parse(index_set_q_v);

      // check if everything is right
      validate_hypercube_edge_index_set(*edge_tree);
      validate_hypercube_quad_index_set(*quad_tree);

      // clean up
      delete mesh;
      delete edge_tree;
      delete quad_tree;
    }

    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  } // hexa_index_tree_test()

  void tetra_index_tree_test() const
  {

    try
    {
      // typedef index set type
      typedef IndexSet<2> IndexSetTypeEV;
      typedef IndexSet<3> IndexSetTypeTV;

      // create mesh
      TetraMesh* mesh = create_big_tetra_mesh_3d();
      SEdgeIndexTree* edge_tree = nullptr;
      TriaIndexTree* tria_tree = nullptr;

      // fetch the triangle-vertex- and the edge-vertex index set
      const IndexSetTypeEV& index_set_e_v = mesh->get_index_set<1,0>();
      const IndexSetTypeTV& index_set_t_v = mesh->get_index_set<2,0>();

      // create index-tree objects
      Index vertex_number = mesh->get_num_entities(0);
      edge_tree = new SEdgeIndexTree(vertex_number);
      tria_tree = new TriaIndexTree(vertex_number);

      // parsing
      edge_tree->parse(index_set_e_v);
      tria_tree->parse(index_set_t_v);

      // checking
      validate_simplex_edge_index_set(*edge_tree);
      validate_simplex_triangle_index_set(*tria_tree);

      // clean up
      delete mesh;
      delete edge_tree;
      delete tria_tree;
    }

    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  } // simplex_index_tree_test()

  void hexa_index_calculator_test() const
  {
    // typedef calculators
    typedef IndexCalculator<Shape::Hypercube<2>, 1> Quad_Edge_Calc;
    typedef IndexCalculator<Shape::Hypercube<3>, 1> Cube_Edge_Calc;
    typedef IndexCalculator<Shape::Hypercube<3>, 2> Cube_Quad_Calc;

    try
    {
      // typedef index set types
      typedef IndexSet<2> IndexSetTypeEV;

      typedef IndexSet<4> IndexSetTypeQV;
      typedef IndexSet<4> IndexSetTypeQE;

      typedef IndexSet<8> IndexSetTypeCV;
      typedef IndexSet<12> IndexSetTypeCE;

      typedef IndexSet<6> IndexSetTypeCQ;

      // create mesh
      HexaMesh* mesh = create_big_tetris_mesh_3d();
      HEdgeIndexTree* edge_tree = nullptr;
      QuadIndexTree* quad_tree = nullptr;

      // fetch the index sets
      const IndexSetTypeEV& index_set_e_v = mesh->get_index_set<1,0>();

      const IndexSetTypeQV& index_set_q_v = mesh->get_index_set<2,0>();
      IndexSetTypeQE& index_set_q_e = mesh->get_index_set<2,1>();

      const IndexSetTypeCV& index_set_c_v = mesh->get_index_set<3,0>();
      IndexSetTypeCE& index_set_c_e = mesh->get_index_set<3,1>();

      IndexSetTypeCQ& index_set_c_q = mesh->get_index_set<3,2>();

      // create index-tree objects
      Index vertex_number = mesh->get_num_entities(0);
      edge_tree = new HEdgeIndexTree(vertex_number);
      quad_tree = new QuadIndexTree(vertex_number);

      // parsing
      edge_tree->parse(index_set_e_v);
      quad_tree->parse(index_set_q_v);

      if(!Quad_Edge_Calc::compute(*edge_tree, index_set_q_v, index_set_q_e))
        throw String("Failed to compute edges-at-quad index set");
      if(!Cube_Edge_Calc::compute(*edge_tree, index_set_c_v, index_set_c_e))
        throw String("Failed to compute edges-at-hexa index set");
      if(!Cube_Quad_Calc::compute(*quad_tree, index_set_c_v, index_set_c_q))
        throw String("Failed to compute quads-at-hexa index set");

      // check if everything is right
      validate_refined_tetris_mesh_3d(*mesh);

      // clean up
      delete mesh;
      delete edge_tree;
      delete quad_tree;
    }

    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  } // hexa_index_calculator_test()

  void tetra_index_calculator_test() const
  {
    typedef IndexCalculator<Shape::Simplex<2>, 1> Tria_Edge_Calc;
    typedef IndexCalculator<Shape::Simplex<3>, 1> Tetra_Edge_Calc;
    typedef IndexCalculator<Shape::Simplex<3>, 2> Tetra_Tria_Calc;

    try
    {

      // typedef index set types
      typedef IndexSet<2> IndexSetTypeEV;

      typedef IndexSet<3> IndexSetTypeTV;
      typedef IndexSet<3> IndexSetTypeTE;

      typedef IndexSet<4> IndexSetTypeSV;
      typedef IndexSet<6> IndexSetTypeSE;

      typedef IndexSet<4> IndexSetTypeST;

      // create mesh
      TetraMesh* mesh = create_really_big_tetra_mesh_3d();
      SEdgeIndexTree* edge_tree = nullptr;
      TriaIndexTree* tria_tree = nullptr;

      // fetch the index sets
      const IndexSetTypeEV& index_set_e_v = mesh->get_index_set<1,0>();

      const IndexSetTypeTV& index_set_t_v = mesh->get_index_set<2,0>();
      IndexSetTypeTE& index_set_t_e = mesh->get_index_set<2,1>();

      const IndexSetTypeSV& index_set_s_v = mesh->get_index_set<3,0>();
      IndexSetTypeSE& index_set_s_e = mesh->get_index_set<3,1>();

      IndexSetTypeST& index_set_s_t = mesh->get_index_set<3,2>();

      // create index-calculator objects
      Index vertex_number = mesh->get_num_entities(0);
      edge_tree = new SEdgeIndexTree(vertex_number);
      tria_tree = new TriaIndexTree(vertex_number);

      // parsing
      edge_tree->parse(index_set_e_v);
      tria_tree->parse(index_set_t_v);

      // compute index sets
      if(!Tria_Edge_Calc::compute(*edge_tree, index_set_t_v, index_set_t_e))
        throw String("Failed to compute edges-at-tria index set");
      if(!Tetra_Edge_Calc::compute(*edge_tree, index_set_s_v, index_set_s_e))
        throw String("Failed to compute edges-at-tetra index set");
      if(!Tetra_Tria_Calc::compute(*tria_tree, index_set_s_v, index_set_s_t))
        throw String("Failed to compute trias-at-tetra index set");

      // check if everything is right
      validate_refined_big_tetra_mesh_3d(*mesh);

      // clean up
      delete mesh;
      delete edge_tree;
      delete tria_tree;
    }

    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  } // tetra_index_calculator_test()

} index_calculator_test;

/**
 * \brief Class for IndexCalculator::compute tests
 *
 * \test Tests calculation of vertex@subshape information from vertex@shape information.
 *
 * \tparam Shape_
 * Shape for the mesh cells
 *
 * \tparam ShapeDim
 * Dimension of the mesh cells, this determines the the vertex@shape information
 *
 * \tparam SubshapeDim
 * Dimension of the subshape for which the vertex@subshape information is generated and checked
 *
 * \author Jordi Paul
 */
template<template<int> class Shape_, int ShapeDim, int SubshapeDim>
class IndexCalculatorVertexTest
: public TestSystem::TaggedTest<Archs::None, Archs::None>
{
  public:
    /// The complete shape type for the mesh cells
    typedef Shape_<ShapeDim> ShapeType;
    /// Type for the subshapes
    typedef Shape_<SubshapeDim> SubshapeType;

    IndexCalculatorVertexTest() :
      TestSystem::TaggedTest<Archs::None, Archs::None>("index_calculator_vertex-test")
      {
      }

    // run the tests
    virtual void run() const override
    {
      typedef ConformalMesh<ShapeType> MeshType;

      typedef IndexTree<SubshapeType> IndexTreeType;
      typedef typename IndexTreeType::IndexVector IndexVectorType;

      // vertex@shape IndexSet
      typedef IndexSet<Shape::FaceTraits<ShapeType,0>::count> VertAtShapeIndexSetType;
      // vertex@subshape IndexSet
      typedef IndexSet<Shape::FaceTraits<SubshapeType,0>::count> VertAtSubshapeIndexSetType;

      // Generate a mesh of specified ShapeType
      RefineFactory< MeshType, UnitCubeFactory> my_factory(1);
      MeshType my_mesh(my_factory);

      // Now parse the mesh's original IndexSet into an IndexTree
      IndexTreeType original_tree(my_mesh.get_num_entities(0));
      original_tree.parse(my_mesh.template get_index_set<SubshapeDim,0>());

      // This is the IndexSet providing the vertex@shape information
      VertAtShapeIndexSetType my_vert_at_cell_index_set = my_mesh.template get_index_set<ShapeDim,0>();

      // Calculate vertex@subshape information from the VertAtShapeIndexSet
      VertAtSubshapeIndexSetType my_vert_at_subshape_index_set;
      IndexCalculator<ShapeType, SubshapeDim>::compute_vertex_subshape
        (my_vert_at_cell_index_set, my_vert_at_subshape_index_set);

      // Temporary object for passing to find()
      IndexVectorType indices_at_edge;

      // Check if every subshape in the new IndexSet was present in the original
      for(Index i(0); i < my_vert_at_subshape_index_set.get_num_entities(); ++i)
      {
        for(Index j(0); j < Index(IndexTreeType::num_indices); ++j)
          indices_at_edge[j] = my_vert_at_subshape_index_set[i][j];

        std::pair<bool, Index> bi = original_tree.find(indices_at_edge);
        TEST_CHECK_MSG(bi.first,"New subshape not found in orginal IndexSet");
      }

      // Parse the new IndexSet into an IndexTree
      IndexTreeType new_tree(my_vert_at_subshape_index_set.get_num_entities());
      new_tree.parse(my_vert_at_subshape_index_set);
      // Original vert@subshape information
      VertAtSubshapeIndexSetType original_vert_at_subshape_index_set(my_mesh.template get_index_set<SubshapeDim,0>());

      // Check if every subshape in the original IndexSet ist present in the new IndexSet
      for(Index i(0); i < original_vert_at_subshape_index_set.get_num_entities(); ++i)
      {
        for(Index j(0); j < Index(IndexTreeType::num_indices); ++j)
          indices_at_edge[j] = original_vert_at_subshape_index_set[i][j];

        std::pair<bool, Index> bi = new_tree.find(indices_at_edge);
        TEST_CHECK_MSG(bi.first,"Original subshape not found in new IndexSet");
      }
    } // run
};

IndexCalculatorVertexTest<Shape::Hypercube, 2, 1> test_hypercube_2_1;
IndexCalculatorVertexTest<Shape::Hypercube, 3, 1> test_hypercube_3_1;
IndexCalculatorVertexTest<Shape::Hypercube, 3, 2> test_hypercube_3_2;

IndexCalculatorVertexTest<Shape::Simplex, 2, 1> test_simplex_2_1;
IndexCalculatorVertexTest<Shape::Simplex, 3, 1> test_simplex_3_1;
IndexCalculatorVertexTest<Shape::Simplex, 3, 2> test_simplex_3_2;
