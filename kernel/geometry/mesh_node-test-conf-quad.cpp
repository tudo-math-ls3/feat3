#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/tetris_quad.hpp>
#include <kernel/geometry/test_aux/validate_neighbours.hpp>
#include <kernel/geometry/mesh_node.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;
using namespace FEAST::Geometry::TestAux;

typedef ConformalMesh<Shape::Quadrilateral> RootMeshType;
typedef MeshPart<RootMeshType> MeshPartType;
typedef RootMeshNode<RootMeshType> RootMeshNodeType;
typedef MeshPartNode<RootMeshType> MeshPartNodeType;

/**
 * \brief Test class for the MeshNode class template.
 *
 * \test Tests the MeshNode class templates.
 *
 * \author Peter Zajac
 */
class MeshNodeTestConfQuad
  : public TestSystem::TaggedTest<Archs::None, Archs::None>
{
public:
  MeshNodeTestConfQuad() :
    TestSystem::TaggedTest<Archs::None, Archs::None>("mesh_node-test-conf-quad")
  {
  }

  virtual void run() const override
  {
    refine_tree_test();
  }

  void refine_tree_test() const
  {
    // create coarse root mesh
    RootMeshType* root_mesh_c = create_tetris_mesh_2d();

    // create coarse quad submesh
    MeshPartType* submesh_quad_c = create_tetris_quad_submesh_2d();

    // create coarse edge submesh
    MeshPartType* submesh_edge_c = create_tetris_edge_submesh_2d();

    // create coarse edge submesh
    MeshPartType* submesh_quad_edge_c = create_tetris_quad_edge_submesh_2d();

    // create coarse cell subset
    MeshPartType* subset_quad_c = create_tetris_quad_cellsubset_2d();

    // create coarse cell subset
    MeshPartType* subset_quad_edge_c = create_tetris_quad_edge_cellsubset_2d();

    /* ********************************************************************* */

    // create a root mesh node
    RootMeshNodeType* root_node_c = new RootMeshNodeType(root_mesh_c);

    // add quad submesh node
    MeshPartNodeType* subnode_quad_c =
      root_node_c->add_mesh_part_node("0", new MeshPartNodeType(submesh_quad_c));

    // add edge submesh node
    //MeshPartNodeType* subnode_edge_c =
      root_node_c->add_mesh_part_node("1", new MeshPartNodeType(submesh_edge_c));

    // add edge sub-submesh node
    //MeshPartNodeType* subnode_quad_edge_c =
      subnode_quad_c->add_mesh_part_node("2", new MeshPartNodeType(submesh_quad_edge_c));

    // add quad cell subset node
    MeshPartNodeType* subsetnode_quad_c =
      root_node_c->add_mesh_part_node("7", new MeshPartNodeType(subset_quad_c));

    // add edge cell subset node
    //MeshPartNodeType* subsetnode_quad_edge_c =
      subsetnode_quad_c->add_mesh_part_node("42", new MeshPartNodeType(subset_quad_edge_c));

    /* ********************************************************************* */

    // refine the mesh tree
    RootMeshNodeType* root_node_f = root_node_c->refine();

    /* ********************************************************************* */

    // fetch quad submesh node
    MeshPartNodeType* subnode_quad_f = root_node_f->find_mesh_part_node("0");
    TEST_CHECK_MSG(subnode_quad_f != nullptr, "failed to fetch quad submesh node");

    // fetch edge submesh node
    MeshPartNodeType* subnode_edge_f = root_node_f->find_mesh_part_node("1");
    TEST_CHECK_MSG(subnode_edge_f != nullptr, "failed to fetch edge submesh node");

    // fetch quad-edge submesh
    MeshPartNodeType* subnode_quad_edge_f = subnode_quad_f->find_mesh_part_node("2");
    TEST_CHECK_MSG(subnode_quad_edge_f != nullptr, "failed to fetch quad-edge submesh node");

    // fetch quad cell subset node
    MeshPartNodeType* subsetnode_quad_f = root_node_f->find_mesh_part_node("7");
    TEST_CHECK_MSG(subsetnode_quad_f != nullptr, "failed to fetch quad cell subset node");

    // fetch quad-edge cell subset node
    MeshPartNodeType* subsetnode_quad_edge_f = subsetnode_quad_f->find_mesh_part_node("42");
    TEST_CHECK_MSG(subsetnode_quad_edge_f != nullptr, "failed to fetch quad-edge cell subset node");

    /* ********************************************************************* */

    try
    {
      // validate refined root mesh
      validate_refined_tetris_mesh_2d(*root_node_f->get_mesh());

      // validate refined quad submesh
      validate_refined_tetris_quad_submesh_2d(*subnode_quad_f->get_mesh());

      // validate refined edge submesh
      validate_refined_tetris_edge_submesh_2d(*subnode_edge_f->get_mesh());

      // validate refined quad-edge submesh
      validate_refined_tetris_quad_edge_submesh_2d(*subnode_quad_edge_f->get_mesh());

      // validate refined quad cell subset
      validate_refined_tetris_quad_cellsubset_2d(*subsetnode_quad_f->get_mesh());

      // validate refined quad-edge cell subset
      validate_refined_tetris_quad_edge_cellsubset_2d(*subsetnode_quad_edge_f->get_mesh());

      // validate neighbour information
      validate_neighbours(*(root_node_f->get_mesh()));

    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }

    // delete refined and coarse root nodes
    delete root_node_f;
    delete root_node_c;
  }

} mesh_node_test_conf_quad;
