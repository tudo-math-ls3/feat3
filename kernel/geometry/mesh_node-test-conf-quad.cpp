#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/tetris_quad.hpp>
#include <kernel/geometry/mesh_node.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Geometry;


typedef StandardConformalMeshNodePolicy<Shape::Quadrilateral> MeshNodePolicy;
typedef MeshNodePolicy::RootMeshType RootMeshType;
typedef MeshNodePolicy::SubMeshType SubMeshType;
typedef RootMeshNode<MeshNodePolicy> RootMeshNodeType;
typedef SubMeshNode<MeshNodePolicy> SubMeshNodeType;
typedef MeshNodePolicy::CellSubSetType CellSubSetType;
typedef CellSubSetNode<MeshNodePolicy> CellSubSetNodeType;

/**
 * \brief Test class for the MeshNode class template.
 *
 * \test Tests the MeshNode class templates.
 *
 * \author Peter Zajac
 */
class MeshNodeTestConfQuad
  : public TestSystem::TaggedTest<Archs::None, Nil>
{
public:
  MeshNodeTestConfQuad() :
    TestSystem::TaggedTest<Archs::None, Nil>("mesh_node-test-conf-quad")
  {
  }

  virtual void run() const
  {
    refine_tree_test();
  }

  void refine_tree_test() const
  {
    // create coarse root mesh
    RootMeshType* root_mesh_c = TestAux::create_tetris_mesh_2d();

    // create coarse quad submesh
    SubMeshType* submesh_quad_c = TestAux::create_tetris_quad_submesh_2d();

    // create coarse edge submesh
    SubMeshType* submesh_edge_c = TestAux::create_tetris_edge_submesh_2d();

    // create coarse edge submesh
    SubMeshType* submesh_quad_edge_c = TestAux::create_tetris_quad_edge_submesh_2d();

    // create coarse cell subset
    CellSubSetType* subset_quad_c = TestAux::create_tetris_quad_cellsubset_2d();

    // create coarse cell subset
    CellSubSetType* subset_quad_edge_c = TestAux::create_tetris_quad_edge_cellsubset_2d();

    /* ********************************************************************* */

    // create a root mesh node
    RootMeshNodeType* root_node_c = new RootMeshNodeType(root_mesh_c);

    // add quad submesh node
    SubMeshNodeType* subnode_quad_c =
      root_node_c->add_submesh_node(0, new SubMeshNodeType(submesh_quad_c));

    // add edge submesh node
    //SubMeshNodeType* subnode_edge_c =
      root_node_c->add_submesh_node(1, new SubMeshNodeType(submesh_edge_c));

    // add edge sub-submesh node
    //SubMeshNodeType* subnode_quad_edge_c =
      subnode_quad_c->add_submesh_node(2, new SubMeshNodeType(submesh_quad_edge_c));

    // add quad cell subset node
    CellSubSetNodeType* subsetnode_quad_c =
      root_node_c->add_subset_node(7, new CellSubSetNodeType(subset_quad_c));

    // add edge cell subset node
    //CellSubSetNodeType* subsetnode_quad_edge_c =
      subsetnode_quad_c->add_subset_node(42, new CellSubSetNodeType(subset_quad_edge_c));

    /* ********************************************************************* */

    // refine the mesh tree
    RootMeshNodeType* root_node_f = root_node_c->refine();

    /* ********************************************************************* */

    // fetch quad submesh node
    SubMeshNodeType* subnode_quad_f = root_node_f->find_submesh_node(0);
    TEST_CHECK_MSG(subnode_quad_f != nullptr, "failed to fetch quad submesh node");

    // fetch edge submesh node
    SubMeshNodeType* subnode_edge_f = root_node_f->find_submesh_node(1);
    TEST_CHECK_MSG(subnode_edge_f != nullptr, "failed to fetch edge submesh node");

    // fetch quad-edge submesh
    SubMeshNodeType* subnode_quad_edge_f = subnode_quad_f->find_submesh_node(2);
    TEST_CHECK_MSG(subnode_quad_edge_f != nullptr, "failed to fetch quad-edge submesh node");

    // fetch quad cell subset node
    CellSubSetNodeType* subsetnode_quad_f = root_node_f->find_subset_node(7);
    TEST_CHECK_MSG(subsetnode_quad_f != nullptr, "failed to fetch quad cell subset node");

    // fetch quad-edge cell subset node
    CellSubSetNodeType* subsetnode_quad_edge_f = subsetnode_quad_f->find_subset_node(42);
    TEST_CHECK_MSG(subsetnode_quad_edge_f != nullptr, "failed to fetch quad-edge cell subset node");

    /* ********************************************************************* */

    try
    {
      // validate refined root mesh
      TestAux::validate_refined_tetris_mesh_2d(*root_node_f->get_mesh());

      // validate refined quad submesh
      TestAux::validate_refined_tetris_quad_submesh_2d(*subnode_quad_f->get_mesh());

      // validate refined edge submesh
      TestAux::validate_refined_tetris_edge_submesh_2d(*subnode_edge_f->get_mesh());

      // validate refined quad-edge submesh
      TestAux::validate_refined_tetris_quad_edge_submesh_2d(*subnode_quad_edge_f->get_mesh());

      // validate refined quad cell subset
      TestAux::validate_refined_tetris_quad_cellsubset_2d(*subsetnode_quad_f->get_subset());

      // validate refined quad-edge cell subset
      TestAux::validate_refined_tetris_quad_edge_cellsubset_2d(*subsetnode_quad_edge_f->get_subset());
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
