// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/tetris_quad.hpp>
#include <kernel/geometry/test_aux/validate_neighbors.hpp>
#include <kernel/geometry/mesh_node.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Geometry;
using namespace FEAT::Geometry::TestAux;

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
  : public TestSystem::UnitTest
{
public:
  MeshNodeTestConfQuad() :
    TestSystem::UnitTest("mesh_node-test-conf-quad")
  {
  }

  virtual ~MeshNodeTestConfQuad()
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
    std::unique_ptr<RootMeshNodeType> root_node_c = RootMeshNodeType::make_unique(std::unique_ptr<RootMeshType>(root_mesh_c));

    // add quad submesh node
    MeshPartNodeType* subnode_quad_c =
      root_node_c->add_mesh_part_node("0", MeshPartNodeType::make_unique(std::unique_ptr<MeshPartType>(submesh_quad_c)));

    // add edge submesh node
    //MeshPartNodeType* subnode_edge_c =
      root_node_c->add_mesh_part_node("1",  MeshPartNodeType::make_unique(std::unique_ptr<MeshPartType>(submesh_edge_c)));

    // add edge sub-submesh node
    //MeshPartNodeType* subnode_quad_edge_c =
      subnode_quad_c->add_mesh_part_node("2",  MeshPartNodeType::make_unique(std::unique_ptr<MeshPartType>(submesh_quad_edge_c)));

    // add quad cell subset node
    MeshPartNodeType* subsetnode_quad_c =
      root_node_c->add_mesh_part_node("7",  MeshPartNodeType::make_unique(std::unique_ptr<MeshPartType>(subset_quad_c)));

    // add edge cell subset node
    //MeshPartNodeType* subsetnode_quad_edge_c =
      subsetnode_quad_c->add_mesh_part_node("42",  MeshPartNodeType::make_unique(std::unique_ptr<MeshPartType>(subset_quad_edge_c)));

    /* ********************************************************************* */

    // refine the mesh tree
    std::unique_ptr<RootMeshNodeType> root_node_f = root_node_c->refine_unique();

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

      // validate neighbor information
      validate_neighbors(*(root_node_f->get_mesh()));

    }
    catch(const String& msg)
    {
      TEST_CHECK_MSG(false, msg);
    }
  }

} mesh_node_test_conf_quad;
