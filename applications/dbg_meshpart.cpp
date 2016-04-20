#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/util/mesh_streamer.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAST;
using namespace FEAST::Geometry;

typedef ConformalMesh<Shape::Simplex<2>> MeshType;
typedef MeshAtlas<MeshType> AtlasType;
typedef RootMeshNode<MeshType> MeshNodeType;

int main(int argc, char** argv)
{
  Runtime::initialise(argc, argv);

  std::cout << "Creating streamer..." << std::endl;
  // create streamer
  MeshStreamer streamer;
  streamer.parse_mesh_file("unit-circle-tria.txt");

  std::cout << "Creating atlas..." << std::endl;
  // create atlas
  AtlasType atlas(streamer);

  std::cout << "Creating base-mesh root node..." << std::endl;
  // create base root node
  MeshNodeType base_node(streamer, &atlas);

  std::cout << "Refining base-mesh root node..." << std::endl;
  MeshNodeType* ref_base_node_tmp = base_node.refine();
  MeshNodeType* ref_base_node = ref_base_node_tmp;

  // create rank graph
  Adjacency::Graph ranks_at_elem(16, 16, 16);
  for(Index i(0); i < 16; ++i)
  {
    ranks_at_elem.get_domain_ptr()[i] = i;
    ranks_at_elem.get_image_idx()[i] = (i / 4);
  }
  ranks_at_elem.get_domain_ptr()[16] = 16;
  std::vector<Index> ranks0, ranks1, ranks2, ranks3;
  ranks0.push_back(1);
  ranks0.push_back(2);
  ranks0.push_back(3);
  ranks1.push_back(0);
  ranks1.push_back(2);
  ranks1.push_back(3);
  ranks2.push_back(0);
  ranks2.push_back(1);
  ranks2.push_back(3);
  ranks3.push_back(0);
  ranks3.push_back(1);
  ranks3.push_back(2);

  // create all partitions
  std::cout << "Creating patches..." << std::endl;
  MeshNodeType* patch_0 = ref_base_node->extract_patch(0, ranks_at_elem, ranks0);
  MeshNodeType* patch_1 = ref_base_node->extract_patch(1, ranks_at_elem, ranks1);
  MeshNodeType* patch_2 = ref_base_node->extract_patch(2, ranks_at_elem, ranks2);
  MeshNodeType* patch_3 = ref_base_node->extract_patch(3, ranks_at_elem, ranks3);

  std::cout << "Creating mesh streamers..." << std::endl;
  auto *streamer_b = ref_base_node->create_mesh_writer();
  auto* streamer_0 = patch_0->create_mesh_writer();
  auto* streamer_1 = patch_1->create_mesh_writer();
  auto* streamer_2 = patch_2->create_mesh_writer();
  auto* streamer_3 = patch_3->create_mesh_writer();

  std::cout << "Writing patches..." << std::endl;
  streamer_b->write_mesh_file("part_base.txt");
  streamer_0->write_mesh_file("patch_0.txt");
  streamer_1->write_mesh_file("patch_1.txt");
  streamer_2->write_mesh_file("patch_2.txt");
  streamer_3->write_mesh_file("patch_3.txt");

  return Runtime::finalise();
}
