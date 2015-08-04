#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/util/mesh_streamer.hpp>

using namespace FEAST;
using namespace FEAST::Geometry;

String get_file_title(const String& filename)
{
  // find last slash
  std::size_t p = filename.find_last_of("\\/");
  if(p == filename.npos)
    p = 0;
  else
    ++p;

  // fine last dot
  std::size_t q = filename.find_last_of(".");
  if(q == filename.npos)
    return filename.substr(p);
  else
    return filename.substr(p, q-p);
}

template<typename MeshType_>
int run(MeshStreamer& my_streamer, const String& filename, const Index lvl_min, const Index lvl_max)
{
  // create atlas
  std::cout << "Creating mesh atlas..." << std::endl;
  MeshAtlas<MeshType_>* atlas = nullptr;
  try
  {
    atlas = new MeshAtlas<MeshType_>(my_streamer);
  }
  catch(std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    return 1;
  }

  // create mesh node
  std::cout << "Creating mesh node..." << std::endl;
  RootMeshNode<MeshType_>* node = nullptr;
  try
  {
    node = new RootMeshNode<MeshType_>(my_streamer, atlas);
//    node->generate_meshpart_topologies();
    node->adapt();
  }
  catch(std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    return 1;
  }

  // get all mesh part names
  std::deque<String> part_names = node->get_mesh_part_names();

  // refine
  for(Index lvl(0); lvl <= lvl_max; ++lvl)
  {
    if(lvl > 0)
    {
      std::cout << "Refining up to level " << lvl << "..." << std::endl;
      auto* old = node;
      node = old->refine();
      delete old;
    }

    if(lvl < lvl_min)
      continue;

    // Create a VTK exporter for our mesh
    FEAST::String vtkname = get_file_title(filename) + "." + stringify(lvl);
    std::cout << "Writing file '" << vtkname << ".vtu'..." << std::endl;
    Geometry::ExportVTK<MeshType_> exporter(*node->get_mesh());

    std::vector<double> vtx_data(node->get_mesh()->get_num_entities(0), 0.0);

    // loop over all mesh parts
    for(auto it = part_names.begin(); it != part_names.end(); ++it)
    {
      // get the mesh part
      MeshPart<MeshType_>* part = node->find_mesh_part(*it);
      if(part == nullptr)
        continue;

      // get the vertex target set
      TargetSet& trg = part->template get_target_set<0>();

      // mark all indexes vertices
      for(Index i(0); i < trg.get_num_entities(); ++i)
        vtx_data[trg[i]] = 1.0;

      // add variable
      exporter.add_scalar_vertex(*it, vtx_data.data());

      // unmark vertices
      for(Index i(0); i < trg.get_num_entities(); ++i)
        vtx_data[trg[i]] = 0.0;
    }

    exporter.write(vtkname);
  }

  delete node;
  delete atlas;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    std::cout << "USAGE: mesh2vtk <meshfile> [<level-min>] [<level-max>]" << std::endl;
    return 0;
  }

  // Input file name, required
  FEAST::String filename(argv[1]);

  // Refinement level
  Index lvl_min(0);
  if((argc > 2) && (!String(argv[2]).parse(lvl_min)))
  {
    std::cerr << "ERROR: Failed to parse '" << argv[2] << "' as refinement level!" << std::endl;
    return 1;
  }
  Index lvl_max(lvl_min);
  if((argc > 3) && (!String(argv[3]).parse(lvl_max)))
  {
    std::cerr << "ERROR: Failed to parse '" << argv[3] << "' as refinement level!" << std::endl;
    return 1;
  }

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<1>> Simplex1Mesh;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>> Simplex2Mesh;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>> Simplex3Mesh;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>> Hypercube2Mesh;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>> Hypercube3Mesh;

  // Create a MeshStreamer and read the mesh file
  MeshStreamer my_streamer;
  std::cout << "Parsing mesh file '" << filename << "'..." << std::endl;
  try
  {
    my_streamer.parse_mesh_file(filename);
  }
  catch(const std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    return 1;
  }

  // This is the raw mesh data my_streamer read from filename
  auto& mesh_data = my_streamer.get_root_mesh_node()->mesh_data;
  // Marker int for the MeshType
  int mesh_type = mesh_data.mesh_type;
  // Marker int for the ShapeType
  int shape_type = mesh_data.shape_type;

  if(mesh_type != mesh_data.mt_conformal)
  {
    std::cerr << "ERROR: Unsupported mesh type; only conformal meshes are supported" << std::endl;
    return 1;
  }

  // Call the run() method of the appropriate wrapper class
  if(shape_type == mesh_data.st_edge)
    return run<Simplex1Mesh>(my_streamer, filename, lvl_min, lvl_max);
  if(shape_type == mesh_data.st_tria)
    return run<Simplex2Mesh>(my_streamer, filename, lvl_min, lvl_max);
  if(shape_type == mesh_data.st_tetra)
    return run<Simplex3Mesh>(my_streamer, filename, lvl_min, lvl_max);
  if(shape_type == mesh_data.st_quad)
    return run<Hypercube2Mesh>(my_streamer, filename, lvl_min, lvl_max);
  if(shape_type == mesh_data.st_hexa)
    return run<Hypercube3Mesh>(my_streamer, filename, lvl_min, lvl_max);

  // If no MeshType from the list was in the file, return 1
  return 1;
}
