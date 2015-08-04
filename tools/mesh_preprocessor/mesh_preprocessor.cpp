#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/util/mesh_streamer.hpp>
#include <kernel/util/simple_arg_parser.hpp>

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
int run(MeshStreamer& my_streamer, const String& filename, const Index lvl_min, const Index lvl_max,
const std::deque<String>& deduct_from_bottom,
const std::deque<String>& deduct_from_top,
const std::deque<String>& deduct_topology)
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

    // Deduct TargetSets from bottom to top for specified MeshParts
    for(auto& it : deduct_from_bottom)
    {
      typename RootMeshNode<MeshType_>::MeshPartType* mesh_part(node->find_mesh_part(it));
      if(mesh_part != nullptr)
      {
        mesh_part->template deduct_target_sets_from_bottom<0>(node->get_mesh()->get_index_set_holder());
        std::cout << "Deducting target sets from bottom to top for MeshPart " << mesh_part->get_identifier() << std::endl;
      }
      else
      {
        std::cerr << "ERROR: MeshPart " << it << " was not found but was specified via --deduct_from_bottom!" << std::endl;
        return 1;
      }
    }

    // Deduct TargetSets from top to bottom for specified MeshParts
    for(auto& it : deduct_from_top)
    {
      typename RootMeshNode<MeshType_>::MeshPartType* mesh_part(node->find_mesh_part(it));
      if(mesh_part != nullptr)
      {
        mesh_part->template deduct_target_sets_from_top<MeshType_::shape_dim>
          (node->get_mesh()->get_index_set_holder());
        std::cout << "Deducting target sets from top to bottom for MeshPart " << mesh_part->get_identifier() << std::endl;
      }
      else
      {
        std::cerr << "ERROR: MeshPart " << it << " was not found but was specified via --deduct_from_top!" << std::endl;
        return 1;
      }
    }

    // Deduct topology for specified MeshParts
    for(auto& it : deduct_topology)
    {
      typename RootMeshNode<MeshType_>::MeshPartType* mesh_part(node->find_mesh_part(it));
      if(mesh_part != nullptr)
      {
        mesh_part->deduct_topology(node->get_mesh()->get_index_set_holder());
        std::cout << "Deducting topology for MeshPart " << mesh_part->get_identifier() << std::endl;
      }
      else
      {
        std::cerr << "ERROR: MeshPart " << it << " was not found but was specified via --deduct_topology!" << std::endl;
        return 1;
      }
    }

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
  // create an argument parser
  SimpleArgParser args(argc, argv);

  // add all supported options
  args.support("deduct_from_bottom");
  args.support("deduct_from_top");
  args.support("deduct_topology");
  args.support("filenames");
  args.support("levels");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
    return 1;
  }

  auto* filenames_pair(args.query("filenames"));
  if(filenames_pair == nullptr)
  {
    std::cerr << "ERROR: No filenames specified!" << std::endl;
    return 1;
  }
  std::deque<String> filenames(filenames_pair->second);

  std::deque<String> deduct_from_bottom;
  auto* deduct_from_bottom_pair(args.query("deduct_from_bottom"));
  if(deduct_from_bottom_pair != nullptr)
    deduct_from_bottom = deduct_from_bottom_pair->second;

  std::deque<String> deduct_from_top;
  auto* deduct_from_top_pair(args.query("deduct_from_top"));
  if(deduct_from_top_pair != nullptr)
    deduct_from_top = deduct_from_top_pair->second;

  std::deque<String> deduct_topology;
  auto* deduct_topology_pair(args.query("deduct_topology"));
  if(deduct_from_bottom_pair != nullptr)
    deduct_topology = deduct_topology_pair->second;

  Index lvl_min(0);
  Index lvl_max(0);
  int level_flag(args.parse("levels", lvl_min, lvl_max));
  if(level_flag < 0)
  {
    std::cerr << "ERROR: Invalid level information " << args.get_arg(-level_flag) << std::endl;
    return 1;
  }
  if(level_flag == 1)
  {
    lvl_max = lvl_min;
    lvl_min = 0;
  }


  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<1>> Simplex1Mesh;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>> Simplex2Mesh;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>> Simplex3Mesh;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>> Hypercube2Mesh;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>> Hypercube3Mesh;

  // Create a MeshStreamer and read the mesh file
  MeshStreamer my_streamer;
//  std::cout << "Parsing mesh files '" << filenames[0] << "'..." << std::endl;
  try
  {
    my_streamer.parse_multiple_files(filenames);
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

  String& filename(filenames[0]);
  // Call the run() method of the appropriate wrapper class
  if(shape_type == mesh_data.st_edge)
    return run<Simplex1Mesh>(my_streamer, filename, lvl_min, lvl_max, deduct_from_bottom, deduct_from_top, deduct_topology);
  if(shape_type == mesh_data.st_tria)
    return run<Simplex2Mesh>(my_streamer, filename, lvl_min, lvl_max, deduct_from_bottom, deduct_from_top, deduct_topology);
  if(shape_type == mesh_data.st_tetra)
    return run<Simplex3Mesh>(my_streamer, filename, lvl_min, lvl_max, deduct_from_bottom, deduct_from_top, deduct_topology);
  if(shape_type == mesh_data.st_quad)
    return run<Hypercube2Mesh>(my_streamer, filename, lvl_min, lvl_max, deduct_from_bottom, deduct_from_top, deduct_topology);
  if(shape_type == mesh_data.st_hexa)
    return run<Hypercube3Mesh>(my_streamer, filename, lvl_min, lvl_max, deduct_from_bottom, deduct_from_top, deduct_topology);

  // If no MeshType from the list was in the file, return 1
  return 1;
}
