#include<kernel/geometry/conformal_factories.hpp>
#include<kernel/geometry/export_vtk.hpp>
#include<kernel/geometry/mesh_streamer_factory.hpp>
#include<kernel/geometry/mesh_node.hpp>
#include<kernel/geometry/mesh_part.hpp>
#include<kernel/util/mesh_streamer.hpp>
#include<kernel/util/simple_arg_parser.hpp>

using namespace FEAST;
using namespace FEAST::Geometry;

/**
 * \brief Wrapper struct for writing a streamed mesh to vtk
 *
 * \tparam MeshType_
 * Type of the mesh, i.e. ConformalMesh<Simplex<2>>
 *
 * \author Jordi Paul
 *
 */
template<typename MeshType_>
struct MeshStreamerApp
{
  /**
   * \brief Routine that does the actual work
   *
   * \param[in] my_streamer
   * MeshStreamer that contains the data from the mesh file.
   *
   * \param[in] level
   * Number of refines.
   */
  static int run(MeshStreamer& my_streamer, const Index level)
  {
    RefineFactory<MeshType_, MeshStreamerFactory> my_factory(level, my_streamer);
    MeshType_ my_mesh(my_factory);

    //typedef typename MeshType_::ShapeType ShapeType;
    //typedef StandardConformalMeshNodePolicy<ShapeType> MeshNodePolicy;

    //typedef typename MeshNodePolicy::RootMeshType MeshType;
    //typedef RootMeshNode<MeshNodePolicy> RootMeshNodeType;

    //RootMeshNodeType my_root_mesh_node(my_streamer);

    //MeshStreamerFactory<MeshPart<MeshType_>> my_mesh_part_factory(my_streamer,"1");
    //MeshPart<MeshType_>  my_mesh_part(my_mesh_part_factory);
    //MeshStreamer::MeshNode* root_mesh_node = my_streamer.get_root_mesh_node();
    //my_streamer.write_mesh_file("streamed_mesh.txt");
    //std::cout << "MeshPart num_entities = [";
    //for(int d(0); d <= my_mesh_part.shape_dim; ++d)
    //  std::cout << " " << my_mesh_part.get_num_entities(d);
    //std::cout << "]" << std::endl;

    String output_name = "mesh_streamer_" + stringify(level) + ".vtk";
    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType_> exporter(my_mesh);
    std::cout << "Writing file " << output_name << std::endl;
    // finally, write the VTK file
    exporter.write(output_name);
    return 0;
  }
};


/**
 * \cond internal
 *
 * Mesh Streamer Application
 *
 */
int main(int argc, char* argv[])
{
  // Creata a parser for command line arguments.
  SimpleArgParser args(argc, argv);

  if( args.check("help") > -1 || args.num_args()==1)
  {
    std::cout << "Mesh Streamer Application usage: " << std::endl;
    std::cout << "Required arguments: --filename [String]: Path to a FEAST mesh file." << std::endl;
    std::cout << "Optional arguments: --level [unsigned int]: Number of refines, defaults to 0." << std::endl;
    exit(1);
  }
  // Specify supported command line switches
  args.support("level");
  args.support("filename");
  args.support("help");
  // Refinement level
  Index level(0);
  // Input file name, required
  FEAST::String filename;
  // Get unsupported command line arguments
  std::deque<std::pair<int,String> > unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
  }

  // Check and parse --filename
  if(args.check("filename") != 1 )
    throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --filename");
  else
  {
    args.parse("filename", filename);
    std::cout << "Reading mesh from file " << filename << std::endl;
  }

  // Check and parse --level
  if(args.check("level") != 1)
    std::cout << "No refinement level specified, defaulting to 0." << std::endl;
  else
  {
    args.parse("level", level);
    std::cout << "Refinement level " << level << std::endl;
  }

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<1>> Simplex1Mesh;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>> Simplex2Mesh;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>> Simplex3Mesh;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>> Hypercube2Mesh;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>> Hypercube3Mesh;

  // Create a MeshStreamer and read the mesh file
  MeshStreamer my_streamer;
  my_streamer.parse_mesh_file(filename);

  // This is the raw mesh data my_streamer read from filename
  auto& mesh_data = my_streamer.get_root_mesh_node()->mesh_data;
  // Marker int for the MeshType
  int mesh_type = mesh_data.mesh_type;
  // Marker int for the ShapeType
  int shape_type = mesh_data.shape_type;

  ASSERT(mesh_type == mesh_data.mt_conformal, "This application only works for conformal meshes!");

  // Call the run() method of the appropriate wrapper class
  if(shape_type == mesh_data.st_edge)
    return MeshStreamerApp<Simplex1Mesh>::run(my_streamer, level);
  if(shape_type == mesh_data.st_tria)
    return MeshStreamerApp<Simplex2Mesh>::run(my_streamer, level);
  if(shape_type == mesh_data.st_tetra)
    return MeshStreamerApp<Simplex3Mesh>::run(my_streamer, level);
  if(shape_type == mesh_data.st_quad)
    return MeshStreamerApp<Hypercube2Mesh>::run(my_streamer, level);
  if(shape_type == mesh_data.st_hexa)
    return MeshStreamerApp<Hypercube3Mesh>::run(my_streamer, level);

  // If no MeshType from the list was in the file, return 1
  return 1;
}
/// \endcond
