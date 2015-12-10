#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/util/mesh_streamer.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <kernel/util/runtime.hpp>

/**
 * Mesh preprocessing tool, mainly copy-pasted from tools/mesh2vtk.
 *
 * Some meshes (or rather, their MeshParts) are given just by a set of points without any higher-dimensional topology.
 * This topology can be inferred by first deducting higher-dimensional target mappings.
 *
 * Assume we have a 3d mesh with a MeshPart representing a 2d boundary of that mesh and that the MeshPart consists of
 * only vertices. We can now deduct a topology of edges and faces on the MeshPart by checking which edges and faces
 * of the original mesh are made up entirely of vertices belonging to the MeshPart.
 * This is "from bottom target deduction".
 *
 * "from top target deduction" is even simpler.
 *
 * Having higher-dimensional entities in a MeshPart changes how it behaves under refinement. Namely, if it has edges,
 * newly inserted vertices on these edge belong to the refined MeshPart.
 *
 * All of this deduction is carried out by routines in the kernel, but deciding at runtime if and when to call them
 * is not very safe. So the assumption is that a streamed mesh has complete information OR that the user really knows
 * what she is doing and uses the kernel routines correctly.
 *
 * For the cases where a mesh does NOT have complete information (i.e. the mesh file was written by some 3rd party
 * application, or converted from another format), this preprocessing tool can do the deriving and writes a new
 * mesh file containing everything needed.
 *
 * \author Jordi Paul
 */

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

      auto* mesh_writer = node->create_mesh_writer();
      FEAST::String preprocessed_name = "preprocessed_"+get_file_title(filename) + "." + stringify(lvl)+".txt";
      std::cout << "Writing mesh " << preprocessed_name << "..." << std::endl;

      mesh_writer->write_mesh_file(preprocessed_name);
      delete mesh_writer;
    }

    delete node;
    if(atlas != nullptr)
      delete atlas;

    return 0;
  }

int main(int argc, char** argv)
{
  FEAST::Runtime::initialise(argc, argv);
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
    FEAST::Runtime::abort();
  }

  auto* filenames_pair(args.query("filenames"));
  if(filenames_pair == nullptr)
  {
    std::cerr << "ERROR: No filenames specified!" << std::endl;
    FEAST::Runtime::abort();
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
  if(deduct_topology_pair != nullptr)
    deduct_topology = deduct_topology_pair->second;

  Index lvl_min(0);
  Index lvl_max(0);
  int level_flag(args.parse("levels", lvl_min, lvl_max));
  if(level_flag < 0)
    throw InternalError("ERROR: Invalid level information "+args.get_arg(-level_flag));

  if(level_flag == 1)
  {
    lvl_max = lvl_min;
    lvl_min = 0;
  }

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<1>, 1, 1, Real> Simplex1Mesh_1d;
  typedef Geometry::ConformalMesh<Shape::Simplex<1>, 2, 2, Real> Simplex1Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Simplex<1>, 3, 3, Real> Simplex1Mesh_3d;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> Simplex2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, 3, Real> Simplex2Mesh_3d;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, Real> Simplex3Mesh;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> Hypercube2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, 3, Real> Hypercube2Mesh_3d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, Real> Hypercube3Mesh;

  // Create a MeshStreamer and read the mesh file
  MeshStreamer my_streamer;

  try
  {
    my_streamer.parse_multiple_files(filenames);
  }
  catch(const std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    FEAST::Runtime::abort();
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
    FEAST::Runtime::abort();
  }

  String& filename(filenames[0]);
  // Call the run() method of the appropriate wrapper class
  int ret(1);
  if(shape_type == mesh_data.st_edge)
  {
    if(mesh_data.coord_per_vertex == 1)
      ret = run<Simplex1Mesh_1d>(my_streamer, filename, lvl_min, lvl_max,
      deduct_from_bottom, deduct_from_top, deduct_topology);
    if(mesh_data.coord_per_vertex == 2)
      ret = run<Simplex1Mesh_2d>(my_streamer, filename, lvl_min, lvl_max,
      deduct_from_bottom, deduct_from_top, deduct_topology);
    if(mesh_data.coord_per_vertex == 3)
      ret = run<Simplex1Mesh_3d>(my_streamer, filename, lvl_min, lvl_max,
      deduct_from_bottom, deduct_from_top, deduct_topology);

    std::cerr << "ERROR: Unsupported world dim " << mesh_data.coord_per_vertex << " for Simplex<1> mesh!"<< std::endl;
    FEAST::Runtime::abort();
  }

  if(shape_type == mesh_data.st_tria)
  {
    if(mesh_data.coord_per_vertex == 2)
      return run<Simplex2Mesh_2d>(my_streamer, filename, lvl_min, lvl_max,
      deduct_from_bottom, deduct_from_top, deduct_topology);
    if(mesh_data.coord_per_vertex == 3)
      return run<Simplex2Mesh_3d>(my_streamer, filename, lvl_min, lvl_max,
      deduct_from_bottom, deduct_from_top, deduct_topology);
    std::cerr << "ERROR: Unsupported world dim " << mesh_data.coord_per_vertex << " for Simplex<2> mesh!"<< std::endl;
    FEAST::Runtime::abort();
  }

  if(shape_type == mesh_data.st_tetra)
    ret = run<Simplex3Mesh>(my_streamer, filename, lvl_min, lvl_max,
    deduct_from_bottom, deduct_from_top, deduct_topology);

  if(shape_type == mesh_data.st_quad)
  {
    if(mesh_data.coord_per_vertex == 2)
      ret = run<Hypercube2Mesh_2d>(my_streamer, filename, lvl_min, lvl_max,
      deduct_from_bottom, deduct_from_top, deduct_topology);
    if(mesh_data.coord_per_vertex == 3)
      ret = run<Hypercube2Mesh_3d>(my_streamer, filename, lvl_min, lvl_max,
      deduct_from_bottom, deduct_from_top, deduct_topology);

    std::cerr << "ERROR: Unsupported world dim " << mesh_data.coord_per_vertex << " for Hypercube<2> mesh!"<< std::endl;
    FEAST::Runtime::abort();
  }

  if(shape_type == mesh_data.st_hexa)
    ret = run<Hypercube3Mesh>(my_streamer, filename, lvl_min, lvl_max, deduct_from_bottom, deduct_from_top, deduct_topology);

  ret = ret | FEAST::Runtime::finalise();
  // If no MeshType from the list was in the file, return 1
  return ret;
}
