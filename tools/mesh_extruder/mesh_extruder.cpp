#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/geometry/mesh_extruder.hpp>
#include <kernel/geometry/partition_set.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

namespace MeshExtruderTool
{
  using namespace FEAT;
  using namespace Geometry;

  typedef ConformalMesh<Shape::Quadrilateral> QuadMesh;
  typedef ConformalMesh<Shape::Hexahedron> HexaMesh;

  typedef QuadMesh::VertexSetType QuadVertexSet;
  typedef HexaMesh::VertexSetType HexaVertexSet;

  typedef QuadMesh::IndexSetHolderType QuadIdxSetHolder;
  typedef HexaMesh::IndexSetHolderType HexaIdxSetHolder;

  typedef MeshPart<QuadMesh> QuadPart;
  typedef MeshPart<HexaMesh> HexaPart;

  typedef QuadPart::TargetSetHolderType QuadTrgSetHolder;
  typedef HexaPart::TargetSetHolderType HexaTrgSetHolder;

  typedef QuadPart::MeshAttributeType QuadAttrib;
  typedef HexaPart::MeshAttributeType HexaAttrib;

  typedef MeshAtlas<QuadMesh> QuadAtlas;
  typedef MeshAtlas<HexaMesh> HexaAtlas;

  typedef RootMeshNode<QuadMesh> QuadNode;
  typedef RootMeshNode<HexaMesh> HexaNode;

  int main(int argc, char** argv)
  {
    if(argc < 8)
    {
      std::cout << "USAGE: " << argv[0] << " <in-file> <out-file> <planes> <z-min> <z-max> <z-min-name> <z-max-name>" << std::endl;
      return 0;
    }

    // input/output filenames
    const String ifilename(argv[1]);
    const String ofilename(argv[2]);

    // number of slices
    Index slices(0);
    if(!String(argv[3]).parse(slices))
    {
      std::cerr << "ERROR: Failed to parse '" << argv[3] << "' as slice count" << std::endl;
      return 1;
    }

    // Z-Min/Max
    Real z_min(0.0), z_max(0.0);
    if(!String(argv[4]).parse(z_min))
    {
      std::cerr << "ERROR: Failed to parse '" << argv[4] << "' as z-min" << std::endl;
      return 1;
    }
    if(!String(argv[5]).parse(z_max))
    {
      std::cerr << "ERROR: Failed to parse '" << argv[5] << "' as z-min" << std::endl;
      return 1;
    }

    // z-min/-max mesh-part names
    const String z_min_name(argv[6]);
    const String z_max_name(argv[7]);

    // validation
    if(z_min + 1E-3 >= z_max)
    {
      std::cerr << "ERROR: z_min must be strictly less than z_max" << std::endl;
      return 1;
    }
    if(slices <= Index(0))
    {
      std::cerr << "ERROR: slice count must be > 0" << std::endl;
      return 1;
    }

    // input atlas and mesh node
    QuadAtlas quad_atlas;
    QuadNode quad_node(nullptr, &quad_atlas);
    PartitionSet quad_part_set;

    // output atlas and mesh node
    HexaAtlas hexa_atlas;
    HexaNode hexa_node(nullptr, &hexa_atlas);
    PartitionSet hexa_part_set;

    std::cout << "Reading input mesh file '" << ifilename << "'..." << std::endl;
    {
      // try to open input file
      std::ifstream ifs(ifilename, std::ios_base::in);
      if(!ifs.is_open() || !ifs.good())
      {
        std::cerr << "ERROR: Failed to open '" << ifilename << "' as input file" << std::endl;
        return 1;
      }

      // create mesh reader for input file
      MeshFileReader mesh_reader(ifs);

      // read root markup
      mesh_reader.read_root_markup();

      // check mesh type
      if(mesh_reader.get_meshtype_string() != "conformal:hypercube:2:2")
      {
        std::cerr << "ERROR: Only 2D quadrilateral input meshes are supported" << std::endl;
        return 1;
      }

      // parse
      mesh_reader.parse(quad_node, quad_atlas, &quad_part_set);
    }

    // extrude
    {
      // create mesh extruder
      MeshExtruder<QuadMesh> mesh_extruder(slices, z_min, z_max, z_min_name, z_max_name);

      // extrude atlas
      std::cout << "Extruding mesh atlas..." << std::endl;
      mesh_extruder.extrude_atlas(hexa_atlas, quad_atlas);

      // extrude root mesh
      std::cout << "Extruding mesh node..." << std::endl;
      mesh_extruder.extrude_root_node(hexa_node, quad_node, &hexa_atlas);

      // extrude partition set
      std::cout << "Extruding partition set..." << std::endl;
      mesh_extruder.extrude_partition_set(hexa_part_set, quad_part_set);
    }

    std::cout << "Writing output mesh file '" << ofilename << "'..." << std::endl;
    {
      // try to open output file
      std::ofstream ofs(ofilename, std::ios_base::out|std::ios_base::trunc);
      if(!ofs.is_open() || !ofs.good())
      {
        std::cerr << "ERROR: Failed to open '" << ofilename << "' as output file" << std::endl;
        return 1;
      }

      // create mesh writer for output file
      MeshFileWriter mesh_writer(ofs);

      // write mesh and atlas
      mesh_writer.write(&hexa_node, &hexa_atlas, &hexa_part_set);
    }

    // okay
    std::cout << "Finished!" << std::endl;
    return 0;
  }
} // namespace MeshExtruder

int main(int argc, char** argv)
{
  try
  {
    return MeshExtruderTool::main(argc, argv);
  }
  catch(std::exception& exc)
  {
    std::cerr << exc.what() << std::endl;
    return 1;
  }
}
