// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/geometry/mesh_extruder.hpp>
#include <kernel/geometry/partition_set.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

namespace MeshExtruderTool
{
  using namespace FEAT;
  using namespace Geometry;

  typedef double CoordType;

  typedef ConformalMesh<Shape::Quadrilateral, 2, CoordType> QuadMesh;
  typedef ConformalMesh<Shape::Hexahedron, 3, CoordType> HexaMesh;

  typedef QuadMesh::VertexSetType QuadVertexSet;
  typedef HexaMesh::VertexSetType HexaVertexSet;

  typedef QuadMesh::IndexSetHolderType QuadIdxSetHolder;
  typedef HexaMesh::IndexSetHolderType HexaIdxSetHolder;

  typedef MeshPart<QuadMesh> QuadPart;
  typedef MeshPart<HexaMesh> HexaPart;

  typedef QuadPart::TargetSetHolderType QuadTrgSetHolder;
  typedef HexaPart::TargetSetHolderType HexaTrgSetHolder;

  typedef QuadPart::AttributeSetType QuadAttrib;
  typedef HexaPart::AttributeSetType HexaAttrib;

  typedef MeshAtlas<QuadMesh> QuadAtlas;
  typedef MeshAtlas<HexaMesh> HexaAtlas;

  typedef RootMeshNode<QuadMesh> QuadNode;
  typedef RootMeshNode<HexaMesh> HexaNode;

  void print_help(const SimpleArgParser& args)
  {
    std::cout << "USAGE: " << args.get_arg(0) << " <Options>\n";
    std::cout << "\n";
    std::cout << "Valid Options:\n";
    std::cout << "--------------\n";
    std::cout << args.get_supported_help() << "\n";
  }

  void main(int argc, char** argv)
  {
    SimpleArgParser args(argc, argv);

    args.support("in", "<input-files>\nMandatory: Specifies the 2D input mesh files.\n");
    args.support("out", "<output-file>\nMandatory: Specifies the 3D output mesh file.\n");
    args.support("z-list", "<z_0> <z_1> ... <z_n>>\nMandatory: Specifies a list of Z-coordinates.\nMutually exclusive with --z-range.\n");
    args.support("z-range", "<z-min> <z-max> <n>\nMandatory: Specifies the Z-coordinate range and the number of slices in Z-direction.\nMutually exclusive with --z-list\n");
    args.support("z-names", "<name-min> <name-max>\nMandatory: Specifies the names for the Z-boundary meshparts.\n");
    args.support("origin", "<x> <y>\nOptional: Specifies the 2D transformation origin. Default to 0, 0, 0.\n");
    args.support("offset", "<x> <y> <z>\nOptional: Specifies the 3D transformation offset. Default to 0, 0, 0.\n");
    args.support("angles", "<yaw> <pitch> <roll>\nOptional: Specifies the 3D transformation angles in revolutions. Default to 0, 0, 0.\n");

    if(argc < 2)
    {
      print_help(args);
      return;
    }

    // query unsupported
    auto unsupp = args.query_unsupported();
    if(!unsupp.empty())
    {
      for(const auto& v : unsupp)
        std::cerr << "ERROR: unknown option '--" << v.second << "'\n";
      std::cerr << "\n";
      print_help(args);
      FEAT::Runtime::abort();
    }

    // query input filenames
    const auto* qq = args.query("in");
    if(qq == nullptr)
    {
      std::cerr << "ERROR: mandatory option '--in <...>' is missing\n";
      print_help(args);
      FEAT::Runtime::abort();
    }
    std::deque<String> ifilenames = qq->second;
    if(ifilenames.empty())
    {
      std::cerr << "ERROR: mandatory option '--in <...>' did not specify any input files\n";
      print_help(args);
      FEAT::Runtime::abort();
    }

    // query output filename
    String ofilename;
    if(args.parse("out", ofilename) != 1)
    {
      std::cerr << "ERROR: mandatory option '--out <...>' is invalid or missing\n";
      print_help(args);
      FEAT::Runtime::abort();
    }

    // we need either --z-range or --z-list
    if((args.check("z-range") < 0) && (args.check("z-list") < 0))
    {
      std::cerr << "ERROR: neither --z-range <...> nor --z-list <...> given\n";
      std::cerr << "You have to specify one of those two options.\n";
      print_help(args);
      FEAT::Runtime::abort();
    }
    if((args.check("z-range") >= 0) && (args.check("z-list") >= 0))
    {
      std::cerr << "ERROR: both --z-range <...> and --z-list <...> given, but they are mutually exclusive\n";
      print_help(args);
      FEAT::Runtime::abort();
    }

    // query Z-ranges
    std::deque<CoordType> z_list;
    if(args.check("z-list") >= 0)
    {
      std::deque<String> sz = args.query("z-list")->second;
      if(sz.size() < std::size_t(2))
      {
        std::cerr << "ERROR: at least 2 coordinates are required for --z-list <...>\n";
        FEAT::Runtime::abort();
      }
      for(const auto& s : sz)
      {
        CoordType x(0);
        if(!s.parse(x))
        {
          std::cerr << "ERROR: failed to parse --z-list coordinate: " << s << "\n";
          FEAT::Runtime::abort();
        }
        z_list.push_back(x);
      }
      // ensure that the z-coordinates are ascending
      for(std::size_t i(0); i+1u < z_list.size(); ++i)
      {
        if(z_list[i] + 1E-5 >= z_list[i+1u])
        {
          std::cerr << "ERROR: --z-list coordinates must be strictly ascending\n";
          FEAT::Runtime::abort();
        }
      }
    }
    else //if(args.check("z-range") >= 0)
    {
      Index slices(0);
      CoordType z_min(0.0), z_max(0.0);

      if(args.parse("z-range", z_min, z_max, slices) != 3)
      {
        std::cerr << "ERROR: mandatory option '--z <...>' is invalid or missing\n";
        print_help(args);
        FEAT::Runtime::abort();
      }
      if(z_min + 1E-3 >= z_max)
      {
        std::cerr << "ERROR: z_min must be strictly less than z_max\n";
        print_help(args);
        FEAT::Runtime::abort();
      }
      if(slices <= Index(0))
      {
        std::cerr << "ERROR: slice count must be > 0\n";
        print_help(args);
        FEAT::Runtime::abort();
      }
      for(Index i(0); i <= slices; ++i)
        z_list.push_back(z_min + (z_max-z_min)*CoordType(i)/CoordType(slices));
    }

    // query meshpart names
    String z_min_name, z_max_name;
    if(args.parse("z-names", z_min_name, z_max_name) != 2)
    {
      std::cerr << "ERROR: mandatory option '--z-names <...>' is invalid or missing\n";
      print_help(args);
      FEAT::Runtime::abort();
    }

    // query origin
    Real origin_x(0.0), origin_y(0.0);
    if(args.parse("origin", origin_x, origin_y) < 0)
    {
      std::cerr << "ERROR: Failed to parse '--origin <x> <y>' option\n";
      print_help(args);
      FEAT::Runtime::abort();
    }

    // query offset
    Real offset_x(0.0), offset_y(0.0), offset_z(0.0);
    if(args.parse("offset", offset_x, offset_y, offset_z) < 0)
    {
      std::cerr << "ERROR: Failed to parse '--offset <x> <y> <z>' option\n";
      print_help(args);
      FEAT::Runtime::abort();
    }

    // query angles
    Real angle_y(0.0), angle_p(0.0), angle_r(0.0);
    if(args.parse("angles", angle_y, angle_p, angle_r) < 0)
    {
      std::cerr << "ERROR: Failed to parse '--angles <yaw> <pitch> <roll>' option\n";
      print_help(args);
      FEAT::Runtime::abort();
    }
    else
    {
      // convert revolutions to radians
      const Real pi2 = Real(2) * Math::pi<Real>();
      angle_y *= pi2;
      angle_p *= pi2;
      angle_r *= pi2;
    }

    // input atlas and mesh node
    QuadAtlas quad_atlas;
    QuadNode quad_node(nullptr, &quad_atlas);
    PartitionSet quad_part_set;

    // output atlas and mesh node
    HexaAtlas hexa_atlas;
    HexaNode hexa_node(nullptr, &hexa_atlas);
    PartitionSet hexa_part_set;

    {
      // create mesh reader for input file
      MeshFileReader mesh_reader;
      std::deque<std::shared_ptr<std::ifstream>> streams;

      // add all input files
      for(std::size_t i(0); i < ifilenames.size(); ++i)
      {
        String filename = ifilenames.at(i);

        // try to open file
        streams.emplace_back(std::make_shared<std::ifstream>(filename, std::ios_base::in));
        std::ifstream& ifs = *streams.back();
        if(!ifs.is_open() || !ifs.good())
        {
          std::cerr << "ERROR: Failed to open mesh file '" << (filename) << "'\n";
          FEAT::Runtime::abort();
        }

        // add to mesh reader
        mesh_reader.add_stream(ifs);
      }

      // read root markup
      mesh_reader.read_root_markup();

      // check mesh type
      if(mesh_reader.get_meshtype_string() != "conformal:hypercube:2:2")
      {
        std::cerr << "ERROR: Only 2D quadrilateral input meshes are supported\n";
        FEAT::Runtime::abort();
      }

      // parse
      mesh_reader.parse(quad_node, quad_atlas, &quad_part_set);
    }

    // extrude
    {
      // create mesh extruder
      //MeshExtruder<QuadMesh> mesh_extruder(slices, z_min, z_max, z_min_name, z_max_name);
      MeshExtruder<QuadMesh> mesh_extruder(z_list, z_min_name, z_max_name);

      // set offset, origin and angles
      mesh_extruder.set_origin(origin_x, origin_y);
      mesh_extruder.set_offset(offset_x, offset_y, offset_z);
      mesh_extruder.set_angles(angle_y, angle_p, angle_r);

      // extrude atlas
      std::cout << "Extruding mesh atlas...\n";
      mesh_extruder.extrude_atlas(hexa_atlas, quad_atlas);

      // extrude root mesh
      std::cout << "Extruding mesh node...\n";
      mesh_extruder.extrude_root_node(hexa_node, quad_node, &hexa_atlas);

      // extrude partition set
      std::cout << "Extruding partition set...\n";
      mesh_extruder.extrude_partition_set(hexa_part_set, quad_part_set);
    }

    std::cout << "Writing output mesh file '" << ofilename << "'...\n";
    {
      // try to open output file
      std::ofstream ofs(ofilename, std::ios_base::out|std::ios_base::trunc);
      if(!ofs.is_open() || !ofs.good())
      {
        std::cerr << "ERROR: Failed to open '" << ofilename << "' as output file\n";
        FEAT::Runtime::abort();
      }

      // create mesh writer for output file
      MeshFileWriter mesh_writer(ofs);

      // write mesh and atlas
      mesh_writer.write(&hexa_node, &hexa_atlas, &hexa_part_set);
    }

    // okay
    std::cout << "Finished!\n";
  }
} // namespace MeshExtruder

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  try
  {
    MeshExtruderTool::main(argc, argv);
  }
  catch(std::exception& exc)
  {
    std::cerr << exc.what() << "\n";
    FEAT::Runtime::abort();
  }
  return 0;
}
