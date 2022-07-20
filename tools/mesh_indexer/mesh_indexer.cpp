// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/parsed_hit_test_factory.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

namespace MeshIndexerTool
{
  using namespace FEAT;

  void display_help()
  {
    std::cout << std::endl;
    std::cout << "mesh-indexer: Creates a FEAT mesh file from vertex coordinate and element index files" << std::endl;
    std::cout << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "--out <meshfile>" << std::endl;
    std::cout << "Specifies the filename of the output mesh file that is to be generated." << std::endl;
    std::cout << std::endl;
    std::cout << "--vtx <filename>" << std::endl;
    std::cout << "Specifies the name of the input text file that contains the vertex coordinates." << std::endl;
    std::cout << "Each line of this file should contain the coordinate tuple of a single vertex" << std::endl;
    std::cout << "and the number of coordinates must match the mesh shape specified by '--shape'." << std::endl;
    std::cout << "Empty lines and lines beginning with the hash character '#' are ignored. " << std::endl;
    std::cout << std::endl;
    std::cout << "--idx <filename>" << std::endl;
    std::cout << "Specifies the name of the input text file that contains the vertices-at-element indices." << std::endl;
    std::cout << "Each line of this file should contain the vertex-index tuple of a single element" << std::endl;
    std::cout << "and the number of indices must match the mesh shape specified by '--shape'." << std::endl;
    std::cout << "Empty lines and lines beginning with the hash character '#' are ignored. " << std::endl;
    std::cout << "Note: The index of the first vertex is 0, not 1." << std::endl;
    std::cout << std::endl;
    std::cout << "--shape <shape>" << std::endl;
    std::cout << "Specifies the shape of the mesh to be generated. Must be one of the following:" << std::endl;
    std::cout << " h2  : Hypercube<2> mesh with 2D coordinates" << std::endl;
    std::cout << " h3  : Hypercube<3> mesh with 3D coordinates" << std::endl;
    std::cout << " h23 : Hypercube<2> mesh with 3D coordinates" << std::endl;
    std::cout << " s2  : Simplex<2> mesh with 2D coordinates" << std::endl;
    std::cout << " s3  : Simplex<3> mesh with 3D coordinates" << std::endl;
    std::cout << " s23 : Simplex<2> mesh with 3D coordinates" << std::endl;
    std::cout << std::endl;
    std::cout << "Optional arguments:" << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "--bnd" << std::endl;
    std::cout << "Specifies that a single mesh part for the entire boundary is to be generated." << std::endl;
    std::cout << "If specified, the generated meshpart will be named 'bnd'." << std::endl;
    std::cout << std::endl;
    std::cout << "--parts <name1> <formula1> [<name2> <formula2> ...]" << std::endl;
    std::cout << "Specifies a set of name-formula argument pairs which are used to generate meshparts" << std::endl;
    std::cout << "by using the Geometry::ParsedHitTestFactory class. The first component of each pair" << std::endl;
    std::cout << "specifies the name for the mesh part, whereas the second component specifies the formula" << std::endl;
    std::cout << "in x,y,z coordinates, which is to be used for the hit test of the mesh part." << std::endl;
    std::cout << "A vertex or edge/face/cell will be contained in the meshpart if the formula evaluates" << std::endl;
    std::cout << "to a positive value in its coordinates or midpoint coordinates, respectively." << std::endl;
    std::cout << "Please note that this option can only be used if FEAT is configured and linked against" << std::endl;
    std::cout << "the 'fparser' third-party library." << std::endl;
    std::cout << std::endl;
    std::cout << "--help" << std::endl;
    std::cout << "Displays this message" << std::endl;
  }

  int parse_vtx(std::vector<Real>& vtx, const String& filename, const int world_dim)
  {
    vtx.reserve(10000u * std::size_t(world_dim));

    // try to open the file
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.is_open() || !ifs.good())
    {
      std::cerr << "ERROR: failed to open vertex file '" << filename << "'" << std::endl;
      return 1;
    }

    std::cout << "Parsing vertex file '" << filename << "'..." << std::endl;

    String line;
    Index line_no(0), count(0);
    std::deque<String> v;
    while(!ifs.eof())
    {
      ++line_no;
      std::getline(ifs, line);
      line.trim();
      if(line.empty())
        continue;
      if(line.front() == '#')
        continue;

      v = line.split_by_whitespaces();
      if(int(v.size()) != world_dim)
      {
        std::cerr << "ERROR: found " << v.size() << " tokens in line " << line_no << ", but expected " << world_dim << std::endl;
        return 1;
      }

      for(std::size_t i(0); i < std::size_t(world_dim); ++i)
      {
        Real t(0.0);
        if(!v[i].parse(t))
        {
          std::cerr << "ERROR: failed to parse vertex coord '" << v[i] << "' in line " << line_no << std::endl;
          return 1;
        }
        vtx.push_back(t);
      }
      ++count;
    }

    ifs.close();
    return 0;
  }

  int parse_idx(std::vector<Index>& idx, const String& filename, int const num_corners, const Index num_vertices)
  {
    idx.reserve(10000u * std::size_t(num_corners));

    // try to open the file
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.is_open() || !ifs.good())
    {
      std::cerr << "ERROR: failed to open index file '" << filename << "'" << std::endl;
      return 1;
    }

    std::cout << "Parsing index file '" << filename << "'..." << std::endl;

    String line;
    Index line_no(0), count(0);
    std::deque<String> v;
    while(!ifs.eof())
    {
      ++line_no;
      std::getline(ifs, line);
      line.trim();
      if(line.empty())
        continue;
      if(line.front() == '#')
        continue;

      v = line.split_by_whitespaces();
      if(int(v.size()) != num_corners)
      {
        std::cerr << "ERROR: found " << v.size() << " tokens in line " << line_no << ", but expected " << num_corners << std::endl;
        return 1;
      }

      for(std::size_t i(0); i < std::size_t(num_corners); ++i)
      {
        Index t(0);
        if(!v[i].parse(t))
        {
          std::cerr << "ERROR: failed to parse index '" << v[i] << "' in line " << line_no << std::endl;
          return 1;
        }
        if(t >= num_vertices)
        {
          std::cerr << "ERROR: vertex index " << t << " out of bounds in line " << line_no << ", expected < " << num_vertices << std::endl;
          return 1;
        }
        idx.push_back(t);
      }
      ++count;
    }

    ifs.close();
    return 0;
  }

  template<typename VtxSet_>
  void build_vertex_set(VtxSet_& vtx_set, const std::vector<Real>& vtx)
  {
    const Index num_verts = vtx_set.get_num_vertices();
    const int num_coords = vtx_set.get_num_coords();

    auto it = vtx.begin();
    for(Index i(0); i < num_verts; ++i)
    {
      for(int j(0); j < num_coords; ++j, ++it)
      {
        vtx_set[i][j] = *it;
      }
    }
  }

  template<typename IdxSet_>
  void build_index_set(IdxSet_& idx_set, const std::vector<Index>& idx)
  {
    const Index num_elems = idx_set.get_num_entities();
    const int num_idx = idx_set.get_num_indices();

    auto it = idx.begin();
    for(Index i(0); i < num_elems; ++i)
    {
      for(int j(0); j < num_idx; ++j, ++it)
      {
        idx_set(i, j) = *it;
      }
    }
  }

#ifdef FEAT_HAVE_FPARSER
  template<typename MeshType_>
  bool build_meshparts(Geometry::RootMeshNode<MeshType_>& mesh_node, const std::deque<String>& pts)
  {
    if((pts.size() % 2u) != 0u)
    {
      std::cerr << "ERROR: invalid number of parameters for option --parts, expected an even count" << std::endl;
      return false;
    }

    auto it = pts.begin();
    auto it_end = pts.end();
    while(it != it_end)
    {
      // get part name and formula
      String name = *it;
      ++it;
      String formula = *it;
      ++it;

      std::cout << "Creating meshpart '" << name << "' by hit-test formula '" << formula << "'..." << std::endl;

      try
      {
        // try to parse the formula
        Geometry::ParsedHitTestFactory<MeshType_, MeshType_::world_dim> factory(*mesh_node.get_mesh());

        // parse formula
        factory.parse(formula);

        // create mesh part
        auto mesh_part = factory.make_unique();

        // add to mesh node
        mesh_node.add_mesh_part(name, std::move(mesh_part));
      }
      catch(std::exception& exc)
      {
        std::cerr << "ERROR: in boundary function formula '" << formula << "'" << std::endl;
        std::cerr << exc.what() << std::endl;
        return false;
      }
    }

    // okay
    return true;
  }
#endif

  template<typename MeshType_>
  int run_shape(SimpleArgParser& args, std::vector<Real>& vtx, std::vector<Index>& idx)
  {
    typedef typename MeshType_::ShapeType ShapeType;
    static constexpr int shape_dim = MeshType_::shape_dim;
    static constexpr int world_dim = MeshType_::world_dim;
    static constexpr int num_corners = Shape::FaceTraits<ShapeType, 0>::count;

    const Index num_verts = Index(vtx.size()) / Index(world_dim);
    const Index num_elems = Index(idx.size()) / Index(num_corners);

    Index num_entities[] = {0, 0, 0, 0};
    num_entities[0] = num_verts;
    num_entities[shape_dim] = num_elems;

    // create a root mesh node
    Geometry::RootMeshNode<MeshType_> mesh_node(std::unique_ptr<MeshType_>(new MeshType_(num_entities)));
    MeshType_& mesh = *mesh_node.get_mesh();

    // build vertex set
    build_vertex_set(mesh.get_vertex_set(), vtx);

    // build vertices-at-element index set
    build_index_set(mesh.template get_index_set<shape_dim, 0>(), idx);

    // deduct topology
    std::cout << "Deducting mesh topology..." << std::endl;
    mesh.deduct_topology_from_top();

    // write final dimensions
    std::cout << "Mesh entity counts:";// << std::endl;
    for(int i(0); i <= shape_dim; ++i)
      std::cout << " " << mesh.get_num_entities(i);
      //std::cout << "Dimension " << i << ": " << mesh.get_num_entities(i) << std::endl;
    std::cout << std::endl;

    // create boundary mesh part
    if(args.check("bnd") >= 0)
    {
      std::cout << "Creating boundary meshpart..." << std::endl;
      Geometry::BoundaryFactory<MeshType_> bnd_factory(mesh);
      mesh_node.add_mesh_part("bnd", bnd_factory.make_unique());
    }

#ifdef FEAT_HAVE_FPARSER
    if(args.check("parts") > 0)
    {
      std::deque<String> parts = args.query("parts")->second;
      if(!build_meshparts(mesh_node, parts))
        return 1;
    }
#else // no FEAT_HAVE_FPARSER
    if(args.check("parts") >= 0)
    {
      std::cerr << "ERROR: you need to compile and link with 'fparser' to build mesh parts" << std::endl;
      return 1;
    }
#endif // FEAT_HAVE_FPARSER

    String out_name;
    if(args.parse("out", out_name) < 1)
    {
      std::cerr << "ERROR: mandatory output mesh file must be specified via --out <filename>" << std::endl;
      return 1;
    }

    // write mesh file
    std::cout << "Writing output mesh file '" << out_name << "'..." << std::endl;

    std::ofstream ofs(out_name, std::ios_base::out);
    if(!ofs.is_open() || !ofs.good())
    {
      std::cerr << "ERROR: Failed to open output file '" << out_name << "'" << std::endl;
      return 1;
    }

    // create mesh file writer and write
    Geometry::MeshFileWriter mesh_writer(ofs);

    mesh_writer.write(&mesh_node);

    // done
    return 0;
  }

  int main(int argc, char** argv)
  {
    typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
    typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, Real> S2M3D;
    typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
    typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
    typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, Real> H2M3D;
    typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

    SimpleArgParser args(argc, argv);

    args.support("help");
    args.support("out");
    args.support("vtx");
    args.support("idx");
    args.support("shape");
    args.support("bnd");
    args.support("parts");

    // need help?
    if((argc < 2) || (args.check("help") > -1))
    {
      display_help();
      return 0;
    }

    String vtx_file, idx_file;
    if(args.parse("vtx", vtx_file) < 1)
    {
      std::cerr << "ERROR: mandatory vertex coordinate input file must be specified via --vtx <filename>" << std::endl;
      display_help();
      return 1;
    }
    if(args.parse("idx", idx_file) < 1)
    {
      std::cerr << "ERROR: mandatory element indices input file must be specified via --idx <filename>" << std::endl;
      display_help();
      return 1;
    }

    String shape_type;
    if(args.parse("shape", shape_type) < 1)
    {
      std::cerr << "ERROR: mandatory shape type option --shape <shapename> is missing" << std::endl;
      display_help();
      return 1;
    }

    int world_dim = 0, num_corners = 0;

    // determine world dimension and shape dimension
    if(shape_type.compare_no_case("h2") == 0)
    {
      world_dim = 2;
      num_corners = 4;
    }
    if(shape_type.compare_no_case("h23") == 0)
    {
      world_dim = 3;
      num_corners = 4;
    }
    if(shape_type.compare_no_case("h3") == 0)
    {
      world_dim = 3;
      num_corners = 8;
    }
    if(shape_type.compare_no_case("s2") == 0)
    {
      world_dim = 2;
      num_corners = 3;
    }
    if(shape_type.compare_no_case("s23") == 0)
    {
      world_dim = 3;
      num_corners = 3;
    }
    if(shape_type.compare_no_case("s3") == 0)
    {
      world_dim = 3;
      num_corners = 4;
    }

    // parse vertex file
    std::vector<Real> vtx;
    int rtn = parse_vtx(vtx, vtx_file, world_dim);
    if(rtn != 0)
      return rtn;

    // compute number of parsed vertices
    const Index num_vertices = Index(vtx.size()) / Index(world_dim);
    std::cout << "Parsed " << num_vertices << " vertices" << std::endl;

    // parse index file
    std::vector<Index> idx;
    rtn = parse_idx(idx, idx_file, num_corners, num_vertices);
    if(rtn != 0)
      return rtn;

    // compute number of parsed elements
    const Index num_elements = Index(idx.size()) / Index(num_corners);
    std::cout << "Parsed " << num_elements << " elements" << std::endl;

    if(shape_type.compare_no_case("h2") == 0)
      return run_shape<H2M2D>(args, vtx, idx);
    if(shape_type.compare_no_case("h23") == 0)
      return run_shape<H2M3D>(args, vtx, idx);
    if(shape_type.compare_no_case("h3") == 0)
      return run_shape<H3M3D>(args, vtx, idx);
    if(shape_type.compare_no_case("s2") == 0)
      return run_shape<S2M2D>(args, vtx, idx);
    if(shape_type.compare_no_case("s23") == 0)
      return run_shape<S2M3D>(args, vtx, idx);
    if(shape_type.compare_no_case("s3") == 0)
      return run_shape<S3M3D>(args, vtx, idx);

    std::cerr << "ERROR: invalid shape type '" << shape_type << "'" << std::endl;
    display_help();
    return 1;
  }
} // namespace MeshIndexerTool

int main(int argc, char** argv)
{
  FEAT::Runtime::initialize(argc, argv);
  int rtn = MeshIndexerTool::main(argc, argv);
  FEAT::Runtime::finalize();
  return rtn;
}
