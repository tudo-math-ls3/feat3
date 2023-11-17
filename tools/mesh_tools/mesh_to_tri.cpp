// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/boundary_factory.hpp>

using namespace FEAT;

static void display_help()
{
  std::cout << std::endl;
  std::cout << "mesh2tri: Converts a mesh from FEAT 3 format to legacy FEATFLOW 1/2 TRI format" << std::endl;
  std::cout << std::endl;
  std::cout << "Mandatory arguments:" << std::endl;
  std::cout << "--------------------" << std::endl;
  std::cout << " --mesh <path to mesh file(s)>" << std::endl;
  std::cout << "Specifies the sequence of input mesh files paths." << std::endl;
  std::cout << std::endl;
  std::cout << "Optional arguments:" << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << " --tri <path to tri file>" << std::endl;
  std::cout << "Specifies the path of the output TRI file." << std::endl;
  std::cout << "If not given, the name of the first input mesh file is used." << std::endl;
  std::cout << std::endl;
  //            ---------1---------2---------3---------4---------5---------6---------7---------8---------9---------0
  //           "123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
  std::cout << " --knpr [<mesh-part-names...>]" << std::endl;
  std::cout << "Specifies how the nodal property array KNPR is to be defined." << std::endl;
  std::cout << "If not given, then the entire nodal property array is formatted to 0." << std::endl;
  std::cout << "If given without any mesh-part names, then the nodal property is set to 0" << std::endl;
  std::cout << "for all interior vertices and to 1 for all boundary vertices." << std::endl;
  std::cout << "If at least one mesh-part name is given, then all vertices belonging to the" << std::endl;
  std::cout << "first mesh-part get the nodal property 1, all vertices belonging to the second" << std::endl;
  std::cout << "mesh-part get the nodal property 2, etc., and all remaining vertices get the" << std::endl;
  std::cout << "nodal property 0. If multiple mesh-parts are to be combined so that all their" << std::endl;
  std::cout << "vertices belong to the same nodal property group, then you can specify that" << std::endl;
  std::cout << "group of mesh-parts by specifying their names in a double-quoted string, e.g." << std::endl;
  std::cout << "  --knpr bnd:l \"bnd:t bnd:b\" bnd:r" << std::endl;
  std::cout << "will yield 3 nodal property groups:" << std::endl;
  std::cout << "bnd:l -> 1, bnd:t -> 2, bnd:b -> 2, bnd:r -> 3, the rest -> 0" << std::endl;
  std::cout << std::endl;
  std::cout << " --level [lvl_max [lvl_min]]" << std::endl;
  std::cout << "Specifies the minimum and maximum refinement levels." << std::endl;
  std::cout << "If not given, only level 0 is processed." << std::endl;
  std::cout << std::endl;
  std::cout << " --help" << std::endl;
  std::cout << "Displays this message" << std::endl;
}

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

//template<typename Mesh_>
void mask_nodal_properties(std::vector<int>& node_prop, const Geometry::TargetSet& trg_set, int prop)
{
  for(Index i(0); i < trg_set.get_num_entities(); ++i)
    node_prop.at(trg_set[i]) = prop;
}

template<typename Mesh_>
int build_nodal_properties(SimpleArgParser& args, std::vector<int>& node_prop, const Geometry::RootMeshNode<Mesh_>& mesh_node)
{
  // nothing to do?
  if(args.check("knpr") < 0)
    return 0;

  // create nodal properties from entire boundary?
  if(args.check("knpr") == 0)
  {
    std::cout << "Building nodal property #1 from entire boundary..." << std::endl;
    Geometry::BoundaryFactory<Mesh_> bnd_factory(*mesh_node.get_mesh());
    Geometry::MeshPart<Mesh_> bnd_part(bnd_factory);
    mask_nodal_properties(node_prop, bnd_part.template get_target_set<0>(), 1);
    return 1;
  }

  // loop over all given mesh parts
  const std::deque<String> snpr = args.query("knpr")->second;
  for(std::size_t k(0); k < snpr.size(); ++k)
  {
    // split
    std::deque<String> part_names = snpr.at(k).split_by_whitespaces();

    std::cout << "Building nodal property #" << k << " from mesh part(s) '" <<
      stringify_join(part_names, "', '") << "'..." << std::endl;

    // loop over all mesh-parts
    for(const String& part_name : part_names)
    {
      const auto* part = mesh_node.find_mesh_part(part_name);
      if(part == nullptr)
      {
        std::cerr << "ERROR: Could not find mesh part named '" << part_name << "'" << std::endl;
        return -1;
      }

      mask_nodal_properties(node_prop, part->template get_target_set<0>(), int(k));
    }
  }

  return int(snpr.size());
}

void write_kvert(std::ostream& ofs, const Geometry::IndexSet<3>& idx_set, Shape::Simplex<2>)
{
  for(Index i(0); i < idx_set.get_num_entities(); ++i)
    ofs << (idx_set(i,0)+1u) << ' ' << (idx_set(i,1)+1u)  << ' ' << (idx_set(i,2)+1u) << '\n';
}

void write_kvert(std::ostream& ofs, const Geometry::IndexSet<4>& idx_set, Shape::Simplex<3>)
{
  for(Index i(0); i < idx_set.get_num_entities(); ++i)
    ofs << (idx_set(i,0)+1u) << ' ' << (idx_set(i,1)+1u)  << ' ' << (idx_set(i,2)+1u) << ' ' << (idx_set(i,3)+1u) << '\n';
}

void write_kvert(std::ostream& ofs, const Geometry::IndexSet<4>& idx_set, Shape::Hypercube<2>)
{
  // swap vertex 2 and 3
  // 2--3     3--2
  // |  |  -> |  |
  // 0--1     0--1
  for(Index i(0); i < idx_set.get_num_entities(); ++i)
    ofs << (idx_set(i,0)+1u) << ' ' << (idx_set(i,1)+1u) << ' ' << (idx_set(i,3)+1u) << ' ' << (idx_set(i,2)+1u) << '\n';
}

void write_kvert(std::ostream& ofs, const Geometry::IndexSet<8>& idx_set, Shape::Hypercube<3>)
{
  // swap vertex 2 and 3, and swap vertex 6 and 7
  for(Index i(0); i < idx_set.get_num_entities(); ++i)
    ofs << (idx_set(i,0)+1u) << ' ' << (idx_set(i,1)+1u) << ' ' << (idx_set(i,3)+1u) << ' ' << (idx_set(i,2)+1u) << ' '
        << (idx_set(i,4)+1u) << ' ' << (idx_set(i,5)+1u) << ' ' << (idx_set(i,7)+1u) << ' ' << (idx_set(i,6)+1u) << '\n';
}

template<typename Mesh_>
bool write_tri(const String& filename, const Geometry::RootMeshNode<Mesh_>& mesh_node, const std::vector<int>& node_prop, const int nbct, const String& mesh_filename)
{
  std::cout << "Writing '" << filename << "'...\n";
  std::ofstream ofs(filename);
  if(!ofs)
  {
    std::cerr << "ERROR: Failed to open '" << filename << "' for writing!" << std::endl;
    return false;
  }

  // write the two leading comment lines
  ofs << "Generated by FEAT 3 mesh3tri tool from '" << mesh_filename << "'\n";
  ofs << "The quick brown fox jumps over the lazy dog\n";

  typedef typename Mesh_::ShapeType ShapeType;

  const Mesh_& mesh = *mesh_node.get_mesh();
  static constexpr int world_dim = Mesh_::world_dim;

  // write entity counts
  if(Mesh_::shape_dim == 2)
  {
    ofs << " " << stringify(mesh.get_num_elements());
    ofs << " " << stringify(mesh.get_num_vertices());
    ofs << " 0";
    ofs << " " << Shape::FaceTraits<ShapeType, 0>::count;
    ofs << " " << nbct;
    ofs << " NEL NVT NMT NVE NBCT\n";
  }
  else // 3D
  {
    ofs << " " << stringify(mesh.get_num_elements());
    ofs << " " << stringify(mesh.get_num_vertices());
    ofs << " " << nbct;
    ofs << " " << Shape::FaceTraits<ShapeType, 0>::count;
    ofs << " " << Shape::FaceTraits<ShapeType, 1>::count;
    ofs << " " << Shape::FaceTraits<ShapeType, 2>::count;
    ofs << " NEL,NVT,NBCT,NVE,NEE,NAE\n";
  }

  // write vertices
  ofs << "DCORVG\n";
  const auto& vtx = mesh.get_vertex_set();
  for(Index i(0); i < vtx.get_num_vertices(); ++i)
  {
    const auto& v = vtx[i];
    ofs << v[0];
    for(int j(1); j < world_dim; ++j)
      ofs << ' ' << v[j];
    ofs << '\n';
  }

  // write vertex-at-element indices
  ofs << "KVERT\n";
  ShapeType shape_type;
  write_kvert(ofs, mesh.template get_index_set<ShapeType::dimension, 0>(), shape_type);

  // write nodal properties
  ofs << "KNPR\n";
  for(std::size_t i(0); i < node_prop.size(); ++i)
    ofs << node_prop[i] << '\n';

  ofs.close();
  return true;
}

template<typename Mesh_>
int run_xml(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader, const String& mesh_filename)
{
  // parse levels
  Index lvl_min(0);
  Index lvl_max(0);
  args.parse("level", lvl_max, lvl_min);

  // our output filename
  String out_name = get_file_title(mesh_filename);
  args.parse("tri", out_name);

  // create an empty atlas and a root mesh node
  auto atlas = Geometry::MeshAtlas<Mesh_>::make_unique();
  auto node = Geometry::RootMeshNode<Mesh_>::make_unique(nullptr, atlas.get());

  // try to parse the mesh file
#ifndef DEBUG
  try
#endif
  {
    std::cout << "Parsing mesh files..." << std::endl;
    // Now parse the mesh file
    mesh_reader.parse(*node, *atlas);
  }
#ifndef DEBUG
  catch(std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    return 1;
  }
#endif

  // adapt coarse mesh
  node->adapt();

  // refine
  for(Index lvl(0); lvl <= lvl_max; ++lvl)
  {
    if(lvl > 0)
    {
      std::cout << "Refining up to level " << lvl << "..." << std::endl;
      node = node->refine_unique();
    }

    if(lvl < lvl_min)
      continue;

    // build nodal properties
    std::vector<int> node_prop(node->get_mesh()->get_num_vertices(), 0);
    int nbct = build_nodal_properties(args, node_prop, *node);
    if(nbct < 0)
      return 1;

    // write TRI file
    String tri_name = out_name;
    if(lvl > 0u)
      (tri_name += ".") += stringify(lvl);
    write_tri(tri_name + ".tri", *node, node_prop, nbct, mesh_filename);
  }

  return 0;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

  SimpleArgParser args(argc, argv);

  // need help?
  if((argc < 2) || (args.check("help") > -1))
  {
    display_help();
    return 0;
  }

  args.support("mesh");
  args.support("level");
  args.support("tri");
  args.support("knpr");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;

    display_help();
    return 1;
  }

  int num_mesh_files = args.check("mesh");
  if(num_mesh_files < 1)
  {
    std::cerr << "ERROR: You have to specify at least one meshfile with --mesh <files...>" << std::endl;
    display_help();
    return 1;
  }

  // get our filename deque
  auto mpars = args.query("mesh");
  XASSERT(mpars != nullptr);
  const std::deque<String>& filenames = mpars->second;

  // create an empty mesh file reader
  Geometry::MeshFileReader mesh_reader;
  mesh_reader.add_mesh_files(filenames);

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  std::cout << "Mesh Type: " << mtype << std::endl;

  if(mtype == "conformal:hypercube:2:2")
    return run_xml<H2M2D>(args, mesh_reader, filenames.front());
  if(mtype == "conformal:hypercube:3:3")
    return run_xml<H3M3D>(args, mesh_reader, filenames.front());
  if(mtype == "conformal:simplex:2:2")
    return run_xml<S2M2D>(args, mesh_reader, filenames.front());
  if(mtype == "conformal:simplex:3:3")
    return run_xml<S3M3D>(args, mesh_reader, filenames.front());

  std::cout << "ERROR: unsupported mesh type!" << std::endl;

  return 1;
}
