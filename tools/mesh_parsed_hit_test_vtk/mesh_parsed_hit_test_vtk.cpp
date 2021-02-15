// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/parsed_hit_test_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/string.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

static void display_help()
{
  std::cout << std::endl;
  std::cout << "mesh2vtk: Converts a mesh from FEAT format to VTK format" << std::endl;
  std::cout << std::endl;
  std::cout << "Mandatory arguments:" << std::endl;
  std::cout << "--------------------" << std::endl;
  std::cout << " --mesh <path to mesh file(s)>" << std::endl;
  std::cout << "Specifies the sequence of input mesh files paths." << std::endl;
  std::cout << std::endl;
  std::cout << " --function <string with function formula>" << std::endl;
  std::cout << "Specifies the hit function for creating a MeshPart." << std::endl;
  std::cout << std::endl;
  std::cout << "Optional arguments:" << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << " --vtk <path to vtk file>" << std::endl;
  std::cout << "Specifies the path of the output VTK file." << std::endl;
  std::cout << "If not given, the name of the first input mesh file is used." << std::endl;
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

template<typename Mesh_>
int run_xml(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader, const String& filename)
{
  static constexpr int world_dim = Mesh_::world_dim;

  // parse levels
  Index lvl_min(0);
  Index lvl_max(0);
  if(args.parse("level", lvl_max, lvl_min) < 2)
    lvl_min = lvl_max;

  // create an empty atlas and a root mesh node
  Geometry::MeshAtlas<Mesh_>* atlas = new Geometry::MeshAtlas<Mesh_>();
  Geometry::RootMeshNode<Mesh_>* node = new Geometry::RootMeshNode<Mesh_>(nullptr, atlas);

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
      auto* old = node;
      node = old->refine();
      delete old;
    }

    if(lvl < lvl_min)
      continue;

    // get our mesh
    Mesh_& mesh = *node->get_mesh();

    //Parse function formula
    String formula(""); //#?# runtime::abord()?
    if (!(args.parse("function", formula) == 1))
      std::cerr << "ERROR: The given formula could not get parsed into fparser!" << std::endl;

    // create a trafo for our mesh
    Trafo::Standard::Mapping<Mesh_> trafo(mesh);

    // Create a VTK exporter for our mesh
    FEAT::String vtkname = filename + "." + stringify(lvl);
    std::cout << "Writing file '" << vtkname << ".vtu'..." << std::endl;
    Geometry::ExportVTK<Mesh_> exporter(mesh);

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // create mesh part
    Geometry::ParsedHitTestFactory<Mesh_, world_dim> hit_test(mesh, formula); // #?#
    MeshPart<Mesh_> mesh_part(hit_test);

    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // get the vertex target set
    TargetSet& trg = mesh_part.template get_target_set<0>();

    std::vector<double> vtx_data(mesh.get_num_entities(0), 0.0);

    // mark all indexes vertices
    for(Index i(0); i < trg.get_num_entities(); ++i)
      vtx_data[trg[i]] = 1.0;

    // add variable
    exporter.add_vertex_scalar("meshpart", vtx_data.data());

    exporter.write(vtkname);
  }

  delete node;
  delete atlas;

  return 0;
}

int run(int argc, char* argv[])
{
  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, Real> S2M3D;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, Real> H1M1D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, Real> H1M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, Real> H1M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, Real> H2M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

  SimpleArgParser args(argc, argv);

  // need help?
  if((argc < 2) || (args.check("help") > -1))
  {
    display_help();
    return 0;
  }

  args.support("mesh");
  args.support("vtk");
  args.support("level");
  args.support("function");

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

  // our VTK output filename
  String vtk_name = get_file_title(filenames.back());
  args.parse("vtk", vtk_name);

  // create an empty mesh file reader
  Geometry::MeshFileReader mesh_reader;
  mesh_reader.add_mesh_files(filenames);

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  std::cout << "Mesh Type: " << mtype << std::endl;

  if(mtype == "conformal:hypercube:1:1")
    return run_xml<H1M1D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:1:2")
    return run_xml<H1M2D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:1:3")
    return run_xml<H1M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:2:2")
    return run_xml<H2M2D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:2:3")
    return run_xml<H2M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:3:3")
    return run_xml<H3M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:simplex:2:2")
    return run_xml<S2M2D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:simplex:2:3")
    return run_xml<S2M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:simplex:3:3")
    return run_xml<S3M3D>(args, mesh_reader, vtk_name);

  std::cout << "ERROR: unsupported mesh type!" << std::endl;

  return 1;
}

int main(int argc, char* argv[])
{
  Runtime::initialize(argc, argv);
  int ret = run(argc, argv);
  Runtime::finalize();
  return ret;
}
