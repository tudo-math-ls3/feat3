// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_svg.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/runtime.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

static void display_help()
{
  std::cout << std::endl;
  std::cout << "mesh2svg: Converts a mesh from FEAT format to SVG format" << std::endl;
  std::cout << std::endl;
  std::cout << "Mandatory arguments:" << std::endl;
  std::cout << "--------------------" << std::endl;
  std::cout << " --mesh <path to mesh file(s)>" << std::endl;
  std::cout << "Specifies the sequence of input mesh files paths." << std::endl;
  std::cout << std::endl;
  std::cout << "Optional arguments:" << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << " --svg <path to svg file>"  << std::endl;
  std::cout << "Specifies the path of the output SVG file." << std::endl;
  std::cout << "If not given, the name of the first input mesh file is used." << std::endl;
  std::cout << std::endl;
  std::cout << " --level [lvl_max [lvl_min]]" << std::endl;
  std::cout << "Specifies the minimum and maximum refinement levels." << std::endl;
  std::cout << "If not given, only level 0 is processed." << std::endl;
  std::cout << std::endl;
  std::cout << " --box [width height]" << std::endl;
  std::cout << "Specifies the bounding box of the figure in millimeters." << std::endl;
  std::cout << "If not given, a bounding box of 100 x 100 millimeters is used." << std::endl;
  std::cout << std::endl;
  std::cout << " --stroke [width]" << std::endl;
  std::cout << "Specifies the stroke width of the edges in millimeters." << std::endl;
  std::cout << "If not given, a stroke width of 0.1 millimeters is used." << std::endl;
  std::cout << std::endl;
  std::cout << " --extra [offset]" << std::endl;
  std::cout << "Specifies the extra offset of the figure in millimeters." << std::endl;
  std::cout << "If not given, an offset of 0.5 millimeters is used." << std::endl;
  std::cout << std::endl;
  std::cout << " --no-adapt" << std::endl;
  std::cout << "Do not adapt mesh after refinement." << std::endl;
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
  Index lvl_min(0), lvl_max(0);
  double box_x(100.0), box_y(100.0), stroke(0.1), extra(0.5);

  bool adapt = (args.check("no-adapt") < 0);

  args.parse("level", lvl_max, lvl_min);
  args.parse("box", box_x, box_y);
  args.parse("stroke", stroke);
  args.parse("extra", extra);

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
  if(adapt)
    node->adapt();

  // choose adapt mode
  const AdaptMode adapt_mode = adapt ? AdaptMode::chart : AdaptMode::none;

  // refine
  for(Index lvl(0); lvl <= lvl_max; ++lvl)
  {
    if(lvl > 0)
    {
      std::cout << "Refining up to level " << lvl << "..." << std::endl;
      node = node->refine_unique(adapt_mode);
    }

    if(lvl < lvl_min)
      continue;

    FEAT::String svgname = filename + "." + stringify(lvl) + ".svg";
    std::cout << "Writing file '" << svgname << "'..." << std::endl;

    Geometry::ExportSVG::write(svgname, *node->get_mesh(), box_x, box_y, stroke, extra);
  }

  return 0;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;

  SimpleArgParser args(argc, argv);

  // need help?
  if(args.check("help") > -1)
  {
    display_help();
    return 0;
  }

  args.support("mesh");
  args.support("level");
  args.support("svg");
  args.support("no-adapt");
  args.support("box");
  args.support("stroke");
  args.support("extra");

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

  // our SVG output filename
  String svg_name = get_file_title(filenames.back());
  args.parse("svg", svg_name);

  // create the mesh file reader
  Geometry::MeshFileReader mesh_reader;
  mesh_reader.add_mesh_files(filenames);

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  std::cout << "Mesh Type: " << mtype << std::endl;

  if(mtype == "conformal:hypercube:2:2")
    return run_xml<H2M2D>(args, mesh_reader, svg_name);
  if(mtype == "conformal:simplex:2:2")
    return run_xml<S2M2D>(args, mesh_reader, svg_name);

  std::cout << "ERROR: unsupported mesh type!" << std::endl;

  return 1;
}
