#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_eps.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

static void display_help()
{
  std::cout << "mesh2eps: Converts a mesh from FEAT format to EPS format" << std::endl;
  std::cout << std::endl;
  std::cout << "Mandatory arguments:" << std::endl;
  std::cout << " --meshfile [path to mesh file]" << std::endl;
  std::cout << "Optional arguments:" << std::endl;
  std::cout << " --chartfile [path to chart file]" <<  std::endl;
  std::cout << " --level [lvl_max lvl_min]" << std::endl;
  std::cout << " --box [width height]" << std::endl;
  std::cout << " --stroke [width]" << std::endl;
  std::cout << " --no-adapt" << std::endl;
  std::cout << " --help: Displays this message" << std::endl;
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
int run_xml(Geometry::MeshFileReader* mesh_reader, Geometry::MeshFileReader* chart_reader,
  const String& filename, Index lvl_min, Index lvl_max, bool adapt, double box_x, double box_y, double stroke)
{
  // create an empty atlas and a root mesh node
  Geometry::MeshAtlas<Mesh_>* atlas = new Geometry::MeshAtlas<Mesh_>();
  Geometry::RootMeshNode<Mesh_>* node = new Geometry::RootMeshNode<Mesh_>(nullptr, atlas);

  XASSERT(mesh_reader != nullptr);

  // try to parse the mesh file
#ifndef DEBUG
  try
#endif
  {
    // Parse the chart file first, if it was given
    if(chart_reader != nullptr)
      chart_reader->parse(*node, *atlas);

    // Now parse the mesh file
    mesh_reader->parse(*node, *atlas);
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
      auto* old = node;
      node = old->refine(adapt_mode);
      delete old;
    }

    if(lvl < lvl_min)
      continue;

    FEAT::String epsname = get_file_title(filename) + "." + stringify(lvl) + ".eps";
    std::cout << "Writing file '" << epsname << "'..." << std::endl;

    Geometry::ExportEPS::write(epsname, *node->get_mesh(), box_x, box_y, stroke);
  }

  delete node;
  delete atlas;

  return 0;
}

int main(int argc, char* argv[])
{
  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> H2M2D;

  Runtime::initialise(argc, argv);

  SimpleArgParser args(argc, argv);

  // need help?
  if(args.check("help") > -1)
  {
    display_help();
    return Runtime::finalise();
  }

  args.support("meshfile");
  args.support("chartfile");
  args.support("level");
  args.support("no-adapt");
  args.support("box");
  args.support("stroke");

  String chart_file_name("");
  String mesh_file_name("");

  Index lvl_min(0);
  Index lvl_max(0);

  double box_x(100.0), box_y(100.0), stroke(0.1);

  Geometry::MeshFileReader* mesh_reader(nullptr);
  Geometry::MeshFileReader* chart_reader(nullptr);

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

  if(args.check("meshfile")!=1)
  {
    std::cerr << "ERROR: You have to specify a mesh with --meshfile" << std::endl;
    display_help();
    return 1;
  }

  args.parse("meshfile", mesh_file_name);
  std::ifstream ifs_mesh(mesh_file_name, std::ios_base::in);

  if(!ifs_mesh.is_open() || !ifs_mesh.good())
  {
    std::cerr << "ERROR: Failed to open mesh file '" << mesh_file_name << "'" << std::endl;
    return 1;
  }
  std::cout << "Reading mesh from " << mesh_file_name << std::endl;

  // create a reader and read the root markup
  mesh_reader = new Geometry::MeshFileReader(ifs_mesh);
  // create input file stream
  mesh_reader->read_root_markup();
  // get mesh type
  const String mtype = mesh_reader->get_meshtype_string();

  // Create a MeshFileReader for the optional chart file
  std::ifstream* ifs_chart(nullptr);
  if(args.check("chartfile")==1)
  {
    args.parse("chartfile", chart_file_name);
    ifs_chart = new std::ifstream(chart_file_name, std::ios_base::in);
    if(!ifs_chart->is_open() || !ifs_chart->good())
    {
      std::cerr << std::endl << "ERROR: Failed to open chart file '" << chart_file_name << "'" << std::endl;
      delete mesh_reader;
      return 1;
    }
    std::cout << "Reading chart file from " << chart_file_name << std::endl;
    chart_reader = new Geometry::MeshFileReader(*ifs_chart);
  }

  if(args.check("level") > 0)
  {
    args.parse("level", lvl_max, lvl_min);
  }

  bool adapt = (args.check("no-adapt") < 0);

  args.parse("box", box_x, box_y);
  args.parse("stroke", stroke);

  int ret(1);
  if(mtype == "conformal:hypercube:2:2")
    ret = run_xml<H2M2D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max, adapt, box_x, box_y, stroke);
  if(mtype == "conformal:simplex:2:2")
    ret = run_xml<S2M2D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max, adapt, box_x, box_y, stroke);

  // Clean up
  if(ifs_chart != nullptr)
    delete ifs_chart;
  if(mesh_reader != nullptr)
    delete mesh_reader;
  if(chart_reader != nullptr)
    delete chart_reader;

  Runtime::finalise();

  return ret;

}
