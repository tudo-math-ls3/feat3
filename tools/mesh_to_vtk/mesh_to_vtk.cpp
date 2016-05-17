#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAST;
using namespace FEAST::Geometry;

static void display_help()
{
  std::cout << "mesh2vtk: Reading meshes from files" << std::endl;
  std::cout << "Mandatory arguments:" << std::endl;
  std::cout << " --meshfile [path to mesh file]" << std::endl;
  std::cout << "Optional arguments:" << std::endl;
  std::cout << " --chartfile [path to chart file]" <<  std::endl;
  std::cout << " --level [lvl_max lvl_min]" << std::endl;
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
int run_xml(Geometry::MeshFileReader* mesh_reader, Geometry::MeshFileReader* chart_reader, const String& filename,
Index lvl_min, Index lvl_max)
{
  // create an empty atlas and a root mesh node
  Geometry::MeshAtlas<Mesh_>* atlas = new Geometry::MeshAtlas<Mesh_>();
  Geometry::RootMeshNode<Mesh_>* node = new Geometry::RootMeshNode<Mesh_>(nullptr, atlas);

  ASSERT_(mesh_reader != nullptr);

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
  node->adapt();

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
    Geometry::ExportVTK<Mesh_> exporter(*node->get_mesh());

    std::vector<double> vtx_data(node->get_mesh()->get_num_entities(0), 0.0);

    // loop over all mesh parts
    for(auto it = part_names.begin(); it != part_names.end(); ++it)
    {
      // get the mesh part
      MeshPart<Mesh_>* part = node->find_mesh_part(*it);
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

    // For every chart in the atlas, compute the distance of every mesh vertex to it
    if(atlas != nullptr)
    {
      const auto& vtx = node->get_mesh()->get_vertex_set();

      typename Mesh_::CoordType* distances
        (new typename Mesh_::CoordType[vtx.get_num_vertices()]);

      for(const auto& it:atlas->get_mesh_chart_map())
      {
        String fieldname("dist_"+it.first);
        for(Index i(0); i < vtx.get_num_vertices(); ++i)
          distances[i] = it.second->dist(vtx[i]);

        exporter.add_scalar_vertex(fieldname, distances);
      }
      delete[] distances;
    }

    exporter.write(vtkname);
  }

  delete node;
  delete atlas;

  return 0;
}

int main(int argc, char* argv[])
{
  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, 3, Real> S2M3D;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, Real> S3M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, 1, Real> H1M1D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, 2, Real> H1M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, 3, Real> H1M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> H2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, 3, Real> H2M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, Real> H3M3D;

  Runtime::initialise(argc, argv);

  SimpleArgParser args(argc, argv);

  args.support("meshfile");
  args.support("chartfile");
  args.support("level");

  String chart_file_name("");
  String mesh_file_name("");

  Index lvl_min(0);
  Index lvl_max(0);

  Geometry::MeshFileReader* mesh_reader(nullptr);
  Geometry::MeshFileReader* chart_reader(nullptr);

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() || args.check("help") > -1)
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;

    display_help();

    Runtime::abort();
    return 1;
  }

  if(args.check("meshfile")!=1)
  {
    display_help();
    throw InternalError(__func__,__FILE__,__LINE__, "You have to specify a mesh with --meshfile");
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
      std::cerr << "ERROR: Failed to open chart file '" << chart_file_name << "'" << std::endl;
      Runtime::abort();
      return 1;
    }
    std::cout << "Reading chart file from " << chart_file_name << std::endl;
    chart_reader = new Geometry::MeshFileReader(*ifs_chart);
  }

  if(args.check("level") > 0)
  {
    args.parse("level", lvl_max, lvl_min);
  }

  int ret(1);
  if(mtype == "conformal:hypercube:1:1")
    ret = run_xml<H1M1D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:1:2")
    ret = run_xml<H1M2D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:1:3")
    ret = run_xml<H1M3D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:2:2")
    ret = run_xml<H2M2D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:2:3")
    ret = run_xml<H2M3D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:3:3")
    ret = run_xml<H3M3D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);
  if(mtype == "conformal:simplex:2:2")
    ret = run_xml<S2M2D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);
  if(mtype == "conformal:simplex:2:3")
    ret = run_xml<S2M3D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);
  if(mtype == "conformal:simplex:3:3")
    ret = run_xml<S3M3D>(mesh_reader, chart_reader, mesh_file_name, lvl_min, lvl_max);

  // Clean up
  delete ifs_chart;
  delete mesh_reader;
  delete chart_reader;

  Runtime::finalise();

  return ret;

}
