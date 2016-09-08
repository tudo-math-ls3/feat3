#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
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
  std::cout << std::endl;
  std::cout << "mesh2vtk: Converts a mesh from FEAT format to VTK format" << std::endl;
  std::cout << std::endl;
  std::cout << "Mandatory arguments:" << std::endl;
  std::cout << " --mesh [path to mesh file(s)]" << std::endl;
  std::cout << std::endl;
  std::cout << "Optional arguments:" << std::endl;
  std::cout << " --vtk [path to vtk file]" << std::endl;
  std::cout << " --level [lvl_max lvl_min]" << std::endl;
  std::cout << " --no-adapt: Do not adapt mesh after refinement" << std::endl;
  std::cout << " --no-dist: Do not compute distance to charts" << std::endl;
  std::cout << " --no-proj: Do not compute projection to charts" << std::endl;
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
int run_xml(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader, const String& filename)
{
  // parse levels
  Index lvl_min(0);
  Index lvl_max(0);
  args.parse("level", lvl_max, lvl_min);

  // check for adaption
  bool adapt = (args.check("no-adapt") < 0);

  // compute distance functions?
  bool calc_dist = (args.check("no-dist") < 0);

  // compute projection fields?
  bool calc_proj = (args.check("no-proj") < 0);

  // create an empty atlas and a root mesh node
  Geometry::MeshAtlas<Mesh_>* atlas = new Geometry::MeshAtlas<Mesh_>();
  Geometry::RootMeshNode<Mesh_>* node = new Geometry::RootMeshNode<Mesh_>(nullptr, atlas);
  Geometry::PartitionSet part_set;

  // try to parse the mesh file
#ifndef DEBUG
  try
#endif
  {
    std::cout << "Parsing mesh files..." << std::endl;
    // Now parse the mesh file
    mesh_reader.parse(*node, *atlas, &part_set);
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

  // get all mesh part names
  std::deque<String> part_names = node->get_mesh_part_names();

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

    // Create a VTK exporter for our mesh
    FEAT::String vtkname = filename + "." + stringify(lvl);
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
      exporter.add_vertex_scalar(*it, vtx_data.data());

      // unmark vertices
      for(Index i(0); i < trg.get_num_entities(); ++i)
        vtx_data[trg[i]] = 0.0;
    }

    // For every chart in the atlas, compute the distance of every mesh vertex to it
    if(calc_dist)
    {
      const auto& vtx = node->get_mesh()->get_vertex_set();

      std::vector<double> distances(vtx.get_num_vertices(), 0.0);

      for(const auto& it:atlas->get_mesh_chart_map())
      {
        for(Index i(0); i < vtx.get_num_vertices(); ++i)
          distances[i] = it.second->signed_dist(vtx[i]);

        exporter.add_vertex_scalar("dist:"+it.first, distances.data());
      }
    }

    // compute projection fields?
    if(calc_proj)
    {
      typedef typename Mesh_::VertexSetType VertexSetType;
      typedef typename VertexSetType::VertexType VertexType;

      const auto& vtx = node->get_mesh()->get_vertex_set();

      std::vector<double> prj_x(vtx.get_num_vertices(), 0.0);
      std::vector<double> prj_y(vtx.get_num_vertices(), 0.0);
      std::vector<double> prj_z(vtx.get_num_vertices(), 0.0);

      Tiny::Vector<double, 3> pt;
      pt[0] = pt[1] = pt[2] = 0.0;
      auto& wpt = pt.template size_cast<Mesh_::world_dim>();

      for(const auto& it : atlas->get_mesh_chart_map())
      {
        const auto& chart = *it.second;

        // skip charts which cannot perform implicit projection
        if(!chart.can_implicit())
          continue;

        for(Index i(0); i < vtx.get_num_vertices(); ++i)
        {
          // prject vertex
          wpt = chart.project(vtx[i]);
          // subtract vertex
          wpt -= vtx[i];
          // copy data
          prj_x[i] = pt[0];
          prj_y[i] = pt[1];
          prj_z[i] = pt[2];
        }

        exporter.add_vertex_vector("proj:" + it.first, prj_x.data(), prj_y.data(), prj_z.data());
      }
    }

    // loop over all partitions
    for(const auto& part : part_set.get_partitions())
    {
      // does this partition refer to the current level?
      if(Index(part.get_level()) != lvl)
        continue;

      // allocate and fill rank vector
      std::vector<double> r(std::size_t(node->get_mesh()->get_num_entities(Mesh_::shape_dim)), -1.0);
      for(Index i(0); i < part.size(); ++i)
      {
        auto it = part.get_patches().image_begin(i);
        auto jt = part.get_patches().image_end(i);
        for(; it != jt; ++it )
          r[*it] = double(i);
      }

      // build the name:
      String name = "partition:";
      if(!part.get_name().empty())
        name += part.get_name() + ":";
      name += stringify(part.size());

      exporter.add_cell_scalar(name, r.data());
    }

    exporter.write(vtkname);
  }

  delete node;
  delete atlas;

  return 0;
}

int run(int argc, char* argv[])
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

  SimpleArgParser args(argc, argv);

  // need help?
  if(args.check("help") > -1)
  {
    display_help();
    return 0;
  }

  args.support("mesh");
  args.support("vtk");
  args.support("level");
  args.support("no-adapt");
  args.support("no-dist");
  args.support("no-proj");

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
  std::deque<std::shared_ptr<std::ifstream>> streams;

  // our VTK output filename
  String vtk_name = get_file_title(filenames.back());
  args.parse("vtk", vtk_name);

  // create an empty mesh file reader
  Geometry::MeshFileReader mesh_reader;

  // open a stream for each filename
  for(auto it = filenames.begin(); it != filenames.end(); ++it)
  {
    // try to open the file
    streams.emplace_back(std::make_shared<std::ifstream>(*it, std::ios_base::in));
    std::ifstream& ifs = *streams.back();
    if(!ifs.is_open() || !ifs.good())
    {
      std::cerr << "ERROR: Failed to open mesh file '" << (*it) << "'" << std::endl;
      return 1;
    }

    // add stream to reader
    std::cout << "Adding '" << (*it) << "' to mesh file reader..." << std::endl;
    mesh_reader.add_stream(ifs);
  }

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  std::cout << "Mesh Type: " << mtype << std::endl;

  int ret(1);
  if(mtype == "conformal:hypercube:1:1")
    ret = run_xml<H1M1D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:1:2")
    ret = run_xml<H1M2D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:1:3")
    ret = run_xml<H1M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:2:2")
    ret = run_xml<H2M2D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:2:3")
    ret = run_xml<H2M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:3:3")
    ret = run_xml<H3M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:simplex:2:2")
    ret = run_xml<S2M2D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:simplex:2:3")
    ret = run_xml<S2M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:simplex:3:3")
    ret = run_xml<S3M3D>(args, mesh_reader, vtk_name);

  return ret;
}

int main(int argc, char* argv[])
{
  Runtime::initialise(argc, argv);
  int ret = run(argc, argv);
  Runtime::finalise();
  return ret;
}
