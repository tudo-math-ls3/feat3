// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/parsed_hit_test_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

static void display_help()
{
  std::cout << "\n";
  std::cout << "mesh2vtk: Converts a mesh from FEAT format to VTK format" << "\n";
  std::cout << "\n";
  std::cout << "Mandatory arguments:" << "\n";
  std::cout << "--------------------" << "\n";
  std::cout << " --mesh <path to mesh file(s)>" << "\n";
  std::cout << "Specifies the sequence of input mesh files paths." << "\n";
  std::cout << "\n";
  std::cout << "Optional arguments:" << "\n";
  std::cout << "-------------------" << "\n";
  std::cout << " --vtk <path to vtk file>" << "\n";
  std::cout << "Specifies the path of the output VTK file." << "\n";
  std::cout << "If not given, the name of the first input mesh file is used." << "\n";
  std::cout << "\n";
  std::cout << " --level [lvl_max [lvl_min]]" << "\n";
  std::cout << "Specifies the minimum and maximum refinement levels." << "\n";
  std::cout << "If not given, only level 0 is processed." << "\n";
  std::cout << "\n";
  std::cout << " --origin <x> [<y> [<z>]]" << "\n";
  std::cout << "Specifies the translation origin." << "\n";
  std::cout << "This translation is subtracted before applying the rotation." << "\n";
  std::cout << "If not given, all origins are set to 0." << "\n";
  std::cout << "\n";
  std::cout << " --angles <angle>" << "\n";
  std::cout << " --angles <yaw> <pitch> <roll>" << "\n";
  std::cout << "Specifies the rotation angle (2D) or the yaw-pitch-roll angles (3D)." << "\n";
  std::cout << "If not given, all angles are set to 0." << "\n";
  std::cout << "Note: all angles are given in revolutions (full turns):" << "\n";
  std::cout << "      1 revolution = 360 degrees = 2pi radians" << "\n";
  std::cout << "\n";
  std::cout << " --offset <x> [<y> [<z>]]" << "\n";
  std::cout << "Specifies the translation offset." << "\n";
  std::cout << "This translation is added after applying the rotation." << "\n";
  std::cout << "If not given, all offsets are set to 0." << "\n";
  std::cout << "\n";
  std::cout << " --no-adapt" << "\n";
  std::cout << "Do not adapt mesh after refinement." << "\n";
  std::cout << "\n";
  std::cout << " --no-dist" << "\n";
  std::cout << "Do not compute distance to charts." << "\n";
  std::cout << "\n";
  std::cout << " --no-proj" << "\n";
  std::cout << "Do not compute projection to charts" << "\n";
  std::cout << "\n";
  std::cout << " --no-volume" << "\n";
  std::cout << "Do not compute cell volumes" << "\n";
  std::cout << "\n";
  std::cout << " --hit-test <name1> <formula1> [<name2> <formula2> ...]" << "\n";
  std::cout << "Specifies a set of name-formula argument pairs which are used to generate meshparts" << "\n";
  std::cout << "by using the Geometry::ParsedHitTestFactory class. The first component of each pair" << "\n";
  std::cout << "specifies the name for the mesh part, whereas the second component specifies the formula" << "\n";
  std::cout << "in x,y,z coordinates, which is to be used for the hit test of the mesh part." << "\n";
  std::cout << "A vertex or edge/face/cell will be contained in the meshpart if the formula evaluates" << "\n";
  std::cout << "to a positive value in its coordinates or midpoint coordinates, respectively." << "\n";
  std::cout << "Please note that this option can only be used if FEAT is configured and linked against" << "\n";
  std::cout << "the 'fparser' third-party library." << "\n";
  std::cout << "\n";
  std::cout << " --help" << "\n";
  std::cout << "Displays this message" << "\n";
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

int parse_trafo(SimpleArgParser& args, Tiny::Vector<Real, 1, 1>& offset, Tiny::Vector<Real, 1, 1>& origin, Tiny::Vector<Real, 1, 1>& /*angles*/)
{
  int noff = args.parse("offset", offset[0]);
  int nori = args.parse("origin", origin[0]);
  if(args.check("angles") >= 0)
  {
    std::cerr << "ERROR: option '--angles' is invalid for 1D meshes" << "\n";
    return -1;
  }

  // parse error?
  if(noff < 0)
  {
    std::cerr << "ERROR: could not parse '--offset <x>' parameters" << "\n";
    return -1;
  }
  if(nori < 0)
  {
    std::cerr << "ERROR: could not parse '--origin <x>' parameters" << "\n";
    return -1;
  }

  // okay
  return ((noff + nori) > 0) ? 1 : 0;
}

int parse_trafo(SimpleArgParser& args, Tiny::Vector<Real, 2, 2>& offset, Tiny::Vector<Real, 2, 2>& origin, Tiny::Vector<Real, 2, 2>& angles)
{
  int noff = args.parse("offset", offset[0], offset[1]);
  int nori = args.parse("origin", origin[0], origin[1]);
  int nang = args.parse("angles", angles[0]);

  // parse error?
  if(noff < 0)
  {
    std::cerr << "ERROR: could not parse '--offset <x> <y>' parameters" << "\n";
    return -1;
  }
  if(nori < 0)
  {
    std::cerr << "ERROR: could not parse '--origin <x> <y>' parameters" << "\n";
    return -1;
  }
  if(nang < 0)
  {
    std::cerr << "ERROR: could not parse '--angles <angle>" << "\n";
    return -1;
  }

  // invalid argument count?
  if((noff > 0) && (noff != 2))
  {
    std::cerr << "ERROR: invalid number of '--offset <x> <y>' parameters" << "\n";
    return -1;
  }
  if((nori > 0) && (nori != 2))
  {
    std::cerr << "ERROR: invalid number of '--origin <x> <y>' parameters" << "\n";
    return -1;
  }
  if((nang > 0) && (nang != 1))
  {
    std::cerr << "ERROR: invalid number of '--angles <angle>' parameters" << "\n";
    return -1;
  }

  // okay
  return ((noff + nori + nang) > 0) ? 1 : 0;
}

int parse_trafo(SimpleArgParser& args, Tiny::Vector<Real, 3, 3>& offset, Tiny::Vector<Real, 3, 3>& origin, Tiny::Vector<Real, 3, 3>& angles)
{
  int noff = args.parse("offset", offset[0], offset[1], offset[2]);
  int nori = args.parse("origin", origin[0], origin[1], origin[2]);
  int nang = args.parse("angles", angles[0], angles[1], angles[2]);

  // parse error?
  if(noff < 0)
  {
    std::cerr << "ERROR: could not parse '--offset <x> <y> <z>' parameters" << "\n";
    return -1;
  }
  if(nori < 0)
  {
    std::cerr << "ERROR: could not parse '--origin <x> <y> <z>' parameters" << "\n";
    return -1;
  }
  if(nang < 0)
  {
    std::cerr << "ERROR: could not parse '--angles <yaw> <pitch> <roll>' parameters" << "\n";
    return -1;
  }

  // invalid argument count?
  if((noff > 0) && (noff != 3))
  {
    std::cerr << "ERROR: invalid number of '--offset <x> <y> <z>' parameters" << "\n";
    return -1;
  }
  if((nori > 0) && (nori != 3))
  {
    std::cerr << "ERROR: invalid number of '--origin <x> <y> <z>' parameters" << "\n";
    return -1;
  }
  if((nang > 0) && (nang != 3))
  {
    std::cerr << "ERROR: invalid number of '--angles <yaw> <pitch> <roll>' parameters" << "\n";
    return -1;
  }

  // okay
  return ((noff + nori + nang) > 0) ? 1 : 0;
}

#ifdef FEAT_HAVE_FPARSER
template<typename MeshType_>
bool build_meshparts(Geometry::RootMeshNode<MeshType_>& mesh_node, const std::map<String, String>& pts)
{
  for(auto it = pts.begin(); it != pts.end(); ++it)
  {
    // get part name and formula
    String name = it->first;
    String formula = it->second;

    std::cout << "Creating meshpart '" << name << "' by hit-test formula '" << formula << "'..." << "\n";

    try
    {
      // try to parse the formula
      Geometry::ParsedHitTestFactory<MeshType_> factory(*mesh_node.get_mesh());

      // parse formula
      factory.parse(formula);

      // make sure to erase any meshpart with that name
      mesh_node.remove_mesh_part(name);

      // create mesh part
      auto mesh_part = factory.make_unique();

      // add to mesh node
      mesh_node.add_mesh_part(name, std::move(mesh_part));
    }
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: in boundary function formula '" << formula << "'" << "\n";
      std::cerr << exc.what() << "\n";
      return false;
    }
  }

  // okay
  return true;
}
#endif

template<typename Mesh_>
int run_xml(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader, const String& filename)
{
  static constexpr int world_dim = Mesh_::world_dim;

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

  // compute cell volumes?
  bool calc_volume = (args.check("no-volume") < 0);

#ifdef FEAT_HAVE_FPARSER
  std::map<String, String> hit_formulae;
  if(args.check("hit-test") > 0)
  {
    std::deque<String> parts = args.query("hit-test")->second;
    if(parts.size() % 2u != 0u)
    {
      std::cerr << "ERROR: invalid number of parameters for option --hit-test, expected an even count" << "\n";
      return 1;
    }

    for(std::size_t i(0); i < parts.size(); i += 2u)
      hit_formulae.emplace(parts[i], parts[i+1u]);
  }
#else // no FEAT_HAVE_FPARSER
  if(args.check("hit-test") >= 0)
  {
    std::cerr << "ERROR: you need to compile and link with 'fparser' to build mesh parts" << "\n";
    return 1;
  }
#endif // FEAT_HAVE_FPARSER

  // create an empty atlas and a root mesh node
  auto atlas = Geometry::MeshAtlas<Mesh_>::make_unique();
  auto node = Geometry::RootMeshNode<Mesh_>::make_unique(nullptr, atlas.get());
  Geometry::PartitionSet part_set;

  // try to parse the mesh file
#ifndef DEBUG
  try
#endif
  {
    std::cout << "Parsing mesh files..." << "\n";
    // Now parse the mesh file
    mesh_reader.parse(*node, *atlas, &part_set);
  }
#ifndef DEBUG
  catch(std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << "\n";
    return 1;
  }
  catch(...)
  {
    std::cerr << "ERROR: unknown exception" << "\n";
    return 1;
  }
#endif

  // adapt coarse mesh
  if(adapt)
    node->adapt();

  // choose adapt mode
  const AdaptMode adapt_mode = adapt ? AdaptMode::chart : AdaptMode::none;

  // parse the transformation (if given)
  Tiny::Vector<Real, world_dim> origin(Real(0)), offset(Real(0)), angles(Real(0));
  int ntrafo = parse_trafo(args, offset, origin, angles);
  if(ntrafo < 0)
    return 1;
  else if(ntrafo > 0)
  {
    // fractions to radians
    angles *= Real(2) * Math::pi<Real>();

    // transform atlas and root mesh
    atlas->transform(origin, angles, offset);
    node->get_mesh()->transform(origin, angles, offset);
  }

  // refine
  for(Index lvl(0); lvl <= lvl_max; ++lvl)
  {
    if(lvl > 0)
    {
      std::cout << "Refining up to level " << lvl << "..." << "\n";
      node = node->refine_unique(adapt_mode);
    }

    if(lvl < lvl_min)
      continue;

    // build hit-test mesh-parts
#ifdef FEAT_HAVE_FPARSER
    if(!build_meshparts(*node, hit_formulae))
      return 1;
#endif

    // get all mesh part names
    std::deque<String> part_names = node->get_mesh_part_names();

    // get our mesh
    Mesh_& mesh = *node->get_mesh();

    // create a trafo for our mesh
    Trafo::Standard::Mapping<Mesh_> trafo(mesh);

    // Create a VTK exporter for our mesh
    FEAT::String vtkname = filename + "." + stringify(lvl);
    std::cout << "Writing file '" << vtkname << ".vtu'..." << "\n";
    Geometry::ExportVTK<Mesh_> exporter(mesh);

    std::vector<double> vtx_data(mesh.get_num_entities(0), 0.0);

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

    // compute cell volumes?
    if(calc_volume)
    {
      std::vector<double> volume(mesh.get_num_elements(), 0.0);

      for(Index i(0); i < mesh.get_num_elements(); ++i)
      {
        volume[i] = trafo.compute_vol(i);
      }

      exporter.add_cell_scalar("volume", volume.data());
    }

    // For every chart in the atlas, compute the distance of every mesh vertex to it
    if(calc_dist)
    {
      const auto& vtx = mesh.get_vertex_set();

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
      const auto& vtx = mesh.get_vertex_set();

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
      std::vector<double> r(std::size_t(mesh.get_num_entities(Mesh_::shape_dim)), -1.0);
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

  return 0;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

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
  args.support("no-adapt");
  args.support("no-dist");
  args.support("no-proj");
  args.support("no-volume");
  args.support("angles");
  args.support("offset");
  args.support("origin");
  args.support("hit-test");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << "\n";

    display_help();
    return 1;
  }

  int num_mesh_files = args.check("mesh");
  if(num_mesh_files < 1)
  {
    std::cerr << "ERROR: You have to specify at least one meshfile with --mesh <files...>" << "\n";
    display_help();
    return 1;
  }

  // get our filename deque
  auto mpars = args.query("mesh");
  XASSERT(mpars != nullptr);
  const std::deque<String>& filenames = mpars->second;

  // our VTK output filename
  String vtk_name = get_file_title(filenames.front());
  args.parse("vtk", vtk_name);

  // create an empty mesh file reader
  Geometry::MeshFileReader mesh_reader;
  mesh_reader.add_mesh_files(filenames);

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  std::cout << "Mesh Type: " << mtype << "\n";

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

  std::cout << "ERROR: unsupported mesh type!" << "\n";

  return 1;
}
