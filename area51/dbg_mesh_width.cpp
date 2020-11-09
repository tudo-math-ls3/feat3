// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

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

  Tiny::Vector<double, 3> ray_3(0.0);
  ray_3[0] = 1.0;

  // parse levels
  Index lvl_min(0);
  Index lvl_max(0);
  args.parse("level", lvl_max, lvl_min);
  args.parse("ray", ray_3[0], ray_3[1], ray_3[2]);

  Tiny::Vector<double, world_dim> ray(ray_3.template size_cast<world_dim>());
  ray.normalize();

  std::cout << "Ray = [ " << stringify(ray[0]);
  for(int i(1); i < world_dim; ++i)
    std::cout << " , " << stringify(ray[i]);
  std::cout << " ]" << std::endl;

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
    mesh_reader.parse(*node, *atlas, nullptr);
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

    // create a trafo for our mesh
    Trafo::Standard::Mapping<Mesh_> trafo(mesh);
    typename decltype(trafo)::template Evaluator<>::Type trafo_eval(trafo);

    // Create a VTK exporter for our mesh
    FEAT::String vtkname = filename + "." + stringify(lvl);
    std::cout << "Writing file '" << vtkname << ".vtu'..." << std::endl;
    Geometry::ExportVTK<Mesh_> exporter(mesh);

    std::vector<double> volume(mesh.get_num_elements(), 0.0);
    std::vector<double> width(mesh.get_num_elements(), 0.0);

    for(Index i(0); i < mesh.get_num_elements(); ++i)
    {
      trafo_eval.prepare(i);
      volume[i] = trafo_eval.volume();
      //volume[i] = Math::pow(trafo_eval.volume(), 1.0 / double(world_dim));
      width[i] = trafo_eval.width_directed(ray);
      trafo_eval.finish();
    }

    exporter.add_cell_scalar("volume", volume.data());
    exporter.add_cell_scalar("width", width.data());

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
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

  SimpleArgParser args(argc, argv);

  // need help?
  if((argc < 2) || (args.check("help") > -1))
  {
    return 0;
  }

  args.support("mesh");
  args.support("vtk");
  args.support("level");
  args.support("ray");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;

    return 1;
  }

  int num_mesh_files = args.check("mesh");
  if(num_mesh_files < 1)
  {
    std::cerr << "ERROR: You have to specify at least one meshfile with --mesh <files...>" << std::endl;
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

  if(mtype == "conformal:hypercube:2:2") return run_xml<H2M2D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:hypercube:3:3") return run_xml<H3M3D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:simplex:2:2") return run_xml<S2M2D>(args, mesh_reader, vtk_name);
  if(mtype == "conformal:simplex:3:3") return run_xml<S3M3D>(args, mesh_reader, vtk_name);

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
