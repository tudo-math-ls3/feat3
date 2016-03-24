#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>

using namespace FEAST;
using namespace FEAST::Geometry;

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
int run_xml(Geometry::MeshFileReader& reader, const String& filename, Index lvl_min, Index lvl_max)
{
  // create an empty atlas and a root mesh node
  Geometry::MeshAtlas<Mesh_>* atlas = new Geometry::MeshAtlas<Mesh_>();
  Geometry::RootMeshNode<Mesh_>* node = new Geometry::RootMeshNode<Mesh_>(nullptr, atlas);

  // try to parse the mesh file
#ifndef DEBUG
  try
#endif
  {
    reader.parse(*node, *atlas);
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

    exporter.write(vtkname);
  }

  delete node;
  delete atlas;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    std::cout << "USAGE: " << argv[0] << " <meshfile> [<level-min>] [<level-max>]" << std::endl;
    return 0;
  }

  // Input file name, required
  FEAST::String filename(argv[1]);

  // Refinement level
  Index lvl_min(0);
  if((argc > 2) && (!String(argv[2]).parse(lvl_min)))
  {
    std::cerr << "ERROR: Failed to parse '" << argv[2] << "' as refinement level!" << std::endl;
    return 1;
  }
  Index lvl_max(lvl_min);
  if((argc > 3) && (!String(argv[3]).parse(lvl_max)))
  {
    std::cerr << "ERROR: Failed to parse '" << argv[3] << "' as refinement level!" << std::endl;
    return 1;
  }

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

  // create input file stream
  std::ifstream ifs(filename, std::ios_base::in);
  if(!ifs.is_open() || !ifs.good())
  {
    std::cerr << "ERROR: Failed to open '" << filename << "'" << std::endl;
    return 1;
  }

  // create a reader and read the root markup
  MeshFileReader reader(ifs);
  reader.read_root_markup();

  // get mesh type
  const String mtype = reader.get_meshtype_string();

  if(mtype == "conformal:hypercube:1:1") return run_xml<H1M1D>(reader, filename, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:1:2") return run_xml<H1M2D>(reader, filename, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:1:3") return run_xml<H1M3D>(reader, filename, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:2:2") return run_xml<H2M2D>(reader, filename, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:2:3") return run_xml<H2M3D>(reader, filename, lvl_min, lvl_max);
  if(mtype == "conformal:hypercube:3:3") return run_xml<H3M3D>(reader, filename, lvl_min, lvl_max);
  if(mtype == "conformal:simplex:2:2")   return run_xml<S2M2D>(reader, filename, lvl_min, lvl_max);
  if(mtype == "conformal:simplex:2:3")   return run_xml<S2M3D>(reader, filename, lvl_min, lvl_max);
  if(mtype == "conformal:simplex:3:3")   return run_xml<S3M3D>(reader, filename, lvl_min, lvl_max);

  // unsupported mesh type
  std::cerr << "ERROR: invalid mesh type '" << mtype << "'" << std::endl;
  return 1;
}
