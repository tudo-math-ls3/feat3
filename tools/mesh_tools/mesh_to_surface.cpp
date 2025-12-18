#include <kernel/base_header.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/cgal.hpp>
#include <kernel/util/simple_arg_parser.hpp>


#ifndef FEAT_HAVE_CGAL
static_assert(false)
#endif





template<typename ShapeType_>
void write_out(FEAT::Geometry::MeshFileReader& reader, const FEAT::String& off_file)
{
  typedef ShapeType_ ShapeType;
  typedef FEAT::Geometry::ConformalMesh<ShapeType> MeshType;
  typedef FEAT::Geometry::RootMeshNode<MeshType> RootMeshNode;
  typedef FEAT::Geometry::MeshAtlas<MeshType> MeshAtlas;
  typedef FEAT::Geometry::MeshPart<MeshType> MeshPartType;

  MeshAtlas mesh_atlas;
  std::unique_ptr<RootMeshNode> root_mesh_node = reader.parse(mesh_atlas);
  MeshType& mesh = *root_mesh_node->get_mesh();
  FEAT::Geometry::FacetFlipper<ShapeType> flipper;
  flipper.reorient(mesh.get_index_set_holder());
  FEAT::Geometry::BoundaryFactory<MeshType> bnd_factory(mesh);
  MeshPartType outer_boundary(bnd_factory);

  auto cgal_wrapper = FEAT::Geometry::cgal_wrapper_from_mesh(mesh, outer_boundary);

  cgal_wrapper.write_off(off_file);
}




int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard guard(argc, argv);

  FEAT::SimpleArgParser args(argc, argv);
  args.support("mesh", "The 3D meshfile to be converted to a off file");
  args.support("out", "The off file to be writen to");

  if(args.check("mesh") <= 0)
  {
    XABORTM("You have to provide option --mesh 'meshfile.xml'!");
  }
  FEAT::String mesh_file(args.query("mesh")->second.front());
  FEAT::String off_file(mesh_file.split_by_charset(".").front());
  if(args.check("out") > 0)
  {
    off_file = args.query("out")->second.front();
  }

  FEAT::Geometry::MeshFileReader reader;
  reader.add_mesh_file(mesh_file);
  reader.read_root_markup();

  XASSERTM(reader.get_shape_dim() == 3, "Mesh is not 3 dimesnional");


  std::cout << "Writing surface of " << mesh_file << " to " << off_file << "... ";

  const auto shapetype = reader.get_shape_type();
  switch(shapetype)
  {
    case FEAT::Geometry::MeshFileReader::ShapeType::simplex:
      write_out<FEAT::Shape::Simplex<3>>(reader, off_file);
      break;
    case FEAT::Geometry::MeshFileReader::ShapeType::hypercube:
      write_out<FEAT::Shape::Hypercube<3>>(reader, off_file);
      break;
    default:
      XABORTM(FEAT::String("Unkonw case"));
  }
  std::cout << "Done\n";

  return 0;
}