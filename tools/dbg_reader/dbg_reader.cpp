#include <kernel/geometry/mesh_reader_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>

using namespace FEAST;

typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;
typedef Geometry::ConformalSubMesh<Shape::Quadrilateral> QuadSubMesh;
typedef Geometry::CellSubSet<Shape::Quadrilateral> QuadCellSet;
typedef Geometry::StandardRefinery<QuadMesh> Refinery;
typedef Geometry::MeshReaderFactory<QuadMesh> MeshReaderFactory;
typedef Geometry::MeshReaderFactory<QuadSubMesh> SubMeshReaderFactory;
typedef Geometry::MeshReaderFactory<QuadCellSet> CellSetReaderFactory;
typedef Geometry::ExportVTK<QuadMesh> Exporter;

template<typename T_, typename CS_>
void filter_vertex_vector(T_ x[], CS_& cell_set, T_ value = T_(0))
{
  Geometry::TargetSet& trg(cell_set.template get_target_set<0>());
  for(Index i(0); i < trg.get_num_entities(); ++i)
    x[trg[i]] = value;
}

int main(int /*argc*/, char** /*argv*/)
{
  // create reader and parse
  MeshReader reader;
  reader.parse_mesh_file("../../data/meshes/"
    //"unit-square-mesh.txt");
    "bench1-mesh.txt");

  // create reader factory
  MeshReaderFactory mesh_factory(reader);

  // create mesh
  QuadMesh mesh(mesh_factory);

  // create cell-set reader factory
  CellSetReaderFactory cell_set_factory(reader, "1"); // 1 -> outer rectangular boundary component

  // create cell-set
  QuadCellSet cell_set(cell_set_factory);

  // create sub-mesh reader factory
  SubMeshReaderFactory sub_mesh_factory(reader, "2"); // 2 -> inner circular boundary component

  // create sub-mesh
  QuadSubMesh sub_mesh(sub_mesh_factory);

  Index num_verts = mesh.get_num_entities(0);
  double* x = new double[num_verts];
  for(Index i(0); i < num_verts; ++i)
    x[i] = 0.0;

  filter_vertex_vector(x, cell_set, +1.0);
  filter_vertex_vector(x, sub_mesh, -1.0);

  // create exporter
  Exporter vtk(mesh);

  vtk.add_scalar_vertex("bla", x);

  delete x;

  // export VTK
  vtk.write("./bench1.vtk");
}
