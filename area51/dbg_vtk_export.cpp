// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "kernel/shape.hpp"
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_part.hpp>

namespace DbgVTKExport
{
  using namespace FEAT;

  using MeshType3D = Geometry::ConformalMesh<Shape::Hexahedron>;
  using MeshPartType3D = Geometry::MeshPart<MeshType3D>;

  using MeshType2D = Geometry::ConformalMesh<Shape::Quadrilateral>;
  using MeshPartType2D = Geometry::MeshPart<MeshType2D>;

  int main()
  {
    {
      Geometry::RefinedUnitCubeFactory<MeshType3D> factory(3);
      MeshType3D mesh(factory);

      // Export mesh faces
      Geometry::ExportVTK<MeshType3D, 2> face_exporter(mesh);
      face_exporter.write("inner_faces");

      // Export boundary
      MeshPartType3D boundary = Geometry::make_boundary_meshpart(mesh);
      Geometry::ExportVTK<MeshType3D, 2> boundary_exporter(mesh, boundary);

      std::vector<double> data(boundary.get_num_entities(2));

      for(std::size_t i(0); i < boundary.get_num_entities(2); i++)
      {
        data[i] = double(3 * i);
      }
      boundary_exporter.add_cell_scalar("Data", data.data());
      boundary_exporter.write("boundary");
    }

    {
      Geometry::RefinedUnitCubeFactory<MeshType2D> factory(3);
      MeshType2D mesh(factory);

      // Export edges
      Geometry::ExportVTK<MeshType2D, 1> edge_exporter(mesh);
      edge_exporter.write("edges");

      // Export boundary
      MeshPartType2D boundary = Geometry::make_boundary_meshpart(mesh);
      Geometry::ExportVTK<MeshType2D, 1> boundary_exporter(mesh, boundary);

      std::vector<double> data(boundary.get_num_entities(1));

      for(std::size_t i(0); i < boundary.get_num_entities(1); i++)
      {
        data[i] = double(3 * i);
      }
      boundary_exporter.add_cell_scalar("Data", data.data());
      boundary_exporter.write("boundary");
    }


    return 0;
  }
}

int main(int /*argc*/, char* /*argv*/[])
{
  DbgVTKExport::main();
}
