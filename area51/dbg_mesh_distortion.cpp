
// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We start our little tutorial with a batch of includes...

// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime

// FEAT-Geometry includes
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>            // for RefinedUnitCubeFactory
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart
#include <kernel/geometry/mesh_distortion.hpp>

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT, so use the namespace here.
using namespace FEAT;

namespace DbgMeshDistortion
{
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  typedef double DataType;

  void main(Index level)
  {
    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);

    MeshType mesh_1(mesh_factory);
    MeshType mesh_2(mesh_factory);
    MeshType mesh_3(mesh_factory);
    //mesh_1.create_permutation(Geometry::PermutationStrategy::lexicographic);

    Geometry::ExportVTK<MeshType> exporter_1(mesh_1);
    Geometry::ExportVTK<MeshType> exporter_2(mesh_2);
    Geometry::ExportVTK<MeshType> exporter_3(mesh_3);
    exporter_1.write("./mesh_orig");

    DataType mesh_size = 0.125;
    DataType dsh = mesh_size/4.0;

    Geometry::MeshDistortion<MeshType> mesh_distortion_1(mesh_1);
    mesh_distortion_1.distort_uniform(dsh);
    exporter_1.write("./uniform_dis");

    Geometry::MeshDistortion<MeshType> mesh_distortion_2(mesh_2);
    mesh_distortion_2.distort_uniform(dsh);
    mesh_distortion_2.distort_shortest_edge_uniform();
    exporter_2.write("./shortest_edge_uniform_dis");

    Geometry::MeshDistortion<MeshType> mesh_distortion_3(mesh_3);
    mesh_distortion_3.distort_uniform(dsh);
    mesh_distortion_3.distort_shortest_edge_local();
    exporter_3.write("./shortest_edge_local_dis");

    std::cout << "Finished!\n" << std::endl;

  } // void main(...)
} // namespace DbgMeshDistortion

// Here's our main function
int main(int argc, char* argv[])
{
  Runtime::initialize(argc, argv);
  Index level(3);

  if(argc > 1)
  {
    int ilevel(0);
    if(!String(argv[argc-1]).parse(ilevel) || (ilevel < 1))
    {
      // Failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[argc-1] << "' as refinement level." << std::endl;
      std::cerr << "Note: The last argument must be a positive integer." << std::endl;
      // Abort our runtime environment
      Runtime::abort();
    }
    // If parsing was successful, use the given information
    level = Index(ilevel);
  }

  DbgMeshDistortion::main(level);

  return Runtime::finalize();
}
