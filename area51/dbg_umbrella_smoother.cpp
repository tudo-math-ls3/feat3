
// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We start our little tutorial with a batch of includes...

// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/runtime.hpp>                         // for Runtime

// FEAT-Geometry includes
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>            // for RefinedUnitCubeFactory
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart
#include <kernel/geometry/umbrella_smoother.hpp>

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler
#include <kernel/geometry/mesh_file_reader.hpp>            // for MeshFileReader
#include <kernel/util/simple_arg_parser.hpp>               // for SimpleArgParser


// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

using namespace FEAT;

namespace DbgUmbrellaSmoother
{
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  typedef double DataType;

  template<typename Shape_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader_1, Geometry::MeshFileReader& mesh_reader_2);

  void main(int argc, char* argv[])
  {
    // create args parser
    SimpleArgParser args(argc, argv);

    // mandatory: mesh
    args.support("mesh", "<filename...>\n"
      "Mandatory: Specifies the filename of the mesh file that are to be parsed.\n");
    // optional: refinement level
    args.support("level", "<n>\n"
      "Optional: Specifies the refinement level. If not given, defaults to 3.\n");

    // check for unsopported options
    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      {
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'\n";
      }
      Runtime::abort();
    }

    // check if mandatory mesh was given
    if(args.check("mesh") <= 0)
    {
      std::cerr << "\nERROR: You have to specify at least one mesh file via '--mesh <files...>'\n";
      Runtime::abort();
    }

    // try to parse mesh name
    String filename(" ");
    if(args.parse("mesh", filename) < 0)
    {
      std::cerr << "\nERROR: Failed to parse '--mesh' parameter\n";
      Runtime::abort();
    }

    std::cout << "Reading Mesh ... \n";
    Geometry::MeshFileReader mesh_reader_1;
    Geometry::MeshFileReader mesh_reader_2;

    // try to read meshfile
    try
    {
      mesh_reader_1.add_mesh_file(filename);
      mesh_reader_2.add_mesh_file(filename);
    }
    catch (const std::exception& exc)
    {
      std::cerr << "\nERROR: " << exc.what() << "\n";
      Runtime::abort();
    }

    // try to read root markup
    try
    {
      mesh_reader_1.read_root_markup();
      mesh_reader_2.read_root_markup();
    }
    catch(const std::exception& exc)
    {
      std::cerr << "\nERROR: " << exc.what() << "\n";
      Runtime::abort();
    }

    const String mesh_type = mesh_reader_1.get_meshtype_string();

    // check if we have a mesh type
    if(mesh_type.empty())
    {
      std::cout << "\n";
      std::cout << "ERROR: Mesh file(s) did not provide a mesh-type!\n";
      Runtime::abort();
    }

    if     (mesh_type == "conformal:hypercube:1:1") // 1D mesh
      run<Shape::Hypercube<1>>(args, mesh_reader_1, mesh_reader_2);
    else if(mesh_type == "conformal:hypercube:2:2") // 2D quadrilateral mesh
      run<Shape::Hypercube<2>>(args, mesh_reader_1, mesh_reader_2);
    else if(mesh_type == "conformal:hypercube:3:3") // 3D hexahedron mesh
      run<Shape::Hypercube<3>>(args, mesh_reader_1, mesh_reader_2);
    else if(mesh_type == "conformal:simplex:2:2")   // 2D triangle mesh
      run<Shape::Simplex<2>>(args, mesh_reader_1, mesh_reader_2);
    else if(mesh_type == "conformal:simplex:3:3")   // 3D tetrahedron mesh
      run<Shape::Simplex<3>>(args, mesh_reader_1, mesh_reader_2);
    else
    {
      std::cout << "\nERROR: unsupported mesh type!\n";
      Runtime::abort();
    }
  }// void main(...)

  template<typename Shape_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader_1, Geometry::MeshFileReader& mesh_reader_2)
  {
    // parse level (if given)
    Index level(3);
    if(args.parse("level", level) < 0)
    {
      std::cerr << "\nERROR: Failed to parse '--level' parameter\n";
      Runtime::abort();
    }

    // create mesh atlas
    Geometry::MeshAtlas<MeshType> mesh_atlas_1;
    Geometry::MeshAtlas<MeshType> mesh_atlas_2;
    typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;

    // create mesh node
    std::unique_ptr<RootMeshNodeType> mesh_node_1;
    std::unique_ptr<RootMeshNodeType> mesh_node_2;

    // parse
    mesh_node_1 = mesh_reader_1.parse(mesh_atlas_1);
    mesh_node_2 = mesh_reader_2.parse(mesh_atlas_2);

    std::cout << "Refining Mesh to Level " << level << "...\n";
    for(Index lvl(0); lvl < level; ++lvl)
    {
      mesh_node_1 = mesh_node_1->refine_unique();
      mesh_node_2 = mesh_node_2->refine_unique();
    }
    MeshType& mesh_1 = *mesh_node_1->get_mesh();
    MeshType& mesh_2 = *mesh_node_2->get_mesh();

    std::cout << "Writing original mesh ... \n";
    Geometry::ExportVTK<MeshType> exporter_1(mesh_1);
    exporter_1.write("./mesh_orig");

    std::cout << "Smoothing original mesh uniformly ... \n";
    Geometry::UmbrellaSmoother<MeshType> mesh_smoother_1(mesh_1);
    mesh_smoother_1.smooth(100, 1e-10, true);

    std::cout << "Writing uniformly smoothed mesh ... \n";
    exporter_1.write("./uniform_smooth");

    std::cout << "Smoothing original mesh scale dependently ... \n";
    Geometry::UmbrellaSmoother<MeshType> mesh_smoother_2(mesh_2);
    mesh_smoother_2.smooth(10000, 1e-10, false);

    std::cout << "Writing scale dependently smoothed mesh ... \n";
    Geometry::ExportVTK<MeshType> exporter_2(mesh_2);
    exporter_2.write("./scale_dependent_smooth");
    std::cout << "Finished!\n\n";
  } // void run<Shape_>(...)
}// namespace DbgUmbrellaSmoother

// Here's our main function
int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  DbgUmbrellaSmoother::main(argc, argv);

  return 0;
}
