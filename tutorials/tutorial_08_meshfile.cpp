// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// \brief FEAT Tutorial 08: MeshFileReader
//
// This file contains a simple Poisson solver that demonstrates the usage of
// the MeshFileReader to parse a 1D/2D/3D mesh from a FEAT 3 XML mesh file.
//
// The PDE to be solved reads:
//
//    -Laplace(u) = 1          in the domain
//             u  = 0          on (parts of) the boundary
//
// For simplicity, the right-hand-side function of this Poisson PDE is set to 1
// and the boundary conditions are homogeneous Dirichlet boundary conditions,
// as this is a useful choice for any possible domain.
//
// The purpose of this tutorial is to demonstrate the usage of the MeshFileReader
// as well as to demonstrate how to write a simple application that works for
// both 2D and 3D domains at runtime.
//
// This tutorial is mostly based on Tutorial 01, although it also uses the
// SimpleArgParser class, which has been presented in Tutorial 04.
//
// ATTENTION:
// ==========
// In contrast to all previous tutorials, this tutorial application can *NOT* be
// executed without any command line arguments in a default setting, as it requires
// the user to specify an input mesh file at the very least.
//
// Note:
// You can find several commonly used mesh-files in the "data/meshes" subdirectory
// of your FEAT 3 root directory. Most notably, you can choose one of the following
// mesh files for the 2D/3D unit-square/cube domain:
// - "unit-square-quad.xml": 2D unit-square mesh (1 quadrilateral element)
// - "unit-square-tria.xml": 2D unit-square mesh (4 triangle elements)
// - "unit-cube-hexa.xml"  : 3D unit-cube mesh (1 hexahedral element)
// - "unit-cube-tetra.xml" : 3D unit-cube mesh (6 tetrahedral elements)
//
//
// Supported Command Line Parameters
// =================================
// This tutorials application supports the following command line parameters:
//
// --mesh <filenames...>
// Mandatory: Specifies the filename(s) of at least one input mesh file.
//
// --level <n>
// Optional: Specifies the mesh refinement level.
// If not given, level 0 (unrefined) is used.
//
// --dbc <names...>
// Optional: Specifies the names of the mesh-parts for Dirichlet boundary conditions.
// If not given, all mesh-parts of the mesh-file(s) are used.
//
// --vtk <name>
// Optional: Specifies the filename of the output VTK file.
// If not given, no VTK file output is written.
//
// \author Peter Zajac
//

// We start our little tutorial with a batch of includes...

// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime
#include <kernel/util/simple_arg_parser.hpp>               // for SimpleArgParser

// FEAT-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_node.hpp>                   // NEW: for RootMeshNode
#include <kernel/geometry/mesh_file_reader.hpp>            // NEW: for MeshFileReader

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>                      // for ConstantFunction

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler
#include <kernel/assembly/domain_assembler_helpers.hpp>    // for Assembly::assemble_***
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for ForceFunctional

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter

// FEAT-Solver includes
#include <kernel/solver/ilu_precond.hpp>                   // for ILUPrecond
#include <kernel/solver/pcg.hpp>                           // for PCG

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial08
{
  // Our LAFEM containers work in main memory.
  typedef Mem::Main MemType;
  // Our data arrays should be double precision.
  typedef double DataType;
  // Use the default index type for indexing.
  typedef Index IndexType;

  // Use the standard dense vector
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;
  // Use the standard CSR matrix format
  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;
  // Use the unit-filter for Dirichlet boundary conditions
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> FilterType;

  // In the previous tutorials, we have defined the shape-type of the mesh as a (compile-time
  // constant) typedef at this point prior to defining the actual application code.
  // However, in this tutorial the shape of the mesh is specified by the input mesh file(s) at
  // runtime, so we have a conflict between the code design used in FEAT, which requires the
  // shape-type at compile-time, and the fact that input mesh files are read at runtime.
  //
  // There are two possible solutions to this problem:
  //
  // 1) We could (just as in all previous tutorials) typedef the shape-type as, say,
  //    Shape::Quadrilateral (or some other type, of course) and then check whether the shape-type
  //    of the input mesh files matches this chosen definition. In case of a mismatch, we would
  //    then simply emit an error message stating the unsupported shape of the input mesh file
  //    and abort the application.
  //
  // 2) We could out-source the actual application code into a separate function template (named e.g.
  //    'run'), which is templatized in the mesh shape-type, and then use an if-else cascade to
  //    call the corresponding shape-type specialization of that function template based on the
  //    shape-type of the input mesh file(s).
  //
  // In this tutorial, we choose the second solution, which is (slightly) more complex to implement
  // than the first one. For this, we have to split the application code over two functions:
  //
  // 1) The 'main' function, which performs most of the basic initialization prior to the actual
  //    mesh object construction. In this tutorial, this boils down to setting up the SimpleArgParser
  //    and creating and initializing the MeshFileReader to obtain the shape-type of the input
  //    mesh files.
  //
  // 2) The 'run' function template, which is templatized in the shape-type that is to be used.
  //    This function then contains the 'actual' application code that starts with our well-known
  //    geometry typedefs (MeshType, TrafoType, etc.) as well as the creation of the mesh object.
  //    The remainder of this function is then more or less equivalent the previous tutorials.


  // Here comes the forward declaration of the "run" function template; it is templatized in the
  // shape-type, which is chosen from the input mesh file(s), and expects the SimpleArgParser and
  // MeshFileReader objects as arguments, which are created in the "main" function.
  // This 'run' function template is implemented right below the following "main" function.
  template<typename Shape_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader);


  // Here's our main function
  void main(int argc, char* argv[])
  {
    // First of all, let's create a SimpleArgParser
    SimpleArgParser args(argc, argv);

    // Now let's add the supported options
    // For the sake of simplicity, the actual arguments are restricted to:
    // 1) Mandatory: The input mesh option:
    args.support("mesh", "<filenames...>\n"
      "Mandatory: Specifies the filenames of the mesh files that are to be parsed.\n");
    // 2) Optional: The refinement level option:
    args.support("level", "<n>\n"
      "Optional: Specifies the refinement level. If not given, defaults to 0.\n");
    // 3) Optional: The names of the Dirichlet BC mesh-parts:
    args.support("dbc", "<meshparts...>\n"
      "Optional: Specifies the names of the Dirichlet boundary mesh-parts.\n"
      "If not given, all mesh-parts in the mesh-file(s) are used.\n");
    // 4) Optional: The VTK output filename option:
    args.support("vtk", "<filename>\n"
      "Optional: Specifies the filename of the output VTU file.\n"
      "If not given, no output VTU file is created.\n");

    // Before checking for unsupported arguments, let's check if the user supplied any arguments
    // at all, as otherwise this application would call Runtime::abort() further below:
    if(args.num_args() <= 1)
    {
      std::cout << std::endl;
      std::cout << "Info: For this tutorial, you need to specify at least one" << std::endl;
      std::cout << "input mesh file via the '--mesh <filenames...>' option." << std::endl;
      std::cout << std::endl;

      std::cout << "Supported Options:" << std::endl;
      std::cout << args.get_supported_help();
      return;
    }

    // As usual, we check for unsupported arguments and abort with an error message if necessary:
    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      // Print unsupported options
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      {
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'" << std::endl;
      }

      // And abort
      Runtime::abort();
    }

    // Now let us check whether all mandatory options have been given.
    // In our case, the only mandatory option is the '--mesh' option, which specifies the path(s)
    // to the mesh file(s) that describe our computational domain.
    if(args.check("mesh") <= 0)
    {
      // This tutorial requires at least one mesh file...
      std::cerr << std::endl << "ERROR: You have to specify at least one mesh file via '--mesh <files...>'" << std::endl;

      // ...so we cannot continue without one.
      Runtime::abort();
    }

    // At this point, you may ask why one should supply more than one mesh file here.
    // The reason for this is that the FEAT mesh file(s) do not only contain a mesh but
    // also possibly other additional important data such as
    // 1) charts for describing the analytic domain boundaries
    // 2) mesh-parts for describing the discrete domain boundary parts
    // 3) partitionings for parallel simulations

    // Now all these parts are usually combined in a single mesh file for the sake of
    // convenience, but it is also possible to split this data among several disjoint
    // files. This is especially interesting if you have several different meshes that
    // discretize the same analytic domain and you want to "outsource" the common
    // analytic domain description into a separate common file.

    // Therefore, we always have to expect that the user does not supply just one filename,
    // but a set of filenames that we have to pass to the MeshFileReader, so thus we are
    // always dealing with a deque of filename strings instead of a single string.

    // Due to the above check we know that the user supplied the '--mesh' option, so
    // we can query the deque of strings that represent the parameters of the option:
    const std::deque<String>& filenames = args.query("mesh")->second;

    // For convenience, we'll print the filenames to the console by utilizing the
    // 'stringify_join' function that will stringify each item in the deque and
    // concatenate them using a separator string:
    std::cout << "Mesh Files: " << stringify_join(filenames, " ") << std::endl;

    // Now we can create the actual reader object. Note that we need only one mesh reader
    // object even if there are multiple files to be read.
    Geometry::MeshFileReader mesh_reader;

    // Next we have to pass the filename deque to the reader's 'add_mesh_files' function.
    // However, if one of the files does not exist, the following function may throw
    // an instance of the 'FEAT::FileNotFound' exception, so we enclose this function
    // call in a try-catch block:
    try
    {
      // try to add the filenames to the mesh reader:
      mesh_reader.add_mesh_files(filenames);
    }
    catch(const std::exception& exc)
    {
      // Something went wrong; probably one of the files could not be opened...
      std::cerr << std::endl << "ERROR: " << exc.what() << std::endl;
      Runtime::abort();
    }

    // Okay, if we arrive there, then the reader has successfully opened all
    // input files. However, the reader did not actually start parsing any of
    // those files yet, so we will initiate the actual parsing process by
    // calling the 'read_root_markup' function. Again, this function may
    // throw various types of exceptions -- most notably Xml::Scanner related
    // exceptions that derive from Xml::Error -- so we also put this call
    // inside a try-catch block:
    try
    {
      // Try to read the root markup: at this point the reader analyses the first
      // line(s) of the input files to determine whether the files are actually
      // FEAT mesh files. If not, then the parser will throw an instance of the
      // 'Xml::SyntaxError' (in case a file is not an XML file) or 'Xml::GrammarError'
      // (in case the file is an XML file but not a FEAT mesh file).
      mesh_reader.read_root_markup();
    }
    catch(const std::exception& exc)
    {
      // Oops...
      std::cerr << std::endl << "ERROR: " << exc.what() << std::endl;
      Runtime::abort();
    }

    // If we arrive here, then the reader has successfully started parsing
    // the root nodes of all input files and is now able to tell us what
    // type of mesh is stored in the input file(s) (assuming that the mesh
    // files provide this information in the root markup, which is not
    // mandatory). So we can now query the mesh type from the reader:
    const String mesh_type = mesh_reader.get_meshtype_string();

    // Let us first check whether the mesh-file(s) provide us with a mesh-type.
    // If this information is missing, then the user probably did not specify all
    // required files -- or the mesh file simply did not specify the mesh-type.
    if(mesh_type.empty())
    {
      std::cout << std::endl;
      std::cout << "ERROR: Mesh file(s) did not provide a mesh-type!" << std::endl;
      std::cout << std::endl;
      std::cout << "Did you supply all required mesh-files?" << std::endl;
      Runtime::abort();
    }

    // And we'll print the mesh type string to the console:
    std::cout << "Mesh Type : " << mesh_type << std::endl;

    // Now comes the interesting part:
    // The "mesh_type" string specifies the mesh-type and therefore the shape-type
    // of the mesh that is stored in the mesh-file(s). At this point, we check all
    // five mesh-types that are supported by this tutorial application using the
    // following if-else cascade and call the corresponding "run" function template
    // specialization for the required shape-type and pass our SimpleArgParser and
    // MeshFileReader objects as parameters to it:
    if     (mesh_type == "conformal:hypercube:1:1") // 1D mesh
      run<Shape::Hypercube<1>>(args, mesh_reader);
    else if(mesh_type == "conformal:hypercube:2:2") // 2D quadrilateral mesh
      run<Shape::Hypercube<2>>(args, mesh_reader);
    else if(mesh_type == "conformal:hypercube:3:3") // 3D hexahedron mesh
      run<Shape::Hypercube<3>>(args, mesh_reader);
    else if(mesh_type == "conformal:simplex:2:2")   // 2D triangle mesh
      run<Shape::Simplex<2>>(args, mesh_reader);
    else if(mesh_type == "conformal:simplex:3:3")   // 3D tetrahedron mesh
      run<Shape::Simplex<3>>(args, mesh_reader);
    else
    {
      // The mesh-type is either invalid or not supported. In fact, there are some
      // other valid mesh-type specifiers (e.g. sub-dimensional meshes), which this
      // (and most other) application does not support by construction.
      std::cout << std::endl << "ERROR: unsupported mesh type!" << std::endl;
      Runtime::abort();
    }

    // And that's it for the main function.
  } // void main(...)


  // And this is where the magic happens: our 'run' function template.
  template<typename Shape_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
  {
    // First, define our shape-type based on the template parameter. At this point,
    // 'ShapeType' may refer to any of those types that were specified in the if-else
    // cascade in the main function above.
    typedef Shape_ ShapeType;

    // Let's print out the shape-type name just for convenience:
    std::cout << "Shape Type: " << ShapeType::name() << std::endl;

    // Now that we know the shape-type that we want to use in this specialization of the
    // 'run' function template, we can continue with the remaining typedefs as usual:

    // Use the unstructured conformal mesh class
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    // Define the corresponding mesh-part type
    typedef Geometry::MeshPart<MeshType> MeshPartType;
    // Use the standard transformation mapping
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    // Use the Lagrange-1 element
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // At this point, we have to introduce a new class that has not been used in any
    // previous tutorial so far: the (mesh) atlas, which is merely a class that manages
    // a named set of charts. A chart is a geometric object that is used for the analytic
    // description of the computational domain, such as e.g. a circle, a Bezier-Spline,
    // a sphere or a surface triangulation. At this point, we could provide a typedef
    // for the mesh atlas type, but we can also skip this, as we require the actual
    // type only for the one following variable declaration.
    // The mesh atlas class is templatized in the mesh type and we do not have to pass
    // anything to its constructor:
    Geometry::MeshAtlas<MeshType> mesh_atlas;

    // In addition, we have to define the (root) mesh-node type (which has also been used in
    // tutorial 06). This class is responsible for managing a mesh as well as a named set of
    // mesh parts defined on that particular mesh. As we are going to need this type several
    // times, we'll use a typedef for convenience here:
    typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;

    // Next, let us create a shared-pointer for the mesh-node, which will be assigned to
    // a new object by the mesh file reader in the next step.
    std::shared_ptr<RootMeshNodeType> mesh_node;

    // Okay, let the parsing begin:
    std::cout << "Parsing mesh files..." << std::endl;

    // At this point, a lot of errors could occur if either one of the files is corrupt/invalid
    // or if the set of different input files is inconsistent/mismatching. Therefore, we enclose
    // the actual call of the 'parse' function in a try-catch block:
    try
    {
      // In its most simple form (precisely: overload), the only input parameter is a reference
      // to the mesh atlas and the function returns a shared-pointer to a newly allocated
      // mesh-node object, which we assign to the previously declared variable:
      mesh_node = mesh_reader.parse(mesh_atlas);
    }
    catch(const std::exception& exc)
    {
      // That's not good...
      std::cerr << std::endl << "ERROR: " << exc.what() << std::endl;
      Runtime::abort();
    }

    // If we come out here, then the user has actually succeeded in choosing a consistent, valid
    // and supported set of input mesh files and the mesh file reader has successfully finished
    // its job. As a first action with our new mesh-node, we'll query the names of all mesh-parts
    // from the mesh-node, so that we can print them to the console.
    // The 'true' parameter for the following function call specifies that we want to query only
    // the names of the mesh-parts that were parsed from the mesh-file and that we are NOT
    // interested in any additional mesh-parts that the management code had to create for internal
    // use (mesh-parts for "internal use" are characterized by a leading underscore in their name).
    std::deque<String> meshpart_names = mesh_node->get_mesh_part_names(true);
    std::cout << "Mesh Parts: " << stringify_join(meshpart_names, " ") << std::endl;

    // In many cases, the mesh-file only contains a relatively "coarse" mesh that first has to be
    // refined a few times to obtain a mesh that is fine enough for a finite element discretization.
    // The main reason for this is that FEAT is a software package that uses geometric multigrid as
    // one of its core components, and geometric multigrid needs a mesh hierarchy.
    // However, we do not want to deal with multigrid in this tutorial, but we probably still need
    // to refine the parsed mesh, so let us check whether the user has specified a desired level
    // and, if so, try to parse it:
    Index level(0);
    if(args.parse("level", level) < 0)
    {
      // That's not meant to happen...
      std::cerr << std::endl << "ERROR: Failed to parse '--level' parameter" << std::endl;

      // and abort our program
      Runtime::abort();
    }

    // Did the user specify a level > 0?
    if(level > Index(0))
      std::cout << "Refining Mesh to Level " << level << "..." << std::endl;

    // Okay, so let's refine the mesh-node up to the desired level:
    for(Index lvl(0); lvl < level; ++lvl)
    {
      // Refining the mesh-node is really easy: we just have to call the 'refine_shared' function
      // of the mesh-node, which gives us a shared_ptr of the refined mesh node.
      mesh_node = mesh_node->refine_shared();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Okay, now that we have refined the mesh-node up to the desired level, we can now obtain a
    // reference to the actual mesh from it by calling the 'get_mesh' function and dereferencing
    // the returned pointer:
    MeshType& mesh = *mesh_node->get_mesh();

    // We'll print out the number of mesh elements just for kicks:
    std::cout << "Number of Elements: " << mesh.get_num_elements() << std::endl;

    // From this point on, the remainder of this tutorial's code is mostly identical to the
    // code from tutorial 01. The only interesting remaining difference is the assembly of the
    // Dirichlet boundary conditions, which comes after the

    // As usual, we can create the trafo and the FE space now:
    TrafoType trafo(mesh);
    SpaceType space(trafo);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Linear system assembly: Matrix + Vectors

    // The matrix/vector assembly is virtually identical to the previous tutorials.

    // Allocate matrix and assemble its structure
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    // Allocate vectors
    VectorType vec_sol = matrix.create_vector_r();
    VectorType vec_rhs = matrix.create_vector_l();

    // Create a domain assembler on all mesh elements
    Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
    domain_assembler.compile_all_elements();

    // Choose a cubature rule
    String cubature_name = "auto-degree:3";

    // First of all, format the matrix entries to zero.
    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // Assemble the Laplace operator:
    Assembly::Common::LaplaceOperator laplace_operator;
    Assembly::assemble_bilinear_operator_matrix_1(
      domain_assembler, matrix, laplace_operator, space, cubature_name);

    // Assemble the right-hand-side function:
    Analytic::Common::ConstantFunction<ShapeType::dimension> one_func(1.0);
    Assembly::assemble_force_function_vector(
      domain_assembler, vec_rhs, one_func, space, cubature_name);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    // The next step is the assembly of the Dirichlet boundary conditions.
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // In most previous tutorials, we have created a single mesh-part for the whole boundary by
    // using the BoundaryFactory. We do not want to do this here, as the mesh-file typically
    // contains the set of mesh-parts that is required for the definition of boundary conditions.
    // Note that in general, a mesh-file may also contain mesh-parts that are not meant for
    // boundary condition assembly but for some other purposes. Therefore, in a real-life
    // application, the user has to specify the list of all mesh-part names that (s)he
    // wants to define boundary conditions on. In this tutorial, the user may do this by
    // specifying the list of boundary mesh-part names as parameters for the '--dbc' option.
    // If this option was not given in the command line arguments, then we will add all
    // mesh-parts of the mesh-file to the unit filter assembler.

    // Did the user specify Dirichlet boundary condition mesh-parts via '--dbc <names..>' ?
    if(args.check("dbc") > 0)
    {
      // Yes, so choose these as our mesh-part names list:
      meshpart_names = args.query("dbc")->second;
    }

    // Loop over all Dirichlet BC mesh-part names
    for(const auto& part_name : meshpart_names)
    {
      // Try to find a mesh-part with that name in the mesh-node:
      const MeshPartType* mesh_part = mesh_node->find_mesh_part(part_name);

      // If the user specified the mesh-part names manually by using the '--dbc' option,
      // then the desired mesh-part may not actually exist (e.g. due to a typo), so we have
      // to check this case here:
      if(mesh_part == nullptr)
      {
        std::cerr << std::endl << "ERROR: no mesh-part named '" << part_name << "' could be found!" << std::endl;
        Runtime::abort();
      }

      // Okay, that mesh-part exists, so add it to the assembler:
      std::cout << "Adding mesh-part '" << part_name << "' to Dirichlet assembler..." << std::endl;
      unit_asm.add_mesh_part(*mesh_part);
    }

    // Let's create our filter and assemble it
    FilterType filter;
    unit_asm.assemble(filter, space);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    // Apply the filter onto the system matrix...
    filter.filter_mat(matrix);
    // ...the right-hand-side vector...
    filter.filter_rhs(vec_rhs);
    // ...and the solution vector.
    filter.filter_sol(vec_sol);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    std::cout << std::endl << "Solving System..." << std::endl;

    // Create an ILU(0) preconditioner
    auto precond = Solver::new_ilu_precond(matrix, filter);

    // Create a PCG solver
    auto solver = Solver::new_pcg(matrix, filter, precond);

    // Enable convergence plot
    solver->set_plot_mode(Solver::PlotMode::iter);
    solver->set_max_iter(1000);

    // Initialize the solver
    solver->init();

    // Solve our linear system
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // Release the solver
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // The filename of our VTK file
    String vtk_name;
    if(args.parse("vtk", vtk_name) > 0)
    {
      std::cout << std::endl;
      std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(mesh);

      // add the solution vector
      exporter.add_vertex_scalar("sol", vec_sol.elements());

      // finally, write the VTK file
      exporter.write(vtk_name);
    }

    // That's it for today.
    std::cout << std::endl << "Finished!" << std::endl;
  } // void run<Shape_>(...)

} // namespace Tutorial08

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialize the runtime
  Runtime::initialize(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #08: MeshFileReader" << std::endl;

  // call the tutorial's main function
  Tutorial08::main(argc, argv);

  // Finalize the runtime
  return Runtime::finalize();
}
