// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// We start our little tutorial with a batch of includes...

// Misc. FEAT includes
#include <kernel/runtime.hpp>                              // for Runtime
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/simple_arg_parser.hpp>               // NEW: for SimpleArgParser
#include <kernel/util/stop_watch.hpp>                      // also for Runtime


// FEAT-Geometry includes
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>            // for RefinedUnitCubeFactory
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange2/element.hpp>              // the Lagrange-2 Element (aka "Q2")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>                      // for SineBubbleFunction, ConstantFunction
#include <kernel/analytic/parsed_function.hpp>             // NEW: for ParsedScalarFunction
#include <kernel/analytic/auto_derive.hpp>                 // NEW: for AutoDerive

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>    // for Assembly::assemble_***
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for LaplaceFunctional / ForceFunctional

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter

// FEAT-Solver includes
#include <kernel/solver/ssor_precond.hpp>                  // for SSORPrecond
#include <kernel/solver/pcg.hpp>                           // for PCG

// FEAT-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_node.hpp>                   // NEW: for RootMeshNode
#include <kernel/geometry/mesh_file_reader.hpp>            // NEW: for MeshFileReader

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial04
{
  // Note:
  // This tutorial works for all shape-types. However, the console output is partially 'hard-wired'
  // for the 2D case to keep the code simple.

  // Our data arrays should be double precision.
  typedef double DataType;
  // Use the default index type for indexing.
  typedef Index IndexType;

  // Use the standard dense vector
  typedef LAFEM::DenseVector<DataType, IndexType> VectorType;
  // Use the standard CSR matrix format
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixType;
  // Use the unit-filter for Dirichlet boundary conditions
  typedef LAFEM::UnitFilter<DataType, IndexType> FilterType;

  template<typename Shape_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader);

  // Here's our main function
  void main(int argc, char* argv[])
  {
    // As a very first step, we create our argument parser that we will use to analyze
    // the command line arguments. All we need to do is to create an instance of the
    // SimpleArgParser class and pass the 'argc' and 'argv' parameters to its constructor.
    // Note that the SimpleArgParser does not modify 'argc' and 'argv' in any way.
    SimpleArgParser args(argc, argv);

    // We will keep track of whether the caller needs help, which is the case when an
    // invalid option has been supplied or '--help' was given as a command line argument.
    // In this case, we will write out some help information to std::cout later.
    bool need_help = false;

    // Now, before we start to interact with our argument parser object, there are a
    // few things to be discussed about the technical details.
    //
    // The argument parser distinguishes between three argument types:
    // Any argument beginning with a double-hyphen ('--') is called an 'option', whereas
    // any other argument is called a 'parameter', which is associated with the last
    // option preceding the parameter (if there is one). Finally, any argument preceding
    // the first option is silently ignored by the parser - this includes the very first
    // argument which is always the path of the application's binary executable.
    //
    // Confused?
    // Assume the caller entered the following at the command line:
    //
    // "./my_application foobar --quiet --level 5 2 --file ./myfile"
    //
    // This call has a total of 8 arguments:
    // 0: "./my_application" is the path of the application's binary.
    // 1: "foobar" is some argument that is ignored by our SimpleArgParser,
    //    because there is no option preceding it
    // 2: "--quiet" is an option without any parameters
    // 3: "--level" is another option with two parameters:
    // 4: "5" is the first parameter for the option "--level"
    // 5: "2" is the second parameter for the option "--level"
    // 6: "--file" is yet another option with one parameter:
    // 7: "./myfile" is the one and only parameter for "--file"
    //
    //
    // We will now tell the argument parser which options this application supports,
    // along with a short description of what the corresponding option is meant to do.
    // Although this step is not mandatory, it is highly recommended for two reasons:
    // 1. By telling the parser all our supported options, we may later on instruct the
    //    parser to check if the caller has supplied any unsupported options, so that
    //    we may print an appropriate error message. Without this, any mistyped option
    //    (e.g. '--levle' instead of '--level') will simply be ignored by the parser
    //    without emitting any warning or error message.
    // 2. Moreover, we can instruct the parser to give us a string containing all
    //    the supported options and their descriptions which we may then print out
    //    to the screen in case that the user requires help.

    // So, let's start adding our supported options by calling the parser's 'support'
    // function: The first argument is the name of the option *without* the leading
    // double-hyphen, whereas the second argument is a short description used for the
    // formatting of the argument list description.
    args.support("help", "\nDisplays this help information.\n");
    args.support("level", "<n>\nSets the mesh refinement level.\n");
    args.support("mesh", "<filenames...>\n"
      "Mandatory: Specifies the filenames of the mesh files that are to be parsed.\n");



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

    // Now that we have added all supported options to the parser, we call the 'query_unsupported'
    // function to check whether the user has supplied any unsupported options in the command line:
    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();

    // If the caller did not specify any unsupported (or mistyped) options, the container returned
    // by the query_unsupported function is empty, so we first check for that:
    if(!unsupported.empty())
    {
      // Okay, we have at least one unsupported option.
      // Each entry of the container contains a pair of an int specifying the index of the
      // faulty command line argument as well as the faulty argument itself as a string.
      // We loop over all unsupported arguments and print them out:
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      {
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'\n";
      }

      // We could abort program execution here, but instead we remember that we have to print
      // the help information containing all supported options later:
      need_help = true;
    }

    // Now let us check whether all mandatory options have been given.
    // In our case, the only mandatory option is the '--mesh' option, which specifies the path(s)
    // to the mesh file(s) that describe our computational domain.
    if(args.check("mesh") <= 0)
    {
      // This tutorial requires at least one mesh file...
      std::cerr << std::endl << "ERROR: You have to specify at least one mesh file via '--mesh <files...>'" << std::endl;

      need_help = true;
    }

    // Check whether the caller has specified the '--help' option:
    need_help = need_help || (args.check("help") >= 0);

    // Print help information?
    if(need_help)
    {
      // Okay, let's display some help.
      std::cout << "\n";
      std::cout << "USAGE: " << args.get_arg(0) << " <options>\n";
      std::cout << "\n";

      // Now we call the 'get_supported_help' function of the argument parser - this
      // will give us a formatted string containing the supported options and their
      // corresponding short descriptions  which we supplied before:
      std::cout << "Supported options:\n";
      std::cout << args.get_supported_help();

      // We abort program execution here:
      return;
    }

    // Okay, if we come out here, then all supplied options (if any) are supported and
    // the caller did not explicitly ask for help by specifying '--help', so we can now
    // start the actual parsing.

    // First of all, we will check for 'basic' options, which do not require any parameters.
    // These options are often used to enable or disable certain functionality and are just
    // interpreted as either 'being given' or 'not being given'.

    // For this, we can use the 'check' function of the parser. The only argument of this
    // function is the name of the option that we want to check for (without the leading
    // double-hyphen), and the function returns an int, which is:
    //  * = -1, if the option was not supplied by the caller
    //  * =  0, if the option was given without any parameters
    //  * = n > 0, if the option was given with n parameters


    // Next, we will parse all options which require parameters.
    // For this, we first define all variables with default values and then we will
    // continue parsing, possibly replacing the initial values by ones parsed from
    // the command line arguments specified by the caller.

    // We have two basic values, which we can parse as option parameters here:

    // Our mesh refinement level
    Index level(3);

    // Now let's start parsing the command line arguments.
    // For this, we call the 'parse' function of the argument parser and supply
    // the name of the desired option as well as the variable(s) to be parsed.
    // This function returns an int which specifies how many parameters have
    // been parsed successfully or which parameter was not parsed, i.e.
    // if the return value is
    // * = 0    , then either the option was not given at all or it was given
    //            but without any parameters.
    // * = n > 0, then the option was given and the first n parameters were
    //            parsed successfully.
    // * = n < 0, if the option was given, but the (-n)-th command line argument
    //            could not be parsed into the variable supplied to the function

    // We first try to parse the option '--level' and then check the return value:
    int iarg_level = args.parse("level", level);
    if(iarg_level < 0)
    {
      // In this case, we have an error, as the corresponding command line
      // argument could not be parsed, so print out an error message:
      std::cerr << "ERROR: Failed to parse '" << args.get_arg(-iarg_level) << "'";
      std::cerr << "as parameter for option '--level'\n";
      std::cerr << "Expected: a non-negative integer\n";

      // and abort our program
      Runtime::abort();
    }

    // Okay, that was the interesting part.
    // The remainder of this tutorial is more or less the same as in Tutorial 01 and 02.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Create Mesh, Boundary, Trafo and Space

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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
    // Use the Lagrange-2 element
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

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

    // Next, let us create a unique pointer for the mesh-node, which will be assigned to
    // a new object by the mesh file reader in the next step.
    std::unique_ptr<RootMeshNodeType> mesh_node;

    // Okay, let the parsing begin:
    std::cout << "Parsing mesh files..." << std::endl;

    // At this point, a lot of errors could occur if either one of the files is corrupt/invalid
    // or if the set of different input files is inconsistent/mismatching. Therefore, we enclose
    // the actual call of the 'parse' function in a try-catch block:
    try
    {
      // In its most simple form (precisely: overload), the only input parameter is a reference
      // to the mesh atlas and the function returns a unique pointer to a newly allocated
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
      // Refining the mesh-node is really easy: we just have to call the 'refine_unique' function
      // of the mesh-node, which gives us a unique_ptr of the refined mesh node.
      mesh_node = mesh_node->refine_unique();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Okay, now that we have refined the mesh-node up to the desired level, we can now obtain a
    // reference to the actual mesh from it by calling the 'get_mesh' function and dereferencing
    // the returned pointer:
    MeshType& mesh = *mesh_node->get_mesh();

    // We'll print out the number of mesh elements just for kicks:
    std::cout << "Number of Elements: " << mesh.get_num_elements() << std::endl;

    std::cout << "Creating Mesh on Level " << level << "...\n";


    // And create the boundary
    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);

    std::cout << "Creating Trafo and Space...\n";

    // Let's create a trafo object now.
    TrafoType trafo(mesh);

    // Create the desire finite element space.
    SpaceType space(trafo);

    //create Stopwatch
    StopWatch clock;


    // create dof-mapping
    std::cout<<"---------------------------------------------------------------\n";
    clock.start();
    Adjacency::Graph dof_graph(Space::DofMappingRenderer::render(space));
    clock.stop();
    std::cout<<"Creating took "<<clock.elapsed()<<" seconds\n";
    std::cout<<"---------------------------------------------------------------\n";
    clock.reset();

    // transposing
    clock.start();
    Adjacency::Graph dof_support(Adjacency::RenderType::transpose, dof_graph);
    clock.stop();
    std::cout<<"Transposing took "<<clock.elapsed()<<" seconds\n";
    std::cout<<"---------------------------------------------------------------\n";
    clock.reset();

    //Render composite dof-mapping/dof-support graph
    clock.start();
    Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify, dof_support, dof_graph);
    clock.stop();
    std::cout<<"Render composite dof-mapping took "<<clock.elapsed()<<" seconds\n";
    std::cout<<"---------------------------------------------------------------\n";
    clock.reset();

    //Sorting
    clock.start();
    dof_adjactor.sort_indices();
    clock.stop();
    std::cout<<"Sorting took "<<clock.elapsed()<<" seconds\n";
    std::cout<<"---------------------------------------------------------------\n";
  }// void run(...)
} // namespace Tutorial04

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialize the runtime
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #04: Parser\n";

  // call the tutorial's main function
  Tutorial04::main(argc, argv);

  return 0;
}
