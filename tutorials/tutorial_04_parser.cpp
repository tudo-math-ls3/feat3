// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// \brief FEAT Tutorial 04: Parser demonstration
//
// This file contains a simple Poisson/Laplace solver for the unit square domain.
//
// The PDE to be solved reads:
//
//    -Laplace(u) = f          in the domain [0,1]x[0,1]
//             u  = g          on the boundary
//
// with runtime user-specified functions 'u', 'f' and 'g'.
// See the section 'The Problem Definition' at the end of this comment block for details.
//
// The purpose of this tutorial is to demonstrate the usage of two parser classes,
// which can be used to define the behavior and underlying problem of this application
// up to a certain extend. These two classes are:
//
// 1. The 'SimpleArgParser' class:
// -------------------------------
// This class is a basic light-weight command line argument parser offered by FEAT.
// Although the features offered by this class are quite limited (thus 'simple'), it
// can be used to parse simple parameter as e.g. mesh refinement levels, stopping
// criterions or VTK output filenames.
//
// 2. The 'ParsedScalarFunction' class:
// ------------------------------------
// This class implements the AnalyticFunction interface (which has been presented in
// tutorial 02) acting as a wrapper around the function parser offered by the 'fparser'
// third-party library. This class is constructed from a string (e.g. "2*x^2-y") at
// runtime and offers the possibility of evaluating the function represented by the
// string in the system assembly and/or post-processing.
//
// In combination with the SimpleArgParser, the ParsedScalarFunction class gives us a
// convenient (yet not impressively efficient) way to specify reference solutions,
// right-hand-side and/or boundary condition functions for our Poisson equation from
// the command line at runtime! Yay!
//
// Important Note #1:
// ==================
// The ParsedScalarFunction class is only defined if FEAT was configured with the
// 'fparser' build-id tag, which enables the use of the corresponding third-party library.
// If you configure without the corresponding token, then this application will compile
// without support for the ParsedScalarFunction class, therefore seriously limiting the
// functionality offered by this tutorial. In this case, the reference solution,
// right-hand-side and boundary condition functions are those coinciding to the
// sine-bubble solution used in Tutorial 01.
//
// Important Note #2:
// ==================
// It is highly recommended to specify function formulae in double-quotation marks
// when passing them as arguments to the program call from the command line, because
// the command line interpreter (bash, csh, cmd) may 'swallow' some characters (e.g.
// the circumflex ^ which acts as the power operator) or may misinterpret them as
// internal operators otherwise, thus resulting in potentially erroneous behavior.
// In short: Type
//                    ./tutorial-04-parser --u "2*x+y*2^3" --g "4-x"
// instead of
//                    ./tutorial-04-parser --u  2*x+y*2^3  --g  4-x
//
// The Problem Definition
// ======================
// As already mentioned above, this tutorial application solves a Poisson PDE with
// a caller-defined right-hand-side function 'f', a Dirichlet boundary condition
// function 'g' as well as a (optional) reference solution function 'u'.
// The actual problem depends on which of those options were actually supplied by
// the caller via the command line, so there are various cases:
//
// 1. The right-hand-side 'f' of the PDE is chosen by the following rules:
//    1.1: If the caller explicitly specified this function by supplying
//           --f <formula>
//         at the command line, then this formula is used for the definition of 'f'.
//    1.2: Else if the caller explicitly specified a reference solution function
//         by supplying
//           --u <formula>
//         then the Laplacian of this reference solution is used for 'f'.
//    1.3: Otherwise the right-hand-side function 'f' is constant zero.
//
// 2. The boundary condition function 'g' is chosen in analogy to 'f':
//    2.1: If the caller explicitly specified this function by supplying
//           --g <formula>
//         at the command line, then this formula is used for the definition of 'g'.
//    2.2: Else if the caller explicitly specified a reference solution function
//         by supplying
//           --u <formula>
//         then the boundary values of this reference solution are used for 'g'.
//    2.3: Otherwise the boundary condition function 'g' is constant zero.
//
// Finally, if the caller does not specify any of the three functions via the
// command line, this tutorial uses the sine-bubble from Tutorial 01 as a fallback,
// i.e. specifying none of the three functions is equivalent to just specifying
//     --u "sin(pi*x)*sin(pi*y)"
//
// Note:
// If no reference solution function is supplied, then this application does not
// (as it can not) provide the computation of L2/H1-errors at the end of the
// simulation.
//
// \author Peter Zajac
//

// We start our little tutorial with a batch of includes...

// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime
#include <kernel/util/simple_arg_parser.hpp>               // NEW: for SimpleArgParser

// FEAT-Geometry includes
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>            // for RefinedUnitCubeFactory
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")

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
#include <kernel/assembly/domain_assembler_helpers.hpp>    // for Assembly::assemble_***
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

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial04
{
  // Note:
  // This tutorial works for all shape-types. However, the console output is partially 'hard-wired'
  // for the 2D case to keep the code simple.

  // Once again, we use quadrilaterals.
  typedef Shape::Quadrilateral ShapeType;
  // Use the unstructured conformal mesh class
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  // Define the corresponding mesh-part type
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  // Use the standard transformation mapping
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  // Use the Lagrange-1 element
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

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
    args.support("plot", "\nDisplay the convergence plot of the linear solver.\n");
    args.support("vtk", "<filename>\nSpecifies the filename for the VTK exporter.\n"
      "If this option is not specified, no VTK file will be written\n");

    // Beside the support for the 'SimpleArgParser' class, which is a build-in feature of FEAT,
    // this tutorial also demonstrates the usage of the 'ParsedScalarFunction' later on.
    // However, this functionality is only available if the 'fparser' library was included,
    // so we need to use an #ifdef here to include the corresponding options.
#ifdef FEAT_HAVE_FPARSER
    // Add our options related to the functionality offered by the fparser library.
    args.support("u", "<formula>\nSpecifies the reference solution u.\n");
    args.support("f", "<formula>\nSpecifies the right-hand-side force function f.\n");
    args.support("g", "<formula>\nSpecifies the Dirichlet boundary condition function g.\n");
#else // no fparser support
    // Specifying functions at runtime is not supported, so let' display an annoying note instead
    std::cout << std::endl;
    std::cout << "Important Note:" << std::endl;
    std::cout << "This application binary has been configured and build without support for" << std::endl;
    std::cout << "the 'fparser' third-party library, which is required for the specification" << std::endl;
    std::cout << "of custom solution, right-hand-side and boundary condition functions." << std::endl;
    std::cout << "To enable this functionality, please re-configure your FEAT build" << std::endl;
    std::cout << "by specifying 'fparser' as an additional part of your build-id." << std::endl;
    std::cout << std::endl;
#endif // FEAT_HAVE_FPARSER

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
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'" << std::endl;
      }

      // We could abort program execution here, but instead we remember that we have to print
      // the help information containing all supported options later:
      need_help = true;
    }

    // Check whether the caller has specified the '--help' option:
    need_help = need_help || (args.check("help") >= 0);

    // Print help information?
    if(need_help)
    {
      // Okay, let's display some help.
      std::cout << std::endl;
      std::cout << "USAGE: " << args.get_arg(0) << " <options>" << std::endl;
      std::cout << std::endl;

      // Now we call the 'get_supported_help' function of the argument parser - this
      // will give us a formatted string containing the supported options and their
      // corresponding short descriptions  which we supplied before:
      std::cout << "Supported options:" << std::endl;
      std::cout << args.get_supported_help();

      // In case that we built with support for the fparser library, we also print
      // some additional information regarding the specification of the functions.
#ifdef FEAT_HAVE_FPARSER
      std::cout << "Remarks regarding function formulae:" << std::endl;
      std::cout << "The <formula> parameters of the options '--u', '--f' and '--g'" << std::endl;
      std::cout << "are expected to be function formulae in the variables ";
      switch(ShapeType::dimension)
      {
      case 1:
        std::cout << "'x'" << std::endl;
        break;
      case 2:
        std::cout << "'x' and 'y'" << std::endl;
        break;
      case 3:
        std::cout << "'x', 'y' and 'z'" << std::endl;
        break;
      }
      std::cout << "As an example, one may specify" << std::endl;
      std::cout << "             --u \"sin(pi*x)*sin(pi*y)\"" << std::endl;
      std::cout << "to define the reference solution to be the sine-bubble." << std::endl;
      std::cout << "For a full list of supported expressions and build-in functions, refer to" << std::endl;
      std::cout << "http://warp.povusers.org/FunctionParser/fparser.html#functionsyntax" << std::endl;
      std::cout << std::endl;
      std::cout << "Important Note:" << std::endl;
      std::cout << "Although it may not be required in all cases, it is highly recommended" << std::endl;
      std::cout << "that you enclose the function formulae in quotation marks as shown in" << std::endl;
      std::cout << "example above. If not quoted, your command line interpreter (e.g. bash," << std::endl;
      std::cout << "csh, cmd) may misinterpret special characters, thus possibly leading to" << std::endl;
      std::cout << "incorrect program behavior." << std::endl;

#endif // FEAT_HAVE_FPARSER

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

    // Currently, the only option which falls into this category is the '--plot' option,
    // which enables the printing of the solver convergence plot. For this, we check
    // whether the result of the 'check' function is non-negative, in which case the
    // option was supplied, possibly with some parameter that we will ignore in this case:
    bool solver_plot = (args.check("plot") >= 0);

    // Next, we will parse all options which require parameters.
    // For this, we first define all variables with default values and then we will
    // continue parsing, possibly replacing the initial values by ones parsed from
    // the command line arguments specified by the caller.

    // We have two basic values, which we can parse as option parameters here:

    // Our mesh refinement level
    Index level(3);
    // The filename of our VTK file
    String vtk_name("");

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
      std::cerr << "as parameter for option '--level'" << std::endl;
      std::cerr << "Expected: a non-negative integer" << std::endl;

      // and abort our program
      Runtime::abort();
    }

    // Next, we check for the '--vtk' option, which specifies the filename of the
    // VTK file to be written. If this option is not given or if it has no parameters,
    // then we will not write a VTK file at all, so we just check whether the return
    // value of the parse function is positive:
    bool want_vtk = (args.parse("vtk", vtk_name) > 0);

    // Now we will parse the formulae for our PDE functions if the fparser library in enabled.

#ifdef FEAT_HAVE_FPARSER
    // Let's initialize our reference solution, the rhs and the dbc function formulae
    // to empty strings
    String formula_u("");
    String formula_f("");
    String formula_g("");

    // And let's parse the corresponding command line arguments:
    bool have_u = (args.parse("u", formula_u) == 1);
    bool have_f = (args.parse("f", formula_f) == 1);
    bool have_g = (args.parse("g", formula_g) == 1);

    // If the caller did not specify any of the above, we use the sine-bubble:
    if(!(have_u || have_f || have_g))
    {
      formula_u = "sin(pi*x)*sin(pi*y)";
      have_u = true;
    }

    // At this point, we have at least one function formula. Now, we need to create
    // three instances of the ParsedScalarFunction class template, which we will pass on
    // to our assembly functions later on. The only template parameter is the dimension
    // of the function to be parsed:
    Analytic::ParsedScalarFunction<ShapeType::dimension> rhs_function; // right-hand-side
    Analytic::ParsedScalarFunction<ShapeType::dimension> dbc_function; // boundary conditions

    // In the case of the reference solution function, we also require the computation of
    // derivatives for the assembly of the right-hand-side (if 'f' is not given explicitly)
    // and for the computation of errors in the post-processing step.
    // Unfortunately, the ParsedScalarFunction cannot compute the derivatives by itself, so
    // we need to put it into an 'AutoDerive' function wrapper - this one will add
    // the numeric computation of derivatives to our ParsedScalarFunction automagically.
    Analytic::AutoDerive<Analytic::ParsedScalarFunction<ShapeType::dimension>> sol_function;

    // We have three ParsedScalarFunction object, but we still need to supply them with
    // our (or the caller's) function formulae:
    if(have_u)
    {
      // We have a formula for our reference solution function 'u', so we simply
      // call the 'parse' function of the ParsedScalarFunction object and supply the
      // formula to it. Note that the ParsedScalarFunction::parse function may throw
      // an exception if the user entered garbage at the command line, so we
      // need to put the parse call into a try-catch block here:
      try
      {
        // Let's give it a try
        sol_function.parse(formula_u);
      }
      catch(const std::exception& exc)
      {
        // Oops...
        std::cerr << "ERROR: Cannot parse expression '" << formula_u << "' as function 'u(x,y)'" << std::endl;
        // Let's also print the exception's message, which might contain useful information
        // regarding the cause of the error
        std::cerr << exc.what() << std::endl;
        Runtime::abort();
      }
    }
    // Okay, now let's repeat the same for 'f' and 'g'
    if(have_f)
    {
      try
      {
        rhs_function.parse(formula_f);
      }
      catch(const std::exception& exc)
      {
        std::cerr << "ERROR: Cannot parse expression '" << formula_f << "' as function 'f(x,y)'" << std::endl;
        std::cerr << exc.what() << std::endl;
        Runtime::abort();
      }
    }
    if(have_g)
    {
      try
      {
        dbc_function.parse(formula_g);
      }
      catch(const std::exception& exc)
      {
        std::cerr << "ERROR: Cannot parse expression '" << formula_g << "' as function 'g(x,y)'" << std::endl;
        std::cerr << exc.what() << std::endl;
        Runtime::abort();
      }
    }

    // Okay, if the arrived here, then all formulae specified by the caller (if any) were
    // parsed successfully, so we'll write our functions to cout for convenience.
    // Note:
    // At this point, we always output 'x,y' as the arguments of u, f and g, even if this
    // tutorial was modified to use a 1D or 3D shape. This is only to keep the following
    // code simple and to avoid a nested if-orgy here.
    std::cout << std::endl;
    std::cout << "Function summary:" << std::endl;
    if(have_u)
      std::cout << "u(x,y) = " << formula_u << std::endl;
    else
      std::cout << "u(x,y) = - unknown -" << std::endl;
    if(have_f)
      std::cout << "f(x,y) = " << formula_f << std::endl;
    else if(have_u)
      std::cout << "f(x,y) = -Laplace(u)" << std::endl;
    else
      std::cout << "f(x,y) = 0" << std::endl;
    if(have_g)
      std::cout << "g(x,y) = " << formula_g << std::endl;
    else if(have_u)
      std::cout << "g(x,y) = u(x,y)" << std::endl;
    else
      std::cout << "g(x,y) = 0" << std::endl;
    std::cout << std::endl;

#else
    // Without the fparser library, we use the sine-bubble as a solution.
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;
    // The following two functions will not be used and are therefore just 'dummies',
    // which are only declared here to avoid even more #ifdef's in the following code.
    Analytic::Common::ConstantFunction<ShapeType::dimension> dbc_function(0.0);
    Analytic::Common::ConstantFunction<ShapeType::dimension> rhs_function(0.0);
    const bool have_u = true;
    const bool have_f = false;
    const bool have_g = false;
#endif

    // Okay, that was the interesting part.
    // The remainder of this tutorial is more or less the same as in Tutorial 01 and 02.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Create Mesh, Boundary, Trafo and Space

    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;

    // Create the mesh
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    // And create the boundary
    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);

    std::cout << "Creating Trafo and Space..." << std::endl;

    // Let's create a trafo object now.
    TrafoType trafo(mesh);

    // Create the desire finite element space.
    SpaceType space(trafo);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly

    std::cout << "Allocating matrix and vectors..." << std::endl;

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
    String cubature_name = "auto-degree:5";

    // First of all, format the matrix entries to zero.
    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // Create the pre-defined Laplace operator:
    Assembly::Common::LaplaceOperator laplace_operator;

    // And assemble that operator
    Assembly::assemble_bilinear_operator_matrix_1(
      domain_assembler, matrix, laplace_operator, space, cubature_name);

    // Now, let us assemble the right-hand-side function.
    if(have_f)
    {
      // We have an rhs function given explicitly, so we use this one.
      // We use the 'assemble_force_function_vector', which will wrap the given analytical function
      // into an Assembly::Common::ForceFunctional for us automatically.
      Assembly::assemble_force_function_vector(domain_assembler, vec_rhs, rhs_function, space, cubature_name);
    }
    else if(have_u)
    {
      // We do not have an explicit rhs function, but we have a solution function,
      // so let's compute the RHS functional from the solution.
      Assembly::Common::LaplaceFunctional<decltype(sol_function)> functional(sol_function);
      Assembly::assemble_linear_functional_vector(domain_assembler, vec_rhs, functional, space, cubature_name);
    }
    // else: We have neither u nor f given, so we leave the RHS vector zero.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    // The next step is the assembly of the inhomogeneous Dirichlet boundary conditions.
    // For this task, we require a Unit-Filter assembler:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // Add the one and only boundary component to the assembler:
    unit_asm.add_mesh_part(boundary);

    // Let's create our filter
    FilterType filter;

    // Do we have a boundary condition function given?
    if(have_g)
    {
      // Yes, so use it for the assembly of the filter.
      unit_asm.assemble(filter, space, dbc_function);
    }
    else if(have_u)
    {
      // We don't have an explicit boundary condition function, but we have a solution function,
      // we will use this one for our boundary conditions.
      unit_asm.assemble(filter, space, sol_function);
    }
    else
    {
      // We have neither a boundary condition function nor a solution function, so we
      // initialize homogeneous Dirichlet boundary conditions.
      unit_asm.assemble(filter, space);
    }

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

    std::cout << "Solving System..." << std::endl;

    // Create a SSOR preconditioner
    auto precond = Solver::new_ssor_precond(matrix, filter);

    // Create a PCG solver
    auto solver = Solver::new_pcg(matrix, filter, precond);

    // Enable convergence plot
    if(solver_plot)
    {
      solver->set_plot_mode(Solver::PlotMode::iter);
    }

    // Initialize the solver
    solver->init();

    // Solve our linear system
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // Release the solver
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    if(have_u)
    {
      std::cout << std::endl;
      std::cout << "Computing errors against reference solution..." << std::endl;

      // Compute the error norms:
      auto error_info = Assembly::integrate_error_function<1>(
        domain_assembler, sol_function, vec_sol, space, cubature_name);

      // Print the error norms to the console
      std::cout << "Error Analysis:" << std::endl;
      std::cout << error_info.print_norms() << std::endl;
    }

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    if(want_vtk)
    {
      std::cout << std::endl;
      std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

      // project solution and right-hand-side vectors
      VectorType vertex_sol, vertex_rhs;
      Assembly::DiscreteVertexProjector::project(vertex_sol, vec_sol, space);
      Assembly::DiscreteVertexProjector::project(vertex_rhs, vec_rhs, space);

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(mesh);

      // add the vertex-projection of our solution and rhs vectors
      exporter.add_vertex_scalar("sol", vertex_sol.elements());
      exporter.add_vertex_scalar("rhs", vertex_rhs.elements());

      // finally, write the VTK file
      exporter.write(vtk_name);
    }

    // That's it for today.
    std::cout << std::endl << "Finished!" << std::endl;
  } // void main(...)
} // namespace Tutorial04

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialize the runtime
  Runtime::initialize(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #04: Parser" << std::endl;

  // call the tutorial's main function
  Tutorial04::main(argc, argv);

  // Finalize the runtime
  return Runtime::finalize();
}
