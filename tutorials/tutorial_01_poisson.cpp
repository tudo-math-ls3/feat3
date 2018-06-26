//
// \brief FEAT Tutorial 01: Poisson solver (TM)
//
// This file contains a simple prototypical Poisson solver for the unit square domain.
//
// The PDE to be solved reads:
//
//    -Laplace(u) = f          in the domain [0,1]x[0,1]
//             u  = 0          on the boundary
//
// The analytical solution u is given as
//
//         u(x,y) = sin(pi*x) * sin(pi*y)
//
// and its corresponding force (right-hand-side) f is
//
//         f(x,y) = 2 * pi^2 * u(x,y)
//
//
// The purpose of this tutorial is to demonstrate the basic program flow of a simple
// scalar stationary linear PDE solver without going to far into the details of what
// magic happens under-the-hood.
//
//
// The basic program flow of this application is as follows:
//
// 1. Define the required types for spatial discretisation and linear algebra.
//
// 2. Create a mesh and a boundary mesh-part by using a factory.
//
// 3. Create a trafo based on the mesh and a finite element space based on the trafo.
//
// 4. Assemble the PDE operator and the force to a matrix and a right-hand-side vector.
//
// 5. Assemble the boundary conditions to a filter.
//
// 6. Solve the resulting linear system using a simple fire-and-forget solver.
//
// 7. Compute the L2- and H1-errors against the analytical reference solution.
//
// 8. Write the result to a VTK file for visual post-processing, if desired.
//
// \author Peter Zajac
//

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

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")
#include <kernel/space/lagrange2/element.hpp>              // the Lagrange-2 Element (aka "Q2")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>                      // for SineBubbleFunction

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for LaplaceFunctional
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/error_computer.hpp>              // for ScalarErrorComputer
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter

// FEAT-Solver includes
#include <kernel/solver/ssor_precond.hpp>                  // for SSORPrecond
#include <kernel/solver/pcg.hpp>                           // for PCG

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT, so use the namespace here.
using namespace FEAT;

// We're opening a new namespace for our tutorial.
// The only reason for this is that some compilers may give us warnings about "shadowing" types
// otherwise -- this would not be dramatic, but somewhat annoying...
namespace Tutorial01
{
  // We start with a set of typedefs, which make up the basic configuration for this
  // tutorial application. The general idea of these typedefs is (1) to avoid typing and
  // (2) specialise FEAT to do, out of the many possibilities, only what we want to do in
  // this tutorial.

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // First, we need to specify what type of mesh and what type of finite element space we want
  // to use in this tutorial application. For this, FEAT uses a "template type nesting" approach,
  // in which have to define the corresponding classes by nesting them.

  // The first type that we need to choose is the "shape type" of the mesh, i.e. the type of
  // cells/elements that our mesh should contain. At this point, we implicitly also choose whether
  // this will be a 1D, 2D or 3D tutorial, depending on the dimension of the chosen shape-type.
  // There are 5 shape-types available in FEAT and for this tutorial we pick quadrilateral elements.

  //typedef Shape::Hypercube<1> ShapeType;  // 1D
  //typedef Shape::Triangle ShapeType;      // 2D, same as Shape::Simplex<2>
  typedef Shape::Quadrilateral ShapeType;   // 2D, same as Shape::Hypercube<2>
  //typedef Shape::Tetrahedron ShapeType;   // 3D, same as Shape::Simplex<3>
  //typedef Shape::Hexahedron ShapeType;    // 3D, same as Shape::Hypercube<3>

  // The next step in the "template type nesting" approach is the choice of a mesh class.
  // We want to employ a simple unstructured conformal mesh, which is implemented by the
  // "Geometry::ConformalMesh" class template. At this point, we have to pass the chosen
  // shape-type as the first template parameter to the class template.
  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  // The next type that we will require later is the "mesh-part" type. Object of this class
  // are used to describe certain mesh regions of some interest, e.g. boundary regions of
  // the mesh that we want to use for the definition of some boundary conditions.
  // The corresponding class is the "Geometry::MeshPart" class template and its one and only
  // template parameter is the mesh type that we have just defined:
  typedef Geometry::MeshPart<MeshType> MeshPartType;

  // The next thing we need is a transformation that our finite element spaces should use.
  // The transformation is responsible for providing the "reference-to-real-cell" mapping
  // functionality that is used by both the finite element space as well as various assembly
  // algorithms. Currently, the only available transformation is the "standard" mapping,
  // which represents a first-order transformation. The corresponding class is the
  // "Trafo::Standard::Mapping" class template and its one and only template parameter is
  // the underlying mesh type that we have defined above:
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  // Finally, as the last step of our "template type nesting" approach is the definition of
  // the actual finite element space(s). In this tutorial we want to stick with the simple
  // standard first-order H1-conforming Q1 (or P1) space, which is implemented by the
  // "Space::Lagrange1::Element" class template. This element family takes the transformation
  // type as the one and only template parameter:

  // Use the Lagrange-1 element (aka "Q1" or "P1"):
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  // Or you could also use the Lagrange-2 element (aka "Q2" or "P2") instead:
  //typedef Space::Lagrange2::Element<TrafoType> SpaceType;


  // Before we continue, let us recapitulate what we have defined so far:
  // Reading the previous "template type nesting" odyssey in a bottom-up manner show us that
  // "SpaceType" is a
  // - Lagrange-1 (Q1/P1) finite element space defined on a
  // - standard transformation mapping defined on a
  // - conformal (unstructured) mesh consisting of
  // - quadrilateral cells/elements

  // At this point, you maybe already have an idea why we use those typedefs:

  // We could simply change the "ShapeType" definition from "Shape::Quadrilateral" to, say,
  // "Shape::Tetrahedron" and "SpaceType" would automatically be switched from a 2D Q1 element
  // to the matching 3D P1 element! This type of modularity is a major ingredient when it
  // comes to writing truly multi-dimensional applications -- or in other (more fancy) words:
  // Write a 2D application and get the 3D version for free! Yay!

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // In the next part, we have to define the linear algebra container types (matrix, vector, filter)
  // that our tutorial will use. At first, it may seem confusing that we have to define a type for
  // vectors via some template typedefs instead of using "the almighty one-size-fits-all vector"
  // class. Again, the reason is modularity: FEAT has much more to offer than just the standard
  // types, although that is all that we want to consider in this tutorial.

  // First, we have to choose the "Memory-Data-Index" type triplet, which specifies the basic types
  // that our linear algebra containers will use for their elements.

  // The first one is the memory type: This "tag class" specifies in which type of memory our
  // matrices and vectors will operate. In this tutorial, we want to stick with the main memory,
  // which is simply the RAM that the CPU has access to. FEAT also supports linear algebra
  // containers which work on GPUs using CUDA, but this will be covered another time.
  typedef Mem::Main MemType;

  // The second type is the data type: This is simply the type of the matrix and vector elements,
  // which is typically the double-precision floating point type aka "double". FEAT also supports
  // other data types as single or even quadruple precision, but we want to avoid that for now.
  typedef double DataType;

  // The third type is the index type: This is used by various containers which also store arrays
  // of indices, such as the row-pointer and column-index arrays of the CSR matrix format.
  // In particular, any (sufficiently large) integer type will do, so we stick to the "Index" type,
  // which corresponds to "unsigned long" by default.
  typedef Index IndexType;


  // Based on the three memory, data and index typedefs, we can now define the vector type.
  // In this tutorial, we want to solve a simple scalar PDE, so we require just a "standard"
  // vector class, which is implemented by the "LAFEM::DenseVector" class template:
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;

  // Furthermore, for the discretised Poisson operator, we require a scalar sparse matrix type.
  // We choose the famous CSR format here, because it is pretty much standard for unstructured FEM:
  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;

  // Finally, we need a filter. Filters are responsible for "enforcing" simple linear constraints
  // such as boundary conditions and are required by the linear solver framework. In this tutorial,
  // we have a scalar PDE with Dirichlet boundary conditions and for this type of problem, we
  // require a so-called "unit-filter":
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> FilterType;

  // That's it for the linear algebra types.

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function.
  void main(Index level)
  {
    // Okay, we have already defined the mesh, trafo and space types, so we can start our
    // actual tutorial code by creating those.

    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;

    // In a "real-life" application, we would either read the mesh from an input file or call
    // some sort of mesh generator library/tool to create a mesh for us, but in this tutorial,
    // we want to stick with a more simple solution, especially to avoid dependencies on external
    // files, libraries or tools. Instead, we will use the "RefinedUnitCubeFactory", which will
    // (as the name suggests) generate a refined unit-square mesh for us.

    // First of all, create a mesh factory object representing a refined unit-square domain
    // and pass the desired refinement level to its constructor:
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);

    // Now create the actual mesh by using that factory:
    MeshType mesh(mesh_factory);

    // Furthermore, we require a mesh-part that represents the boundary of the mesh for the
    // assembly of Dirichlet boundary conditions later on. Again, this would typically come
    // from an external file or a mesh generator, but we stick to a more simple solution.
    // We will utilise the "BoundaryFactory", which will create a boundary mesh-part for a
    // given mesh. Note that this BoundaryFactory class works for any given input mesh and
    // not only for meshes created by the RefinedUnitCubeFactory.

    // Now let's create a boundary factory for our mesh.
    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);

    // And create the boundary mesh-part by using the factory.
    MeshPartType boundary(boundary_factory);

    // Voila, we now have a mesh and a corresponding boundary mesh-part.


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Trafo and Finite Element Space initialisation

    std::cout << "Creating Trafo and Space..." << std::endl;

    // We have already defined the types of the transformation and finite element space,
    // so we just need to create the objects.

    // Let's create a trafo object now. Its only parameter is the mesh that it is defined on.
    TrafoType trafo(mesh);

    // Create the desire finite element space. Its only parameter is the trafo that it is defined on.
    SpaceType space(trafo);


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly

    std::cout << "Allocating matrix and vectors..." << std::endl;

    // Now we need to perform the symbolic matrix assembly, i.e., the computation of the non-zero
    // sparsity pattern and the allocation of the internal matrix arrays. For this, we first create
    // an empty matrix and then call the SymbolicAssembler to create the matrix sparsity pattern:
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    // Note: The "1" at the end of the function name "assemble_matrix_std1" indicates that there
    // is only one finite element space involved in the assembly, i.e. the test- and trial-
    // spaces are identical. The symbolic assembler also provides a function named
    // "assemble_matrix_std2", which takes two (possibly different) finite element spaces as
    // test- and trial-spaces, but we're not going to use it for now. However, this is
    // required for more complex PDEs like, e.g., the Stokes equations.

    // Now that the matrix structure is assembled, we can easily create two vector objects
    // for the solution and the right-hand-side vector by the create_vector_l/r member function
    // of the matrix object. This member function creates a vector of the corresponding
    // type and dimension, so that the vector can acts as a left or right multiplicand
    // in a matrix-vector multiply operation. Our matrix is square, so it does not matter
    // whether we call 'create_vector_l' or 'create_vector_r' here:
    VectorType vec_sol = matrix.create_vector_r();
    VectorType vec_rhs = matrix.create_vector_l();

    // Okay, we now have a matrix and two vectors. The internal arrays of those object are
    // allocated, but their data arrays are still uninitialised (automatic initialisation is
    // not performed due to performance reasons). Before we can continue with the numerical
    // assembly, we first need to reset the matrix and data arrays, i.e. set all entries
    // to zero, which is done by calling the "format" function. This is required because the
    // assembly methods for the matrix and the right-hand-side work in an "additive" way, i.e.
    // the assembled operators/functionals are added onto the given matrix/vector, so we
    // have to start with a null-matrix/-vector.

    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)

    // Before we start assembling anything, we need create a cubature factory, which will be
    // used to generate cubature rules for integration.
    // There are several cubature rules available in FEAT, and a complete list can be queried
    // by compiling and executing the 'cub-list' tool in the 'tools' directory.
    // We choose the 'auto-degree:n' rule here, which will automatically choose an appropriate
    // cubature rule that can integrate polynomials up to degree n exactly. Some other possible
    // choices are also provided here, but commented out.
    Cubature::DynamicFactory cubature_factory(
      "auto-degree:5"          // automatic cubature rule for 5th degree polynomials
    //"gauss-legendre:3"       // 3x3 Gauss-Legendre rule
    //"gauss-lobatto:4"        // 4x4 Gauss-Lobatto rule
    //"newton-cotes-open:5"    // 5x5 'open' Newton-Cotes rule
    //"trapezoidal"            // trapezoidal rule (works for all shape types)
    //"barycentre"             // barycentre rule (not recommended due to insufficient order)
    );

    std::cout << "Assembling system matrix..." << std::endl;

    // We want to assemble the 2D "-Laplace" operator (-u_xx -u_yy).
    // In this tutorial, we use a pre-defined class for the implementation of this operator.
    // The definition of other and more complex PDE operators is discussed in another tutorial.
    Assembly::Common::LaplaceOperator laplace_operator;

    // Next, we call the bilinear operator assembler to assemble the operator into a matrix.
    // In analogy to the SymbolicAssembler class, the corresponding assemble function of the
    // BilinearOperatorAssembler class also has a "1" suffix indicating that there is only one
    // finite element space involved:
    Assembly::BilinearOperatorAssembler::assemble_matrix1(
      matrix,           // the matrix that receives the assembled operator
      laplace_operator, // the operator that is to be assembled
      space,            // the finite element space in use
      cubature_factory  // the cubature factory to be used for integration
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: RHS (linear forms)

    std::cout << "Assembling right-hand-side vector..." << std::endl;

    // The assembly of right-hand-side vector follows pretty much the same generic structure.
    // We use the opportunity to explain how to prescribe systems with analytically-known solutions.

    // In this tutorial, we first choose a reference solution for our PDE and then assemble the
    // corresponding right-hand-side for our Poisson problem. We choose the "sine-bubble" function,
    // which is pre-defined as a 'common function' in the "analytic/common.hpp" header,
    // so we can use it here. The SineBubbleFunction is implemented for 1D, 2D and 3D, so we
    // need to specify the desired dimension, which we can obtain from the ShapeType definition:
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;

    // Next, we need a linear functional that can be applied onto test functions for our
    // right-hand-side vector. The sine-bubble is an eigenfunction of the Laplace operator,
    // so the corresponding right-hand-side function is the solution multiplied by 2*pi^2.
    // We could exploit this by simply using our solution function for the right-hand-side
    // assembly and passing the constant factor as a multiplier for our assembly method,
    // but we will instead use the pre-defined LaplaceFunctional wrapper, which will
    // compute the right-hand-side force for any given solution function based on its
    // second derivatives. The LaplaceFunctional requires the type of the solution function
    // as its one and only template parameter, so we use the decltype specifier here:
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);

    // Now we can call the LinearFunctionalAssembler class to assemble our linear
    // functional into a vector.
    Assembly::LinearFunctionalAssembler::assemble_vector(
      vec_rhs,          // the vector that receives the assembled functional
      force_functional, // the functional that is to be assembled
      space,            // the finite element space in use
      cubature_factory  // the cubature factory to be used for integration
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions..." << std::endl;

    // The next step is the assembly of the homogeneous Dirichlet boundary conditions.
    // For this task, we require a Unit-Filter assembler:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // Now we need to add all boundary mesh parts to the assembler on which we want to prescribe
    // the boundary conditions. In this tutorial, we have only one boundary object which describes
    // the whole domain's boundary, so add it:
    unit_asm.add_mesh_part(boundary);

    // Now, we need to assemble a unit-filter representing homogeneous Dirichlet BCs.
    // This is done by calling the 'assemble' function:
    FilterType filter;
    unit_asm.assemble(filter, space);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions..." << std::endl;

    // We have assembled the boundary conditions, but the linear system does not know about that
    // yet. So we need to apply the filter onto the system matrix and both vectors now.

    // Apply the filter onto the system matrix...
    filter.filter_mat(matrix);

    // ...the right-hand-side vector...
    filter.filter_rhs(vec_rhs);

    // ...and the solution vector.
    filter.filter_sol(vec_sol);

    // Now we have set up the linear system representing our discretised Poisson PDE, including
    // the homogeneous Dirichlet boundary conditions.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    // For this tutorial, we stick to a simple PCG-SSOR solver.

    std::cout << "Solving linear system..." << std::endl;

    // First, we need to create a SSOR preconditioner. Most of the solver implementations
    // have corresponding 'Solver::new_***' functions, which take care of the nasty type deduction
    // process. The returned object is a 'std::shared_ptr<SolverType>', where 'SolverType' is
    // an instance of the actual solver class template with the corresponding template arguments
    // that are deducted from the input arguments of the following function call:
    auto precond = Solver::new_ssor_precond(matrix, filter);

    // Now we create a PCG solver and pass the preconditioner as an additional argument:
    auto solver = Solver::new_pcg(matrix, filter, precond);

    // PCG is an iterative solver, so we want to enable plotting of the convergence process,
    // as otherwise our solver would not print anything to the console:
    solver->set_plot_mode(Solver::PlotMode::iter);

    // Next, we need to initialise the solver. During this call, the solver and all of its
    // sub-solvers and preconditioners allocate required temporary vectors, perform factorisation
    // and all the other stuff that our solvers need to do before they can actually start
    // solving anything.
    solver->init();

    // Solve our linear system; for this, we pass our solver object as well as the initial solution
    // vector, the right-hand-side vector, the matrix and the filter defining the linear system
    // that we intend to solve to the 'Solver::solve' function:
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // Once we do not require the solver anymore, we have to release it. This is done my calling
    // the 'done' member function, which is the counterpart of the 'init' member function, i.e.
    // this will release all temporary vectors and factorisations.
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    // We have a discrete solution now, so have to do something with it.
    // Since this is a simple benchmark problem, we know the analytical solution of our PDE, so
    // we may compute the L2- and H1-errors of our discrete solution against it.

    std::cout << "Computing errors against reference solution..." << std::endl;

    // The class responsible for this is the 'ScalarErrorComputer' assembly class template.
    // The one and only template parameter is the maximum desired error derivative norm, i.e.
    // setting the parameter to
    //   = 0 will compute only the H0- (aka L2-) error
    //   = 1 will compute both the H0- and H1-errors
    //   = 2 will compute the H0-, H1- and H2-errors

    // We have already created the 'sine_bubble' object representing our analytical solution for
    // the assembly of the right-hand-side vector, so we may reuse it for the computation now:

    Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
      vec_sol,          // the coefficient vector of the discrete solution
      sol_function,     // the analytic function object, declared for RHS assembly
      space,            // the finite element space
      cubature_factory  // and the cubature factory used for integration
    );

    // The returned ScalarErrorInfo object contains all computed error norms,
    // so we may print the errors by simply pushing the object to cout:
    std::cout << errors << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we also write the discrete solution to a VTK file, so that it can be admired in Paraview...

    // First of all, build the filename string
    String vtk_name(String("./tutorial-01-poisson-lvl") + stringify(level));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

    // Next, project our solution into the vertices. This is not necessary for Q1, but in case that
    // someone chose to use the Q1~ element instead of Q1, this will be necessary.

    // First, declare a vector that will receive the vertex projection of our solution.
    // We will also project and write out our right-hand-side, just for fun...
    VectorType vertex_sol, vertex_rhs;

    // And use the DiscreteVertexProjector class to do the dirty work:
    Assembly::DiscreteVertexProjector::project(
      vertex_sol,   // the vector that receives the projection
      vec_sol,      // the vector to be projected
      space         // the finite element space in use
    );

    // And the same for the right-hand-side:
    Assembly::DiscreteVertexProjector::project(vertex_rhs, vec_rhs, space);

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // add the vertex-projection of our solution and rhs vectors
    exporter.add_vertex_scalar("sol", vertex_sol.elements());
    exporter.add_vertex_scalar("rhs", vertex_rhs.elements());

    // finally, write the VTK file
    exporter.write(vtk_name);

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
  } // void main(...)
} // namespace Tutorial01

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialise the FEAT runtime environment:
  Runtime::initialise(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #01: Poisson" << std::endl;

  // Specify the desired mesh refinement level, defaulted to 3.
  // Note that FEAT uses its own "Index" type rather than a wild mixture of int, uint, long
  // and such.
  Index level(3);

  // Now let's see if we have command line parameters: This tutorial supports passing
  // the refinement level as a command line parameter, to investigate the behaviour of the L2/H1
  // errors of the discrete solution.
  if(argc > 1)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[argc-1]).parse(ilevel) || (ilevel < 1))
    {
      // Failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[argc-1] << "' as refinement level." << std::endl;
      std::cerr << "Note: The last argument must be a positive integer." << std::endl;
      // Abort our runtime environment
      Runtime::abort();
    }
    // If parsing was successful, use the given information and notify the user
    level = Index(ilevel);
    std::cout << "Refinement level: " << level << std::endl;
  }
  else
  {
    // No command line parameter given, so inform the user that defaults are used
    std::cout << "Refinement level (default): " << level << std::endl;
  }

  // call the tutorial's main function
  Tutorial01::main(level);

  // And finally, finalise our runtime environment. This function returns the 'EXIT_SUCCESS' return code,
  // so we can simply return this as the result of our main function to indicate a successful run.
  return Runtime::finalise();
}
