// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
// u(x,y) = (exp(-(2x - 1)^2) - exp(-1))*(exp(-(2y - 1)^2) - exp(-1)) / (1 - exp(-1))^2
//
// and its corresponding force (right-hand-side) f.
//
//
// The purpose of this tutorial is to demonstrate the basic program flow of a simple
// scalar stationary linear PDE solver without going to far into the details of what
// magic happens under-the-hood.
//
//
// The basic program flow of this application is as follows:
//
// 1. Define the required types for spatial discretization and linear algebra.
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
#include <kernel/runtime.hpp>                              // for Runtime
#include <kernel/util/string.hpp>                          // for String

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
#include <kernel/space/lagrange3/element.hpp>              // the Lagrange-3 Element (aka "Q3")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>                      // for ExpBubbleFunction

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler
#include <kernel/assembly/domain_assembler_helpers.hpp>    // for Assembly::assemble_***
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for LaplaceFunctional
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
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
  // (2) specialize FEAT to do, out of the many possibilities, only what we want to do in
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

  // Use the Lagrange-1 element (aka "Q1" or "P1", depending on whether ShapeType is a
  // simplex or a hypercube shape type):
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  // Or you could also use the Lagrange-2 element (aka "Q2" or "P2") instead:
  //typedef Space::Lagrange2::Element<TrafoType> SpaceType;

  // Or you could also use the Lagrange-3 element (aka "Q3" or "P3") instead:
  //typedef Space::Lagrange3::Element<TrafoType> SpaceType;


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

  // First, we have to choose the "Data-Index" type duo, which specifies the basic types
  // that our linear algebra containers will use for their elements.

  // The first type is the data type: This is simply the type of the matrix and vector elements,
  // which is typically the double-precision floating point type aka "double". FEAT also supports
  // other data types as single or even quadruple precision, but we want to avoid that for now.
  typedef double DataType;

  // The second type is the index type: This is used by various containers which also store arrays
  // of indices, such as the row-pointer and column-index arrays of the CSR matrix format.
  // In particular, any (sufficiently large) integer type will do, so we stick to the "Index" type,
  // which corresponds to "unsigned long" by default.
  typedef Index IndexType;


  // Based on the two data and index typedefs, we can now define the vector type.
  // In this tutorial, we want to solve a simple scalar PDE, so we require just a "standard"
  // vector class, which is implemented by the "LAFEM::DenseVector" class template:
  typedef LAFEM::DenseVector<DataType, IndexType> VectorType;

  // Furthermore, for the discretized Poisson operator, we require a scalar sparse matrix type.
  // We choose the famous CSR format here, because it is pretty much standard for unstructured FEM:
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixType;

  // Finally, we need a filter. Filters are responsible for "enforcing" simple linear constraints
  // such as boundary conditions and are required by the linear solver framework. In this tutorial,
  // we have a scalar PDE with Dirichlet boundary conditions and for this type of problem, we
  // require a so-called "unit-filter":
  typedef LAFEM::UnitFilter<DataType, IndexType> FilterType;

  // That's it for the linear algebra types.

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function.
  void main(Index level)
  {
    // Okay, we have already defined the mesh, trafo and space types, so we can start our
    // actual tutorial code by creating those.

    std::cout << "Creating Mesh on Level " << level << "...\n";

    // In a "real-life" application, we would either read the mesh from an input file or call
    // some sort of mesh generator library/tool to create a mesh for us, but in this tutorial,
    // we want to stick with a more simple solution, especially to avoid dependencies on external
    // files, libraries or tools. Instead, we will use the "RefinedUnitCubeFactory", which will
    // (as the name suggests) generate a refined unit-square mesh for us.

    // First of all, create a mesh factory object representing a refined unit-square domain
    // and pass the desired refinement level to its constructor. The only purpose of this
    // "factory" object is to create a mesh for us and we will not need it anymore after that.
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);

    // Now create the actual mesh by using that factory:
    MeshType mesh(mesh_factory);

    // Furthermore, we require a mesh-part that represents the boundary of the mesh for the
    // assembly of Dirichlet boundary conditions later on. Again, this would typically come
    // from an external file or a mesh generator, but we stick to a more simple solution.
    // We will utilize the "BoundaryFactory", which will create a boundary mesh-part for a
    // given mesh. Note that this BoundaryFactory class works for any given input mesh and
    // not only for meshes created by the RefinedUnitCubeFactory.

    // Now let's create a boundary factory for our mesh.
    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);

    // And create the boundary mesh-part by using the factory.
    MeshPartType boundary(boundary_factory);

    // Voila, we now have a mesh and a corresponding boundary mesh-part.


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Trafo and Finite Element Space initialization

    std::cout << "Creating Trafo and Space...\n";

    // We have already defined the types of the transformation and finite element space,
    // so we just need to create the objects.

    // Let's create a trafo object now. Its only parameter is the mesh that it is defined on.
    // At this point, it is important to mention that the trafo object will keep an internal reference
    // to the mesh object that is passed to its constructor here, so the mesh object has to be alive
    // as long at the trafo object exists. This is not a problem here, because both the mesh and
    // the trafo object reside on the stack and are automatically cleaned up upon the exit of this
    // main function, but it may require some attention when you organize the mesh and trafo objects
    // on the heap, possibly as member variables of some class object as e.g. in tutorial 05.
    TrafoType trafo(mesh);

    // Create the desire finite element space. Its only parameter is the trafo that it is defined on.
    // Again, the space keeps a reference to the trafo object throughout its whole lifetime.
    SpaceType space(trafo);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly

    std::cout << "Allocating matrix and vectors...\n";

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
    // allocated, but their data arrays are still uninitialized (automatic initialization is
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

    // In the next step we are going to assemble a bilinear form into the system matrix and a linear
    // form into the right hand side vector. These (bi)linear forms consist of domain integrals and
    // therefore we need to create a 'domain assembler' here, which will be responsible for the
    // assembly of these forms and we will also use this domain assembler to perform an a posteriori
    // error analysis by computing the H0- and H1-errors against a given reference solution.
    // We create a domain assembler for our trafo, which internally also contains a reference to
    // our mesh. Note that the domain assembler keeps the trafo reference internally.
    Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);

    // Next, we need to tell the assembler on which elements of the mesh we want to assemble.
    // We want to assemble on the whole mesh, which is the usual case, so we just have to call
    // the 'compile_all_elements' function here:
    domain_assembler.compile_all_elements();

    // Before we start assembling anything, we need to chose a cubature rule, which will be
    // used for the integration of the domain integrals.
    // There are several cubature rules available in FEAT, and a complete list can be queried
    // by compiling and executing the 'cub-list' tool in the 'tools' directory.
    // We choose the 'auto-degree:n' rule here, which will automatically choose an appropriate
    // cubature rule that can integrate polynomials up to degree n exactly. Some other possible
    // choices are also provided here, but commented out.
    String cubature_name =
      "auto-degree:5";          // automatic cubature rule for 5th degree polynomials
    //"gauss-legendre:3";       // 3x3 Gauss-Legendre rule
    //"gauss-lobatto:4";        // 4x4 Gauss-Lobatto rule
    //"newton-cotes-open:5";    // 5x5 'open' Newton-Cotes rule
    //"trapezoidal";            // trapezoidal rule (works for all shape types)
    //"barycentre";             // barycentre rule (not recommended due to insufficient order)

    std::cout << "Assembling system matrix...\n";

    // We want to assemble the 2D "-Laplace" operator (-u_xx -u_yy).
    // In this tutorial, we use a pre-defined class for the implementation of this operator.
    // The definition of other and more complex PDE operators is discussed in another tutorial.
    Assembly::Common::LaplaceOperator laplace_operator;

    // Now that we have our operator, we can use the domain assembler to assemble it by using the
    // 'assemble_bilinear_operator_matrix_1' helper function. It may seem confusing that this
    // is a global function and not a member function of the domain assembler in the first place,
    // but this is for technical (dependency) reasons, so you just have to get used to it.
    // In analogy to the SymbolicAssembler class, the corresponding assemble function also has
    // a "1" suffix indicating that there is only one finite element space involved in the process:
    Assembly::assemble_bilinear_operator_matrix_1(
      domain_assembler, matrix, laplace_operator, space, cubature_name);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: RHS (linear forms)

    std::cout << "Assembling right-hand-side vector...\n";

    // The assembly of right-hand-side vector follows pretty much the same generic structure.
    // We use the opportunity to explain how to prescribe systems with analytically-known solutions.

    // In this tutorial, we first choose a reference solution for our PDE and then assemble the
    // corresponding right-hand-side for our Poisson problem. We choose the "exp-bubble" function,
    // which is pre-defined as a 'common function' in the "analytic/common.hpp" header,
    // so we can use it here. The ExpBubbleFunction is implemented for 1D, 2D and 3D, so we
    // need to specify the desired dimension, which we can obtain from the ShapeType definition:
    Analytic::Common::ExpBubbleFunction<ShapeType::dimension> sol_function;

    // Next, we have to define a functional for the right hand side of our PDE. Normally, one
    // would define the right-hand-side function 'f' explicitly here, but since this is a simple
    // benchmark problem we already know the analytical solution, so we can define the right-hand-
    // side functional from it automatically. For this, we use the pre-defined LaplaceFunctional
    // wrapper, which computes the right-hand-side force for any given solution function based on
    // its second derivatives. The LaplaceFunctional requires the type of the solution function
    // as its one and only template parameter, so we use the 'decltype' specifier here, which
    // returns the class type of its argument, so that we do not have to write that out again:
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);

    // Now we can use the domain assembler to assemble the linear form for us by using another
    // helper function, just in analogy to the bilinear form assembly:
    Assembly::assemble_linear_functional_vector(
      domain_assembler, vec_rhs, force_functional, space, cubature_name);


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions...\n";

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

    std::cout << "Imposing boundary conditions...\n";

    // We have assembled the boundary conditions, but the linear system does not know about that
    // yet. So we need to apply the filter onto the system matrix and both vectors now.

    // Apply the filter onto the system matrix...
    filter.filter_mat(matrix);

    // ...the right-hand-side vector...
    filter.filter_rhs(vec_rhs);

    // ...and the solution vector.
    filter.filter_sol(vec_sol);

    // Now we have set up the linear system representing our discretized Poisson PDE, including
    // the homogeneous Dirichlet boundary conditions.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    // For this tutorial, we stick to a simple PCG-SSOR solver.

    std::cout << "Solving linear system...\n";

    // First, we need to create a SSOR preconditioner. Most of the solver implementations
    // have corresponding 'Solver::new_***' functions, which take care of the nasty type deduction
    // process. The returned object is a 'std::shared_ptr<SolverType>', where 'SolverType' is
    // an instance of the actual solver class template with the corresponding template arguments
    // that are deducted from the input arguments of the following function call:
    auto precond = Solver::new_ssor_precond(PreferredBackend::generic, matrix, filter);

    // Now we create a PCG solver and pass the preconditioner as an additional argument:
    auto solver = Solver::new_pcg(matrix, filter, precond);

    // PCG is an iterative solver, so we want to enable plotting of the convergence process,
    // as otherwise our solver would not print anything to the console:
    solver->set_plot_mode(Solver::PlotMode::iter);

    // Set the maximum number of iterations to 1000:
    solver->set_max_iter(1000);

    // Next, we need to initialize the solver. During this call, the solver and all of its
    // sub-solvers and preconditioners allocate required temporary vectors, perform factorization
    // and all the other stuff that our solvers need to do before they can actually start
    // solving anything.
    solver->init();

    // Solve our linear system; for this, we pass our solver object as well as the initial solution
    // vector, the right-hand-side vector, the matrix and the filter defining the linear system
    // that we intend to solve to the 'Solver::solve' function:
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // Once we do not require the solver anymore, we have to release it. This is done by calling
    // the 'done' member function, which is the counterpart of the 'init' member function, i.e.
    // this will release all temporary vectors and factorizations.
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing H0/H1-Errors

    // We have a discrete solution now, so have to do something with it.
    // Since this is a simple benchmark problem, we know the analytical solution of our PDE, so
    // we may compute the H0- and H1-errors of our discrete solution against it.

    std::cout << "\n";
    std::cout << "Computing errors against reference solution...\n";

    // We have to use the domain assembler for the actual error computation and we are going to
    // use another helper function here to do the job for us. This helper function returns a
    // structure that is actually just an instance of the 'Assembly::FunctionIntegralInfo' class
    // template, however, the template arguments are not trivial to determine and so we simply
    // use the 'auto' keyword here to let the compile determine the correct type.
    // Note that this helper function requires an integer template parameter which describes the
    // maximum desired derivative to be used in the norm integral computation, i.e.
    //   = 0 will compute only the H0- (aka L2-), L1- and Lmax-errors
    //   = 1 will compute both the H0-, L1-, Lmax- and H1-errors
    //   = 2 will compute the H0-, L1-, Lmax-1, H1- and H2-errors
    auto error_info = Assembly::integrate_error_function<1>(
      domain_assembler,      // the domain assembler
      sol_function,          // the analytic reference solution function
      vec_sol,               // the coefficient vector of the FEM solution
      space,                 // the finite element space
      cubature_name          // the cubature name used for integration
    );

    // The error_info structure, which is an instance of the Assembly::FunctionIntegralInfo class
    // template defined in <kernel/assembly/function_integral_jobs.hpp>, contains all the computed
    // norms and integrals as member functions, which we could access here individually. However,
    // we simply want to use the 'print_norms' function, which returns a formatted multi-line
    // string containing all the error norm values and
    std::cout << "Error Analysis:\n";
    std::cout << error_info.print_norms() << "\n";

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we also write the discrete solution to a VTK file, so that it can be admired in Paraview...

    // First of all, build the filename string
    String vtk_name(String("./tutorial-01-poisson-lvl") + stringify(level));

    std::cout << "\n";
    std::cout << "Writing VTK file '" << vtk_name << ".vtu'...\n";

    // Next, project our solution into the vertices. This is not necessary for Q1, but in case that someone
    // chose to use the Rannacher-Turek or Crouzeix-Raviart element instead of Q1, this *will* be necessary.

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
    std::cout << "Finished!\n";
  } // void main(...)
} // namespace Tutorial01

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialize the FEAT runtime environment by
  // creating a scope guard object, which will take care of releasing the runtime once we're done:
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #01: Poisson\n";

  // Specify the desired mesh refinement level, defaulted to 3.
  // Note that FEAT uses its own "Index" type rather than a wild mixture of int, uint, long
  // and such.
  Index level(3);

  // Now let's see if we have command line parameters: This tutorial supports passing
  // the refinement level as a command line parameter, to investigate the behavior of the L2/H1
  // errors of the discrete solution.
  if(argc > 1)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[argc-1]).parse(ilevel) || (ilevel < 1))
    {
      // Failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[argc-1] << "' as refinement level.\n";
      std::cerr << "Note: The last argument must be a positive integer.\n";
      // Abort our runtime environment
      Runtime::abort();
    }
    // If parsing was successful, use the given information and notify the user
    level = Index(ilevel);
    std::cout << "Refinement level: " << level << "\n";
  }
  else
  {
    // No command line parameter given, so inform the user that defaults are used
    std::cout << "Refinement level (default): " << level << "\n";
  }

  // call the tutorial's main function
  Tutorial01::main(level);

  // And finally, return the exit code 0 to indicate a successful run.
  return 0;
}
