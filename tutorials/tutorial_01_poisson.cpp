//
// \brief FEAST Tutorial 01: Poisson solver (TM)
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
// The purpose of this tutorial is to demonstate the basic program flow of a simple
// scalar stationary linear PDE solver without going to far into the details of what
// magic happens under-the-hood.
//
//
// The basic program flow of this application is as follows:
//
// 1. Create a mesh and a boundary description by using a factory.
//
// 2. Create a trafo based on the mesh and a finite element space based on the trafo.
//
// 3. Assemble the PDE operator and the force to a matrix and a right-hand-side vector.
//
// 4. Assemble the boundary conditions to a filter.
//
// 5. Solve the resulting linear system using a simple fire-and-forget solver.
//
// 6. Compute the L2- and H1-errors against the analytical reference solution.
//
// 7. Write the result to a VTK file for visual post-processing, if desired.
//
// \author Peter Zajac
//

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We start our little tutorial with a batch of includes...

// Misc. FEAST includes
#include <kernel/util/string.hpp>                          // for String

// FEAST-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/cell_sub_set.hpp>                // for CellSubSet
#include <kernel/geometry/conformal_factories.hpp>         // for RefinedUnitCubeFactor
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK

// FEAST-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAST-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")
//#include <kernel/space/rannacher_turek/element.hpp>      // the Rannacher-Turek Element (aka "Q1~")

// FEAST-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAST-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/dirichlet_assembler.hpp>         // for DirichletAssembler
#include <kernel/assembly/error_computer.hpp>              // for L2/H1-error computation
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for ForceFunctional
#include <kernel/assembly/common_functions.hpp>            // for SineBubbleFunction

// FEAST-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter

// FEAST-LAFEM provisional solver includes
#include <kernel/lafem/preconditioner.hpp>                 // for NonePreconditioner
#include <kernel/lafem/bicgstab.hpp>                       // for BiCGStab

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAST
using namespace FEAST;

// We're opening a new namespace for our tutorial.
namespace Tutorial01
{
  // We start with a set of typedefs, which make up the basic configuration for this
  // tutorial application. The general idea of these typedefs is (1) to avoid typing and
  // (2) specialise FEAST to do, out of the many possibilities, only what we want to do in
  // this tutorial.

  // First, we define the 'shape type' of the mesh that we want to use.
  // We use a quadrilateral mesh for this tutorial, so typedef the quadrilateral shape.
  typedef Shape::Quadrilateral ShapeType;

  // Next, we need to define the memory type, the data type and the algorithm type for
  // the linear algebra containers and operations. Details about these types will be
  // discussed in another tutorial.

  // We want double precision.
  typedef double DataType;

  // Moreover, we use main memory (aka "RAM") for our containers. FEAST supports GPUs and other
  // architectures, but we do not want to make use of anything that does not reside in host memory.
  typedef Mem::Main MemType;

  // And we'll use the Generic algorithm implementation for now. Similar to the last typedef,
  // this essentially means that we turn off any fancyness in the implementation, which is
  // appropriate for the scope of this tutorial -- we do not need GPUs or anything like that.
  typedef Algo::Generic AlgoType;

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function
  int main(Index level)
  {
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Geometry typedefs

    // First, we need to define the type of mesh we want to use.
    // For the purpose of this tutorial, we employ a ConformalMesh, whose first template parameter is
    // the shape-type that we want to use, and which we typedef'ed to Quadrilateral above.
    // All other parameters are optional and we're not interested in using anything else than the
    // standard parameters for now.
    typedef Geometry::ConformalMesh<ShapeType> MeshType;

    // Moreover, we need to define a boundary type, which is required for the assembly of boundary
    // conditions. In this example, we employ the CellSubSet class template for this purpose, which is
    // perfectly appropriate for the pure Dirichlet case we want to solve. Note that this class template
    // is also templatised in the shape that we're using.
    typedef Geometry::CellSubSet<ShapeType> BoundaryType;

    // In this tutorial, we will generate the mesh and its boundary information by using a so-called
    // factory instead of reading them from a file. This is a convenient shortcut that FEAST provides,
    // the more general case is covered in advanced tutorials.
    // For our purpose here, we need to define a RefinedUnitCubeFactory, which will (as the name suggests)
    // generate a refined unit-square mesh for us.
    typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactoryType;

    // Furthermore, we need a factory for the cell sub-set representing our domain's boundary.
    // The BoundaryFactory class template will generate such a cell sub-set representing the whole
    // boundary of the domain, so we also typedef it here.
    typedef Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Geometry initialisation

    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;

    // First of all, create a mesh factory object representing a refined unit-square domain according
    // to the refinement scheme specified by the above type definitions.
    MeshFactoryType mesh_factory(level);

    // Create a mesh by using our factory.
    MeshType mesh(mesh_factory);

    std::cout << "Creating Boundary..." << std::endl;

    // Now let's create a boundary factory for our mesh.
    BoundaryFactoryType boundary_factory(mesh);

    // And create the boundary by using the factory.
    BoundaryType boundary(boundary_factory);

    // Voila, we now have a mesh and a corresponding boundary cell sub-set.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Trafo initialisation

    // The next thing we need is a transformation that our finite element spaces should use.
    // Here, we use the standard transformation, which we'll typedef here.
    // The Mapping class template takes the mesh type as a template parameter.
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    std::cout << "Creating Trafo..." << std::endl;

    // Let's create a trafo object now. Its only parameter is the mesh that it is defined on.
    TrafoType trafo(mesh);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Finite Element Space initialisation

    // Now we need to define the finite element space which shall be used for our Poisson problem.
    // All finite element spaces are parameterised (templated) by the transformation type.

    // Use the Lagrange-1 element (aka "Q1"):
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // Use the Rannacher-Turek element (aka "Q1~"):
    //typedef Space::RannacherTurek::Element<TrafoType> SpaceType;

    std::cout << "Creating Space..." << std::endl;

    // Create the desire finite element space. Its only parameter is the trafo that it is defined on.
    SpaceType space(trafo);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Linear System type definitions

    // The next step is defining the types of vectors and matrices which we want to use.

    // Based on the two data and memory typedefs, we now define the vector type.
    // Since we have a scalar problem, we require a dense vector, which is simply FEAST's name for
    // contiguous storage.
    typedef LAFEM::DenseVector<MemType, DataType> VectorType;

    // Furthermore, for the discretised Poisson operator, we require a sparse matrix type.
    // We choose the CSR format here, because it is pretty much standard.
    typedef LAFEM::SparseMatrixCSR<MemType, DataType> MatrixType;

    // Moreover, we need a boundary condition filter. Filters can be used to apply boundary conditions
    // in matrices and vectors. Since we are using Dirichlet boundary conditions for our problem, we
    // need a so-called "unit filter".
    typedef LAFEM::UnitFilter<MemType, DataType> FilterType;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly
    std::cout << "Allocating and initialising vectors and matrix..." << std::endl;

    // First, allocate two vectors: the solution vector and the right-hand-side vector.
    // The lengths of these vectors generally coincide with the number of DOFs of the space:
    VectorType vec_sol(space.get_num_dofs());
    VectorType vec_rhs(space.get_num_dofs());

    // Now we need to perform the symbolic matrix assembly, i.e., the computation of the nonzero
    // pattern and the allocation of storage for the three CSR arrays.
    // The class which performs this task is the SymbolicMatrixAssembler class.
    MatrixType matrix;
    Assembly::SymbolicMatrixAssembler<>::assemble1(matrix, space);
    // Note: The "1" at the end of the function name "assemble1" indicates, that there
    // is only one finite element space involved in the assembly, i.e. the test- and trial-
    // spaces are identical. The symbolic and numeric operator assemblers also contain function named
    // "assemble2", which takes two (possibly different) finite element spaces as
    // test- and trial-spaces, but we're not going to use it for now. However, this is
    // required for more complex PDEs like, e.g., the Stokes equations.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)

    // Before we start assembling anything, we need create a cubature factory, which will be
    // used to generate cubature rules for integration.
    // There are several cubature rules available in FEAST, and a complete list can be queried
    // by compiling and executing the 'cub-list' tool in the 'tools' directory.
    // We choose the 'auto-degree:n' rule here, which will automatically choose an appropriate
    // cubature rule  that can integrate polynomials up to degree n exactly. Some other possible
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

    // The assembly routines require the matrix to filled with zeros, so we use a LAFEM routine to
    // zero out everything
    matrix.format();

    // We want to assemble the 2D "-Laplace" operator -u_xx -u_yy.
    // Note that so far, we did not consider the notion of operators at all.
    // For this operator, FEAST provides a pre-defined bilinear operator class, of which we first need
    // to create an instance. Check the implementation of this factory, and additional tutorials, on how
    // to implement more complex operators.
    Assembly::Common::LaplaceOperator laplace_operator;

    // Next, we call the bilinear operator assembler to assemble the operator into a matrix.
    Assembly::BilinearOperatorAssembler::assemble_matrix1(
        matrix,           // the matrix that receives the assembled operator
        laplace_operator, // the operator that is to be assembled
        space,            // the finite element space in use
        cubature_factory  // the cubature factory to be used for integration
        );

    std::cout << "Assembling right-hand-side vector..." << std::endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: RHS (linear forms)

    // The assembly follows pretty much the same generic structure. We use the opportunity to
    // explain how to prescribe systems with analytically-known solutions.

    // Initialise the right-hand-side vector entries to zero.
    vec_rhs.format();

    // First, we need an analytic function representing the right-hand-side, or more generally
    // speaking, any force.
    // In our example, the analytical solution of the PDE is the sine-bubble,
    // which is an eigenfunction of the Laplace operator, so the corresponding
    // right-hand-side is the solution multiplied by 2*pi^2. We exploit this by we'll simply
    // using our solution function for the right-hand-side assembly and pass
    // the constant factor as a multiplier for our assembly method. More general linear form
    // asemblies are covered in later tutorials.

    // The sine-bubble function is pre-defined as a 'common function' in the assembly
    // "common_functions.hpp" header, so we'll use it here.
    typedef Assembly::Common::SineBubbleFunction SolFunction;
    SolFunction sol_function;

    // Next, we need a linear functional that can be applied onto test functions.
    // We utilise the Common::ForceFunctional class template to convert our
    // sine-bubble function into a functional.
    Assembly::Common::ForceFunctional<SolFunction> force_functional(sol_function);

    // Now we can call the LinearFunctionalAssembler class to assemble our linear
    // functional into a vector.
    Assembly::LinearFunctionalAssembler::assemble_vector(
        vec_rhs,          // the vector that receives the assembled functional
        force_functional, // the functional that is to be assembled
        space,            // the finite element space in use
        cubature_factory, // the cubature factory to be used for integration
        2.0 * Math::sqr(Math::pi<DataType>()) // this is our factor 2*pi^2
        );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Solution initialisation.

    // Clear the initial solution vector (because its contents are undefined up to now).
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions..." << std::endl;

    // The next step is the assembly of the homogeneous Dirichlet boundary conditions.
    // For this task, we require a Dirichler assembler:
    Assembly::DirichletAssembler<SpaceType> dirichlet_asm(space);

    // Now we need to add all boundary parts to the assembler on which we want to prescribe
    // the boundary conditions. In this tutorial, we have only one boundary object which describes
    // the whole domain's boundary, so add it:
    dirichlet_asm.add_cell_set(boundary);

    // Now, we need to assemble a unit-filter representing homogeneous Dirichlet BCs.
    // This is done by first creating a filter of the corresponding size and then calling the 'assemble'
    // function:
    FilterType filter(space.get_num_dofs());
    dirichlet_asm.assemble(filter);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions..." << std::endl;

    // We have assembled the boundary conditions, but the linear system does not know about that
    // yet. So we need to apply the filter onto the system matrix and both vectors now.

    // Apply the filter onto the system matrix...
    filter.filter_mat<AlgoType>(matrix);

    // ...the right-hand-side vector...
    filter.filter_rhs<AlgoType>(vec_rhs);

    // ...and the solution vector.
    filter.filter_sol<AlgoType>(vec_sol);

    // Now we have set up the linear system representing our discretised Poisson PDE, including
    // the homogene Dirichlet boundary conditions.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    // Now we should configure our fancy ScaRC2 solver, but unfortunately, we're not that far...
    // So we stick to a provisional BiCGStab solver for now...

    // Create a dummy preconditioner
    LAFEM::NonePreconditioner<AlgoType, MatrixType, VectorType> precond;

    std::cout << "Solving linear system..." << std::endl;

    // Fire up the BiCGStab solver
    LAFEM::BiCGStab<AlgoType>::value(
      vec_sol,    // the initial solution vector
      matrix,     // the system matrix
      vec_rhs,    // the right-hand-side vector
      precond,    // the dummy preconditioner
      100,        // maximum number of iterations
      1E-8        // relative tolerance to achieve
      );

    // Unless the BiCGStab died, we should now have our discrete solution.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    // We have a discrete solution now, so have to do something with it.
    // Since this is a simple benchmark problem, we know the analytical solution of our PDE, so
    // we may compute the L2- and H1-errors of our discrete solution against it.

    std::cout << "Computing errors against reference solution..." << std::endl;

    // The classes responsible for this are the 'ScalarErrorComputerXY' assembly classes.
    // We have already created the 'sine_bubble' object representing our analytical solution for
    // the assembly of the right-hand-side vector, so we may reuse it for the computation now:

    DataType l2_error = Assembly::ScalarErrorComputerL2::compute(
      vec_sol,          // the coefficient vector of the discrete solution
      sol_function,     // the analytic function object, declared for RHS assembly
      space,            // the finite element space
      cubature_factory  // and the cubature factory used for integration
      );

    // In analogy, compute the H1-error:
    DataType h1_error = Assembly::ScalarErrorComputerH1::compute(vec_sol, sol_function, space, cubature_factory);

    // Print the errors to cout
    std::cout << "L2-Error: " << scientify(l2_error, 12) << std::endl;
    std::cout << "H1-Error: " << scientify(h1_error, 12) << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we also write the discrete solution to a VTK file, so that it can be admired in Paraview...

    // First of all, build the filename string
    String vtk_name(String("./tutorial-01-poisson-lvl") + stringify(level) + ".vtk");

    std::cout << "Writing VTK file '" << vtk_name << "'..." << std::endl;

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
    exporter.add_scalar_vertex("solution", vertex_sol.elements());
    exporter.add_scalar_vertex("rhs", vertex_rhs.elements());

    // finally, write the VTK file
    exporter.write(vtk_name);

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
    return 0;
  } // int main(...)
} // namespace Tutorial01

// Here's our main function
int main(int argc, char* argv[])
{
  // Print a welcome message
  std::cout << "Welcome to FEAST's tutorial #01: Poisson" << std::endl;

  // Specify the desired mesh refinement level, defaulted to 3.
  // Note that FEAST uses its own "Index" type rather than a wild mixture of int, uint, long
  // and such.
  Index level(3);

  // First of all, let's see if we have command line parameters: This tutorial supports passing
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
      return 1;
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
  return Tutorial01::main(level);
}
