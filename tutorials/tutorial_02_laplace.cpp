//
// \brief FEAT Tutorial 02: Laplace solver (TM)
//
// This file contains a simple prototypical Laplace solver for the unit square domain.
//
// The PDE to be solved reads:
//
//    -Laplace(u) = 0          in the domain [0,1]x[0,1]
//             u  = g          on the boundary
//
// The analytical solution u is given as
//
//         u(x,y) = (x - 1/2)^2 - (y - 1/2)^2
//
// The purpose of this tutorial is to demonstate the implementation of a custom
// AnalyticalFunction class for the solution function 'u' and the boundary condition
// function 'g'.
//
// The basic program flow of this application coincides with the one from the
// 'tutorial_01_poisson' example, therefore the documentation of the code that has
// already been discussed in that tutorial is shortened in this one.
//
// \author Peter Zajac
//

// We start our little tutorial with a batch of includes...

// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime

// FEAT-Geometry includes
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/conformal_factories.hpp>         // for RefinedUnitCubeFactor
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/function.hpp>                    // NEW: for Analytic::Function

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/error_computer.hpp>              // for L2/H1-error computation
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator

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
namespace Tutorial02
{
  // Once again, we use quadrilaterals.
  typedef Shape::Quadrilateral ShapeType;
  // We want double precision.
  typedef double DataType;
  // Use the default index type.
  typedef Index IndexType;
  // Moreover, we use main memory (aka "RAM") for our containers.
  typedef Mem::Main MemType;

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // As a first step, we want to implement an 'Analytic Function' class, which will represent both our
  // analytical solution 'u' (for post-processing) as well as the boundary condition function 'g'.
  // For this, we derive our own solution function:
  class PringlesFunction :
    public Analytic::Function
  {
  public:
    // For analytic functions, we first need to provide information abou the domain and the image
    // sets. The domain of any analytic function in our context is the R^n, so the only thing we
    // have to specify is the domain dimension, which equals 2 as we are solving a 2D PDE:
    static constexpr int domain_dim = 2;

    // Moreover, we need to classify the image set of the function. Here, we have the choice
    // between two variants: either we have a 'scalar' function, i.e. a function that maps into
    // R, or we have a vector field. In our case, this function is a scalar one, so we specify
    // the following typedef:
    typedef Analytic::Image::Scalar ImageType;

    // Next, we need to specify what this function implementation is capable of computing. An object
    // of this class will not only be responsible for computing the function values of our function
    // 'u', but also its derivatives -- namely its gradient and its hessian.
    // To inform the assembler about what computations we can perform, we specify the following
    // the bool values, which specify whether the corresponding information can be computed:

    // We are capable of computing the function values of 'u'.
    static constexpr bool can_value = true;
    // Moreover, we can also compute the gradient of 'u'. This will be required for the computation
    // of the H1-error during the post-processing step.
    static constexpr bool can_grad = true;
    // Moreover, we can also compute the hessian of 'u', although we will not need it in this
    // tutorial.
    static constexpr bool can_hess = true;

    // Up to now, we have only declared our capabilities to the assembler, but we still need to
    // implement our analytical function formula somewhere. For this purpose, we need to implement
    // a so-called 'Evaluator'. The Evaluator is a nested class template, whose template parameter
    // is a so-called 'evaluation traits' class containing various useful typedefs that we will
    // require in a moment.

    // Declare our evaluator and derive it from the base-class version.
    template<typename Traits_>
    class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
    {
    public:
      // We will 'extract' the required data types from the traits class:

      // First, it contains a typedef for the currently used data-type, which will coincide with
      // our global DataType typedef.
      typedef typename Traits_::DataType DataType;

      // Moreover, we have a typedef that specifies the type of the point in which we need
      // to evaluate our function. This is always an instance of the Tiny::Vector class template.
      typedef typename Traits_::PointType PointType;

      // The next thing is the function-value type, which coincides with 'DataType' for scalar functions:
      typedef typename Traits_::ValueType ValueType;

      // Then we have the function-gradient type, which is a Tiny::Vector.
      typedef typename Traits_::GradientType GradientType;

      // And finally, the function-hessian type, which is a Tiny::Matrix.
      typedef typename Traits_::HessianType HessianType;

      // Now we require a *mandatory* constructor that takes a const reference to our
      // analytic function object as its one and only parameter:
      explicit Evaluator(const PringlesFunction&)
      {
      }

      // At the beginning of our PringlesFunction class definition, we told the assembler that we are
      // capable of computing function values, so we have to provide a function for this job.
      // This function is called 'value'; its first parameter is a reference to the 'value'
      // object that we have to fill and its second parameter is the point in which we evaluate:
      void value(ValueType& val, const PointType& point) const
      {
        // We can now return the value of our function
        //  u(x,y) = (x - 1/2)^2 - (y - 1/2)^2
        val = Math::sqr(point[0] - DataType(0.5)) - Math::sqr(point[1] - DataType(0.5));
      }

      // The next function that we need to supply is the gradient evaluation function:
      void gradient(GradientType& grad, const PointType& point) const
      {
        // Set the X-derivative of our function:
        // dx u(x,y) =  2*x - 1
        grad[0] =  DataType(2) * point[0] - DataType(1);
        // Set the Y-derivative of our function:
        // dy u(x,y) = -2*y + 1
        grad[1] = -DataType(2) * point[1] + DataType(1);
      }

      // And finally, one evaluation function for the hessian:
      void hessian(HessianType& hess, const PointType& /*point*/) const
      {
        // The mixed derivatives dx dy u are zero
        hess[0][1] = hess[1][0] = DataType(0);
        // The second XX-derivate: dxx u(x,y) =  2
        hess[0][0] =  DataType(2);
        // The second YY-derivate: dyy u(x,y) = -2
        hess[1][1] = -DataType(2);
      }
    }; // class PringlesFunction::Evaluator<...>
  }; // class PringlesFunction

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our main function
  void main(Index level)
  {
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Create Mesh and Boundary

    // Define the mesh type
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    // Define the boundary type
    typedef Geometry::MeshPart<MeshType> BoundaryType;
    // Define the mesh factory type
    typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactoryType;
    // And define the boundary factory type
    typedef Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;

    // Create the mesh
    MeshFactoryType mesh_factory(level);
    MeshType mesh(mesh_factory);

    std::cout << "Creating Boundary..." << std::endl;

    // And create the boundary
    BoundaryFactoryType boundary_factory(mesh);
    BoundaryType boundary(boundary_factory);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Create Trafo and Space

    // Define the trafo
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    std::cout << "Creating Trafo..." << std::endl;

    // Let's create a trafo object now.
    TrafoType trafo(mesh);

    // Use the Lagrange-1 element (aka "Q1"):
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    std::cout << "Creating Space..." << std::endl;

    // Create the desire finite element space.
    SpaceType space(trafo);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Allocate linear system and perform symbolic assembly

    // Define the vector type
    typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;

    // Define the matrix type
    typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> MatrixType;

    // Define the filter type
    typedef LAFEM::UnitFilter<MemType, DataType, IndexType> FilterType;

    std::cout << "Allocating matrix and vectors..." << std::endl;

    // Allocate matrix and assemble its structure
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    // Allocate vectors
    VectorType vec_sol = matrix.create_vector_r();
    VectorType vec_rhs = matrix.create_vector_l();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Perform numerical matrix assembly

    // Create a cubature factory
    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    std::cout << "Assembling system matrix..." << std::endl;

    // First of all, format the matrix entries to zero.
    matrix.format();

    // Create the pre-defined Laplace operator:
    Assembly::Common::LaplaceOperator laplace_operator;

    // And assemble that operator
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, laplace_operator, space, cubature_factory);

    // Note that the force functional is zero for the Laplace equation, therefore we do not
    // have to assemble the right-hand-side vector.

    // Format the right-hand-side vector entries to zero.
    vec_rhs.format();

    // Finally, clear the initial solution vector.
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions..." << std::endl;

    // The next step is the assembly of the inhomogeneous Dirichlet boundary conditions.
    // For this task, we require a Unit-Filter assembler:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // Add the one and only boundary component to the assembler:
    unit_asm.add_mesh_part(boundary);

    // To assemble inhomogeneous Dirichlet boundary conditions, we need to tell the assembler the
    // function 'g' that defines the boundary values. In this tutorial, we simply pass an object
    // that represents our analytical solution 'u' for this task, as it will give the correct
    // boundary values, of course. So create an instance of our analytical function now:
    PringlesFunction sol_function;

    // And assemble a unit-filter representing inhomogene Dirichlet BCs; This is done by calling
    // the 'assemble' function, to which we pass our boundary value function object as the third
    // argument to the function:
    FilterType filter;
    unit_asm.assemble(filter, space, sol_function);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions..." << std::endl;

    // Apply the filter onto the system matrix...
    filter.filter_mat(matrix);
    // ...the right-hand-side vector...
    filter.filter_rhs(vec_rhs);
    // ...and the solution vector.
    filter.filter_sol(vec_sol);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    std::cout << "Solving linear system..." << std::endl;

    // Create a SSOR preconditioner
    auto precond = Solver::new_ssor_precond(matrix, filter);

    // Create a PCG solver
    auto solver = Solver::new_pcg(matrix, filter, precond);

    // Enable convergence plot
    solver->set_plot_mode(Solver::PlotMode::iter);

    // Initialise the solver
    solver->init();

    // Solve our linear system
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // Release the solver
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    std::cout << "Computing errors against reference solution..." << std::endl;

    // Compute and print the H0-/H1-errors
    Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
      vec_sol, sol_function, space, cubature_factory);

    std::cout << errors << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we also write the discrete solution to a VTK file, so that it can be admired in Paraview...

    // First of all, build the filename string
    String vtk_name(String("./tutorial-02-laplace-lvl") + stringify(level));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

    // project solution and right-hand-side vectors
    VectorType vertex_sol, vertex_rhs;
    Assembly::DiscreteVertexProjector::project(vertex_sol, vec_sol, space);
    Assembly::DiscreteVertexProjector::project(vertex_rhs, vec_rhs, space);

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // add the vertex-projection of our solution and rhs vectors
    exporter.add_vertex_scalar("solution", vertex_sol.elements());
    exporter.add_vertex_scalar("rhs", vertex_rhs.elements());

    // finally, write the VTK file
    exporter.write(vtk_name);

    // That's it for today.
    std::cout << "Finished!" << std::endl;
  } // int main(...)
} // namespace Tutorial02

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialise the runtime
  Runtime::initialise(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #02: Laplace" << std::endl;

  // The desired mesh refinement level, defaulted to 3
  Index level(3);

  // First of all, let's see if we have command line parameters
  if(argc > 1)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[argc-1]).parse(ilevel) || (ilevel < 1))
    {
      // failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[argc-1] << "' as refinement level." << std::endl;
      std::cerr << "Note: The last argument must be a positive integer." << std::endl;
      Runtime::abort();
    }
    // parse successful
    level = Index(ilevel);
  }

  // call the tutorial's main function
  Tutorial02::main(level);

  // Finalise the runtime
  return Runtime::finalise();
}
