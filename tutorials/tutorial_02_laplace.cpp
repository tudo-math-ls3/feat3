// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
// The purpose of this tutorial is to demonstrate the implementation of a custom
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

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/function.hpp>                    // NEW: for Analytic::Function

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler
#include <kernel/assembly/domain_assembler_helpers.hpp>    // for Assembly::assemble_***
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
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
  // Note:
  // This tutorial works only for 2D shapes, i.e. quadrilaterals or triangles.
  // The reason for this is that the implementation of the 'SaddleFunction' class below
  // is restricted to 2D for the sake of keeping the code simple. However, the remainder of this
  // tutorial code works for any shape type.

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

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // As a first step, we want to implement an 'Analytic Function' class, which will represent both our
  // analytical solution 'u' (for post-processing) as well as the boundary condition function 'g'.
  // For this, we derive our own solution function by defining a class that derives from the
  // 'Analytic::Function' base-class:
  class SaddleFunction :
    public Analytic::Function
  {
  public:
    // For analytic functions, we first need to provide information about the domain and the image
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
      explicit Evaluator(const SaddleFunction&)
      {
      }

      // At the beginning of our SaddleFunction class definition, we told the assembler that we
      // are capable of computing function values, so we have to provide a function for this job.
      // This function is called 'value', its one and only parameter is the point in which we want
      // to evaluate and it returns the function value of type 'ValueType', which we have declared
      // by a typedef a few lines above:
      ValueType value(const PointType& point) const
      {
        // We can now return the value of our function
        //  u(x,y) = (x - 1/2)^2 - (y - 1/2)^2
        return Math::sqr(point[0] - DataType(0.5)) - Math::sqr(point[1] - DataType(0.5));
      }

      // The next function that we need to supply is the gradient evaluation function, which also
      // receives the evaluation point and returns the gradient:
      GradientType gradient(const PointType& point) const
      {
        // Create a gradient object (which is a vector):
        GradientType grad;
        // Set the X-derivative of our function:
        // dx u(x,y) =  2*x - 1
        grad[0] =  DataType(2) * point[0] - DataType(1);
        // Set the Y-derivative of our function:
        // dy u(x,y) = -2*y + 1
        grad[1] = -DataType(2) * point[1] + DataType(1);
        // And finally return the gradient:
        return grad;
      }

      // And finally, one evaluation function for the hessian:
      HessianType hessian(const PointType& /*point*/) const
      {
        // Create a hessian object (which is a matrix):
        HessianType hess;
        // The mixed derivatives dx dy u are zero
        hess[0][1] = hess[1][0] = DataType(0);
        // The second XX-derivate: dxx u(x,y) =  2
        hess[0][0] =  DataType(2);
        // The second YY-derivate: dyy u(x,y) = -2
        hess[1][1] = -DataType(2);
        return hess;
      }
    }; // class SaddleFunction::Evaluator<...>
  }; // class SaddleFunction

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our main function
  void main(Index level)
  {
    // Create Mesh, Boundary, Trafo and Space

    std::cout << "Creating Mesh on Level " << level << "..." << "\n";

    // Create the mesh
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    // And create the boundary
    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);

    // Create Mesh, Boundary, Trafo and Space

    std::cout << "Creating Trafo and Space..." << "\n";

    // Let's create a trafo object now.
    TrafoType trafo(mesh);

    // Create the desire finite element space.
    SpaceType space(trafo);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly

    std::cout << "Allocating matrix and vectors..." << "\n";

    // Allocate matrix and assemble its structure
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    // Allocate vectors
    VectorType vec_sol = matrix.create_vector_r();
    VectorType vec_rhs = matrix.create_vector_l();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Perform numerical matrix assembly

    // Create a domain assembler on all mesh elements
    Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
    domain_assembler.compile_all_elements();

    // Choose a cubature rule
    String cubature_name = "auto-degree:5";

    std::cout << "Assembling system matrix..." << "\n";

    // First of all, format the matrix and vector entries to zero.
    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // Create the pre-defined Laplace operator:
    Assembly::Common::LaplaceOperator laplace_operator;

    // And assemble that operator
    Assembly::assemble_bilinear_operator_matrix_1(
      domain_assembler, matrix, laplace_operator, space, cubature_name);

    // Note that the force functional is zero for the Laplace equation, therefore we do not
    // have to assemble the right-hand-side vector.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions..." << "\n";

    // The next step is the assembly of the inhomogeneous Dirichlet boundary conditions.
    // For this task, we require a Unit-Filter assembler:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // Add the one and only boundary component to the assembler:
    unit_asm.add_mesh_part(boundary);

    // To assemble inhomogeneous Dirichlet boundary conditions, we need to tell the assembler the
    // function 'g' that defines the boundary values. In this tutorial, we simply pass an object
    // that represents our analytical solution 'u' for this task, as it will give the correct
    // boundary values, of course. So create an instance of our analytical function now:
    SaddleFunction sol_function;

    // And assemble a unit-filter representing inhomogeneous Dirichlet BCs; This is done by calling
    // the 'assemble' function, to which we pass our boundary value function object as the third
    // argument to the function:
    FilterType filter;
    unit_asm.assemble(filter, space, sol_function);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions..." << "\n";

    // Apply the filter onto the system matrix...
    filter.filter_mat(matrix);
    // ...the right-hand-side vector...
    filter.filter_rhs(vec_rhs);
    // ...and the solution vector.
    filter.filter_sol(vec_sol);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    std::cout << "Solving linear system..." << "\n";

    // Create a SSOR preconditioner
    auto precond = Solver::new_ssor_precond(PreferredBackend::generic, matrix, filter);

    // Create a PCG solver
    auto solver = Solver::new_pcg(matrix, filter, precond);

    // Enable convergence plot
    solver->set_plot_mode(Solver::PlotMode::iter);

    // Initialize the solver
    solver->init();

    // Solve our linear system
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // Release the solver
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    std::cout << "Computing errors against reference solution..." << "\n";

    // Compute the error norms:
    auto error_info = Assembly::integrate_error_function<1>(
      domain_assembler, sol_function, vec_sol, space, cubature_name);

    // Print the error norms to the console
    std::cout << "Error Analysis:" << "\n";
    std::cout << error_info.print_norms() << "\n";


    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we also write the discrete solution to a VTK file, so that it can be admired in Paraview...

    // First of all, build the filename string
    String vtk_name(String("./tutorial-02-laplace-lvl") + stringify(level));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << "\n";

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

    // That's it for today.
    std::cout << "Finished!" << "\n";
  } // void main(...)
} // namespace Tutorial02

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialize the runtime
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #02: Laplace" << "\n";

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
      std::cerr << "ERROR: Failed to parse '" << argv[argc-1] << "' as refinement level." << "\n";
      std::cerr << "Note: The last argument must be a positive integer." << "\n";
      Runtime::abort();
    }
    // parse successful
    level = Index(ilevel);
  }

  // call the tutorial's main function
  Tutorial02::main(level);

  return 0;
}
