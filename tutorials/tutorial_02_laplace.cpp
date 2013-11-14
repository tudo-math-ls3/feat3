//
// \brief FEAST Tutorial 02: Laplace solver (TM)
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
#include <kernel/assembly/analytic_function.hpp>           // NEW: for AnalyticFunction

// FEAST-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter

// FEAST-LAFEM provisional solver includes
#include <kernel/lafem/preconditioner.hpp>                 // for NonePreconditioner
#include <kernel/lafem/bicgstab.hpp>                       // for BiCGStab


// We are using FEAST
using namespace FEAST;

// Once again, we use quadrilaterals.
typedef Shape::Quadrilateral ShapeType;
// We want double precision.
typedef double DataType;
// Moreover, we use main memory (aka "RAM") for our containers.
typedef Mem::Main MemType;
// And we'll use the Generic algorithm implementation for now.
typedef Algo::Generic AlgoType;

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// As a first step, we want to implement an 'AnalyticFunction' class, which will represent both our
// analytical solution 'u' (for post-processing) as well as the boundary condition function 'g'.
// For this, we derive our own solution function:
class PringlesFunction :
  public Assembly::AnalyticFunction
{
public:
  // For analytic functions, we first need to provide information about what this function is
  // capable of computing. An object of this class will not only be responsible for computing the
  // function values of our function 'u', but also its derivatives -- namely its gradient and its
  // hessian.
  // To inform the assembler about what computations we can perform, we define an enumeration,
  // which is just a C++ trick to get compile-time constants.
  // See the documentation of the Assembly::AnalyticFunction class for the names of these constants.
  enum
  {
    // We are capable of computing the function values of 'u'.
    can_value = 1,
    // Moreover, we can also compute the gradient of 'u'. This will be required for the computation
    // of the H1-error during the post-processing step.
    can_grad = 1,
    // Moreover, we can also compute the hessian of 'u', although we will not need it in this
    // tutorial.
    can_hess = 1
  };

  // For the computation of these values, gradients and hessians we require some data from the
  // underlying transformation, namely the real point coordinates. These wishes are managed by
  // so-called 'configurations'. The following class template receives a configuration telling
  // us what the assembler expects us to provide. We're not interested in that class for now,
  // so we simply pass this argument to the base-class version.

  // In our case, we always need the same configuration:
  template<typename AnalyticConfig_>
  struct ConfigTraits :
    public Assembly::AnalyticFunction::ConfigTraits<AnalyticConfig_>
  {
    // This config traits class needs to contain a struct named 'TrafoConfig', which collects
    // our wishes for the transformation. We derive our trafo-config class from the basic version:
    // Define the required trafo configuration:
    struct TrafoConfig :
      public Trafo::ConfigBase
    {
      enum
      {
        // At this point we express our wishes. The only thing that we require from the
        // transformation are the 'real' point coordinates (i.e. 'x' and 'y'), so that we may
        // evaluate our analytic function 'u' and/or its derivatives in these coordinates.
        need_img_point = 1
      };
    };
  }; // class PringlesFunction::ConfigTraits<...>

  // Up to now, we have only declared our capabilities and expressed our wishes to the assembler,
  // but we still need to implement our analytical function formula somewhere. For this purpose,
  // we need to implement a so-called 'Evaluator'. The Evaluator is a nested class template,
  // whose template parameter is a so-called 'evaluation traits' class containing various useful
  // typedefs that we will require in a moment.

  // Declare our evaluator and derive it from the base-class version.
  template<typename EvalTraits_>
  class Evaluator :
    public Assembly::AnalyticFunction::Evaluator<EvalTraits_>
  {
  public:
    // First off, we require a mandatory constructor that takes a const reference to our
    // analytic function object as its one and only parameter:
    explicit Evaluator(const PringlesFunction&)
    {
    }

    // The evaluation traits class 'EvalTraits_' contains a handful of typedefs which we require.
    // The first type is the trafo-evaluator, which defines the mapping from the reference cell
    // to the real cell.
    typedef typename EvalTraits_::TrafoEvaluator TrafoEvaluator;

    // Next comes the trafo data, which contains all required data from the transformation
    // in the current cubature point. This includes, e.g., the coordinates of the current point,
    // the jacobian matrix, the hessian and other useful stuff that we're not interested in
    // for this simple function.
    typedef typename EvalTraits_::TrafoData TrafoData;

    // Moreover, it contains a typedef for the currently used data-type, which will coincide with
    // our global DataType typedef.
    typedef typename EvalTraits_::DataType DataType;

    // The next thing is the function-value type, which coincides with 'DataType' for scalar functions:
    typedef typename EvalTraits_::ValueType ValueType;

    // Then we have the function-gradient type, which is a Tiny::Vector.
    typedef typename EvalTraits_::GradientType GradientType;

    // And finally, the function-hessian type, which is a Tiny::Matrix.
    typedef typename EvalTraits_::HessianType HessianType;

    // Each evaluator is equipped with a 'prepare'-'finish' function pair, which is called by
    // the assembler for each cell during the assembly process.

    // We first provide our prepare function. This function is called by the assembler each
    // time it starts assembling on a new cell of the mesh, so the evaluator has the chance
    // to pre-compute data required for its evaluation.
    // The only parameter to this function is a const reference to the trafo-evaluator:
    void prepare(const TrafoEvaluator& /* trafo_eval */)
    {
      // Our function does not require any initialisation...
    }

    // Moreover, once the assembler is finished with one particular cell of the mesh, it calls
    // the evaluator's finish function, giving it the chance to clean up any auxiliary data
    // created by the prepare function:
    void finish()
    {
      // We did not initialise anything, so we don't have to clean up either...
    }

    // At the beginning of our PringlesFunction class definition, we told the assembler that we are
    // capable of computing function values, so we have to provide a function for this job.
    // This function is called 'value'; its only parameter is a const reference to the trafo data
    // object and its return type is the function value type:
    ValueType value(const TrafoData& tau) const
    {
      // The 'tau' parameter contains the coordinates of our evaluation point. These are stored
      // in the member named 'img_point', which is simply a tuple of scalar coordinates.
      // We can now return the value of our function
      //  u(x,y) = (x - 1/2)^2 - (y - 1/2)^2
      return Math::sqr(tau.img_point[0] - DataType(0.5)) - Math::sqr(tau.img_point[1] - DataType(0.5));
    }

    // The next function that we need to supply is the gradient evaluation function:
    GradientType gradient(const TrafoData& tau) const
    {
      // create an auxiliary gradient:
      GradientType grad;
      // Set the X-derivative of our function: dx u(x,y) = 2*x - 1
      grad(0) =  DataType(2) * tau.img_point[0] - DataType(1);
      // Set the Y-derivative of our function: dy u(x,y) = -2*y + 1
      grad(1) = -DataType(2) * tau.img_point[1] + DataType(1);
      // and return the gradient:
      return grad;
    }

    // And finally, one evaluation function for the hessian:
    HessianType hessian(const TrafoData& tau) const
    {
      // create an auxiliary hessian:
      HessianType hess;
      // The mixed derivatives dx dy u are zero
      hess(0,1) = hess(1,0) = DataType(0);
      // The second XX-derivate: dxx u(x,y) =  2
      hess(0,0) =  DataType(2);
      // The second YY-derivate: dyy u(x,y) = -2
      hess(1,1) = -DataType(2);
      // and return the hessian:
      return hess;
    }
  }; // class PringlesFunction::Evaluator<...>
}; // class PringlesFunction

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// Here's our main function
int main(int argc, char* argv[])
{
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
      return 1;
    }
    // parse successful
    level = Index(ilevel);
  }

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Create Mesh and Boundary

  // Define the mesh type
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  // Define the boundary type
  typedef Geometry::CellSubSet<ShapeType> BoundaryType;
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

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Allocate linear system and perform symbolic assembly

  // Define the vector type
  typedef LAFEM::DenseVector<MemType, DataType> VectorType;

  // Define the matrix type
  typedef LAFEM::SparseMatrixCSR<MemType, DataType> MatrixType;

  // Define the filter type
  typedef LAFEM::UnitFilter<MemType, DataType> FilterType;

  std::cout << "Allocating vectors and matrix..." << std::endl;

  // Allocate vectors
  VectorType vec_sol(space.get_num_dofs());
  VectorType vec_rhs(space.get_num_dofs());

  // Allocate matrix and assemble its structure
  MatrixType matrix;
  Assembly::SymbolicMatrixAssembler<>::assemble1(matrix, space);

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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

  // Please note that the force functional is zero for the Laplace equation, therefore we do not
  // have to assemble the right-hand-side vector.

  // Format the right-hand-side vector entries to zero.
  vec_rhs.format();

  // Finally, clear the initial solution vector.
  vec_sol.format();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Boundary Condition assembly

  std::cout << "Assembling boundary conditions..." << std::endl;

  // The next step is the assembly of the inhomogeneous Dirichlet boundary conditions.
  // For this task, we require a Dirichler assembler:
  Assembly::DirichletAssembler<SpaceType> dirichlet_asm(space);

  // Add the one and only boundary component to the assembler:
  dirichlet_asm.add_cell_set(boundary);

  // To assemble inhomogeneous Dirichlet boundary conditions, we need to tell the assembler the
  // function 'g' that defines the boundary values. In this tutorial, we simply pass an object
  // that represents our analytical solution 'u' for this task, as it will give the correct
  // boundary values, of course. So create an instance of our analytical function now:
  PringlesFunction sol_function;

  // And assemble a unit-filter representing inhomogene Dirichlet BCs; This is done by calling
  // the 'assemble' function, to which we pass our boundary value function object as the second
  // argument to the function:
  FilterType filter(space.get_num_dofs());
  dirichlet_asm.assemble(filter, sol_function);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Boundary Condition imposition

  std::cout << "Imposing boundary conditions..." << std::endl;

  // Apply the filter onto the system matrix...
  filter.filter_mat<AlgoType>(matrix);
  // ...the right-hand-side vector...
  filter.filter_rhs<AlgoType>(vec_rhs);
  // ...and the solution vector.
  filter.filter_sol<AlgoType>(vec_sol);

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Solver set-up

  // Create a dummy preconditioner
  LAFEM::NonePreconditioner<AlgoType, MatrixType, VectorType> precond;

  std::cout << "Solving linear system..." << std::endl;

  // Fire up the BiCGStab solver
  LAFEM::BiCGStab<AlgoType>::value(vec_sol, matrix, vec_rhs, precond, 100, 1E-8);

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Post-Processing: Computing L2/H1-Errors

  std::cout << "Computing errors against reference solution..." << std::endl;

  // We have already created an object representing our solution function for the assembly of the
  // boundary conditions, so we may re-use it for the error computation now.

  // Compute the L2-error:
  DataType l2_error = Assembly::ScalarErrorComputerL2::compute(vec_sol, sol_function, space, cubature_factory);

  // Compute the H1-error:
  DataType h1_error = Assembly::ScalarErrorComputerH1::compute(vec_sol, sol_function, space, cubature_factory);

  // And prints the errors to cout
  std::cout << "L2-Error: " << scientify(l2_error, 12) << std::endl;
  std::cout << "H1-Error: " << scientify(h1_error, 12) << std::endl;

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Post-Processing: Export to VTK file

  // Finally, we also write the discrete solution to a VTK file, so that it can be admired in Paraview...

  // First of all, build the filename string
  String vtk_name(String("./tutorial-02-laplace-lvl") + stringify(level) + ".vtk");

  std::cout << "Writing VTK file '" << vtk_name << "'..." << std::endl;

  // project solution and right-hand-side vectors
  VectorType vertex_sol, vertex_rhs;
  Assembly::DiscreteVertexProjector::project(vertex_sol, vec_sol, space);
  Assembly::DiscreteVertexProjector::project(vertex_rhs, vec_rhs, space);

  // Create a VTK exporter for our mesh
  Geometry::ExportVTK<MeshType> exporter(mesh);

  // add the vertex-projection of our solution and rhs vectors
  exporter.add_scalar_vertex("solution", vertex_sol.elements());
  exporter.add_scalar_vertex("rhs", vertex_rhs.elements());

  // finally, write the VTK file
  exporter.write(vtk_name);

  // That's it for today.
  std::cout << "Finished!" << std::endl;
  return 0;
}
