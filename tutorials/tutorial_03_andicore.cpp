// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// \brief FEAT Tutorial 03: Andicore solver (TM)
//
// This file contains a simple prototypical Anisotropic Diffusion-Convection-Reaction
// solver for the unit square domain.
//
// The PDE to be solved reads:
//
//    -div(A*grad(u)) + dot(b,grad(u)) + c*u = f          in the domain [0,1]x[0,1]
//                                         u = 0          on the boundary
//
// where
//  - A is a (constant) positive-definite matrix in R^{2x2}
//  - b is a (constant) vector in R^2
//  - c is a (constant) non-negative scalar
//
// The differential operator written as a weak bilinear operator is given as:
//
//  a(phi,psi) = < A*grad(phi),grad(psi) > + < b, grad(phi) > * psi + c*phi*psi
//
// The analytical solution u is given as
//
//         u(x,y) = sin(pi*x) * sin(pi*y)
//
// The purpose of this tutorial is to demonstrate
//
// 1. The implementation of more complex differential operators including custom data.
//
// 2. The implementation of a custom functional for the right-hand-side assembly.
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
#include <kernel/analytic/common.hpp>                      // for SineBubbleFunction

// FEAT-Assembly includes
#include <kernel/assembly/linear_functional.hpp>           // NEW: for LinearOperator
#include <kernel/assembly/bilinear_operator.hpp>           // NEW: for BilinearOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler
#include <kernel/assembly/domain_assembler_helpers.hpp>    // for Assembly::assemble_***
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter

// FEAT-Solver includes
#include <kernel/solver/jacobi_precond.hpp>                  // for JacobiPrecond
#include <kernel/solver/bicgstab.hpp>                      // for BiCGStab

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial03
{
  // Note:
  // This tutorial works only for 2D shapes, i.e. quadrilaterals or triangles.
  // The reason for this is that the implementation of the 'AndicoreData' class below as well as
  // its setup at the beginning of the 'main' function is restricted to 2D for the sake of keeping
  // the code simple. However, the remainder of this tutorial code works for any shape type.

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

  // This tutorial covers an "anisotropic diffusion-convection-reaction" PDE, which is slightly more
  // complex than the Poisson example, so we'll need to define a few custom classes for the assembly
  // before we can start with our actual main function.

  // The first thing we require is a class that will contain the data for our PDE, namely the
  // matrix A, the vector b and the scalar c:
  class AndicoreData
  {
  public:
    // First of all, we need a R^{2x2} matrix. The class for this job is the "Tiny" matrix:
    Tiny::Matrix<Real, 2, 2> a;

    // Note that 'Real' is a FEAT typedef, which equals 'double'.

    // Next, the R^2 vector b. For this one we have the "Tiny" vector:
    Tiny::Vector<Real, 2> b;

    // And the scalar c:
    Real c;
  }; // class AndicoreData

  // Now we need to define a custom bilinear operator for our PDE. There is no pre-defined "andicore"
  // operator in the kernel, so we have to write one by hand. Each class that represents a simple
  // integral-based bilinear form is derived from the Assembly::BilinearOperator class:
  class AndicoreOperator :
    public Assembly::BilinearOperator
  {
  protected:
    // Our bilinear operator requires the data that is encapsulated in the 'AndicoreData' class, so
    // we require a (const) reference to such an object:
    const AndicoreData& _data;

  public:
    // And a simple initialization constructor for our operator:
    explicit AndicoreOperator(const AndicoreData& data) :
      _data(data)
    {
    }

    // The first thing that is required from a bilinear operator are three so-called 'config tags',
    // which tell the assembly routines which information from the trafo, the test- and the
    // trial-spaces are required for this operator.

    // We start off with the trafo configuration, which is a constexpr TrafoTags member variable
    // named 'trafo_config'. Our operator does not require any information from the trafo,
    // so we simply set it to TrafoTags::none:
    static constexpr TrafoTags trafo_config = TrafoTags::none;

    // Test- and trial-space configurations are specified in the same way as trafo configurations.
    // Our bilinear operator requires both test-/trial-function values and gradients, so we
    // need to combine the corresponding value with the OR operator:
    static constexpr SpaceTags test_config = SpaceTags::value | SpaceTags::grad;

    // And we specify the very same thing for the trial space:
    static constexpr SpaceTags trial_config = SpaceTags::value | SpaceTags::grad;

    // Up to now, we have a reference to our data and have expressed our wishes to the assembler,
    // but we still need to implement the actual operator somewhere. For this purpose, we need to
    // implement a so-called 'evaluator'. The evaluator is a class template, whose template is a
    // so-called 'assembly traits' class containing various information that we will discuss in
    // a moment.

    // Declare the evaluator class template and derive it from the base-class version:
    template<typename AsmTraits_>
    class Evaluator :
      public Assembly::BilinearOperator::Evaluator<AsmTraits_>
    {
    public:
      // The assembly traits class 'AsmTraits_' contains a handful of typedefs which we require.
      // The first type is the trafo data, which contains all required data from the transformation
      // in the current cubature point. This includes e.g. the coordinates of the current point,
      // the jacobian matrix, the hessian and other useful stuff that we're not interested in
      // for this operator.
      typedef typename AsmTraits_::TrafoData TrafoData;

      // The next type is the trial basis function data type. This is the type of the object that
      // contains the function value and gradient of the trial function, which we have requested
      // from the assembler in the 'TrialConfig' class (typedef) before:
      typedef typename AsmTraits_::TrialBasisData TrialBasisData;

      // In analogy, there is a test basis function data type which we require:
      typedef typename AsmTraits_::TestBasisData TestBasisData;

      // Moreover, the assembly traits contains the type of the trafo evaluator.
      typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;

      // Next, we also need the scalar DataType typedef.
      typedef typename AsmTraits_::DataType DataType;

      // Finally, we need the ValueType typedef, which is identical to DataType for scalar operators.
      typedef DataType ValueType;

    protected:
      // The first thing that we require is again a const reference to our data object:
      const AndicoreData& _data;

    public:
      // And now the constructor: An instance of this evaluator class is created by the assembler,
      // which passes a const reference to the corresponding operator to the evaluator's constructor.
      // By this approach, we can get the reference to our data object:
      explicit Evaluator(const AndicoreOperator& operat) :
        _data(operat._data)
      {
      }

      // As all other evaluator classes (e.g. the one from the analytic function in
      // 'tutorial_02_laplace'), also the bilinear operator evaluator has a 'prepare-finish'
      // function pair. However, as we do not require them, we may simply skip implementing them,
      // as the base-class we derived from already contains a 'do-nothing' implementation for those.

      // Finally, we can provide the actual evaluation function, which performs the actual work:
      // This function returns a (scalar) value and its parameters are const references to...
      ValueType eval(
        const TrialBasisData& phi,  // ...the trial basis data object
        const TestBasisData& psi    // ...the test basis data object
        ) const
      {
        // And this is where the magic happens.
        // First, we have to take a look at our weak bilinear operator:
        //  a(phi,psi) = dot(A*grad(phi),grad(psi)) + dot(b, grad(phi)) * psi + c*phi*psi

        // The operator is composed of three additive terms:
        // 1. the diffusion  :  dot(A*grad(phi),grad(psi))
        // 2. the convection :  dot(b, grad(phi) * psi
        // 3. the reaction   :  c * phi * psi

        // We will now combine these three terms into one statement:
        return
          // First, we take the diffusion term, which is a matrix-based scalar product of the trial-
          // and test-function gradients. We have stored the matrix 'A' as a Tiny::Matrix in our
          // AndicoreData object, and the basis function gradients are of type Tiny::Vector, so we
          // can call the matrix's member function to compute the scalar product for us:
          _data.a.scalar_product(phi.grad, psi.grad) +

          // Next, the convection term: This a dot-product of our vector 'b' and the gradient of
          // the trial function multiplied by the test function value:
          dot(_data.b, phi.grad) * psi.value +

          // Finally, the reaction term, which is simply a product of the scalar 'c' and the trial-
          // and test-function values:
          _data.c * phi.value * psi.value;
      }
    }; // class AndicoreOperator::Evaluator<...>
  }; // class AndicoreOperator


  // The next class we need to define is a corresponding functional for the right-hand-side, which
  // matches our analytical solution. For this task, we write the functional as a template which
  // accepts a class implementing the AnalyticalSolution interface as a template parameter, see the
  // 'tutorial_02_laplace' example for a documentation of this interface and its implementation.
  // This class template will combine the derivatives of the analytical solution and the data of
  // our operator into a corresponding linear functional for the right-hand side, which derives
  // from the LinearFunctional base class defined in the assembly subsystem:
  template<typename SolFunction_>
  class AndicoreFunctional :
    public Assembly::LinearFunctional
  {
  protected:
    // The contents of a class implementing the LinearFunctional interface are similar to the
    // BilinearOperator interface that we have implemented above for our operator.

    // First, we require a reference to our data object:
    const AndicoreData& _data;

    // Moreover, we require a reference to the solution function object for the evaluation
    // of the right-hand-side functional.
    const SolFunction_& _sol_func;

  public:
    // And again, we need a corresponding constructor:
    explicit AndicoreFunctional(const AndicoreData& data, const SolFunction_& sol_func) :
      _data(data),
      _sol_func(sol_func)
    {
    }

    // In analogy to the bilinear operator, we first need to provide a configuration for the trafo
    // and one for the test space (there exists no trial space for functionals).
    // We will need to evaluate our solution function and for this, we must tell the trafo
    // that we require image point coordinates, so let's state this:
    static constexpr TrafoTags trafo_config = TrafoTags::img_point;

    // Now for the test-space configuration. For a linear functional, we usually require test
    // function values for the evaluation, so define a corresponding test-space config:
    static constexpr SpaceTags test_config = SpaceTags::value;

    // In analogy to operators, functionals are evaluated using 'evaluator' class templates as well,
    // so let's declare one of those and derive it from the base-class evaluator:
    template<typename AsmTraits_>
    class Evaluator :
      public Assembly::LinearFunctional::Evaluator<AsmTraits_>
    {
    public:
      // The 'AsmTraits_' class template, which is passed to this evaluator template by the
      // assembler, is the same as the one passed on for bilinear operators, so we can again
      // query our required types from it. As functionals do not operate on trial basis functions,
      // we only require the trafo data and the test basis data types:
      typedef typename AsmTraits_::TrafoData TrafoData;
      typedef typename AsmTraits_::TestBasisData TestBasisData;

      // Moreover, we need another type that we have not used before: the 'trafo evaluator', which
      // is again an evaluator class similar to the one we are currently building up.
      typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;

      // Next, we also need the DataType typedef.
      typedef typename AsmTraits_::DataType DataType;

      // Now comes the interesting part:
      // We want to evaluate the analytical solution, and for this, we need to get its corresponding
      // evaluator, which is again a class template similar to the one that you are currently reading,
      // see the 'SaddleFunction' from the 'tutorial_02_laplace' example.
      // However, its template parameter is *not* our assembly traits type 'AsmTraits_', but another
      // type called 'analytic evaluation traits', which we need to provide here.
      // Luckily, there is another template class in the 'Analytic' namespace which takes care
      // of all the necessary typedef'ing. The only two template parameters are the data type which
      // we want to use internally and the function class itself:
      typedef Analytic::EvalTraits<DataType, SolFunction_> AnalyticEvalTraits;

      // The return value of the analytic function, which can be scalar or vector-valued in general,
      // determines the ValueType for this functional, too. Since we are using a scalar function
      // in this tutorial, ValueType is again identical to DataType here:
      typedef typename AnalyticEvalTraits::ValueType ValueType;

      // With that type, we can now define the type of the solution function evaluator:
      typedef typename SolFunction_::template Evaluator<AnalyticEvalTraits> SolEvaluator;

    protected:
      // As usual, a reference to our data object:
      const AndicoreData& _data;

      // Now, declare one of those solution function evaluators:
      SolEvaluator _sol_eval;

      // Finally, declare a variable that will store the value of the force in the current cubature point:
      ValueType _force_value;

    public:
      // And a matching constructor, which receives a reference to the functional object as its
      // only parameter. This constructor will first get a reference to our custom data object and
      // will create the solution function evaluator using its mandatory constructor:
      explicit Evaluator(const AndicoreFunctional& functional) :
        _data(functional._data),
        _sol_eval(functional._sol_func),
        _force_value(DataType(0))
      {
      }

      // Each type of evaluator has a 'prepare' and a 'finish' member function which are called
      // for each cell of the mesh that the assembler iterates over. However, as this simple
      // functional does not require any type of additional initialization, we do not need to
      // implement these function, as our base-class provides empty implementations for us.

      // Now we come to learn a new function, which is present in all bilinear operator and
      // linear functional evaluator classes. The next function is called for each cubature point
      // before the actual evaluation operator is called. At this point, we have the chance to
      // pre-compute any data that depends on the current cubature point, but not on the test-function.
      void set_point(const TrafoData& tau)
      {
        // In analogy to the bilinear operator, this functional is given by three additive terms:
        // 1. the diffusion  : -dot(A,hess(u)) * psi
        // 2. the convection :  dot(b,grad(u)) * psi
        // 3. the reaction   :      c *    u   * psi

        // For this task, we need to call the 'value', 'gradient' and 'hessian' functions of the
        // solution function evaluator (see the SaddleFunction in Tutorial 02), so we need to
        // obtain the corresponding types of the return values from the AnalyticEvalTraits class,
        // which we have defined above. All we have to do now is to call the corresponding functions
        // of the function's evaluator class and supply the image point of our trafo data as the
        // evaluation point:
        typename AnalyticEvalTraits::ValueType    value = _sol_eval.value(tau.img_point);
        typename AnalyticEvalTraits::GradientType grad  = _sol_eval.gradient(tau.img_point);
        typename AnalyticEvalTraits::HessianType  hess  = _sol_eval.hessian(tau.img_point);

        // Finally, we can combine these terms to pre-compute the value of our force
        // functional, so that we just need to multiply by 'psi' later on:
        _force_value =
          // First, the diffusion. This is the matrix-dot product of 'A' and the hessian of 'u':
          - dot(_data.a, hess)
          // Next, the convection. This is the vector-dot-product of 'b' and the gradient of 'u':
          + dot(_data.b, grad)
          // Finally, the reaction. This is a simple product:
          + _data.c * value;
      }

      // Once again, we implement our evaluation function, which is quite similar to the one
      // for bilinear operators, with the only difference being that there is no trial function.
      // It takes a const reference to the test-function data object as a parameter
      // and returns the evaluation of the functional in the current cubature point:
      ValueType eval(const TestBasisData& psi) const
      {
        // As we have already computed the value of our force functional in the "set_point" function
        // above, we just need to multiply by the function value of 'psi':
        return _force_value * psi.value;
      }
    }; // class AndicoreFunctional<...>::Evaluator<...>
  }; // class AndicoreFunctional<...>

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function
  void main(Index level)
  {
    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Initial data setup

    // First of all, we create and initialize the data object for our PDE:
    AndicoreData andicore_data;

    // Choose a (non-symmetric) positive definite matrix A for the diffusion
    andicore_data.a(0,0) =  1.2;
    andicore_data.a(0,1) = -0.7;
    andicore_data.a(1,0) = -0.4;
    andicore_data.a(1,1) =  0.9;

    // Choose some vector for the convection
    andicore_data.b(0) =  0.4;
    andicore_data.b(1) = -0.2;

    // Choose a non-negative value for the reaction
    andicore_data.c = 0.3;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

    // First of all, format the matrix entries to zero.
    matrix.format();

    // At this point, we need to create an object representing our custom bilinear operator:
    AndicoreOperator andicore_operator(andicore_data);

    // And assemble that operator
    Assembly::assemble_bilinear_operator_matrix_1(
      domain_assembler, matrix, andicore_operator, space, cubature_name);

    std::cout << "Assembling right-hand-side vector..." << "\n";

    // Format the right-hand-side vector entries to zero.
    vec_rhs.format();

    // The next task is the assembly of the right-hand-side vector.
    // As our functional is parameterized in the solution function, we first need to create an
    // object representing that solution function. We choose the standard sine-bubble:
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;

    // Now, we utilize our custom AndicoreFunctional to obtain a functional that combines
    // the desired analytical solution with our andicore PDE operator:
    AndicoreFunctional<decltype(sol_function)> force_functional(andicore_data, sol_function);

    // And assemble our linear functional:
    Assembly::assemble_linear_functional_vector(
      domain_assembler, vec_rhs, force_functional, space, cubature_name);

    // Finally, clear the initial solution vector.
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions..." << "\n";

    // In this example, we assemble homogeneous Dirichlet boundary conditions:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // Add our only boundary component:
    unit_asm.add_mesh_part(boundary);

    // And assemble the filter
    FilterType filter;
    unit_asm.assemble(filter, space);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions..." << "\n";

    // Apply the filter onto the system matrix...
    filter.filter_mat(matrix);

    // ...the right-hand-side vector...
    filter.filter_rhs(vec_rhs);

    // ...and the solution vector.
    filter.filter_sol(vec_sol);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    std::cout << "Solving linear system..." << "\n";

    // Create a SPAI preconditioner
    auto precond = Solver::new_jacobi_precond(matrix, filter);

    // Create a preconditioned BiCGStab solver
    auto solver = Solver::new_bicgstab(matrix, filter, precond);

    // Enable convergence plot
    solver->set_plot_mode(Solver::PlotMode::iter);

    // Initialize the solver
    solver->init();

    // Solve our linear system
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // And release the solver
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    std::cout << "Computing errors against reference solution..." << "\n";

    // Compute the error norms:
    auto error_info = Assembly::integrate_error_function<1>(
      domain_assembler, sol_function, vec_sol, space, cubature_name);

    // Print the error norms to the console
    std::cout << "Error Analysis:" << "\n";
    std::cout << error_info.print_norms() << "\n";

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we also write the discrete solution to a VTK file, so that it can be admired in Paraview...

    // First of all, build the filename string
    String vtk_name(String("./tutorial-03-andicore-lvl") + stringify(level));

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
} // namespace Tutorial03

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialize the runtime
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #03: Andicore" << "\n";

  // The desired mesh refinement level.
  Index level(3);

  // First of all, let's see if we have command line parameters.
  if(argc > 1)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[argc-1]).parse(ilevel) || (ilevel < 1))
    {
      // failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[argc-1] << "' as refinement level." << "\n";
      std::cerr << "Note: The last argument must be a positive integer." << "\n";
      // Abort our runtime environment
      Runtime::abort();
    }
    // parse successful
    level = Index(ilevel);
  }

  // call the tutorial's main function
  Tutorial03::main(level);

  return 0;
}
