//
// \brief FEAST Tutorial 03: Andicore solver (TM)
//
// This file contains a simple prototypical Anisotropic Diffusion-Convection-Reaction
// solver for the unit square domain.
//
// The PDE to be solved reads:
//
//    -grad(A*grad(u)) + dot(b,grad(u)) + c*u = f         in the domain [0,1]x[0,1]
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

// Misc. FEAST includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime

// FEAST-Geometry includes
#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/conformal_factories.hpp>         // for RefinedUnitCubeFactor
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart

// FEAST-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAST-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")

// FEAST-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAST-Assembly includes
#include <kernel/assembly/linear_functional.hpp>           // NEW: for LinearOperator
#include <kernel/assembly/bilinear_operator.hpp>           // NEW: for BilinearOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/error_computer.hpp>              // for L2/H1-error computation
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/common_functions.hpp>            // for SineBubbleFunction

// FEAST-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter

// FEAST-Solver includes
#include <kernel/solver/spai_precond.hpp>                  // for SPAIPrrecond
#include <kernel/solver/bicgstab.hpp>                      // for BiCGStab

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAST
using namespace FEAST;

// We're opening a new namespace for our tutorial.
namespace Tutorial03
{
  // Once again, we use quadrilaterals.
  typedef Shape::Quadrilateral ShapeType;
  // We want double precision.
  typedef double DataType;
  // Moreover, we use main memory (aka "RAM") for our containers.
  typedef Mem::Main MemType;

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // This tutorial covers an "anisotropic diffusion-convection-reaction" PDE, which is slightly more
  // complex than the poisson example, so we'll need to define a few custom classes for the assembly
  // before we can start with our actual main function.

  // The first thing we require is a class that will contain the data for our PDE, namely the
  // matrix A, the vector b and the scalar c:
  class AndicoreData
  {
  public:
    // First of all, we need a R^{2x2} matrix. The class for this job is the "Tiny" matrix:
    Tiny::Matrix<Real, 2, 2> a;

    // Note that 'Real' is a FEAST typedef, which equals 'double'.

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
    // And a simple initialisation constructor for our operator:
    explicit AndicoreOperator(const AndicoreData& data) :
      _data(data)
    {
    }

    // The first thing that is required from a bilinear operator are three so-called 'config tag'
    // classes, which tell the assembly routines which information from the trafo, the test- and the
    // trial-spaces are required for this operator.

    // We start off with the trafo configuration. Our operator does not require any information
    // from the trafo, so we can use the pre-defined 'base' config from the 'kernel/trafo/base.hpp'
    // header.
    typedef Trafo::ConfigBase TrafoConfig;

    // For the space configurations, it's not that easy. We have to write our own test- and
    // trial-space configurations, as our bilinear operator requires both test/trial function
    // values and gradients.
    // For this, we derive our 'TestConfig' class/struct from the 'Space::ConfigBase' class...
    struct TestConfig :
      public Space::ConfigBase
    {
      // We can now express our wishes.

      // Our bilinear operator requires test function values for the convection and rection terms:
      static constexpr bool need_value = true;

      // Moreover, we require the test function gradients for the diffusion term:
      static constexpr bool need_grad = true;

      // And that's it, all the other possibilities in 'kernel/trafo/base.hpp' are defaulted to false.
    };

    // Now we should do the same for the trial-space configuration 'TrialConfig'. However, as it
    // is identical to the test space configuration, we'll simply use a typedef here...
    typedef TestConfig TrialConfig;

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

      // Finally, we also need the DataType typedef.
      typedef typename AsmTraits_::DataType DataType;

      // As all other evaluator classes (e.g. the one from the analytic function in
      // 'tutorial_02_laplace'), also the bilinear operator evaluator has a 'prepare-finish'
      // function pair. However, as we do not require them, we may simply skip implementing them,
      // as the base-class we derived from already contains a 'do-nothing' implementation for those.

      // Finally, we can provide the actual evaluation operator, which performs the actual work:
      // This operator returns a scalar value and its parameters are const references to...
      DataType operator()(
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
    // In this case, we may require some data from the transformation, as the solution function
    // will require it. So the first task is to find out what the solution function needs from
    // the transformation. For this, however, we first have to tell the solution function what we
    // need from it -- see the 'ConfigTraits' class template of the 'PringlesFunction' class
    // of the 'tutorial_02_laplace' example code.

    // So we define a 'configuration' for our analytic solution function, which derives from the
    // basic configuration for analytic functions defined in 'kernel/trafo/base.hpp':
    struct AnalyticConfig :
      public Trafo::AnalyticConfigBase
    {
      // We need solution function values for the reaction term
      static constexpr bool need_value = true;
      // We need solution gradients for the convection term
      static constexpr bool need_grad = true;
      // We need solution hessians for the diffusion term
      static constexpr bool need_hess = true;
    };

    // The trafo requirements of this functional are identical to the ones requested by the
    // analytical solution function, so we simply make a typedef for its trafo configuration.
    typedef typename SolFunction_::template ConfigTraits<AnalyticConfig>::TrafoConfig TrafoConfig;

    // Now for the test-space configuration. For a linear functional, we usually require test
    // function values for the evaluation, so build a corresponding test-space config:
    struct TestConfig :
      public Space::ConfigBase
    {
      // Give us test function values
      static constexpr bool need_value = true;
    };

    // In analogy to operators, functionals are evaluated using 'evaluator' class templates as well,
    // so let's declare one of those and derive it from the base-class evaluator:
    template<typename AsmTraits_>
    class Evaluator :
      public Assembly::LinearFunctional::Evaluator<AsmTraits_>
    {
    protected:
      // As usual, a reference to our data object:
      const AndicoreData& _data;

      // The 'AsmTraits_' class template, which is passed to this evaluator template by the
      // assembler, is the same as the one passed on for bilinear operators, so we can again
      // query our required types from it. As functionals do not operate on trial basis functions,
      // we only require the trafo data and the test basis data types:
      typedef typename AsmTraits_::TrafoData TrafoData;
      typedef typename AsmTraits_::TestBasisData TestBasisData;

      // Moreover, we need another type that we have not used before: the 'trafo evaluator', which
      // is again an evaluator class similar to the one we are currently building up.
      typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;

      // Finally, we also need the DataType typedef.
      typedef typename AsmTraits_::DataType DataType;

      // Now comes the interesting part:
      // We want to evaluate the analytical solution, and for this, we need to get its corresponding
      // evaluator, which is again a class template similar to the one that you are currently reading,
      // see the 'PringlesFunction' from the 'tutorial_02_laplace' example.
      // However, its template parameter is *not* our assembly traits type 'AsmTraits_', but another
      // type called 'analytic evaluation traits', which can be queried from our assembly traits.
      // For convenience, we'll make a typedef of this traits class:
      typedef typename AsmTraits_::AnalyticEvalTraits AnalyticEvalTraits;

      // With that type, we can now define the type of the solution function evaluator:
      typedef typename SolFunction_::template Evaluator<AnalyticEvalTraits> SolEvaluator;

      // Now, declare one of those solution function evaluators:
      SolEvaluator _sol_eval;

      // Finally, declare a variable that will store the value of the force in the current cubature point:
      DataType _force_value;

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
      // for each cell of the mesh that the assembler iterates over.  Although uur functional
      // itself does not require this initialisation, the evaluator of our analytical solution
      // may do, so we have to manually call it here:
      void prepare(const TrafoEvaluator& trafo_eval)
      {
        // prepare the analytical solution evaluator:
        _sol_eval.prepare(trafo_eval);
      }

      // Also, the solution function evaluator may wish to finish itself...
      void finish()
      {
        _sol_eval.finish();
      }

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

        // At this point, we can already combine these terms to pre-compute the value of our force
        // functional, so that we just need to multiply by 'psi' lateron:
        _force_value =
          // First, the diffusion. This is the matrix-dot product of 'A' and the hessian of 'u':
          - dot(_data.a, _sol_eval.hessian(tau))
          // Next, the convection. This is the vector-dot-product of 'b' and the gradient of 'u':
          + dot(_data.b, _sol_eval.gradient(tau))
          // Finally, the reaction. This is a simple product:
          + _data.c * _sol_eval.value(tau);
      }

      // Once again, we implement our evaluation operator, which is quite similar to the one
      // for bilinear operators, with the only difference being that there is no trial function.
      // It takes a const reference to the test-function data object as a parameter
      // and returns the evaluation of the functional in the current cubature point:
      DataType operator()(const TestBasisData& psi) const
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

    // First of all, we create and initialiase the data object for our PDE:
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

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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

    // Allocate matrix
    MatrixType matrix;
    Assembly::SymbolicMatrixAssembler<>::assemble1(matrix, space);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Perform numerical matrix assembly

    // Create a cubature factory
    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    std::cout << "Assembling system matrix..." << std::endl;

    // First of all, format the matrix entries to zero.
    matrix.format();

    // At this point, we need to create an object representing our custom bilinear operator:
    AndicoreOperator andicore_operator(andicore_data);

    // Next, we call the bilinear operator assembler to assemble the operator into a matrix.
    Assembly::BilinearOperatorAssembler::assemble_matrix1( matrix, andicore_operator, space, cubature_factory);

    std::cout << "Assembling right-hand-side vector..." << std::endl;

    // Format the right-hand-side vector entries to zero.
    vec_rhs.format();

    // The next task is the assembly of the right-hand-side vector.
    // As our functional is parameterised in the solution function, we first need to create an
    // object representing that solution function. We choose the standard sine-bubble:
    typedef Assembly::Common::SineBubbleFunction SolFunction;
    SolFunction sol_function;

    // Now, we utilise our custom AndicoreFunctional to obtain a functional that combines
    // the desired analytical solution with our andicore PDE operator:
    AndicoreFunctional<Assembly::Common::SineBubbleFunction> force_functional(andicore_data, sol_function);

    // Now we can call the LinearFunctionalAssembler class to assemble our linear
    // functional into a vector.
    Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, force_functional, space, cubature_factory);

    // Finally, clear the initial solution vector.
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions..." << std::endl;

    // In this example, we assemble homogene Dirichlet boundary conditions:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // Add our only boundary component:
    unit_asm.add_mesh_part(boundary);

    // And assemble the filter
    FilterType filter;
    unit_asm.assemble(filter, space);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions..." << std::endl;

    // Apply the filter onto the system matrix...
    filter.filter_mat(matrix);

    // ...the right-hand-side vector...
    filter.filter_rhs(vec_rhs);

    // ...and the solution vector.
    filter.filter_sol(vec_sol);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    std::cout << "Solving linear system..." << std::endl;

    // Create a SPAI preconditioner
    auto precond = Solver::new_spai_precond(matrix, filter, matrix.layout());

    // Create a preconditioned BiCGStab solver
    auto solver = Solver::new_bicgstab(matrix, filter, precond);

    // Enable convergence plot
    solver->set_plot(true);

    // Initialise the solver
    solver->init();

    // Solve our linear system
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // And release the solver
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    std::cout << "Computing errors against reference solution..." << std::endl;

    // We have already created an object representing our solution function for the assembly of the
    // right-hand-side vector, so we may re-use it for the error computation now.

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
    String vtk_name(String("./tutorial-03-andicore-lvl") + stringify(level));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

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
  } // int main(...)
} // namespace Tutorial03

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialise the FEAST runtime environment:
  Runtime::initialise(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAST's tutorial #03: Andicore" << std::endl;

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
      std::cerr << "ERROR: Failed to parse '" << argv[argc-1] << "' as refinement level." << std::endl;
      std::cerr << "Note: The last argument must be a positive integer." << std::endl;
      // Abort our runtime environment
      Runtime::abort();
    }
    // parse successful
    level = Index(ilevel);
  }

  // call the tutorial's main function
  Tutorial03::main(level);

  // Finalise the runtime
  return Runtime::finalise();
}
