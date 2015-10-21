//
// \brief FEAST Tutorial x1: Tutorial 03 equipped with a V-cycle multigrid solver
//
// This file contains a clone of tutorial 03 (aka "andicore"), enhanced with a simple V-cycle
// multigrid solver. Also included is the necessary data management to maintain a hierarchy of
// assembled operators and right-hand-sides.
//
// See tutorial-01 through 03 for the actual problem being solved.
//
// \author Peter Zajac
// \author Dominik Goeddeke
//

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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

// FEAST-Analytic includs
#include <kernel/analytic/common.hpp>                      // for SineBubbleFunction

// FEAST-Assembly includes
#include <kernel/assembly/linear_functional.hpp>           // for LinearOperator
#include <kernel/assembly/bilinear_operator.hpp>           // for BilinearOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/error_computer.hpp>              // for L2/H1-error computation
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/grid_transfer.hpp>               // for GridTransfer

// FEAST-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter

// FEAST-Solver includes
#include <kernel/solver/jacobi_precond.hpp>                // for JacobiPrecond
#include <kernel/solver/bicgstab.hpp>                      // for BiCGStab

// we need std::vector as a container for our level hierarchy
#include <vector>

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAST
using namespace FEAST;

// We're opening a new namespace for our tutorial.
namespace TutorialX1
{
  // Once again, we use quadrilaterals.
  typedef Shape::Quadrilateral ShapeType;
  // We want double precision.
  typedef double DataType;
  // Moreover, we use main memory (aka "RAM") for our containers.
  typedef Mem::Main MemType;

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // This tutorial covers an "anisotropic diffusion-convection-reaction" PDE, which is slightly more
  // complex than the Poisson example, so we'll need to define a few custom classes for the assembly
  // (aka problem definition) before we can start with our actual main function.

  // The first thing we require is a class that contains the data for our PDE, namely the diffusion,
  // reaction and convection coefficients that will eventually be used to assemble the matrix A, the
  // vector b and the scalar c for our diffusion-convection-reaction problem
  //    -grad(A*grad(u)) + dot(b,grad(u)) + c*u = f         in the domain [0,1]x[0,1]
  //                                         u = 0          on the boundary
  class AndicoreData
  {
  public:
    // First of all, we need a R^{2x2} matrix. The class for this job is the "Tiny" matrix:
    Tiny::Matrix<DataType, 2, 2> a;

    // Next, the R^2 vector b. For this one we have the "Tiny" vector:
    Tiny::Vector<DataType, 2> b;

    // And the scalar c:
    DataType c;
  }; // class AndicoreData

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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
      // This enum is just a trick to define compile-time constant values...
      enum
      {
        // We can now express our wishes.

        // Our bilinear operator requires test function values for the convection and rection terms:
        need_value = 1,

        // Moreover, we require the test function gradients for the diffusion term:
        need_grad = 1
      };
      // And that's it, all the other possibilities in 'kernel/trafo/base.hpp' are defaulted to zero.
    }; // struct TestConfig

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

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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
    // that we require image point coordinates, so let's state this by deriving our TrafoConfig
    // from the 'Trafo::ConfigBase' and then expressing our wishes:
    struct TrafoConfig :
      public Trafo::ConfigBase
    {
      // We (more precise: our solution function) need image point coordinates
      static constexpr bool need_img_point = true;
    };

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
      // type called 'analytic evaluation traits', which we need to provide here.
      // Luckily, there is another template class in te 'Analytic' namespace which takes care
      // of all the necessary typedefing. The only two template parameters are the data type which
      // we want to use internally and the function class itself:
      typedef Analytic::EvalTraits<DataType, SolFunction_> AnalyticEvalTraits;

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
      // for each cell of the mesh that the assembler iterates over. However, as this simple
      // functional does not require any type of additional initialisation, we do not need to
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
        // solution function evalutuar (see the PringlesFunction in Tutorial 02), so we first need
        // to create the objects where the functions can write their values to.
        // We do this by using the types defined in the AnalyticEvalTraits class which we have
        // defined above:
        typename AnalyticEvalTraits::ValueType value(DataType(0));   // function value
        typename AnalyticEvalTraits::GradientType grad(DataType(0)); // function gradient
        typename AnalyticEvalTraits::HessianType hess(DataType(0));  // function hessian

        // Now we call the solution function's evaluator to compute these values. For this, we
        // provide the value object as well as the image point from our trafo data as parameters:
        _sol_eval.value   (value, tau.img_point);
        _sol_eval.gradient(grad , tau.img_point);
        _sol_eval.hessian (hess , tau.img_point);

        // Finally, we can combine these terms to pre-compute the value of our force
        // functional, so that we just need to multiply by 'psi' lateron:
        _force_value =
          // First, the diffusion. This is the matrix-dot product of 'A' and the hessian of 'u':
          - dot(_data.a, hess)
          // Next, the convection. This is the vector-dot-product of 'b' and the gradient of 'u':
          + dot(_data.b, grad)
          // Finally, the reaction. This is a simple product:
          + _data.c * value;
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

  //
  // Now we need a class to keep track of multigrid levels.
  //
  // See the main() below on how this class is actually supposed to be used.
  //
  class SystemLevel
  {
  public:
    // Define the mesh type
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    // Define the boundary type
    typedef Geometry::MeshPart<MeshType> BoundaryType;
    // Define the trafo...
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    // and the space, as usual
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // Define the vector type
    typedef LAFEM::DenseVector<MemType, DataType> VectorType;
    // Define the matrix type
    typedef LAFEM::SparseMatrixCSR<MemType, DataType> MatrixType;
    // Define the filter type
    typedef LAFEM::UnitFilter<MemType, DataType> FilterType;

    // Define the smoother; this is the base-class for all solvers and preconditioners
    typedef Solver::SolverBase<VectorType> SmootherType;

    // Define everything that is actually needed:
    MeshType mesh;
    BoundaryType boundary;

    // The trafo-space pair
    TrafoType trafo;
    SpaceType space;

    // Our system matrix, four vectors and the boundary condition filter
    MatrixType mat_sys;
    VectorType vec_sol;
    VectorType vec_rhs;
    VectorType vec_def;
    VectorType vec_cor;
    FilterType filter;

    // Prolongation and restriction matrices
    MatrixType mat_prol;
    MatrixType mat_rest;

    // A pointer for our smoother or coarse mesh preconditioner
    std::shared_ptr<SmootherType> smoother;

  public:
    // The CTOR creates the mesh and the boundary from the factories, afterwards initialises
    // the trafo and the space and finally allocates the vectors and the filter.
    explicit SystemLevel(Geometry::Factory<MeshType>& mesh_factory, Geometry::Factory<BoundaryType>& boundary_factory) :
      mesh(mesh_factory),
      boundary(boundary_factory),
      trafo(mesh),
      space(trafo),
      vec_sol(space.get_num_dofs()),
      vec_rhs(space.get_num_dofs()),
      vec_def(space.get_num_dofs()),
      vec_cor(space.get_num_dofs()),
      filter(space.get_num_dofs()),
      smoother(nullptr)
    {
    }

    // as does the DTOR
    virtual ~SystemLevel()
    {
    }

    // now things get interesting. This aux. routine does all the magic to switch from
    // a coarser to a finer level, including the assembly of a transfer operator as a
    // rectangular matrix.
    SystemLevel* refine() const
    {
      Geometry::StandardRefinery<MeshType> mesh_refinery(mesh);
      Geometry::StandardRefinery<BoundaryType> boundary_refinery(boundary, mesh);
      SystemLevel* fine = new SystemLevel(mesh_refinery, boundary_refinery);
      fine->assemble_prolrest(*this);
      return fine;
    }

    void assemble_system(const AndicoreData& data)
    {
      // create matrix structure
      Assembly::SymbolicMatrixAssembler<>::assemble1(mat_sys, space);

      // Create a cubature factory
      Cubature::DynamicFactory cubature_factory("auto-degree:5");

      // First of all, format the matrix entries to zero.
      mat_sys.format();

      // At this point, we need to create an object representing our custom bilinear operator:
      AndicoreOperator andicore_operator(data);

      // Next, we call the bilinear operator assembler to assemble the operator into a matrix.
      Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_sys, andicore_operator, space, cubature_factory);

      // In this example, we assemble homogene Dirichlet boundary conditions:
      Assembly::UnitFilterAssembler<MeshType> unit_asm;

      // Add our only boundary component:
      unit_asm.add_mesh_part(boundary);

      // And assemble the filter
      unit_asm.assemble(filter, space);

      // apply filter onto matrix
      filter.filter_mat(mat_sys);
    }

    void assemble_rhs_sol(const AndicoreData& data)
    {
      vec_rhs.format();

      // The next task is the assembly of the right-hand-side vector.
      typedef Analytic::Common::SineBubbleFunction<2> SolFunction;
      SolFunction sol_function;

      // Create a cubature factory
      Cubature::DynamicFactory cubature_factory("auto-degree:5");

      // Now, we utilise our custom AndicoreFunctional to obtain a functional that combines
      // the desired analytical solution with our andicore PDE operator:
      AndicoreFunctional<SolFunction> force_functional(data, sol_function);

      // Now we can call the LinearFunctionalAssembler class to assemble our linear
      // functional into a vector.
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, force_functional, space, cubature_factory);

      // Finally, clear the initial solution vector.
      vec_sol.format();

      // apply filter
      filter.filter_rhs(vec_rhs);
      filter.filter_sol(vec_sol);
    }

    // This function assembles the prolongation and restriction matrices for the grid transfer between
    // the next coarser level and this level.
    void assemble_prolrest(const SystemLevel& coarse)
    {
      // Symbolic matrix assembly of prolongation matrix
      Assembly::SymbolicMatrixAssembler<Assembly::Stencil::StandardRefinement>::assemble(mat_prol, space, coarse.space);
      // Allocate a weight vector and format it to zero
      VectorType weight(space.get_num_dofs(), DataType(0));

      Cubature::DynamicFactory cubature_factory("auto-degree:5");

      // Now assemble the prolongation matrix and its corresponding weight vector
      Assembly::GridTransfer::assemble_prolongation(mat_prol, weight, space, coarse.space, cubature_factory);

      // Scale the prolongation matrix rows by the inverted weights
      weight.component_invert(weight);
      mat_prol.scale_rows(mat_prol, weight);

      // Finally, transpose to obtain the restriction matrix
      mat_rest.transpose(mat_prol);
    }

    // You guess: sets the smoother for this level
    void set_smoother(std::shared_ptr<SmootherType> new_smoother)
    {
      smoother = new_smoother;
    }

    // Returns the smoother pointer
    std::shared_ptr<SmootherType> get_smoother()
    {
      return smoother;
    }

    // This function computes the current defect on this level based on the rhs and sol vectors,
    // and returns the defect's norm.
    DataType compute_defect()
    {
      mat_sys.apply(vec_def, vec_sol, vec_rhs, -DataType(1));
      filter.filter_def(vec_def);
      return vec_def.norm2();
    }

    // This function applies the smoother for a desired number of Richardson steps.
    // Note: Upon entry, this function silently assumes that the defect vector has been calculated.
    // Note: Upon exit, the defect vector is up-to-date.
    void smooth(Index numsteps)
    {
      for(Index i(0); i < numsteps; ++i)
      {
        // apply smoother
        smoother->apply(vec_cor, vec_def);
        // update solution vector
        vec_sol.axpy(vec_cor, vec_sol);
        // compute new defect vector
        compute_defect();
      }
    }

    // Restricts the defect on this level onto the rhs of the coarse level.
    // Moreover, this function also formats the coarse sol vector to zero and copies the
    // rhs into the defect vector.
    void restrict_def(SystemLevel& coarse) const
    {
      // Restrict defect into coarse rhs
      mat_rest.apply(coarse.vec_rhs, vec_def);
      // Apply the defect filter (Yes, the defect and not the rhs filter!)
      coarse.filter.filter_def(coarse.vec_rhs);
      // Format sol vector
      coarse.vec_sol.format();
      // Copy rhs to defect
      coarse.vec_def.copy(coarse.vec_rhs);
    }

    // Prolongates the solution of the coarse level onto this level and updates this
    // level's solution vector.
    void prolongate_sol(const SystemLevel& coarse)
    {
      // Prolongate coarse sol into correction
      mat_prol.apply(vec_cor, coarse.vec_sol);
      // Apply correction filter
      filter.filter_cor(vec_cor);
      // Update our sol vector.
      vec_sol.axpy(vec_cor, vec_sol);
    }

    // This function computes the L2/H1-errors of the solution vector against the analytic solution.
    // Note: this makes only sense on the finest level.
    void compute_errors() const
    {
      Analytic::Common::SineBubbleFunction<2> sol_function;

      Cubature::DynamicFactory cubature_factory("auto-degree:5");

      // Compute and print the H0-/H1-errors
      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
        vec_sol, sol_function, space, cubature_factory);

      std::cout << errors << std::endl;
    }

    // Writes the rhs and solution vector into a VTK file.
    void write_vtk(String filename) const
    {
      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(mesh);

      // add the vertex-projection of our solution and rhs vectors
      exporter.add_scalar_vertex("sol", vec_sol.elements());
      exporter.add_scalar_vertex("rhs", vec_rhs.elements());

      // finally, write the VTK file
      exporter.write(filename);
    }
  };

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function
  void main(Index lvl_max, Index lvl_min)
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

    // our level vector
    std::vector<SystemLevel*> levels;

    // create coarse level (this is a bit tricky)
    {
      Geometry::RefinedUnitCubeFactory<SystemLevel::MeshType> unit_cube_factory(lvl_min);
      SystemLevel::MeshType unit_cube_mesh(unit_cube_factory);
      Geometry::BoundaryFactory<SystemLevel::MeshType> boundary_factory(unit_cube_mesh);
      levels.push_back(new SystemLevel(unit_cube_factory, boundary_factory));
    }

    // create other levels
    for(Index i(lvl_min); i < lvl_max; ++i)
      levels.push_back(levels.back()->refine());

    // assemble systems levels
    for(auto it(levels.begin()); it != levels.end(); ++it)
    {
      (*it)->assemble_system(andicore_data);

      // create a damped Jacobi smoother
      (*it)->set_smoother(Solver::new_jacobi_precond((*it)->mat_sys, (*it)->filter, DataType(0.7)));
    }

    // now let's create a coarse-grid solver
    // we use a BiCGStab solver and utilise the Jacobi smoother on the coarse grid as a preconditioner
    // Note: We need to create the coarse grid solver on the heap because it references data which
    // is owned by our level objects. These are going to be deleted at the end of this function and,
    // if we created the solver on the stack, it would have orphaned references afterwards.
    auto coarse_solver = Solver::new_bicgstab(levels.front()->mat_sys, levels.front()->filter, levels.front()->smoother);

    // configure the coarse grid solver
    coarse_solver->set_tol_rel(1E-5);
    coarse_solver->set_max_iter(50);

    // let's initialise the coarse grid solver
    coarse_solver->init();

    // assemble rhs/sol on finest level
    levels.back()->assemble_rhs_sol(andicore_data);

    // compute initial defect
    DataType def_init = levels.back()->compute_defect();
    std::cout << "Iteration   0 | Defect: " << stringify_fp_sci(def_init) << std::endl;

    // the main multigrid loop
    for(Index iter(1); iter <= 20; ++iter)
    {
      // restriction loop
      for(std::size_t j(levels.size()-1); j > 0; --j)
      {
        // apply pre-smoother
        levels.at(j)->smooth(4);
        // restrict defect
        levels.at(j)->restrict_def(*levels.at(j-1));
      }

      // apply the coarse grid solver and check its status
      Solver::Status cgs_status = coarse_solver->apply(levels.front()->vec_sol, levels.front()->vec_rhs);
      if(!Solver::status_success(cgs_status))
      {
        std::cerr << "ERROR: Coarse grid solver broke down!" << std::endl;
        break;
      }

      // prolongation loop
      for(std::size_t j(1); j < levels.size(); ++j)
      {
        // prolongate coarse grid solution
        levels.at(j)->prolongate_sol(*levels.at(j-1));
        // compute new defect (This can be skipped if the post-smoother is deactivated)
        levels.at(j)->compute_defect();
        // apply the post-smoother
        levels.at(j)->smooth(4);
      }

      // compute new defect
      DataType def = levels.back()->compute_defect();
      std::cout << "Iteration " << stringify(iter).pad_front(3) << " | Defect: " << stringify_fp_sci(def) << std::endl;

      // check for relative stopping criterion
      if(def < 1E-8 * def_init)
        break;
    }

    // release and delete our coarse-grid solver
    coarse_solver->done();
    coarse_solver.reset();

    // compute errors
    levels.back()->compute_errors();
    levels.back()->write_vtk("tutorial_x1_andicore_mg");

    // clean up
    while(!levels.empty())
    {
      delete levels.back();
      levels.pop_back();
    }
  }
} // namespace TutorialX1

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialise runtime
  Runtime::initialise(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAST's tutorial #X1: Andicore-MG" << std::endl;

  // The desired mesh refinement level.
  Index level_min(2), level_max(6);

  // First of all, let's see if we have command line parameters.
  if(argc > 1)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[1]).parse(ilevel) || (ilevel < 1))
    {
      // failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[1] << "' as maximum refinement level." << std::endl;
      std::cerr << "Note: The first argument must be a positive integer." << std::endl;
      Runtime::abort();
    }
    // parse successful
    level_max = Index(ilevel);
  }
  if(argc > 2)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[2]).parse(ilevel) || (ilevel < 1))
    {
      // failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[2] << "' as minimum refinement level." << std::endl;
      std::cerr << "Note: The second argument must be a positive integer." << std::endl;
      Runtime::abort();
    }
    // parse successful
    level_min = Index(ilevel);
    if(level_max < level_min)
    {
      // minimum level greater than maximum level
      std::cerr << "ERROR: Minimum level " << level_min << " is greather than the maximum level " << level_max << std::endl;
      Runtime::abort();
    }
  }

  // print selected levels
  std::cout << "Level-Min: " << level_min << std::endl;
  std::cout << "Level-Max: " << level_max << std::endl;

  // call the tutorial's main function
  TutorialX1::main(level_max, level_min);

  // Finalise runtime
  return Runtime::finalise();
}
