// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// \brief FEAT Tutorial 05: Geometric Multigrid Solver (TM)
//
// This file contains a simple prototypical multigrid Poisson solver for the unit square domain.
//
// The purpose of this tutorial is to demonstrate how to set up a simple geometrical multigrid
// solver for a linear PDE. In this example, the PDE to be solved is exactly the same as
// in Tutorial 01.
//
// The basic program flow of this application is as follows:
//
// 1. Define a class that acts as a container for the data that is required on each level
//    on our multigrid mesh hierarchy.
//
// 2. Create a coarse mesh level and successively refine it to obtain a mesh hierarchy.
//
// 3. Assemble the system matrix and system filter for each level of the hierarchy.
//
// 4. Assemble the grid transfer operators for each level pair of the hierarchy.
//
// 5. Assemble the right-hand-side vector and set up the initial solution vector.
//
// 6. Set up the multigrid solver
//
// 7. Solve the system using multigrid.
//
// 8. Compute the L2- and H1-errors against the analytical reference solution.
//
// 9. Write the result to a VTK file for visual post-processing.
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
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/error_computer.hpp>              // for L2/H1-error computation
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for ForceFunctional
#include <kernel/assembly/grid_transfer.hpp>               // NEW: for GridTransfer

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/transfer.hpp>                       // NEW: for Transfer

// FEAT-Solver includes
#include <kernel/solver/pcg.hpp>                           // for PCG
#include <kernel/solver/richardson.hpp>                    // NEW: for Richardson
#include <kernel/solver/jacobi_precond.hpp>                // NEW: for JacobiPrecond
#include <kernel/solver/multigrid.hpp>                     // NEW: for MultiGrid

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial05
{
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

  // The next one is new: this is the class that is responsible for the grid transfer.
  // In contrast to the other LAFEM containers that we have typedefed before, the 'Transfer'
  // class template requires a matrix type instead of the usual mem-data-index type triplet.
  typedef LAFEM::Transfer<MatrixType> TransferType;

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // In contrast to the previous tutorials, we need to store various data on multiple levels and
  // not just one the finest one, because we want to use geometric multigrid as a solver.
  // In the most simple case, which is what we consider in this tutorial, we need to store the
  // following data for each level in our multigrid mesh hierarchy:
  //
  // 1) the mesh, trafo and FE space objects
  // 2) the system matrix and the system filter
  // 3) grid transfer matrices (except for the coarsest level)
  //
  // This simple Level class will act as a "container" that encapsulates these required objects
  // for each level. We declare all members as "public" and do not implement any member functions
  // (except for the constructor) to keep the program flow as simple as possible.
  class Level
  {
  public:
    // The mesh for this level
    MeshType mesh;
    // The trafo for this level
    TrafoType trafo;
    // The space for this level
    SpaceType space;

    // The system matrix for this level
    MatrixType matrix;
    // The system filter for this level
    FilterType filter;

    // The grid-transfer operators for this level
    TransferType transfer;

    // This constructor will create a mesh, trafo and space based on the mesh factory.
    // Note that the other member variables of this level (matrix, filter, etc.) are
    // not initialized here - this will be performed in the main function below.
    explicit Level(Geometry::Factory<MeshType>& mesh_factory) :
      mesh(mesh_factory),
      trafo(mesh),
      space(trafo)
    {
    }
  }; // class Level

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function
  void main(const Index level_max, const Index level_min)
  {
    // We need a container that will store all our Level objects.
    // For the sake of simplicity, we use a simple std::deque of shared pointers for this task:
    std::deque<std::shared_ptr<Level>> levels;

    // Now we need to create our level hierarchy. For this, we will use the
    // RefinedUnitCubeFactory to generate the coarse mesh on the desired level,
    // and then successively refine the mesh by using the StandardRefinery.

    // At this point we have to choose in which order we want to store our level hierarchy:
    // coarse-to-fine or fine-to-coarse?
    // The answer is: fine-to-coarse, i.e. the first entry in our level deque stores the
    // finest mesh level, whereas the last entry in our deque stores the coarsest mesh level.
    // This may seem unintuitive at first, but there is indeed a good reason to do so when
    // it comes to setting up a parallel multigrid on hierarchic partitions. This goes far
    // beyond the scope of this tutorial, so we merely stick to this convention for consistency,
    // as it is used by all the "real-world" parallel applications. So keep in mind:
    // the finest level is at the front of the deque, the coarsest level at its back.

    // First of all, let's create the coarse level
    {
      std::cout << "Creating coarse level " << level_min << "..." << std::endl;

      // Create the coarse mesh factory
      Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level_min);

      // Allocate a new level using the factory and append it to our level deque:
      levels.push_front(std::make_shared<Level>(mesh_factory));
    }

    // Now, let's refine up to the maximum level
    for(Index ilevel(level_min); ilevel < level_max; ++ilevel)
    {
      std::cout << "Refining to level " << (ilevel+1) << "..." << std::endl;

      // The StandardRefinery is another mesh factory class, which creates a new mesh by
      // applying the "standard-regular-refinement" algorithm onto a given input mesh.
      // Create a StandardRefinery object and pass the last mesh to its constructor.
      Geometry::StandardRefinery<MeshType> mesh_factory(levels.front()->mesh);

      // Allocate the next level using the mesh factory and append it to our deque.
      // Remember that we have to push the new level to the front of the deque to achieve
      // the desired "fine-to-coarse" ordering of the levels:
      levels.push_front(std::make_shared<Level>(mesh_factory));
    }

    // At this point, the basic setup (mesh, trafo and space creation) of all levels is complete.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    // In the next step, we will assemble the system matrices and filters for each level.

    // Create a cubature factory; this one can be used for all assembly steps on all levels.
    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    // Now loop over all levels:
    for(Index ilevel(level_min); ilevel <= level_max; ++ilevel)
    {
      std::cout << "Assembling Matrix and Filter for level " << ilevel << "..." << std::endl;

      // Get a reference to the corresponding level for easier member access
      Level& lvl = *levels.at(std::size_t(level_max - ilevel));

      // Assemble the Laplace matrix:
      Assembly::SymbolicAssembler::assemble_matrix_std1(lvl.matrix, lvl.space);
      Assembly::Common::LaplaceOperator laplace_operator;
      lvl.matrix.format();
      Assembly::BilinearOperatorAssembler::assemble_matrix1(lvl.matrix, laplace_operator, lvl.space, cubature_factory);

      // Create the boundary object for the boundary condition assembly
      Geometry::BoundaryFactory<MeshType> boundary_factory(lvl.mesh);
      MeshPartType boundary(boundary_factory);

      // Assemble the unit filter for homogeneous Dirichlet boundary conditions
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      unit_asm.add_mesh_part(boundary);
      unit_asm.assemble(lvl.filter, lvl.space);
    } // end of level loop

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    // Next, we need to assemble the grid transfer (prolongation and restriction) matrices for
    // each pair of consecutive levels.

    // Now loop over all level pairs:
    for(Index ilevel(level_min); ilevel < level_max; ++ilevel)
    {
      std::cout << "Assembling Grid Transfer for level " << ilevel << "..." << std::endl;

      // Get references to the corresponding coarse and fine levels:
      Level& lvl_fine   = *levels.at(std::size_t(level_max - ilevel - 1));
      Level& lvl_coarse = *levels.at(std::size_t(level_max - ilevel));

      // Now let's get the references to the internal prolongation and restriction matrices
      // of the transfer operator:
      MatrixType& mat_prol = lvl_fine.transfer.get_mat_prol();
      MatrixType& mat_rest = lvl_fine.transfer.get_mat_rest();

      // We need to assemble the prolongation matrix, which is used by the multigrid
      // solver to project correction vectors from the coarse mesh onto the current mesh.
      // The assembly of prolongation matrices is quite similar to the assembly of
      // operator matrices (like the Laplace matrix above), i.e. we first need to
      // assemble the matrix structure and then the matrix content.

      // Assemble the prolongation matrix structure:
      Assembly::SymbolicAssembler::assemble_matrix_2lvl(
        mat_prol,           // the prolongation matrix that is to be assembled
        lvl_fine.space,     // the fine-mesh space
        lvl_coarse.space    // the coarse-mesh space
      );

      // As always, format the matrix:
      mat_prol.format();

      // Assemble the contents of the prolongation matrix:
      Assembly::GridTransfer::assemble_prolongation_direct(
        mat_prol,           // the prolongation matrix that is to be assembled
        lvl_fine.space,     // the fine-mesh space
        lvl_coarse.space,   // the coarse-mesh space
        cubature_factory    // the cubature factory to be used for integration
      );

      // That's it for our prolongation matrix.

      // We also need the restriction matrix, which is used by multigrid to project
      // defect vectors from the fine mesh to the coarse mesh.
      // Fortunately, this task is easy, because the restriction matrix is
      // always identical to the transpose of the prolongation matrix:
      mat_rest = mat_prol.transpose();
    } // end of level loop

    // At this point, all levels of our level hierarchy are fully assembled.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    // We have assembled the system matrices and filter, but we still need to set up the
    // right-hand-side vector as well as an initial solution vector. Note that these vectors
    // are required only on the finest level, so there is no loop involved here.

    // Get a reference to the finest level at the front of the deque:
    Level& lvl_fine = *levels.front();

    // Create solution and rhs vectors for the finest level
    VectorType vec_sol = lvl_fine.matrix.create_vector_r();
    VectorType vec_rhs = lvl_fine.matrix.create_vector_r();

    // Format both vectors
    vec_sol.format();
    vec_rhs.format();

    std::cout << "Assembling right-hand-side vector..." << std::endl;

    // Now assemble the right-hand-side vector; this is virtually identical to tutorial 01.

    // Our desired reference solution function:
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;

    // Create a corresponding right-hand-side functional:
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);

    // And assemble the rhs vector:
    Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, force_functional, lvl_fine.space, cubature_factory);

    // As always, we have to apply the boundary condition filter onto our vectors:
    lvl_fine.filter.filter_sol(vec_sol);
    lvl_fine.filter.filter_rhs(vec_rhs);

    // Our system is fully assembled now.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Multigrid Solver set-up

    std::cout << "Setting up Multigrid solver..." << std::endl;

    // All matrices, filters and vectors have been assembled, so we can start to set up our geometric
    // multigrid solver. Naturally, this task is more complex than the setup of a simple PCG-SSOR
    // solver as in tutorial 01, so let's take a look at the "road map":

    // The solver we want to set up is a "standard" geometric multigrid solver, which
    // 1) uses a simple CG solver as the coarse grid solver
    // 2) uses 4 steps of damped Jacobi as a pre-smoother and post-smoother

    // The first object that we require is a "MultiGridHierarchy" object, which is responsible
    // for managing the level hierarchy used by a multigrid solver. Note that this hierarchy
    // object is not a solver/preconditioner -- we will create that one later on.

    // The MultiGridHierarchy class templates has three template parameters and its only
    // constructor parameter is the size of the level hierarchy, i.e. the total number of
    // levels that we want to create:
    auto multigrid_hierarchy = std::make_shared<Solver::MultiGridHierarchy<
      MatrixType,          // the system matrix type
      FilterType,          // the system filter type
      TransferType         // the transfer operator type
      >>( levels.size() ); // the number of levels in the hierarchy

    // Now we need to fill this empty hierarchy object with life, i.e. we have to attach
    // all our matrices and filters to it. Moreover, we also need to create the corresponding
    // coarse grid solver and smoother objects.

    // For each level above the coarse level, we have to create a smoother and attach it to the
    // hierarchy. So let' loop over all levels except for the coarse-most one in *descending* order:
    for(std::size_t i(0); (i+1) < levels.size(); ++i)
    {
      // Get a reference to the corresponding level
      Level& lvl = *levels.at(i);

      // Create a Jacobi preconditioner for the smoother
      auto jacobi = Solver::new_jacobi_precond(lvl.matrix, lvl.filter);

      // Create a Richardson solver for the smoother
      auto smoother = Solver::new_richardson(lvl.matrix, lvl.filter, DataType(0.8), jacobi);

      // Set both the minimum and maximum number of iterations to 4; this will ensure that
      // the smoother always performs exactly 4 iterations. Moreover, this will disable
      // the convergence control, so that no unnecessary defect norm computations are
      // performed.
      smoother->set_max_iter(4);
      smoother->set_min_iter(4);

      // Finally, attach our system matrix, system filter, grid transfer matrices and
      // our smoother to the multigrid hierarchy. This is done by calling the following
      // member function:
      multigrid_hierarchy->push_level(
        lvl.matrix,     // the system matrix for this level
        lvl.filter,     // the system filter for this level
        lvl.transfer,   // the transfer operator for this level
        smoother,       // the pre-smoother
        smoother,       // the post-smoother
        smoother        // the peak-smoother
      );

      // Note:
      // The 'peak-smoother' is the smoother that is applied in an "inner peak"
      // step within a F- or W-cycle. See the Multigrid page in the documentation
      // for more details about the different cycles and their smoother calls.
    }

    // Finally, we have to set up the coarse level:
    {
      // Get a reference to the coarsest level at the back of the deque:
      Level& lvl = *levels.back();

      // We use a simple (unpreconditioned) CG solver as a coarse-grid solver, so let's create one:
      auto coarse_solver = Solver::new_pcg(lvl.matrix, lvl.filter);

      // At this point, we could configure the coarse grid solver, i.e. set tolerances and maximum
      // allowed iterations, but we'll just leave it at its default configuration here.

      // Now we need to attach this solver as well as the system matrix and filter to
      // our multigrid hierarchy. This is done by calling the following member function:
      multigrid_hierarchy->push_level(
        lvl.matrix,       // the coarse-level system matrix
        lvl.filter,       // the coarse-level system filter
        coarse_solver     // the coarse-level solver
      );
    }

    // That's it for the coarse level.

    // Next, we need to create a multigrid preconditioner for our hierarchy. This task
    // is quite simple, as we can use one of the convenience functions to obtain a
    // MultiGrid solver object using the hierarchy we have just set up. At this point,
    // we can also choose which multigrid cycle we want to use. We have the choice
    // between V-, F- and W-cycle, but we stick to the simple V-cycle in this example:
    auto multigrid = Solver::new_multigrid(
      multigrid_hierarchy,          // the multigrid hierarchy object
      Solver::MultiGridCycle::V     // the desired multigrid cycle
    );

    // Alright, now we have a multigrid preconditioner object.

    // However, the 'multigrid' object we have set up so far is just a *single cycle*, i.e.
    // it is merely a preconditioner, so we still need to create a "real solver" around that.
    // As we want to have a "classical multigrid", we just need to use this multigrid object
    // as a preconditioner for a Richardson solver:
    auto solver = Solver::new_richardson(lvl_fine.matrix, lvl_fine.filter, DataType(1), multigrid);

    // Change the solver's plot name (that's what is written to the console) to "Multigrid".
    // This is just for the sake of aesthetics: Without this step, our solver would present
    // itself as "Richardson" ;)
    solver->set_plot_name("Multigrid");

    // And that's it - now we have a "real" multigrid solver ready to go!

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    std::cout << std::endl << "Solving..." << std::endl;

    // Finally, there is one important point to keep in mind when using a multigrid solver:
    // The "multigrid" object, which we have created above, does not initialize its smoothers
    // and coarse grid solver. We have to initialize the multigrid hierarchy by hand *before*
    // initialising the remaining solver tree by calling its init() function.
    // This may seem inconvenient, but this is the only way to go when building more complex
    // multigrid solver configurations like e.g. the ScaRC solver, so we have to live with that.

    // So, initialize the multigrid hierarchy first:
    multigrid_hierarchy->init();

    // Now we can initialize our solver:
    solver->init();

    // As always, enable the convergence plot:
    solver->set_plot_mode(Solver::PlotMode::iter);

    // Solve our linear system
    Solver::solve(*solver, vec_sol, vec_rhs, lvl_fine.matrix, lvl_fine.filter);

    // Release our solver
    solver->done();

    // Release our multigrid hierarchy; this also needs to be done separately *after*
    // the remaining solver tree was released in the call directly above.
    multigrid_hierarchy->done();

    // The remainder of this tutorial is virtually identical to the previous ones.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    std::cout << std::endl << "Computing errors against reference solution..." << std::endl;

    // Compute and print errors
    Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
      vec_sol, sol_function, lvl_fine.space, cubature_factory);

    std::cout << errors << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Build the filename string
    String vtk_name(String("./tutorial-05-multigrid-lvl") + stringify(level_max));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(lvl_fine.mesh);

    // add our solution and rhs vectors
    exporter.add_vertex_scalar("sol", vec_sol.elements());
    exporter.add_vertex_scalar("rhs", vec_rhs.elements());

    // finally, write the VTK file
    exporter.write(vtk_name);

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
  } // void main(...)
} // namespace Tutorial05

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialize our runtime environment
  Runtime::initialize(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #05: Geometric Multigrid" << std::endl;

  Index level_max(3), level_min(1);

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
  }

  // Perform level sanity checks
  if(level_max < level_min)
  {
    // minimum level greater than maximum level
    std::cerr << "ERROR: Minimum level " << level_min << " is greater than the maximum level " << level_max << std::endl;
    Runtime::abort();
  }
  if(level_max == level_min)
  {
    // minimum and maximum levels are equal
    std::cout << "WARNING: Minimum and maximum levels are equal" << std::endl;
  }

  // print selected levels
  std::cout << "Level-Min: " << level_min << std::endl;
  std::cout << "Level-Max: " << level_max << std::endl;

  // call the tutorial's main function
  Tutorial05::main(level_max, level_min);

  // Finalize our runtime environment
  return Runtime::finalize();
}
