// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// \brief FEAT Tutorial 09: Adaptive Mesh
//
// This file contains a simple Poisson solver that demonstrates the usage of the
// AdaptiveMesh class for storing and producing partially refined meshes.
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
// The purpose of this tutorial is to demonstrate the usage of the AdaptiveMesh
// class and to show its integration with other parts of FEAT3.
// It mimics 'tutorial_01_poisson' and the program flow is largely the same.
//
// \author Markus Muegge
//

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We start our little tutorial with a batch of includes...

// Misc. FEAT includes
#include "kernel/geometry/template_sets.hpp"
#include "kernel/geometry/templates/schneiders_data.hpp"
#include <kernel/runtime.hpp>                                    // for Runtime
#include <kernel/util/string.hpp>                                // for String

// FEAT-Geometry includes
#include <kernel/geometry/boundary_factory.hpp>                  // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>                    // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>                  // for RefinedUnitCubeFactory
#include <kernel/geometry/export_vtk.hpp>                        // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                         // for MeshPart
#include <kernel/geometry/subdivision_levels.hpp>                // NEW: for SubdivisionLevels
#include <kernel/geometry/adaptive_mesh.hpp>                     // NEW: for AdaptiveMesh
#include <kernel/geometry/adaptive_mesh_layer.hpp>               // NEW: for AdaptiveMeshLayer

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>                     // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange1/element.hpp>                    // the Lagrange-1 Element (aka "Q1")
#include <kernel/space/lagrange2/element.hpp>                    // the Lagrange-2 Element (aka "Q2")
#include <kernel/space/lagrange3/element.hpp>                    // the Lagrange-3 Element (aka "Q3")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>                   // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>                            // for ExpBubbleFunction

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>                // for SymbolicAssembler
#include <kernel/assembly/domain_assembler.hpp>                  // for DomainAssembler
#include <kernel/assembly/domain_assembler_helpers.hpp>          // for Assembly::assemble_***
#include <kernel/assembly/common_operators.hpp>                  // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>                // for LaplaceFunctional
#include <kernel/assembly/unit_filter_assembler.hpp>             // for UnitFilterAssembler
#include <kernel/assembly/discrete_projector.hpp>                // for DiscreteVertexProjector

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                         // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>                    // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                          // for UnitFilter

// FEAT-Solver includes
#include <kernel/solver/ssor_precond.hpp>                        // for SSORPrecond
#include <kernel/solver/pcg.hpp>                                 // for PCG

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial09
{
  // Once again, we use quadrilaterals
  typedef Shape::Quadrilateral ShapeType;

  // The adaptive mesh requires a foundation mesh to build on.
  // We use the unstructured conformal mesh for this purpose.
  typedef Geometry::ConformalMesh<ShapeType> FoundationMeshType;

  // We also define a mesh part on the foundation mesh
  typedef Geometry::MeshPart<FoundationMeshType> FoundationMeshPartType;

  // The refinement of the adaptive mesh is guided by a template set. The
  // template set prescribes how to refine entities of the mesh while ensuring
  // the mesh stays conformal.
  typedef Geometry::SchneidersTemplates TemplateSetType;

  // Next is the focus of this tutorial, the AdaptiveMesh.
  // It is parameterized by the chosen template set and the mesh shape.
  // The AdaptiveMesh template ensures the chosen template set is compatible with the mesh shape.
  typedef Geometry::AdaptiveMesh<TemplateSetType, ShapeType> AdaptiveMeshType;

  // As we will see later, one of the differences between the ConformalMesh and
  // the AdaptiveMesh is that the AdaptiveMesh stores multiple mesh layers at the
  // same time. This is done to ensure the partial refinement can be easily
  // changed, for example in the case of moving geometry. For compatability with
  // other features the AdaptiveMeshLayer allows focusing on a single layer of
  // the AdaptiveMesh.
  typedef Geometry::AdaptiveMeshLayer<AdaptiveMeshType> AdaptiveMeshLayerType;

  // The AdaptiveMeshLayer provides the same interface as other FEAT mesh classes.
  // We can thus define a MeshPart for a layer of the adaptive mesh, ...
  typedef Geometry::MeshPart<AdaptiveMeshLayerType> MeshPartType;

  // as well as a TrafoType, ...
  typedef Trafo::Standard::Mapping<AdaptiveMeshLayerType> TrafoType;

  // and a finite element space.
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

  // Here's our tutorial's main function.
  void main(Index level)
  {
    // The AdaptiveMesh requires a mesh to partially (or fully) refine further.
    // We start by creating this foundation mesh.

    std::cout << "Creating foundation mesh... " << level << "...\n";

    // Create the foundation mesh factory.
    Geometry::RefinedUnitCubeFactory<FoundationMeshType> mesh_factory(level);

    // Create the foundation mesh using the factory.
    FoundationMeshType foundation_mesh(mesh_factory);


    // Let's create a boundary factory for our mesh.
    Geometry::BoundaryFactory<FoundationMeshType> boundary_factory(foundation_mesh);

    // And create the boundary mesh-part by using the factory.
    FoundationMeshPartType foundation_boundary(boundary_factory);

    std::cout << "Creating adaptive mesh...\n";

    // With the foundation mesh we can now create an adaptive mesh.
    // Note that no work is done in the constructor of the adaptive mesh.
    // All refinement is done in the AdaptiveMesh::adapt method call later.
    auto adaptive_mesh = std::make_shared<AdaptiveMeshType>(foundation_mesh);

    // The refinement of the adaptive mesh is guided by assigning subdivision
    // levels to the vertices of the foundation mesh. Each template set chooses
    // how to interpret the subdivision levels. For our chosen
    // SchneidersTemplates each subdivision level corresponds to decreasing the
    // mesh size by a factor of 3 around the vertex. The SubdivisionLevels class
    // serves as a collection of subdivision levels.
    Geometry::SubdivisionLevels sdls(foundation_mesh.get_num_vertices());

    // In this tutorial we assign subdivision levels based on distance to the domain center.
    FoundationMeshType::VertexSetType& vertex_set = foundation_mesh.get_vertex_set();

    FoundationMeshType::VertexType center {0.5, 0.5};
    for(Index i(0); i < vertex_set.get_num_vertices(); i++)
    {
      Real distance = (vertex_set[i] - center).norm_euclid();

      if(distance > 0.5)
      {
        sdls[i] = 0;
      }
      else if(distance < 0.25)
      {
        sdls[i] = 2;
      }
      else
      {
        sdls[i] = 1;
      }
    }

    // Now we can refine the adaptive mesh. We say we adapt it to the new
    // subdivision levels. We have to choose how the adaptive mesh treats
    // entities of the foundation mesh that do not get further refined. We have
    // chosen ImportBehaviour::All here, which causes all entities of the
    // foundation mesh to get added to the adaptive mesh. This makes sure the
    // domain stays fully covered by the mesh and makes working with the mesh
    // easier. We could have chosen ImportBehaviour::RequiredOnly, in which case
    // only entities that get further refined (and their lower dimensional
    // sub-entities) would get added to the adaptive mesh.
    adaptive_mesh->adapt(sdls, AdaptiveMeshType::ImportBehaviour::All);

    // At this point we have an adpative mesh with 3 layers. Layer 0 is a copy
    // of (parts of) the foundation mesh. Layers 1 and 2 are the further refined
    // layers specified by the subdivision levels.
    std::cout << "Created adapive mesh with " << adaptive_mesh->num_layers() << " layers\n";

    // For compatability with the rest of FEAT3, we can focus on a single layer
    // of the mesh using the Geometry::AdaptiveMeshLayer class.
    AdaptiveMeshLayerType layer(adaptive_mesh, Geometry::Layer {2});

    // As an aside. The AdaptiveMesh mirrors much of the ConformalMesh
    // interface with an additional layer parameter added. Geometry::Layer is a
    // newtype to keep callsites readable. We can, for example, retrieve the
    // amount of edges on layer 2 with:
    // adaptive_mesh.get_num_entities(Geometry::Layer {2}, 1).
    // Without the newtype-wrapper it would be difficult to determine which
    // parameter is the layer and which is the dimension.

    // From here on we can treat the layer like any othe mesh type.
    // We can, for example, export it.
    String vtk_name(String("./tutorial-09-poisson-lvl") + stringify(level) + "-adapt-2");
    Geometry::ExportVTK<AdaptiveMeshLayerType> exporter(layer);
    exporter.write(vtk_name);

    // If we, by whatever means, decide that we need only the left side of the
    // domain refined, we can create new subdivision levels for that
    // refinement...
    for(Index i(0); i < vertex_set.get_num_vertices(); i++)
    {
      Real distance = (vertex_set[i] - center).norm_euclid();

      if(vertex_set[i][0] > 0.5)
      {
        sdls[i] = 0;
        continue;
      }

      if(distance > 0.5)
      {
        sdls[i] = 0;
      }
      else if(distance < 0.25)
      {
        sdls[i] = 2;
      }
      else
      {
        sdls[i] = 1;
      }
    }

    // and re-adapt the adaptive mesh. The mesh will perform the minimal set of
    // actions to conform to the new subdivision levels. It will reuse as much of
    // the existing refinement as possible.
    adaptive_mesh->adapt(sdls, AdaptiveMeshType::ImportBehaviour::All);

    // To solve the Poisson problem on the mesh layer, we also need a boundary
    // mesh part. We can project mesh parts from the foundation mesh onto the
    // adaptive mesh. The projection will result in a mesh part covering the same
    // part of the domain, just using the refined mesh entities of the chosen
    // layer.
    MeshPartType adaptive_boundary =
      adaptive_mesh->project_meshpart<AdaptiveMeshLayerType>(Geometry::Layer {2}, foundation_boundary);

    // Of course, the mesh layer behaves just like a normal mesh, so we could
    // just have used the BoundaryFactory again. But the projection is useful for
    // MeshParts read in from other sources.


    // For completeness, if we don't want to use the adaptive mehs layer, we
    // can export a layer of the adaptive mesh to a normal ConformalMesh using:
    FoundationMeshType exported_mesh = adaptive_mesh->to_conformal_mesh(Geometry::Layer{2});

    // There is a convenience function we want to use the adaptive mesh as a "mesh generator".
    // AdaptiveMesh::create_refined_mesh takes a foundation mesh and any callable that assigns subdivision levels to the mesh's vertices and produces to corresponding ConformalMesh
    //
    FoundationMeshType created_mesh = AdaptiveMeshType::create_refined_mesh(foundation_mesh, [&](Index i) {
      Real distance = (vertex_set[i] - center).norm_euclid();

      if(vertex_set[i][0] > 0.5)
      {
        return 0;
      }

      if(distance > 0.5)
      {
        return 0;
      }
      else if(distance < 0.25)
      {
        return 2;
      }
      else
      {
        return 1;
      }
    });

    // exported_mesh and created_mesh are now two (independent) identical meshes.


    // Voila, we now have a mesh and a corresponding boundary mesh-part.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Trafo and Finite Element Space initialization

    std::cout << "Creating Trafo and Space...\n";

    // We have already defined the types of the transformation and finite element space,
    // so we just need to create the objects.

    TrafoType trafo(layer);
    SpaceType space(trafo);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly

    std::cout << "Allocating matrix and vectors...\n";

    // Allocate matrix and assemble its structure:
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);

    // Allocate vectors:
    VectorType vec_sol = matrix.create_vector_r();
    VectorType vec_rhs = matrix.create_vector_l();

    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)

    // Create a domain assembler on all mesh elements:
    Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
    domain_assembler.compile_all_elements();

    // Choose a cubature rule:
    String cubature_name = "auto-degree:5";

    std::cout << "Assembling system matrix...\n";

    // Create the pre-defined Laplace operator:
    Assembly::Common::LaplaceOperator laplace_operator;

    // And assemble that operator:
    Assembly::assemble_bilinear_operator_matrix_1(
      domain_assembler, matrix, laplace_operator, space, cubature_name);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: RHS (linear forms)

    std::cout << "Assembling right-hand-side vector...\n";

    // Create reference solution:
    Analytic::Common::ExpBubbleFunction<ShapeType::dimension> sol_function;

    // Define right-hand side functional:
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);

    // Assemble right hand side:
    Assembly::assemble_linear_functional_vector(
      domain_assembler, vec_rhs, force_functional, space, cubature_name);


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions...\n";

    // The next step is the assembly of the homogeneous Dirichlet boundary conditions.
    // For this task, we require a Unit-Filter assembler:
    Assembly::UnitFilterAssembler<AdaptiveMeshLayerType> unit_asm;

    // Add the boundary to the assembler:
    unit_asm.add_mesh_part(adaptive_boundary);

    // Now, we need to assemble a unit-filter representing homogeneous Dirichlet BCs.
    // This is done by calling the 'assemble' function:
    FilterType filter;
    unit_asm.assemble(filter, space);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions...\n";


    // Apply the filter onto the system matrix...
    filter.filter_mat(matrix);

    // ...the right-hand-side vector...
    filter.filter_rhs(vec_rhs);

    // ...and the solution vector.
    filter.filter_sol(vec_sol);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    // For this tutorial, we stick to a simple PCG-SSOR solver.

    std::cout << "Solving linear system...\n";

    // Create SSOR preconditioner
    auto precond = Solver::new_ssor_precond(PreferredBackend::generic, matrix, filter);

    // Create a PCG solver
    auto solver = Solver::new_pcg(matrix, filter, precond);

    // Enable convergence plot
    solver->set_plot_mode(Solver::PlotMode::iter);

    // Set the maximum number of iterations to 1000:
    solver->set_max_iter(1000);

    // Initialize the solver
    solver->init();

    // Solve linear system
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // Release the solver
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing H0/H1-Errors

    // We have a discrete solution now, so have to do something with it.
    // Since this is a simple benchmark problem, we know the analytical solution of our PDE, so
    // we may compute the H0- and H1-errors of our discrete solution against it.

    std::cout << "\n";
    std::cout << "Computing errors against reference solution...\n";

    // Compute error norms
    auto error_info = Assembly::integrate_error_function<1>(
      domain_assembler,      // the domain assembler
      sol_function,          // the analytic reference solution function
      vec_sol,               // the coefficient vector of the FEM solution
      space,                 // the finite element space
      cubature_name          // the cubature name used for integration
    );

    // Print error norms to the console
    std::cout << "Error Analysis:\n";
    std::cout << error_info.print_norms() << "\n";

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we also write the discrete solution to a VTK file, so that it can be admired in Paraview...

    // First of all, build the filename string
    vtk_name = (String("./tutorial-09-adaptivemesh-lvl") + stringify(level) + "adapt-left");

    std::cout << "\n";
    std::cout << "Writing VTK file '" << vtk_name << ".vtu'...\n";

    // Project solution and right-hand-side vectors
    VectorType vertex_sol;
    VectorType vertex_rhs;

    Assembly::DiscreteVertexProjector::project(vertex_sol, vec_sol, space);
    Assembly::DiscreteVertexProjector::project(vertex_rhs, vec_rhs, space);

    // Create a VTK exporter for our mesh. The mesh layer has changed since we
    // created the first exporter above, so we require a new one.
    Geometry::ExportVTK<AdaptiveMeshLayerType> result_exporter(layer);

    // add the vertex-projection of our solution and rhs vectors
    result_exporter.add_vertex_scalar("sol", vertex_sol.elements());
    result_exporter.add_vertex_scalar("rhs", vertex_rhs.elements());

    // finally, write the VTK file
    result_exporter.write(vtk_name);

    // That's all, folks.
    std::cout << "Finished!\n";
  } // void main(...)
} // namespace Tutorial09

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
  Tutorial09::main(level);

  // And finally, return the exit code 0 to indicate a successful run.
  return 0;
}
