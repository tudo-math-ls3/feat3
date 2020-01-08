// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// \brief FEAT Tutorial 07: Stokes solver (TM)
//
// This file contains a simple Stokes solver for a parabolic Poiseuille-Flow on the
// unit-square domain.
//
// The PDE to be solved reads:
//
//    -Laplace(v) + Gradient(P) = 0         in the domain [0,1]x[0,1]
//        -Div(v)               = 0
//
// with pure Dirichlet boundary conditions for the velocity on the whole boundary
// as well as vanishing pressure integral mean.
//
// The analytical velocity field v = (v_1, v_2) is given as
//
//         v_1(x,y) = 4*y*(1 - y)
//         v_2(x,y) = 0
//
// The analytical pressure is given as
//
//           p(x,y) = 4 - 8*x
//
//
// The main purpose of this tutorial is to demonstrate the type definitions and usage of
// "meta-containers", i.e. "vectors of vectors" and "matrices of matrices". In contrast to
// some other finite element software packages, we do not assemble "block operators" of
// PDE systems into a single scalar CSR matrix, but instead build "meta-matrices", which
// allow us to work the combined PDE operator as a whole or by using its parts separately.
//
// The program flow of this application roughly coincides to the one from Tutorial 01 (as far
// as Poisson and Stokes can coincide, that is).
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
#include <kernel/space/lagrange2/element.hpp>              // the Lagrange-2 Element (aka "Q2")
#include <kernel/space/discontinuous/element.hpp>          // the Discontinuous Element (aka "P1dc")

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>                      // for ParProfileFunction

// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/mean_filter_assembler.hpp>       // for MeanFilterAssembler
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/gpdv_assembler.hpp>              // for GradPresDivVeloAssembler
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/velocity_analyser.hpp>           // NEW: for VelocityAnalyser

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/dense_vector_blocked.hpp>           // NEW: for DenseVectorBlocked
#include <kernel/lafem/tuple_vector.hpp>                   // NEW: for TupleVector
#include <kernel/lafem/sparse_matrix_bcsr.hpp>             // NEW: for SparseMatrixBCSR
#include <kernel/lafem/saddle_point_matrix.hpp>            // NEW: for SaddlePointMatrix
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/mean_filter.hpp>                    // NEW: for MeanFilter
#include <kernel/lafem/tuple_filter.hpp>                   // NEW: for TupleFilter

// FEAT-Solver includes
#include <kernel/solver/vanka.hpp>                         // NEW: for Vanka preconditioner
#include <kernel/solver/fgmres.hpp>                        // NEW: for FMGRES solver

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT, so use the namespace here.
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial07
{
  // Once again, we use quadrilaterals.
  typedef Shape::Quadrilateral ShapeType;
  // Use the unstructured conformal mesh class
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  // Define the corresponding mesh-part type
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  // Use the standard transformation mapping
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  // Use the Lagrange-2 element for velocity
  typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;

  // Use the P1dc element for pressure
  typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;

  // Finally, store the dimension of the shape as a constant here for convenience,
  // as we will require this several times in a moment:
  static constexpr int dim = ShapeType::dimension;

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Linear System type definitions

  // Our LAFEM containers work in main memory.
  typedef Mem::Main MemType;
  // Our data arrays should be double precision.
  typedef double DataType;
  // Use the default index type for indexing.
  typedef Index IndexType;

  // At this point, we want to define the vector type(s) that we want to use.
  // In contrast to all previous tutorials, which were scalar PDEs, we are now dealing with
  // a system of PDEs and therefore with more than one vector component. More precisely, our
  // desired solution vector is composed of a vector field for the velocity and a scalar
  // function for the pressure. In this case, we require several (different) vector types,
  // which we have to define here.

  // For vector fields, there is a specialised vector class available, namely the DenseVectorBlocked.
  // In addition to the usual Memory-Data-Index type triplet, which we have already used before,
  // we also have to specify the dimension of the vector field as the fourth template parameter.
  typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> VectorVeloType;

  // Next, we need the vector type for the pressure, which is a scalar function, so we can
  // use the already well-known DenseVector here:
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorPresType;

  // Finally, we need to "concatenate" the velocity and pressure vector types to one so-called
  // "meta-vector" type. A "meta-vector" is, generally speaking, a "vector of vectors".
  // There are two class templates which can be used to define such meta-vectors and the one
  // that we will be using here is the "TupleVector", which can concatenate two (or more)
  // vectors of an arbitrary vector type, as long as the Memory-Data-Index type triplets are
  // the same. We call this new type the "system vector type", indicating that this vector
  // type contains all the sub-vector of our full PDE system, which is the velocity-pressure
  // vector pair in this case.
  typedef LAFEM::TupleVector<VectorVeloType, VectorPresType> SystemVectorType;


  // Okay, we have all the required vector types. Now, we need to define the matrix types.

  // The Stokes operator consists of several terms, which correspond to different test-/trial-
  // space pair combinations, namely:

  //      Operator   Matrix  Trial      Test
  //   -Laplace(u) : A       velocity / velocity
  //       Grad(p) : B       pressure / velocity
  //        Div(u) : D       velocity / pressure

  // So in summary, we have the well-known "saddle-point" structure of the system:
  //
  //       [ A  B ] * [ v ] = [ f ]
  //       [ D  0 ]   [ p ]   [ 0 ]

  // As the velocity is a vector field and the pressure is a scalar, we will be dealing with
  // 'blocked' matrices instead of 'scalar' ones as in the previous tutorials.
  // The corresponding matrix type is the 'SparseMatrixBCSR' (note the 'B'), and it expects
  // the dimensions of the blocks as the two additional template parameters:
  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, dim> MatrixTypeA;
  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim,   1> MatrixTypeB;
  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType,   1, dim> MatrixTypeD;

  // Finally, we have to build the 'full matrix', which is a saddle-point matrix consisting
  // of the three above defined matrix blocks.
  typedef LAFEM::SaddlePointMatrix<MatrixTypeA, MatrixTypeB, MatrixTypeD> SystemMatrixType;

  // That's it for the matrix type. Note that although we are dealing with more complex
  // matrix and vector types here, these types still share some of the basic functionality
  // of the simple containers that we have been using before. For example, we can still
  // create compatible vectors by calling the matrix's "create_vector_l/r" functions or
  // format vectors by calling their "format" function.


  // As the last step, we have to define the boundary condition filters for both the
  // velocity and pressure spaces. This is quite similar to the definition of the vector
  // type above, i.e. we first choose the filters for the velocity and pressure spaces and
  // then "concatenate" them by using a "tuple" container.

  // For the velocity, we prescribe in-/out-/no-flow boundary conditions -- or in other
  // word: Dirichlet boundary conditions on the whole boundary. Now since the velocity is
  // a vector field and not a scalar function, we now also require a "blocked" version of
  // the unit-filter, which is the "UnitFilterBlocked":
  typedef LAFEM::UnitFilterBlocked<MemType, DataType, IndexType, dim> FilterVeloType;

  // As there is no Neumann boundary region for the velocity, the pressure is only defined
  // up to an additive constant (update your Stokes solution theory knowledge if this is
  // not clear to you), so we need to fix that. One elegant possibility to circumvent this
  // problem is to enforce that the pressure will have a mean value of zero. For this, we
  // employ a "MeanFilter", which is defined as follows:
  typedef LAFEM::MeanFilter<MemType, DataType, IndexType> FilterPresType;

  // Now in analogy to the vector types, we need to combine the velocity and pressure filters
  // by concatenating them with a "TupleFilter":
  typedef LAFEM::TupleFilter<FilterVeloType, FilterPresType> SystemFilterType;

  // Okay, that's it for the typedef orgy ;)

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function
  void main(Index level)
  {
    // As usual, create a mesh and the corresponding boundary mesh-part first.

    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;

    // Create the mesh
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    // Create the boundary mesh-part
    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Trafo and Finite Element Space initialisation

    std::cout << "Creating Trafo and Spaces..." << std::endl;

    // Create the trafo
    TrafoType trafo(mesh);

    // In the case of Stokes, there are two finite element spaces involved here:
    SpaceVeloType space_velo(trafo); // velocity space
    SpacePresType space_pres(trafo); // pressure space

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly
    std::cout << "Allocating matrix and vectors..." << std::endl;

    // As usual, we create first an empty matrix.
    SystemMatrixType matrix;

    // Our system matrix is a saddle-point matrix which consists of the three sub-matrices
    // named "A", "B" and "D". We can access all of these by calling the "block_a/b/d()"
    // functions, which give us a reference to the corresponding sub-matrix object.
    // As we need to assemble the sub-matrices individually, we begin with the assembly
    // of the sparsity pattern of the sub-matrix A, which is no different from any
    // of the previous tutorials:
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix.block_a(), space_velo);

    // Next, we have to assemble the sparsity pattern of the matrices B and D.
    // Note that "D" is virtually the same as "B^T", but we do not want to exploit this
    // here. In contrast to the sub-matrix "A", the sub-matrices "B" and "D" have different
    // test- and trial-spaces. For this kind of sparsity pattern, we call another assembly
    // function which expects two finite element space objects as the test (first) and
    // trial (second) space arguments:
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix.block_b(), space_velo, space_pres);
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix.block_d(), space_pres, space_velo);

    // That's it for the symbolic assembly of the matrix sparsity patterns.
    // In analogy to all previous tutorials, we can now use our matrix to create compatible
    // vectors for us:
    SystemVectorType vec_sol = matrix.create_vector_r();
    SystemVectorType vec_rhs = matrix.create_vector_l();

    // Finally, format the matrix and the vectors for the upcoming numerical assembly.
    matrix.format();
    vec_rhs.format();
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)

    // First of all, let's create a cubature factory:
    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    // Now we need to assemble the three sub-matrices individually.

    // Let us first assemble the sub-matrix A. The velocity diffusion operator of the Stokes
    // equations can be either the "full" deformation tensor (aka "Du:Dv") or its simplified
    // form, the gradient tensor (aka Laplace). We pick the gradient tensor here, but you
    // may as well use the deformation tensor here:
    Assembly::Common::LaplaceOperatorBlocked<dim> operator_a; // gradient tensor
    //Assembly::Common::DuDvOperatorBlocked<dim>  operator_a; // deformation tensor

    // In analogy to all previous (scalar) tutorials, we now use the BilinearOperatorAssembler
    // to do the dirty work for us. The only notable difference is that we have to call the
    // "assemble_block_matrix1" function instead of its scalar counterpart "assemble_matrix1":
    Assembly::BilinearOperatorAssembler::assemble_block_matrix1(
      matrix.block_a(), operator_a, space_velo, cubature_factory);

    // Next, we need to assemble the sub-matrices B and D for the pressure gradient and
    // velocity divergence. The good news is that this operator pair is required so often
    // that someone (me, actually) decided to write a custom assembler for this task.
    // So instead of setting up some operators, we simply use the "GradPresDivVeloAssembler"
    // for that:
    Assembly::GradPresDivVeloAssembler::assemble(
      matrix.block_b(),  // pressure gradient sub-matrix B
      matrix.block_d(),  // velocity divergence sub-matrix D
      space_velo,        // velocity space
      space_pres,        // pressure space
      cubature_factory   // cubature factory
    );

    // Note that the right-rand-side of our Stokes equations is zero, so we don't have to assemble
    // the contents of the right-rand-side vector either.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    std::cout << "Assembling boundary conditions..." << std::endl;

    // Now we have to assemble the boundary condition filters, which is also performed by
    // assembling the individual parts of the system filter. So, first of all, let us create
    // a new empty filter:
    SystemFilterType filter;

    // We have Dirichlet boundary conditions for the velocity, so we first require a unit-filter
    // assembler:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // As usual, add the whole boundary to the assembler:
    unit_asm.add_mesh_part(boundary);

    // Now we need to chose the boundary condition function. We want to have a parabolic profile
    // flux, so we can use the pre-defined "ParProfileVector" function, which gives us a parabolic
    // profile flux orthogonal to the line segment (x0,y0)--(x1,y1) with maximum magnitude 'v-max':
    Analytic::Common::ParProfileVector profile_function(
      0.0, // x0
      0.0, // y0
      0.0, // x1
      1.0, // y1
      1.0  // v-max
    );

    // We can now assemble the unit-filter for the velocity component of the system filter.
    // As the system filter is a TupleFilter, we can access its first (i.e. velocity) component
    // by calling the 'at' member function template and passing the desired index as a template
    // parameter:
    unit_asm.assemble(
      filter.at<0>(), // velocity unit-filter component of system filter
      space_velo,              // velocity space
      profile_function         // parabolic profile function
    );

    // That's it for the velocity filter, so let us continue with the pressure filter.
    // As mentioned above, we require a "MeanFilter" for the pressure, which can be assembled
    // by the MeanFilterAssembler. Note that there is no special setup required for this assembler,
    // so we can call the "assemble" function directly:
    Assembly::MeanFilterAssembler::assemble(
      filter.at<1>(),  // pressure mean-filter component of system filter
      space_pres,               // pressure space
      cubature_factory          // cubature factory
    );

    // Okay, now we have a mean-filter for the pressure and therefore the assembly of the
    // system filter is complete.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    std::cout << "Imposing boundary conditions..." << std::endl;

    // As usual, apply the filter onto the initial solution and right-hand-side vectors:
    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    // Moreover, we also want to filter the system matrix, just as we did in all previous
    // tutorials (except for tutorial 06). Unfortunately, this is not as easy as previously,
    // because the TupleFilter has no "filter_mat" function -- and this is for a reason.

    // The problem is that there is no general rule on whether and, if so, how filters have
    // to be applied on matrices. This gets even more complicated when you are dealing with
    // systems of PDEs with different (possibly chained) filters for the various components.

    // In this tutorial, we are dealing with a unit-filter/mean-filter pair for a Stokes system.
    // In this case, we have to apply the 'filter_mat' function of the velocity unit-filter
    // into the "A" sub-matrix (in analogy to the previous tutorials):
    filter.at<0>().filter_mat(matrix.block_a());

    // Moreover, we have also have to apply the unit-filter onto the sub-matrix "B", this
    // time calling the "filter_offdiag_row_mat" function to zero our all rows of B, which
    // correspond to test function DOFs affected by Dirichlet boundary conditions:
    filter.at<0>().filter_offdiag_row_mat(matrix.block_b());

    // The mean-filter component of the pressure does not provide a means to filter matrices,
    // thus the "D" sub-matrix remains unfiltered -- this is, by the way, the reason why the
    // commonly (i.e. in literature) used identity "D = B^T" is not really correct from a
    // practical point of view.

    // The Stokes is system is fully assembled now, so we can continue with the linear solver.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    std::cout << "Solving linear system..." << std::endl;

    // We have to set up a linear solver now. As we a dealing with a non-scalar system of the form
    // of a saddle-point system here, our choice of solvers and preconditioners is significantly
    // reduced, as we can only pick solvers which...
    // 1) work for any type of matrix structure (aka "black-box" solvers) or
    // 2) are special solvers for saddle-point systems.

    // The solvers that fall into category 1) are the simple Richardson iteration, any sort
    // of Multigrid as well as all Krylov subspace solvers (PCG, PCR, FMGRES, BiCGStab).

    // Category 2) currently boils down to the "Schur-complement" preconditioner, which splits the
    // saddle-point system into two separate systems, as well as the monolithic "Vanka", which is
    // a sophisticated Block-Gauss-Seidel/Jacobi method that is often used as a multigrid smoother.

    // In this tutorial, we want to keep it as simple as possible, so we stick to a FMGRES(64)
    // solver with a Vanka preconditioner. We could also try BiCGStab as an outer solver, but this
    // one would only work for very small levels, as the solver is simply too "weak" for this
    // type of problem :(

    // So, let us first create a Vanka preconditioner. The first two arguments are the saddle-point
    // system matrix and the system filter. The third mandatory argument is the desired variant of
    // the Vanka, which is represented by the "Solver::VankaType" enumeration. In total, there are
    // eight different Vanka variants available and (unfortunately) it requires some background
    // knowledge and experience to choose a suitable variant for a given problem. In our case, namely
    // a Stokes system with a Q2/P1dc discretisation, the "blocked full multiplicative" variant will
    // do just fine, as long as we apply some under-relaxation of, say, 0.8:
    auto vanka = Solver::new_vanka(
      matrix,                             // the system matrix
      filter,                             // the system filter
      Solver::VankaType::block_full_mult, // the desired Vanka variant
      0.8                                 // the relaxation parameter
    );

    // Next, we can create a FGMRES(64) solver and use the Vanka as a preconditioner.
    // In contrast to PCG and BiCGStab, the FGMRES solver requires a few additional parameters:
    auto solver = Solver::new_fgmres(
      matrix, // the system matrix
      filter, // the system filter
      64,     // the Krylov subspace dimension
      0.0,    // inner pseudo-residual scale factor
      vanka   // the preconditioner
    );

    // Set the plot mode of the solver:
    solver->set_plot_mode(Solver::PlotMode::all);

    // As FGMRES-Vanka is not a really powerful solver configuration, we have to increase the
    // maximum number of allowed iterations. 1000 should be sufficient for up to level 6:
    solver->set_max_iter(1000);

    // Also, reduce the required tolerance to 10^{-5} to speed things up a bit:
    solver->set_tol_rel(1E-5);

    // Okay, our solver is set up, so initialise it now:
    solver->init();

    // Solve our linear system:
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // And release our solver:
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Analyse velocity field

    // At this point, we may compute the L2/H1-errors of the velocity and pressure against
    // reference solutions as in the previous tutorials. However, we want to do something else
    // for a change here, so let us analyse the velocity field, i.e. compute various quantities
    // from our field: divergence, vorticity and its H0/H1-norms.

    std::cout << std::endl << "Performing velocity field analysis..." << std::endl;

    // The class that performs this analysis is the "VelocityAnalyser" and it returns a
    // "VelicityInfo" object of the appropriate datatype and dimension:
    Assembly::VelocityInfo<DataType, dim> velo_info = Assembly::VelocityAnalyser::compute(
      vec_sol.at<0>(), // the velocity field
      space_velo,               // the velocity space
      cubature_factory          // a cubature factory
    );

    // Just as the "ScalarErrorInfo", that is returned by the "ScalarErrorComputer", we can now
    // simply push the velo_info object to cout:
    std::cout << velo_info << std::endl;

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // First of all, build the filename string
    String vtk_name(String("./tutorial-07-stokes-lvl") + stringify(level));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // First, let's add the velocity field to our exporter. As the velocity is a vector field
    // and not a scalar function, we have to call the "add_vertex_vector" function instead of
    // "add_vertex_scalar". Note that we can directly pass the Q2 solution vector to the function
    // without projecting it explicitly to the vertices:
    exporter.add_vertex_vector("velocity", vec_sol.at<0>());

    // Our pressure is a scalar function defined in the P1dc space. Unfortunately, we cannot simply
    // pass the pressure component of the solution vector to the 'add_scalar_cell' function, but we
    // need to project it onto the cells before. This works quite similar to the "vertex projection"
    // that has been demonstrated in the previous tutorials:
    VectorPresType cell_pres;
    Assembly::DiscreteCellProjector::project(
      cell_pres,                // the vector that receives the cell-projection of the pressure
      vec_sol.at<1>(), // the pressure-vector to be projection
      space_pres,               // the pressure space
      cubature_factory          // a cubature factory for the projection.
    );

    // Now we can add the cell-projected pressure to our VTK exporter:
    exporter.add_cell_scalar("pressure", cell_pres.elements());

    // finally, write the VTK file
    exporter.write(vtk_name);

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
  } // int main(...)
} // namespace Tutorial07

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialise the FEAT runtime environment:
  Runtime::initialise(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #07: Stokes" << std::endl;

  // Specify the desired mesh refinement level, defaulted to 3.
  // Note that FEAT uses its own "Index" type rather than a wild mixture of int, uint, long
  // and such.
  Index level(3);

  // Now let's see if we have command line parameters: This tutorial supports passing
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
      // Abort our runtime environment
      Runtime::abort();
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
  Tutorial07::main(level);

  // And finally, finalise our runtime environment. This function returns the 'EXIT_SUCCESS' return code,
  // so we can simply return this as the result of our main function to indicate a successful run.
  return Runtime::finalise();
}
