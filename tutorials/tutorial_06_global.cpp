// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// \brief FEAT Tutorial 06: Global Linear Algebra
//
// This file contains a simple prototypical parallel PCG-Schwarz-ILU solver for the
// Poisson equation on the unit-square domain.
//
// The purpose of this tutorial is to demonstrate how to set up a global linear algebra system
// that is used for the solution of a simple scalar PDE on a partitioned domain in an MPI-based
// parallel environment. This tutorial does NOT demonstrate how to use a partitioner to obtain
// a domain decomposition of a domain, but merely focuses on the steps that have to be performed
// once a decomposition of the domain is available. In this tutorial, the decomposition itself
// is created explicitly by a helper class.
//
// In this example, the PDE to be solved is exactly the same as in Tutorial 01.
//
// This tutorial is only useful when compiled with MPI support and executed via mpirun/mpiexec
// with more than one process. Note that all the classes used in this tutorial also work in
// non-MPI builds, however, this tutorial will then boil down to a simple sequential PCG-ILU
// solver -- or in other words: the code will be mathematically equivalent to tutorial 01,
// with the exception of using ILU instead of SSOR as a preconditioner.
//
// The basic program flow of this application is as follows:
//
//  1. Define the basic spatial discretisation types, as usual.
//
//  2. Define the local and global linear algebra container classes.
//
//  3. Create a distributed communicator for inter-process communication.
//
//  4. Create a mesh-node that represents the patch (sub-domain) of this process
//     along with neighbourhood information required for communication.
//
//  5. Create a trafo and a finite element space for the patch of this process.
//
//  6. Assemble a gate that defines the global (parallel) linear algebra operations.
//
//  7. Create a global matrix and some global vectors.
//
//  8. Assemble the (global) system matrix and the (global) right-hand-side vector.
//
//  9. Assemble the (global) boundary condition filter.
//
// 10. Set up a PCG-Schwarz-ILU solver using a "type-1" Schwarz matrix.
//
// 11. Solve the linear system.
//
// 12. Compute the (global) L2- and H1-errors against the analytical reference solution.
//
// 13. Write the result to a partitioned VTU file set for visual post-processing.
//
// \author Peter Zajac
//

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We start our little tutorial with a batch of includes...

// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime
#include <kernel/util/dist.hpp>                            // NEW: for Dist::Comm

// FEAT-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart
#include <kernel/geometry/mesh_node.hpp>                   // NEW: for RootMeshNode, MeshNode
#include <kernel/geometry/unit_cube_patch_generator.hpp>   // NEW: for UnitCubePatchGenerator

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")

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
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for LaplaceFunctional
#include <kernel/assembly/mirror_assembler.hpp>            // NEW: for MirrorAssembler

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/vector_mirror.hpp>                  // NEW: for VectorMirror

// FEAT-Global includes
#include <kernel/global/gate.hpp>                          // NEW: for Global::Gate
#include <kernel/global/filter.hpp>                        // NEW: for Global::Filter
#include <kernel/global/matrix.hpp>                        // NEW: for Global::Matrix
#include <kernel/global/vector.hpp>                        // NEW: for Global::Vector

// FEAT-Solver includes
#include <kernel/solver/pcg.hpp>                           // for PCG
#include <kernel/solver/schwarz_precond.hpp>               // NEW: for SchwarzPrecond
#include <kernel/solver/ilu_precond.hpp>                   // NEW: for ILUPrecond

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

// We're opening a new namespace for our tutorial.
namespace Tutorial06
{
  // Note:
  // This tutorial works only for Hypercube meshes, as there exists no patch generator for
  // Simplex meshes (yet). If you want to switch from quadrilaterals (2D) to edges (1D) or
  // hexahedra (3D), then you need to manually adjust the communicator size sanity check at
  // the beginning of the 'main' function to check for powers of 2 in the 1D case or
  // powers of 8 in the 3D case in addition to changing the ShapeType definition below.

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


  // In addition to the usual types above, we also need to define another important mesh class
  // for this tutorial: the "mesh-node" type:
  // A mesh-node is a management object, which organises a mesh as well as a set of mesh-parts
  // that refer to the corresponding mesh along with some other data that we are not interested in
  // here. We did not use mesh-nodes in the previous tutorials, as there was only just one mesh
  // and a single mesh-part for the whole domain boundary, which we have managed by ourselves.

  // The class that we require here is a RootMeshNode and its only parameter is the mesh type:
  typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;


  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  // Linear System type definitions

  // Our LAFEM containers work in main memory.
  typedef Mem::Main MemType;
  // Our data arrays should be double precision.
  typedef double DataType;
  // Use the default index type for indexing.
  typedef Index IndexType;

  // Now we have arrived at the first main part of this tutorial: the definition of a global
  // (i.e. distributed) linear algebra system. This involves several additional steps
  // in comparison to the previous (purely sequential) tutorials.

  // As usual, we first define the basic matrix, vector and filter types.
  // We will name them "Local***Type" as there will also be global counterparts,
  // which use the local types as internal containers. In this context, the term "local"
  // refers to "local to this process".

  // Our local matrix type: a standard CSR matrix
  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> LocalMatrixType;

  // Our local vector type: the usual dense vector
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> LocalVectorType;

  // Our local filter type: the unit filter for Dirichlet boundary conditions
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> LocalFilterType;

  // The first new type that we require is a "vector mirror".
  // A mirror is one of the core components of a distributed linear algebra system
  // in the context of domain decomposition, which is used to gather and scatter
  // the data that is shared between several processes over their "halos".
  // Note that we will not use the mirrors directly once we have assembled them,
  // however, we still need to define the types for the definition of our global
  // containers, which will then use the mirrors internally.

  // The vector mirror takes the usual memory, data and index types as template parameters:
  typedef LAFEM::VectorMirror<MemType, DataType, IndexType> VectorMirrorType;

  // Now comes the second core component of a global linear algebra system: the "gate".
  // A gate is responsible for providing basic parallel synchronisation and communication
  // functionality, which is used by the global linear algebra containers that we will
  // define in a moment. We will not use any of the gate's functionality directly after
  // its initial creation, as we will simply work with global matrices and vectors, which
  // internally use the gate's functionality. Note that each gate object is associated
  // with a single finite element space object - so in a more complex application with
  // several different finite element spaces defined on possibly more than just a
  // single mesh (level), there would be one gate per FE space per mesh.

  // The gate class template is defined in the "Global" namespace and it takes two
  // template arguments: the local vector and the vector mirror types:
  typedef Global::Gate<LocalVectorType, VectorMirrorType> GateType;

  // Furthermore, we need to define the global (i.e. distributed) linear algebra types,
  // i.e. global matrix, vector and filter types.

  // First, we define the global matrix type.

  // The template arguments of the global matrix class template are the local matrix type,
  // which is used as an internal container, as well as two vector mirrors, which correspond
  // to the test- and trial-spaces of the matrix, which are identical in our case:
  typedef Global::Matrix<LocalMatrixType, VectorMirrorType, VectorMirrorType> GlobalMatrixType;

  // Next we need a compatible global vector. In analogy to the matrix, the template arguments
  // are the local vector type and the vector mirror type:
  typedef Global::Vector<LocalVectorType, VectorMirrorType> GlobalVectorType;

  // Finally, we also need a global filter, whose template arguments are the local filter
  // as well as the vector mirror type:
  typedef Global::Filter<LocalFilterType, VectorMirrorType> GlobalFilterType;

  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Here's our tutorial's main function
  void main(Index level)
  {
    // The very first thing that we need is a so-called "distributed communicator", which
    // represents a communication interface to the set of all processes that have been started
    // for this tutorial. We won't be using the communicator object for any inter-process
    // communication directly in this tutorial, but this communicator will be used "under-the-hood"
    // by all objects that require global communication.
    // The corresponding class is the Dist::Comm class, which is effectively just a wrapper around
    // a MPI_Comm object. We create a new comm object representing the "world" communicator:
    Dist::Comm comm = Dist::Comm::world();

    // The first difference in this tutorial is that we do not print console output directly to
    // std::cout, but use the communicator's "print" function instead, which ensures that only
    // one single process (in the communicator) prints the output to the console:
    comm.print("Welcome to FEAT's tutorial #06: Global Linear Algebra");

    // Note: The 'print' function automatically appends a line-break after each message,
    // so there's no "\n" at the end --  we only append or prepend line-breaks explicitly
    // if we want to insert empty lines in the output, like in the next message.

    // We also print the number of processes running this job for information:
    comm.print("\nNumber of processes: " + stringify(comm.size()));

    // At this point, we check the number of processes. As this is a simple tutorial without
    // any sophisticated partitioning algorithm, we are stuck to process counts which are
    // powers of 4 (in 2D). Also, we do not allow more than 64 processes, just to simplify
    // the condition of the following if-statement -- the tutorial code itself works even
    // for process counts of 256, 1024, etc.
    if((comm.size() != 1) && (comm.size() != 4) && (comm.size() != 16) && (comm.size() != 64))
    {
      // We pass std::cerr as the first parameter to the 'print' function here:
      comm.print(std::cerr, "ERROR: You must run this tutorial with 1, 4, 16 or 64 processes!");
      Runtime::abort();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Geometry initialisation

    // Now comes the first really interesting part:
    // We need a decomposition (partitioning) of our computational domain.
    //
    // In a "real world" application, we would read in a mesh from a file, apply some sort of
    // partitioning algorithm onto it, and split the domain into "patches" (aka "sub-domains").
    // However, as this is not a tutorial about partitioning, but about global linear algebra and
    // global linear solvers, we stick to an explicit partitioning of the unit-square domain here.
    //
    // The good news is:
    // 1. There exists a Geometry class which can generate such a decomposed domain along with
    //    all the required information that we need -- quite similar to the factories that we
    //    have been using in the previous tutorials.
    // 2. The objects that are returned by this generator, namely a so-called mesh-node, are of
    //    the very same type as the objects that a "real" partitioner would give us.
    //    In other words: If you would replace the following code for the creation of a mesh-node
    //    by a "real" partitioning approach, the remaining code of this tutorial would still
    //    stay the same!
    //
    // The bad news is:
    // 1. This generator class works only for hypercubes (Edges, Quadrilaterals and Hexahedra).
    // 2. The number of processes must be a power of 2 (1D), 4 (2D) or 8 (3D).
    //
    // Assume that we are running this tutorial with 4 or 16 processes, then the generator will
    // decompose the unit-square into 4/16 patches, where each patch contains exactly 1 element:
    //
    //                   bnd:1                                       bnd:1
    //           X###################X                       X###################X
    //           #         |         #                       # 12 | 13 | 14 | 15 #
    //           #    2    |    3    #                       #----+----+----+----#
    //           #         |         #                       #  8 |  9 | 10 | 11 #
    //     bnd:2 #---------+---------# bnd:3           bnd:2 #----+----+----+----# bnd:3
    //           #         |         #                       #  4 |  5 |  6 |  7 #
    //           #    0    |    1    #                       #----+----+----+----#
    //           #         |         #                       #  0 |  1 |  2 |  3 #
    //           X###################X                       X###################X
    //                   bnd:0                                       bnd:0
    //
    // The '#' represent the four boundary edges of our unit-square domain named 'bnd:0' to
    // 'bnd:3' and the 'X' represent the corner vertices of our domain, each of which is
    // part of the two adjacent boundary mesh-parts.
    // We will require these names for the assembly of boundary conditions later on.
    // Note: in 3D there would be six boundary faces named from 'bnd:0' to 'bnd:5'.

    // In contrast to the UnitCubeFactory and BoundaryFactory classes that we have used in the
    // previous tutorials to obtain our mesh, the generator that we want to use will create
    // a mesh-node for us. This mesh-node contains the actual mesh itself, which we will
    // reference later on, as well as the boundary mesh-parts and so-called halos, which we
    // will discuss in a minute.

    // Normally (i.e. in a "real world" application), mesh nodes are managed by a domain controller
    // class, which also takes care of allocating and deleting mesh nodes on the heap.
    // As we do this on foot in this tutorial, we will use a std::shared_ptr for convenience:
    std::shared_ptr<RootMeshNodeType> root_mesh_node;

    // The generator class will not only give us a mesh-node representing our patch of the domain,
    // but it will also tell us the ranks of all processes that manage our neighbour patches.
    // We will require these ranks to set up the communication "mirrors" that are required for
    // the global simulation. For this, we need to create a std::vector of ints, which will be
    // filled with our neighbour process ranks:
    std::vector<int> neighbour_ranks;

    // Now we can call our generator to obtain our patch mesh-node as well as our neighbour ranks.
    // Moreover, the create function returns an index that corresponds to the refinement level
    // of the global unit-square domain - we will require this for the further refinement below.
    Index lvl = Geometry::UnitCubePatchGenerator<MeshType>::create(
      comm.rank(),          // input:  the rank of this process
      comm.size(),          // input:  the total number of processes
      root_mesh_node,       // output: the root-mesh-node shared pointer
      neighbour_ranks);     // output: the neighbour ranks vector

    // At this point, we have our root mesh node as well as the neighbour ranks.
    // Just for fun, we let each process write its neighbour ranks to the console.
    // However, if we just write to std::cout, the output may be scrambled by a poor
    // MPI implementation. To circumvent this, we first build a String object that
    // contains the message to be printed on each process:
    String msg = "Neighbours of process " + stringify(comm.rank()) + ":";
    for(int i : neighbour_ranks)
      msg += " " + stringify(i);

    // Now we pass this message to the 'allprint' function of the communicator.
    // This will print the messages of all processes in rank-ascending order to std::cout,
    // where each output line is prefixed by the rank of the process that generated the
    // corresponding line:
    comm.allprint(msg);

    // Note that the 'allprint' function is somewhat "expensive" because it requires
    // communication between all processes and performs sequential execution;
    // so this function should be only used in exceptional cases such as debugging.
    // Also note that the usual 'print' function does not require communication,
    // so there is no problem with sequentialisation there.

    // As mentioned before, the generator returned a level index that represents the
    // refinement level of the mesh-node that it has generated. We want to print that out:
    comm.print("\nBase Mesh Level: " + stringify(lvl));

    // Now let's refine the mesh up to the level that was passed as a parameter to this function,
    // assuming that it is greater than the base-mesh level the generator gave us:
    if(lvl < level)
    {
      comm.print("Refining Mesh to Level " + stringify(level) + "...");

      for(; lvl < level; ++lvl)
      {
        root_mesh_node = std::shared_ptr<RootMeshNodeType>(root_mesh_node->refine());
      }
    }

    // Now we have a mesh node that represents a single patch of a partitioned unit-square
    // domain on the desired refinement level (or on a higher level) with all the information
    // that we require to set up a global linear algebra system and a corresponding solver.
    //
    // As already mentioned before:
    // Even if we had used some sort of "read-mesh-from-file-and-apply-partitioner" approach
    // instead of the "hand-made" mesh generator class, we would eventually reach the point
    // where the partitioner had given us a root-mesh-node representing our patch.
    //
    // In other words: from this point on, the following code is (almost) totally independent
    // of the way we obtained the root-mesh-node for our patch.

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Trafo and Space initialisation

    comm.print("Creating Trafo and Space...");

    // The creation of a transformation and a finite element space is identical to the
    // previous tutorials, we just need to get the mesh representing our patch of the domain
    // from the mesh-node:
    MeshType& mesh = *root_mesh_node->get_mesh();

    // Create the trafo
    TrafoType trafo(mesh);

    // Create the desired finite element space
    SpaceType space(trafo);

    // Note that both the transformation and the finite element space are defined merely on
    // the local patch, i.e. the trafo and the space do not know about the fact that we are
    // running a parallel simulation. All global interaction and communication is performed
    // only on the linear algebra level, which are going to set up in the next few steps.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Before we can start assembling matrices or vectors, we first need to assemble the gate.

    comm.print("Assembling Gate...");

    // Create the gate object and pass our communicator to it:
    GateType gate(comm);

    // The generator that we have used to obtain a mesh-node for our patch has also given us
    // a vector of the process ranks, which contains the neighbour patches of this process.
    // We need to create a mirror for each of our neighbours, which will then be used for
    // the communication with that neighbour - so let's loop over all neighbour ranks:
    for(auto it = neighbour_ranks.begin(); it != neighbour_ranks.end(); ++it)
    {
      // Get the rank of our neighbour process:
      const int neighbour_rank = *it;

      // As already mentioned before, our mesh-node does not only contain the mesh itself,
      // but also other important stuff, such as the "halos", which describe the overlap of
      // neighboured patches. For the assembly of the mirror, we need to get the halo mesh-part
      // from the mesh-node first:
      const MeshPartType* neighbour_halo = root_mesh_node->get_halo(neighbour_rank);

      // Ensure that we have a halo for this neighbour rank:
      XASSERTM(neighbour_halo != nullptr, "Failed to retrieve neighbour halo!");

      // Now that we have the halo mesh-part, we can create and assemble the corresponding mirror:
      VectorMirrorType neighbour_mirror;

      // Call the MirrorAssembler to do the dirty work for us:
      Assembly::MirrorAssembler::assemble_mirror(
        neighbour_mirror,   // the mirror that is to be assembled
        space,              // the FE space for which we want to assemble the mirror
        *neighbour_halo     // the halo mesh-part that the mirror is to be assembled on
      );

      // Once the mirror is assembled, we give it over to our gate.
      gate.push(
        neighbour_rank,               // the process rank of the neighbour
        std::move(neighbour_mirror)   // the mirror for the neighbour
      );

      // continue with next neighbour rank
    }

    // At this point, our gate contains all the mirrors required for the communication with our
    // neighbours. However, we are not done yet, as the gate needs to perform further internal
    // initialisation. This initialisation is performed by calling the 'compile' function,
    // which requires a local vector for some internal computations. Unfortunately, we do not
    // have a (local) matrix yet, whose create_vector_l/r function could create a vector for us,
    // so we need to create a local vector on foot this time. For DenseVector objects, we can
    // create a vector of the correct size by passing the number of (local) DOFs in the FE space
    // to its constructor:
    gate.compile(LocalVectorType(space.get_num_dofs()));

    // Alright, our gate is ready to go!

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly
    comm.print("Allocating matrix and vectors...");

    // Now that we have a gate, we can finally create our (global and local) matrices, vectors
    // and filters. As we want to set up a parallel solver, we create global matrices and vectors.

    // Before continue, we first need to discuss some details of the concept that the global
    // linear algebra implementation in FEAT is based on.
    // The global matrix, vector and filter containers, whose types we have already defined above,
    // are merely "wrappers" that combine the corresponding internal local (i.e. LAFEM) containers
    // with the global gate(s), which we have just set up. The combination of a local container
    // along with the global gate offers us the possibility of applying global (i.e. parallel/
    // distributed) operations such as matrix-vector-products, vector axpys or norm computations.
    // Each of these global operations consists of a local operation possibly combined with
    // parallel communication. The communication is taken care of by the gate, so the only thing
    // that we still need to provide is the assembly of the corresponding local containers.
    // So in summary, to assemble a "global" matrix we only need to assemble its internal
    // "local" matrix that is defined on the patch of this process -- the gate will take care
    // of everything else. For vectors and filters, the process is similar, i.e. we also need
    // to assemble the corresponding local vectors and filters, but depending on the vector or
    // filter type, some manual synchronisation may be required, as we will see further below.

    // Create a global matrix: the first two mandatory arguments are the gates for the test- and
    // trial spaces, all further arguments (in our case: none) are passed through to the internal
    // local matrix object constructor:
    GlobalMatrixType matrix(&gate, &gate);

    // Each of the global objects internally contains the corresponding local object
    // (in our case: a CSR matrix),  which we can obtain by calling the 'local' member function:
    LocalMatrixType& matrix_local = matrix.local();

    // Now let us assemble the local matrix structure, just as in any previous tutorial:
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_local, space);

    // As already mentioned, the global matrix object is mainly a wrapper around the local
    // matrix object, whose structure we have just assembled. In consequence, we have thus
    // indirectly defined the matrix structure of the global matrix object!

    // Just as usual, we can now use the global matrix to create compatible global vectors by
    // calling the create_vector_r/l function:
    GlobalVectorType vec_sol = matrix.create_vector_r();
    GlobalVectorType vec_rhs = matrix.create_vector_l();

    // In analogy to the global matrix, each global vector contains an internal local vector object
    // (in our case: a dense vector) that we can reference by calling the 'local()' member function:
    LocalVectorType& vec_sol_local = vec_sol.local();
    LocalVectorType& vec_rhs_local = vec_rhs.local();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)

    comm.print("Assembling system matrix...");

    // Create a cubature factory:
    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    // The assembly of a global (distributed) system matrix is quite easy:
    // All that we have to do is to assemble the internal local matrix, which is done just the
    // same way as in the non-parallel case that we have been treating in all previous tutorials!

    // Format the local matrix:
    matrix_local.format();

    // Create the -Laplace(u) operator:
    Assembly::Common::LaplaceOperator laplace_operator;

    // And assemble the numerical matrix content:
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_local, laplace_operator, space, cubature_factory);

    // This is already it: we do not need to do anything else for the global matrix assembly.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: RHS (linear forms)

    comm.print("Assembling right-hand-side vector...");

    // The assembly of the right-hand-side vector is similar to the assembly of the system matrix:
    // We first need to assemble the internal local vector just as in the non-parallel case.
    // However, there is one additional step to be done here right after the local assembly...

    // Initialise the right-hand-side vector entries to zero.
    vec_rhs_local.format();

    // Again, we use the sine-bubble as a reference solution:
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;

    // Create a force functional for our solution:
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);

    // And assemble the local rhs vector:
    Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs_local, force_functional, space, cubature_factory);

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // Now comes the small but significant difference:
    // In contrast to global matrices, global rhs vectors *must* be synchronised manually after the
    // local vector assembly!!! This task is easy but also easily overlooked, as we just have to
    // call the 'sync_0' member function of the global vector object:
    vec_rhs.sync_0();

    // But Beware:
    // If you forget to synchronise the global RHS vector, you will not experience any crashes
    // or solver breakdowns, which would indicate that something is wrong, but you will get
    // wrong results -- and by "wrong" I do not mean "inaccurate" -- I mean *GARBAGE* !

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Solution initialisation.

    // Clear the initial solution vector:
    vec_sol.format();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    comm.print("Assembling boundary conditions...");

    // The assembly of global filters highly depends of the type of the filter.
    // In the case of the unit-filter, which we use for Dirichlet boundary conditions,
    // we only have to assemble the local unit-filter as in the non-parallel case.

    // As usual, we create an instance of the global filter:
    GlobalFilterType filter;

    // And then fetch the internal local filter object:
    LocalFilterType& filter_local = filter.local();

    // Now, we have to assemble the filter.
    // In the case of the unit-filter, which is responsible for Dirichlet boundary conditions,
    // the assembly is quite simple: We just have to assemble the corresponding local filter.
    // So we need a unit-filter assembler:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // In the previous tutorials, we have created a single mesh-part representing the domain
    // boundary using the BoundaryFactory. We could do just the same here, but that would be
    // a waste of compute-time, because the patch generator, which has created the mesh-node
    // for us, also created boundary mesh-parts, which are already contained in the mesh-node.
    // The important difference is that the generator created four mesh-parts, named from
    // "bnd:0" to "bnd:3", representing the four boundary edges of the unit-square domain.
    // Note that if you have changed the ShapeType to hexahedra, the corresponding domain
    // has six boundary faces named "bnd:0" to "bnd:5".

    // To keep this part of the code dimension-independent, we loop up to 2*dimension:
    for(int ibnd(0); ibnd < 2*ShapeType::dimension; ++ibnd)
    {
      // Build the name of the boundary mesh-part and call the mesh-node's "find_mesh_part"
      // function to get a pointer to the mesh-part:
      MeshPartType* bnd_mesh_part = root_mesh_node->find_mesh_part("bnd:" + stringify(ibnd));

      // Do you still remember the domain-decomposition ASCII-art from the beginning of this
      // tutorial? If not, then scroll back up to refreshen your memory and return back here.

      // In the domain-decomposition world, there is one important point to keep in mind:
      // Although our (full) domain is adjacent to all four boundary edges, a single patch
      // may be adjacent only to a sub-set of these boundaries. In fact, if the patch owned
      // by this process is an inner patch (e.g. patch #5 in the right ASCII art above),
      // then this patch is not adjacent to any of the domain's boundaries.

      // If the pointer returned by the "find_mesh_part" function is a nullptr, then either
      // the (full) domain does not contain a mesh-part with that name or the (full) domain
      // contains such a mesh-part, but the patch of this process is not adjacent to it.
      // The first case would be an error (which we do not check here), but the latter case
      // is perfectly valid. Fortunately, this is not a problem for our unit-filter assembly,
      // as we just need to skip this boundary part if our process' patch is not adjacent to it:
      if(bnd_mesh_part != nullptr)
      {
        // Okay, the patch managed by this process is adjacent to this boundary part,
        // so let's add it to our unit-filter assembler:
        unit_asm.add_mesh_part(*bnd_mesh_part);
      }
    }

    // At this point, we have added all boundary mesh-part of our patch to the assembler.
    // Note: If our patch is an inner patch, then we have added no boundary mesh-parts to
    // the assembler, but that's fine as the unit-filter assembler can work with that, as
    // in this case the unit-filter will be an empty filter.

    // Finally, assemble the local filter:
    unit_asm.assemble(filter_local, space);

    // That's it for the unit-filter, i.e. we do not have to do anything else for the global
    // unit-filter.

    // Note: For other filter types (e.g. the mean-filter or the slip-filter), the assembly of the
    // global filter may require more/other steps than just the assembly of the local filter,
    // but we won't go into further details in this tutorial.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    comm.print("Imposing boundary conditions...");

    // Now we impose the boundary conditions into our linear system by applying the global filter
    // onto the global right-hand-side and initial solution vectors:
    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    // In contrast to our previous tutorials, we do *not* filter the global system matrix.
    // The reason is that we require an unfiltered matrix in the next step.
    // (By the way: the global filter object has no 'filter_mat' function, so you would get a
    // compiler error if you tried anyway...)

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Schwarz matrix creation

    comm.print("Setting up Schwarz preconditioner matrix...");

    // In this tutorial, we want to set up an additive overlapping Schwarz preconditioner and for
    // this, we need the so-called "local type-1 matrix". If you don't know what that means, you
    // should consider some literature on Schwarz/ScaRC preconditioners in the FEM context.

    // The easiest way to obtain a local type-1 matrix is to call the "convert_to_1" function of
    // the *unfiltered* global matrix, which will return a new local matrix object, which we will
    // call the "Schwarz" matrix from now on:
    LocalMatrixType schwarz_matrix = matrix.convert_to_1();

    // At this point, we have to apply our filter onto the Schwarz matrix. Note that the Schwarz
    // matrix is a local (and not a global) matrix object, so we need to use the local filter here:
    filter_local.filter_mat(schwarz_matrix);

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Solver set-up

    comm.print("\nSolving linear system...");

    // So, we have been setting up the global linear system for quite a while now and the good news
    // is that setting up a global (i.e. parallel) linear solver is quite trivial now!
    // We just need to decide what type of mixed local/global solver we want and then we can link
    // the individual nodes of the solver tree just as in the simple sequential tutorials.

    // The solver we want to set up in this tutorial is a
    // 1. global PCG preconditioned by
    // 2. an additive overlapping Schwarz preconditioner using
    // 3. an ILU(0) preconditioner as a local "block-solver"

    // As usual, we have to set up the solver in a bottom-up manner, so let's start with
    // with the local ILU preconditioner, which is applied onto the local type-1 matrix
    // (aka Schwarz matrix) and using the local unit-filter:
    auto precond = Solver::new_ilu_precond(schwarz_matrix, filter_local);

    // This ILU is used as a local solver for the (global) Schwarz preconditioner,
    // which is using the global filter:
    auto schwarz = Solver::new_schwarz_precond(precond, filter);

    // And finally, we use the Schwarz preconditioner in our (global) PCG solver:
    auto solver = Solver::new_pcg(matrix, filter, schwarz);

    // Enable the convergence plot.
    solver->set_plot_mode(Solver::PlotMode::iter);

    // The rest of the solver step is identical to the previous tutorials:

    // Initialise the solver.
    solver->init();

    // Solve our linear system:
    Solver::solve(*solver, vec_sol, vec_rhs, matrix, filter);

    // Release the solver:
    solver->done();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Computing L2/H1-Errors

    comm.print("\nComputing errors against reference solution...");

    // As usual, we want to compute the L2- and H1-errors against our reference solution.
    // For this, we first compute the errors on our patch by analysing our local solution vector:
    Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
      vec_sol_local, sol_function, space, cubature_factory);

    // And then we need to synchronise the errors over our communicator to sum up the errors of
    // each patch to obtain the errors over the whole domain:
    errors.synchronise(comm);

    // And let's print the errors to the console; we need to use the "format_string" function here,
    // as the "print" function accepts only String objects as input:
    comm.print(errors.format_string());

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we want to export our solution to a (P)VTU file set.

    // In a parallel simulation, each process will write a separate VTU file, which contains the
    // data that is defined on the patch of the corresponding process. Moreover, one process
    // writes a single additional PVTU file, which can be read by ParaView to visualise
    // the whole domain that consists of all patches.

    // Build the VTK filename; we also append the number of processes to the filename:
    String vtk_name = String("./tutorial-06-global-lvl") + stringify(level) + "-n" + stringify(comm.size());

    comm.print("Writing VTK file '" + vtk_name + ".pvtu'...");

    // Create a VTK exporter for our patch mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // Add the vertex-projection of our (local) solution and rhs vectors
    exporter.add_vertex_scalar("sol", vec_sol_local.elements());
    exporter.add_vertex_scalar("rhs", vec_rhs_local.elements());

    // Finally, write the VTK files by calling the "write" function of the exporter and pass the
    // communicator as a second argument:
    exporter.write(vtk_name, comm);

    // Note: Do not forget the 'comm' argument in the call above as otherwise each process will
    // try to write to the same VTK file, resulting in garbage due to race conditions...

    // That's all, folks.
    comm.print("Finished!");
  } // void main(...)
} // namespace Tutorial06

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialise the FEAT runtime environment:
  Runtime::initialise(argc, argv);

  // Specify the desired mesh refinement level, defaulted to 5.
  Index level(5);

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
  }

  // call the tutorial's main function
  Tutorial06::main(level);

  // And finally, finalise our runtime environment.
  return Runtime::finalise();
}
