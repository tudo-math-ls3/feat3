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
#include <kernel/space/lagrange2/element.hpp>              // the Lagrange-2 Element (aka "Q2")
#include <kernel/space/lagrange3/element.hpp>              // the Lagrange-3 Element (aka "Q3")

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
#include <kernel/lafem/sparse_matrix_factory.hpp>          // for SparseMatrixFactory

//not needed, i think...
// FEAT-Global includes
// #include <kernel/global/gate.hpp>                          // NEW: for Global::Gate
// #include <kernel/global/filter.hpp>                        // NEW: for Global::Filter
// #include <kernel/global/matrix.hpp>                        // NEW: for Global::Matrix
// #include <kernel/global/vector.hpp>                        // NEW: for Global::Vector

// FEAT-Solver includes
#include <kernel/solver/umfpack.hpp>                          // for umfpack
// #include <kernel/solver/pcg.hpp>                           // for PCG
// #include <kernel/solver/schwarz_precond.hpp>               // NEW: for SchwarzPrecond
// #include <kernel/solver/ilu_precond.hpp>                   // NEW: for ILUPrecond

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



using namespace FEAT;

namespace FETI{
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

  // Our local matrix type: a standard CSR matrix
  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> LocalMatrixType;

  // Our local vector type: the usual dense vector
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> LocalVectorType;

  // Our local filter type: the unit filter for Dirichlet boundary conditions
  typedef LAFEM::UnitFilter<MemType, DataType, IndexType> LocalFilterType;

  // The vector mirror takes the usual memory, data and index types as template parameters:
  typedef LAFEM::VectorMirror<MemType, DataType, IndexType> VectorMirrorType;

  //The two umfpack types we will use:
  typedef Solver::Umfpack RegularUmfpack;
  typedef Solver::UmfpackMean FloatingUmfpack;
  //see synch_vec regarding memory...
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> BufferMain;

  //This class handles the local information on each subdomain:
  class LocalDomain
  {
  public:
    IndexType domain_rank;
    IndexType domain_size;
    IndexType neighbour_maxsize;
    // gate_rank describes only relevant neighbours
    std::vector<int> gate_ranks;
    std::vector<VectorMirrorType> gate_mirrors;
    std::vector<bool> gate_floating;

    // size of gate_ranks
    IndexType gate_ranks_size =0;

    //This vector holds the signs for the B operation determined by a simple rule:
    // If the neighbour has a higher rank number, it is +1
    // if the neighbour has a lower one, it is -1
    std::vector<int> gate_signs;
    //for error calculation
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function_backup;

    // mirror dofs is a vector of size size(gate_ranks) +1 saving the first indice in BR_columns for the connected mirror
    // and the overall size in the last entry
    std::vector<IndexType> mirror_dofs;
    LocalFilterType filter_local;
    LocalMatrixType matrix_local;
    std::unique_ptr<RootMeshNodeType> root_mesh_node;


    MeshType* mesh = nullptr;
    std::shared_ptr<TrafoType> trafo;
    std::shared_ptr<SpaceType> space_backup;
    LocalVectorType vec_sol_local;
    LocalVectorType vec_rhs_local;
    LocalVectorType R_vector;
    //This members are for the CG-algorithm
    std::vector<LocalVectorType> residuum;
    std::vector<LocalVectorType> residuum_copy;
    //residuum_opt is the projected residuum
    std::vector<LocalVectorType> residuum_opt;
    //w is the preconditioned residuum
    std::vector<LocalVectorType> w;
    std::vector<LocalVectorType> w_opt;
    std::vector<LocalVectorType> s_vector;
    std::vector<LocalVectorType> Fs_vector;
    LocalVectorType vec_sol_buffer;
    LocalVectorType vec_rhs_buffer;

    //This matrix will contain the rows of the matrix BR connected to the halo of this domain. As only entries connected to the edge of this subdomain will be non zero
    // this can be assembled with communication between neighbours only
    LocalMatrixType BR_columns;
    //there is no reason to handle lambda_0 different to the residuum...
    std::vector<LocalVectorType> lambda_0;
    IndexType lambda_overall_size = 0;
    DataType alpha = 0;
    bool floating;
    // we use two pointers to umfpack, one of which will always be a nullptr
    // why is this shared? Because its shared in its hpp...
    std::shared_ptr<RegularUmfpack> reg_umf;
    std::shared_ptr<FloatingUmfpack> float_umf;



    //default constructor
    LocalDomain():
    domain_rank(IndexType(0)),
    domain_size(IndexType(0)),
    floating(true),
    reg_umf(nullptr),
    float_umf(nullptr)
    {}
    //we use the constructor to initalize and assemble everything for which no communication is needed:
    explicit LocalDomain(IndexType rank, IndexType size, IndexType level) :
    domain_rank(rank),
    domain_size(size),
    floating(true),
    reg_umf(nullptr),
    float_umf(nullptr)
    {
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Geometry initialisation
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

      //we wont need neighbour_ranks itself, so just use it as local variable...
      std::vector<int> neighbour_ranks;

      IndexType lvl = Geometry::UnitCubePatchGenerator<MeshType>::create_unique(
         int(domain_rank),          // input:  the rank of this process
         int(domain_size),          // input:  the total number of processes
         root_mesh_node,       // output: the root-mesh-node shared pointer
         neighbour_ranks );      // output: the neighbour ranks vector



      // Now let's refine the mesh up to the level that was passed as a parameter to this function,
    // assuming that it is greater than the base-mesh level the generator gave us:
    if(lvl < level)
    {
      std::cout << "Refining Mesh to Level " + stringify(level) + "..." << std::endl;

      for(; lvl < level; ++lvl)
      {
        root_mesh_node = root_mesh_node->refine_unique();
      }
    }

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Trafo and Space initialisation
    mesh = root_mesh_node->get_mesh();


    // Create the trafo
    trafo = std::make_shared<TrafoType>(*mesh);
    // Create the desired finite element space
    //SpaceType space(trafo);
    space_backup = std::make_shared<SpaceType>(*trafo);

    for(auto it = neighbour_ranks.begin(); it != neighbour_ranks.end(); ++it)
    {
      // Get the rank of our neighbour process:
      const int neighbour_rank = (*it);

      // As already mentioned before, our mesh-node does not only contain the mesh itself,
      // but also other important stuff, such as the "halos", which describe the overlap of
      // neighboured patches. For the assembly of the mirror, we need to get the halo mesh-part
      // from the mesh-node first:
      const MeshPartType* neighbour_halo = root_mesh_node->get_halo(neighbour_rank);

      // Ensure that we have a halo for this neighbour rank:
      XASSERTM(neighbour_halo != nullptr, "Failed to retrieve neighbour halo!");

      //Ensure the halo does not contain only one dof, this is for now only for 2D meshes, for 3D meshes edges also have to be checked:
      // Question: I want to check how many elements are in vertices, ie 1 dimensional points. If i want to check if there is only one edge
      // this can be taken out... of course we would then have more communication needed
      if(neighbour_halo->get_num_entities(0) == IndexType(1))
      {
        continue;
      }

      // Now that we have the halo mesh-part, we can create and assemble the corresponding mirror:
      VectorMirrorType neighbour_mirror;

      // Call the MirrorAssembler to do the dirty work for us:
      Assembly::MirrorAssembler::assemble_mirror(
        neighbour_mirror,   // the mirror that is to be assembled
        *space_backup,              // the FE space for which we want to assemble the mirror
        *neighbour_halo     // the halo mesh-part that the mirror is to be assembled on
      );

      //first we initialize the residuum, s and Fs with the size of the mirror:
      residuum.push_back(LocalVectorType(neighbour_mirror.num_indices(), 0.));
      residuum_copy.push_back(LocalVectorType(neighbour_mirror.num_indices(), 0.));
      residuum_opt.push_back(LocalVectorType(neighbour_mirror.num_indices(), 0.));
      w.push_back(LocalVectorType(neighbour_mirror.num_indices(), 0.));
      w_opt.push_back(LocalVectorType(neighbour_mirror.num_indices(), 0.));
      s_vector.push_back(LocalVectorType(neighbour_mirror.num_indices(), 0.));
      Fs_vector.push_back(LocalVectorType(neighbour_mirror.num_indices(), 0.));
      lambda_0.push_back(LocalVectorType(neighbour_mirror.num_indices(), 0.));
      // Once the mirror is assembled, we give it to our pseudo-gate.
      gate_ranks.push_back(neighbour_rank);
      gate_mirrors.push_back(std::move(neighbour_mirror));
      gate_floating.push_back(false);
//       ctag_dof_size.push_back(neighbour_halo->get_num_entities(0));
      if(neighbour_rank > int(domain_rank))
      {
        gate_signs.push_back(1);
      }
      else
      {
        gate_signs.push_back(-1);
      }

      // continue with next neighbour rank
    }
    //initialize mirror_dofs with the right size
    gate_ranks_size = gate_ranks.size();
    mirror_dofs = std::vector<IndexType>(gate_ranks_size+1);
    //first we save the maxsize on our own knowlegde
    neighbour_maxsize = gate_ranks_size;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Symbolic linear system assembly

    // Now let us assemble the local matrix structure, just as in any previous tutorial:
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_local, *space_backup);

    vec_sol_local = matrix_local.create_vector_r();
    vec_rhs_local = matrix_local.create_vector_l();
    vec_sol_buffer = matrix_local.create_vector_r();
    vec_rhs_buffer = matrix_local.create_vector_l();
    //initalize R_vector to the right size, regardless if domain is floating
    R_vector = matrix_local.create_vector_l();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Numerical assembly: matrix (bilinear forms)
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
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_local, laplace_operator, *space_backup, cubature_factory);

    // Initialize the right-hand-side vector entries to zero.
    vec_rhs_local.format();
    vec_sol_local.format();

    // Again, we use the sine-bubble as a reference solution:
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;

    // Create a force functional for our solution:
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);

    // And assemble the local rhs vector:
    Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs_local, force_functional, *space_backup, cubature_factory);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly

    // Now, we have to assemble the filter.
    // In the case of the unit-filter, which is responsible for Dirichlet boundary conditions,
    // the assembly is quite simple: We just have to assemble the corresponding local filter.
    // So we need a unit-filter assembler:
    Assembly::UnitFilterAssembler<MeshType> unit_asm;

    // To keep this part of the code dimension-independent, we loop up to 2*dimension:
    for(int ibnd(0); ibnd < 2*ShapeType::dimension; ++ibnd)
    {
      // Build the name of the boundary mesh-part and call the mesh-node's "find_mesh_part"
      // function to get a pointer to the mesh-part:
      MeshPartType* bnd_mesh_part = root_mesh_node->find_mesh_part("bnd:" + stringify(ibnd));
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
        floating = false;
      }
    }

    // At this point, we have added all boundary mesh-part of our patch to the assembler.
    // Note: If our patch is an inner patch, then we have added no boundary mesh-parts to
    // the assembler, but that's fine as the unit-filter assembler can work with that, as
    // in this case the unit-filter will be an empty filter.

    // Finally, assemble the local filter:
    unit_asm.assemble(filter_local, *space_backup);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition imposition

    // Apply the filter onto the system matrix...
    filter_local.filter_mat(matrix_local);
    // Now we impose the boundary conditions into our linear system by applying the local filter
    // onto the local right-hand-side and initial solution vectors:

    filter_local.filter_rhs(vec_rhs_local);
    filter_local.filter_sol(vec_sol_local);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Generating the R vector

    generate_R(floating);

    if(floating)
    {
      float_umf = Solver::new_umfpack_mean(matrix_local, R_vector);
      float_umf->init();
    }
    else
    {
      reg_umf = Solver::new_umfpack(matrix_local);
      reg_umf->init();
    }





    //and here ends the local initialisation...
    }

    //has to be called before deleting this...
    void free_umfpack()
    {
      if(reg_umf != nullptr)
        reg_umf->done();
      if(float_umf != nullptr)
        float_umf->done();
    }

    //for now we just generate a constant one vector...
    void generate_R(bool flo)
    {
      if(flo)
      {
        IndexType size = vec_rhs_local.size();
        R_vector = LocalVectorType(size, DataType(1));
      }
    }

    IndexType find_index(const IndexType neighbour) const
    {
      for(IndexType i(0) ; i < gate_ranks_size ; ++i)
      {
        if(gate_ranks[i] == int(neighbour))
          return i;
      }
      XABORTM("Can not find a neighbour!");
      //the return value is always out of the range of domains... so this will generate an error
      return domain_size + 1;
    }

    void solve(LocalVectorType &solution, LocalVectorType &rhs)
    {
      if(floating)
      {
        float_umf->apply(solution, rhs);
      }
      else
      {
        reg_umf->apply(solution, rhs);
      }
    }

    /**This function assembles the BR matrix for a all domains and the lambda_0.
     * If the domain is not floating, it only gathers the BR of floating neighbours.
     * If the edge is not connected to a floating neighbour, all entries will be simply zero.
     * \param[in] doms, a std::vector of all local Domains
     * \param[in] floating_num, the number of floating_domains
     * \param[in] floating_inverse, a vector which maps the rank number of a floating domain to its "floating" number
     *
     */
    //spezialising for the case, that we now the R_vector is const 1 vector
    // bad thing: we need information about whether the neighbour is floating or not...
    void assemble_BR_columns(std::vector<std::shared_ptr<LocalDomain>> &doms)
    {
      //we will number our Matrix according to the order of the vector mirrors:
      // This matrix will have the following form:
      /*
       *
       * 1 -1 0 0 0 0 0 0 0 0   <- Vector aligned to the first lagrange multi connected to the first mirror in gate_mirrors
       * 1 -1 0 0 0 0 0 0 0 0   <- Vector aligned to the second ....
       *-1  0 0 0 1 0 0 0 0 0   <- Vector aligned to the first ... connected to the second mirror in gate_mirrors
       *-1  0 0 0 1 0 0 0 0 0
       * ^        ^
       * |        |
       * |        Number of the neighbour which the second mirror in gate_mirrors is connected
       * |
       * Number of our domain rank...
       */
      //this should be in the constructor of LocalDomain/ already done in FETIGate
      IndexType num_of_rows(0);
      for(IndexType i(0); i < gate_ranks_size; ++i)
      {
        mirror_dofs[i] = num_of_rows;
        num_of_rows += gate_mirrors[i].num_indices();
      }
      mirror_dofs[gate_ranks_size] = num_of_rows;
      lambda_overall_size = num_of_rows;
      //for now we just use SparseMatrixFactory to assemble the BR_columns matrix
      LAFEM::SparseMatrixFactory<DataType, IndexType> factory(num_of_rows, doms.size());
      if(floating)
      {
        for(IndexType i(0); i < gate_ranks_size; ++i)
        {
          for(IndexType j(mirror_dofs[i]); j < mirror_dofs[i+1]; ++j)
            factory.add(j, domain_rank, gate_signs[i]);
        }
      }

      for(IndexType i(0); i < gate_ranks_size; ++i)
      {
        if(doms[IndexType(gate_ranks[i])]->floating)
        {
          this->gate_floating[i] = true;
          for(IndexType j(mirror_dofs[i]); j < mirror_dofs[i+1]; ++j)
            factory.add(j, IndexType(gate_ranks[i]), -gate_signs[i]);
        }
      }

      BR_columns = factory.make_csr();

    }

    //this extracts the entries of a vector of size lambda_overall_size connected to a specific edge into a LocalVectorType
    // and gives back a rvalue

    LocalVectorType extract_egde_values(const LocalVectorType &vector, const IndexType dom_index) const
    {
      const IndexType start_size = mirror_dofs[dom_index];
      const IndexType local_size = mirror_dofs[dom_index +1] - start_size;
      LocalVectorType tmp(local_size);
      tmp.format();
      for(IndexType i(0); i < local_size; ++i)
      {
        tmp(i, vector(start_size +i));
      }
      return tmp;
    }

    /** The output vector has to have the right sizes
     * also this formats the vector, so any value inside will be lost...
     */
    void extract_residuum(LocalVectorType & output, std::vector<LocalVectorType>& res)
    {
      XASSERTM(output.size() == this->vec_rhs_buffer.size(), "Output vector does not have the right size.");
      output.format();
      for(IndexType i(0); i < this->gate_ranks_size; ++i)
      {
        this->gate_mirrors[i].scatter_axpy(output, res[i], this->gate_signs[i]);
      }

    }

    /**
     * This function appllies B onto an output vector of LocalVectorTypes by taking values in vec_sol_buffer... in general this will be internel data types like residuum etc.
     * Of course output has to have the correct size!
    */
    void apply_B_matrix(std::vector<LocalVectorType> &output, std::vector<std::shared_ptr<LocalDomain>> &domain_vec,
                        DataType sign = 1.)
    {
      XASSERTM(output.size() <= this->gate_ranks_size, "Output vector does not have the right size!");
      //Schleife 1 ( hier direkt ein reicive call
      std::vector<LocalVectorType> dom_buffer_1;
      std::vector<LocalVectorType> dom_buffer_2;
      for(IndexType i(0); i < this->gate_ranks_size; ++i)
      {
        LocalVectorType dom_buffer_loc = this->gate_mirrors[i].create_buffer(this->vec_sol_buffer);
        dom_buffer_loc.format();
        this->gate_mirrors[i].gather(dom_buffer_loc, this->vec_sol_buffer);
        dom_buffer_1.push_back(std::move(dom_buffer_loc));
      }

      //hier eigentlich der send call
      for(IndexType i(0); i < this->gate_ranks_size; ++i)
      {
        IndexType rel_index = domain_vec[IndexType(gate_ranks[i])]->find_index(this->domain_rank);
        LocalVectorType dom_buffer_loc = domain_vec[IndexType(gate_ranks[i])]->gate_mirrors[rel_index].create_buffer(domain_vec[IndexType(gate_ranks[i])]->vec_sol_buffer);
        dom_buffer_loc.format();
        //gather on neighbour domain
        domain_vec[IndexType(gate_ranks[i])]->gate_mirrors[rel_index].gather(dom_buffer_loc, domain_vec[IndexType(gate_ranks[i])]->vec_sol_buffer);
        //now reicive data
        dom_buffer_2.push_back(std::move(dom_buffer_loc));

      }

      //Schleife 2
      //now add both with regard to the gate sign, would be nice, if the if else case could be evaluated at compile time...
      for(IndexType i(0); i < this->gate_ranks_size; ++i)
      {
        XASSERTM(this->gate_signs[i] != 0, "Gate signs are not correctly initialized!");
        if(this->gate_signs[i]*sign > 0)
        {
          output[i].axpy(dom_buffer_2[i], dom_buffer_1[i], -1.);
        }
        else
        {
          output[i].axpy(dom_buffer_1[i], dom_buffer_2[i], -1.);
        }
      }
    }

    /**
     * This function applies K^{-1}B^T on a vector of LocalVectorTypes and saves it into the output_vec
     * Output vector will nearly always will be vec_sol_buffer, so giving a reference as parameter is a little bit wasted...
     */
    void apply_KBT_on_vector(LocalVectorType & output_vec, std::vector<LocalVectorType> & input)
    {
      this->extract_residuum(this->vec_rhs_buffer, input);
      this->solve(output_vec, this->vec_rhs_buffer);
    }

    //This Applies K * B^T on a residuum like structure, meaning a vector of LocalVectorType
    // First we exstract the values through our mirrors, this handles a function and then apply matrix_local
    void apply_precond(LocalVectorType & output_vec, std::vector<LocalVectorType>& input)
    {
      this->extract_residuum(this->vec_rhs_buffer, input);
      this->matrix_local.apply(output_vec, this->vec_rhs_buffer);
    }

    //This initializes s_0 = h_opt
    void init_s_vector()
    {
      XASSERTM(this->s_vector.size() == this->w_opt.size(), "Sizes do not match!");
      for(IndexType i(0); i < this->w_opt.size(); ++i)
      {
        this->s_vector[i].copy(this->w_opt[i]);
      }
    }

    template<typename T>
    T gather_s_dot()
    {
      T lokal_s{0};
      for(IndexType i(0); i < gate_ranks_size ; ++i)
      {
        lokal_s += T(s_vector[i].dot(Fs_vector[i]));
      }
      return lokal_s;
    }

    //this returns the lokal dot product of residuum_opt and w_opt
    template<typename T>
    T gather_residuum_dot()
    {
      T lokal_r{0};
      for(IndexType i(0); i < gate_ranks_size ; ++i)
      {
        lokal_r += T(residuum_opt[i].dot(w_opt[i]));
      }
      return lokal_r;
    }

    template<typename T>
    T gather_local_res_dot()
    {
      T local_res{0};
      for(IndexType i(0); i < gate_ranks_size ; ++i)
      {
        local_res += T(residuum_opt[i].dot(residuum_opt[i]));
      }
      return local_res;
    }

    /**
     * This function add scalar*input onto the designated vector target for vectors of type
     * vector<LocalVectorType>
     */
    void vector_update(std::vector<LocalVectorType> &target, std::vector<LocalVectorType> &input, DataType scalar)
    {
      DataType tmp = 0;
      for(IndexType i(0); i < gate_ranks_size ; ++i)
      {
        XASSERTM(target[i].size() == input[i].size(), "Sizes do not match!");
        for(IndexType j(0); j < target[i].size(); ++j)
        {
          tmp = target[i](j);
          target[i](j, tmp + scalar*input[i](j));
        }
      }
    }

    /**
     * This function add input onto the designated vector target times scalar for vectors of type
     * vector<LocalVectorType>
     */
    void vector_update_alt(std::vector<LocalVectorType> &target, std::vector<LocalVectorType> &input, DataType scalar)
    {
      DataType tmp = 0;
      for(IndexType i(0); i < gate_ranks_size ; ++i)
      {
        XASSERTM(target[i].size() == input[i].size(), "Sizes do not match!");
        for(IndexType j(0); j < target[i].size(); ++j)
        {
          tmp = target[i](j);
          target[i](j, scalar*tmp + input[i](j));
        }
      }
    }

    //copies residuum into a placeholder
    void residuum_start_copy()
    {
      for(IndexType i(0); i < gate_ranks_size; ++i)
      {
        residuum_copy[i].copy(residuum[i]);
      }
    }

    //calculate the solution after lambda and alpha have been determined
    void calc_solution()
    {
      //setup right_side
      this->extract_residuum(this->vec_rhs_buffer, this->lambda_0);
      LocalVectorType new_rhs(vec_rhs_local.size(), 0.);
      new_rhs.axpy(this->vec_rhs_buffer, this->vec_rhs_local);
      //calculate solution without kernel
      this->solve(this->vec_sol_buffer, new_rhs);
      //add kernel
      if(this->floating)
      {
        this->vec_sol_local.axpy(this->R_vector, this->vec_sol_buffer, - this->alpha);
      }
      else // or not
      {
        //I think i have to do this this way to keep the values of the filter
        vec_sol_local.axpy(vec_sol_local, vec_sol_buffer);
      }
    }

    //small hack to calculate errors...
    std::array<DataType,2> calc_errors()
    {
      Cubature::DynamicFactory cubature_factory("auto-degree:5");

      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
      this->vec_sol_local, this->sol_function_backup, *space_backup, cubature_factory);
       return std::array<DataType, 2>{{errors.norm_h0, errors.norm_h1}};
    }

    /**This creates a DenseVector/BufferMain of size n = 2*(neighbour_maxsize+1) representing a row of the global Q matrix, thereby
     * the first n entries represent the Indicies of the non zero entries.
     * Because we will not always have full size vectors these Indicies are encoded the following way:
     * First, always the non diagonal(so Index != domain_num) are saved in ascending order.
     * Then the Diagonal Index is saved, indicating the last vector entry
     *
     * !!!! I could also split those in the real parallel Code... see gather
    */
    BufferMain get_Q_row()
    {
      IndexType n = 2*(this->neighbour_maxsize +1);
      BufferMain buff(n);
      if(!this->floating)
      {
        buff(0, double(this->domain_rank));
        buff(n/2, 1.);
        return buff;
      }
      IndexType counter = 0;
      for(IndexType i(0); i < gate_ranks_size; ++i)
      {
        IndexType neighbour = IndexType(this->gate_ranks[i]);
        if(!gate_floating[i])
        {
          ++counter;
          continue;
        }
        DataType scalar = DataType(this->mirror_dofs[i+1] - this->mirror_dofs[i]);
        buff(i-counter, DataType(neighbour));
        buff(i-counter+(n/2), -scalar);
      }
      buff(gate_ranks_size-counter, DataType(this->domain_rank));
      buff(gate_ranks_size-counter+(n/2), DataType(this->mirror_dofs.back()));
      return buff; //no mobe, because buff is a local object...
    }

  };// class LocalDomain

  //wrapper um umfpack
  struct RegUmf
  {
    RegularUmfpack Q;

    RegUmf(LocalMatrixType const& loc)
    : Q{loc}
    {
      Q.init();
    }
    ~RegUmf()
    {
      Q.done();
    }
    void apply(LocalVectorType &l, LocalVectorType const &r)
    {
      Q.apply(l,r);
    }
  }; //class RegUmf

  /**
   * This is a small class which handles the application of the projection P = (BR)(Q^{-1})(BR)^T
   *
   */
  //now optimized for R_vector = (1,....1) on every domain
  class Projector
  {
  private:
    std::shared_ptr<RegUmf> Q_factorized = nullptr;
    LocalMatrixType Q_matrix;
    LocalVectorType right_side;
    LocalVectorType right_side_buffer;
  public:
    Projector(){}
    Projector( std::vector<std::shared_ptr<LocalDomain>>& domain_vec, IndexType neighbour_maxsize )
    {
      IndexType dom_size = IndexType(domain_vec.size());
      LAFEM::SparseMatrixFactory<DataType, IndexType> factory(dom_size, dom_size);
      //post reicives
      std::vector<BufferMain> recv_bufs;
      std::vector<BufferMain> recv_reqs;
      recv_bufs.resize(dom_size);
      recv_reqs.resize(dom_size);

      for(IndexType i(0); i < dom_size; ++i)
      {
        //this should be BufferMain...
        recv_bufs.at(i) = BufferMain(2*(neighbour_maxsize+1));//, LAFEM::Pinning::disabled);

        //post reicive(we emmulate this by init another buffer
        recv_reqs.at(i) = BufferMain(2*(neighbour_maxsize+1));
      }
      //post sends
      std::vector<BufferMain> send_bufs;
      std::vector<BufferMain> send_reqs;
      send_bufs.resize(dom_size);
      send_reqs.resize(dom_size);
      for(IndexType i(0); i < dom_size; ++i)
      {
        //get the row of i-th domain and save this into a buffer
        // this has to be in main memeory!
        send_bufs.at(i) = domain_vec[i]->get_Q_row();
        send_reqs.at(i) = BufferMain(send_bufs[i].size());

        //post send request
        send_reqs[i].copy(send_bufs[i]);
      }

      //here now the wait() functionality... we just copy the vectors... this should be a broadcast...
      for(IndexType i(0); i <dom_size; ++i)
      {
        //this of course would be data transmission
        recv_reqs[i].copy(send_reqs[i]);

        recv_bufs[i].copy(recv_reqs[i]);
      }

      for(IndexType i(0); i < dom_size ; ++i)
      {
        // we calculate the scalarproduct on the domain itself and between its neighbours in one loop
        // in parallel: each domain constructs its own row...
        for(IndexType j(0); j < neighbour_maxsize +1; ++j)
        {
          factory.add(i, IndexType(recv_bufs.at(i)(j)), recv_bufs.at(i)(j+neighbour_maxsize+1));
          if(IndexType(recv_bufs.at(i)(j)) == i)
            break;
        }
      }
      Q_matrix = factory.make_csr();
      right_side_buffer = Q_matrix.create_vector_r();
      right_side = Q_matrix.create_vector_l();
      right_side_buffer.format();
      right_side.format();
      Q_factorized = std::make_shared<RegUmf>(Q_matrix);
    }
    Projector(const Projector&) = delete;
    Projector & operator=(const Projector&) = delete;
    //rule of 6.... i should implement move constructor and assignment...

    void assemble_lambda_0(std::vector<std::shared_ptr<LocalDomain>> &domain_vec)
    {
      for(IndexType i(0); i < domain_vec.size() ; ++i)
      {
        //negativ sign because ive done this in my matlab code too... have to think about this, but if a \in range(R) also -a should be...
        // so this should also be a admissible start vector...
        if(domain_vec[i]->floating)
            right_side_buffer(i, -domain_vec[i]->vec_rhs_local.dot(domain_vec[i]->R_vector));
        else
          right_side_buffer(i, 0.);
      }
      Q_factorized->apply(right_side, right_side_buffer);
      //now apply each BR to right_side
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        //std::cout << it->BR_columns.columns() << "\n";
        LocalVectorType tmp((*it)->lambda_overall_size);
        tmp.format();
        (*it)->BR_columns.apply(tmp, right_side);
        for(IndexType i(0); i < (*it)->gate_ranks_size; ++i)
        {
          (*it)->lambda_0[i] = (*it)->extract_egde_values(tmp, i);
        }
      }
    }


    void project_residuum(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      //this should not be needed... for now we format these vectors
      right_side.format();
      right_side_buffer.format();

      for(IndexType loc_index(0); loc_index < domain_vec.size() ; ++loc_index)
      {
        if(!domain_vec[loc_index]->floating)
        {
          right_side_buffer(loc_index, 0.);
          continue;
        }
        DataType scalar{0};
        std::shared_ptr<LocalDomain> loc = domain_vec[loc_index];
        //this can be handled faster if we generelly assume that the R_vector is (1,....1)
        for(IndexType i(0); i < loc->gate_ranks_size ; ++i)
        {
          LocalVectorType tmp_R = loc->gate_mirrors[i].create_buffer(loc->R_vector);
          loc->gate_mirrors[i].gather(tmp_R, loc->R_vector);
          DataType val = 0;
          for(IndexType j(0); j < tmp_R.size(); ++j)
          {
            val = tmp_R(j);
            tmp_R(j, loc->gate_signs[i]*val);
          }
          scalar += tmp_R.dot(loc->residuum[i]);
        }
        right_side_buffer(loc_index, scalar);
      }
      Q_factorized->apply(right_side, right_side_buffer);

     //now distribute right side to each domain... as we distributed BR beforehand, this can be applied locally
     for(auto it = domain_vec.begin(); it != domain_vec.end() ; ++it)
     {
       //first we alocate a temp vector
       LocalVectorType tmp((*it)->lambda_overall_size);
       tmp.format();
       (*it)->BR_columns.apply(tmp, right_side);

       //now add the values to the corresponding edges and therefore the residuum vektor
       //this should be a move operation...
       //subtract from the original residuum with axpy
       for(IndexType i(0); i < (*it)->gate_ranks_size; ++i)
       {
         (*it)->residuum_opt[i].axpy((*it)->extract_egde_values(tmp, i), (*it)->residuum[i], -1.);
         //(*it)->residuum_opt[i].copy((*it)->residuum[i]);
       }
     }
    }

    void project_w(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      //this should not be needed... for now we format these vectors
      right_side.format();
      right_side_buffer.format();

      for(IndexType loc_index(0); loc_index < domain_vec.size() ; ++loc_index)
      {
        if(!domain_vec[loc_index]->floating)
        {
          right_side_buffer(loc_index, 0.);
          continue;
        }
        DataType scalar{0};
        std::shared_ptr<LocalDomain> loc = domain_vec[loc_index];
        //this can be handled faster if we generelly assume that the R_vector is (1,....1)
        for(IndexType i(0); i < loc->gate_ranks_size ; ++i)
        {
          LocalVectorType tmp_R = loc->gate_mirrors[i].create_buffer(loc->R_vector);
          loc->gate_mirrors[i].gather(tmp_R, loc->R_vector);
          DataType val = 0;
          for(IndexType j(0); j < tmp_R.size(); ++j)
          {
            val = tmp_R(j);
            tmp_R(j, loc->gate_signs[i]*val);
          }
          scalar += tmp_R.dot(loc->w[i]);
        }
        right_side_buffer(loc_index, scalar);
      }
      Q_factorized->apply(right_side, right_side_buffer);

     //now distribute right side to each domain... as we distributed BR beforehand, this can be applied locally
     for(auto it = domain_vec.begin(); it != domain_vec.end() ; ++it)
     {
       //first we alocate a temp vector
       LocalVectorType tmp((*it)->lambda_overall_size);
       tmp.format();
       (*it)->BR_columns.apply(tmp, right_side);

       //now add the values to the corresponding edges and therefore the residuum vektor
       //this should be a move operation...
       for(IndexType i(0); i < (*it)->gate_ranks_size; ++i)
       {
         (*it)->w_opt[i].axpy((*it)->extract_egde_values(tmp, i), (*it)->w[i], -1.);
         //(*it)->w_opt[i].copy((*it)->w[i]);
       }
     }
    }

    LocalVectorType calc_alpha(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      //this should not be needed... for now we format these vectors
      right_side.format();
      right_side_buffer.format();

      for(IndexType loc_index(0); loc_index < domain_vec.size() ; ++loc_index)
      {
        if(!domain_vec[loc_index]->floating)
        {
          right_side_buffer(loc_index, 0.);
          continue;
        }
        DataType scalar{0};
        std::shared_ptr<LocalDomain> loc = domain_vec[loc_index];
        //this can be handled faster if we generelly assume that the R_vector is (1,....1)
        for(IndexType i(0); i < loc->gate_ranks_size ; ++i)
        {
          LocalVectorType tmp_R = loc->gate_mirrors[i].create_buffer(loc->R_vector);
          loc->gate_mirrors[i].gather(tmp_R, loc->R_vector);
          DataType val = 0;
          for(IndexType j(0); j < tmp_R.size(); ++j)
          {
            val = tmp_R(j);
            tmp_R(j, loc->gate_signs[i]*val);
          }
          scalar += tmp_R.dot(loc->residuum_copy[i]);
        }
        right_side_buffer(loc_index, scalar);
      }
      Q_factorized->apply(right_side, right_side_buffer);


      return std::move(right_side); //std::move needed because right_side is a permanent object....
      //Question: does std::move erase the number of entries?? I think not, so this should be fine

    }

  };//class Projector


    //small function to apply pre_conditioning on residuum
    void precond(std::vector<std::shared_ptr<LocalDomain>> domain_vec)
    {
      //first apply K
      for(auto it = domain_vec.begin(); it != domain_vec.end() ; ++it)
      {
        (*it)->apply_precond((*it)->vec_sol_buffer, (*it)->residuum_opt);
      }
      //apply B, save this into w
      for(auto it = domain_vec.begin(); it != domain_vec.end() ; ++it)
      {
        (*it)->apply_B_matrix((*it)->w, domain_vec);
      }

    }

    //this functions are a bit tricky, as we can not simply take the dot product of each domain and add them, because we
    // then count some of them multiple times...
    // in 2d exactly 2 times... so just half the values...
    DataType gather_residuum_dot(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      DataType global_r{0.};
      //iterate over all LocalDomains while calling the local gathering
      //in parallel this could be optimised by grouping...
      for(auto it = domain_vec.begin(); it != domain_vec.end() ; ++it)
      {
        global_r += (*it)->gather_residuum_dot<DataType>();
      }
      return global_r/2.;
    }

    DataType gather_s_dot(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      DataType global_s{0.};
      //iterate over all LocalDomains while calling the local gathering
      //in parallel this could be optimised by grouping...
      for(auto it = domain_vec.begin(); it != domain_vec.end() ; ++it)
      {
        global_s += (*it)->gather_s_dot<DataType>();
      }
      return global_s/2.;
    }

    DataType norm_res(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      DataType global_res{0.};
      //iterate over all LocalDomains while calling the local gathering
      //in parallel this could be optimised by grouping...
      for(auto it = domain_vec.begin(); it != domain_vec.end() ; ++it)
      {
        global_res += (*it)->gather_local_res_dot<DataType>();
      }
      return std::sqrt(global_res/2.);
    }

    void calc_Fs(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
       {
         (*it)->apply_KBT_on_vector((*it)->vec_sol_buffer, (*it)->s_vector);
       }
       for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
       {
         (*it)->apply_B_matrix((*it)->Fs_vector, domain_vec);
       }
    }

    void calc_Fs_from_lambda(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        (*it)->apply_KBT_on_vector((*it)->vec_sol_buffer, (*it)->lambda_0);
      }

      //now we need to extract the values from vec_sol_buffer. This needs communication between neighbour mirrors
      for(auto it = domain_vec.begin(); it != domain_vec.end() ; ++it)
      {
        (*it)->apply_B_matrix((*it)->Fs_vector, domain_vec);
      }
    }

    void update_lambda(std::vector<std::shared_ptr<LocalDomain>>& domain_vec, DataType scalar)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        (*it)->vector_update((*it)->lambda_0, (*it)->s_vector, scalar);
      }
    }

    void update_residuum(std::vector<std::shared_ptr<LocalDomain>>& domain_vec, DataType scalar)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        (*it)->vector_update((*it)->residuum, (*it)->Fs_vector, - scalar);
      }
    }

    void update_s_vector(std::vector<std::shared_ptr<LocalDomain>>& domain_vec, DataType scalar)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        (*it)->vector_update_alt((*it)->s_vector, (*it)->w_opt, scalar);
      }
    }
    void update_residuum_copy(std::vector<std::shared_ptr<LocalDomain>>& domain_vec, DataType scalar)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        (*it)->vector_update_alt((*it)->residuum_copy, (*it)->Fs_vector, scalar);
      }
    }

    void calc_solution(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        (*it)->calc_solution();
      }
    }

    std::array<DataType,2> print_errors(std::vector<std::shared_ptr<LocalDomain>>& domain_vec)
    {
      DataType h0 = 0.;
      DataType h1 = 0.;
      std::array<DataType, 2> tmp;
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        tmp = (*it)->calc_errors();

        h0 += Math::sqr(tmp[0]);
        h1 += Math::sqr(tmp[1]);
      }
      h0 = Math::sqrt(h0);
      h1 = Math::sqrt(h1);
      //std::cout << "H0-Error: " << h0 << "\n";
      //std::cout << "H1-Error: " << h1 << "\n";
      return std::array<DataType,2>{{h0, h1}};
    }

    //test_functions
    void test_Fs(std::vector<std::shared_ptr<LocalDomain>>& domain_vec, DataType tol)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        for(IndexType i(0); i < (*it)->gate_ranks_size; ++i)
        {
          IndexType this_index = domain_vec[IndexType((*it)->gate_ranks[i])]->find_index((*it)->domain_rank);
          LocalVectorType & first = (*it)->Fs_vector[i];
          LocalVectorType & second = domain_vec[IndexType((*it)->gate_ranks[i])]->Fs_vector[this_index];

          XASSERTM(first.size() == second.size(), "Sizes do not match");
          for(IndexType j(0); j<first.size(); ++j)
          {
            if(std::abs(first(j)-second(j)) >= tol)
            {
              std::cout << "Error occured for Fs between domain " << (*it)->domain_rank << " and domain " << (*it)->gate_ranks[i] << "\n";
              std::cout << "First vector: " << first << "\n";
              std::cout << "Second vector: " << second << "\n";
              XABORTM(" Error occured in FS!");
            }
          }
        }
      }
    }

    void test_residuum(std::vector<std::shared_ptr<LocalDomain>>& domain_vec, DataType tol)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        for(IndexType i(0); i < (*it)->gate_ranks_size; ++i)
        {
          IndexType this_index = domain_vec[IndexType((*it)->gate_ranks[i])]->find_index((*it)->domain_rank);
          LocalVectorType & first = (*it)->residuum[i];
          LocalVectorType & second = domain_vec[IndexType((*it)->gate_ranks[i])]->residuum[this_index];

          XASSERTM(first.size() == second.size(), "Sizes do not match");
          for(IndexType j(0); j<first.size(); ++j)
          {
            if(std::abs(first(j)-second(j)) >= tol)
            {
              std::cout << "Error occured for residuum between domain " << (*it)->domain_rank << " and domain " << (*it)->gate_ranks[i] << "\n";
              std::cout << "First vector: " << first << "\n";
              std::cout << "Second vector: " << second << "\n";
              XABORTM(" Error occured in residuum!");
            }
          }
        }
      }
    }
    void test_lambda(std::vector<std::shared_ptr<LocalDomain>>& domain_vec, DataType tol)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        for(IndexType i(0); i < (*it)->gate_ranks_size; ++i)
        {
          IndexType this_index = domain_vec[IndexType((*it)->gate_ranks[i])]->find_index((*it)->domain_rank);
          LocalVectorType & first = (*it)->lambda_0[i];
          LocalVectorType & second = domain_vec[IndexType((*it)->gate_ranks[i])]->lambda_0[this_index];

          XASSERTM(first.size() == second.size(), "Sizes do not match");
          for(IndexType j(0); j<first.size(); ++j)
          {
            if(std::abs(first(j)-second(j)) >= tol)
            {
              std::cout << "Error occured for lambda between domain " << (*it)->domain_rank << " and domain " << (*it)->gate_ranks[i] << "\n";
              std::cout << "First vector: " << first << "\n";
              std::cout << "Second vector: " << second << "\n";
              XABORTM(" Error occured in lambda");
            }
          }
        }
      }
    }
    void test_s_vector(std::vector<std::shared_ptr<LocalDomain>>& domain_vec, DataType tol)
    {
      for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      {
        for(IndexType i(0); i < (*it)->gate_ranks_size; ++i)
        {
          IndexType this_index = domain_vec[IndexType((*it)->gate_ranks[i])]->find_index((*it)->domain_rank);
          LocalVectorType & first = (*it)->s_vector[i];
          LocalVectorType & second = domain_vec[IndexType((*it)->gate_ranks[i])]->s_vector[this_index];

          XASSERTM(first.size() == second.size(), "Sizes do not match");
          for(IndexType j(0); j<first.size(); ++j)
          {
            if(std::abs(first(j)-second(j)) >= tol)
            {
              std::cout << "Error occured for s between domain " << (*it)->domain_rank << " and domain " << (*it)->gate_ranks[i] << "\n";
              std::cout << "First vector: " << first << "\n";
              std::cout << "Second vector: " << second << "\n";
              XABORTM(" Error occured in s");
            }
          }
        }
      }
    }

  std::array<DataType, 2> main(IndexType level, IndexType domain_num)
  {
    std::vector<std::shared_ptr<LocalDomain>> domain_vec(domain_num);
    DataType tol = 1e-15;
    for(IndexType i(0); i < domain_num; ++i)
    {
      domain_vec[i] = std::make_shared<LocalDomain>(i, domain_num, level);
    }

    //communicate max number of neighbours
    IndexType neighbour_maxsize = 0;
    for(auto it = domain_vec.begin(); it != domain_vec.end(); ++it)
    {
       IndexType local_maxsize = (*it)->neighbour_maxsize;
       if(neighbour_maxsize < local_maxsize)
       {
            neighbour_maxsize = local_maxsize;
       }
    }
    for(auto it = domain_vec.begin(); it != domain_vec.end(); ++it)
    {
      (*it)->neighbour_maxsize = neighbour_maxsize;
    }

    std::cout << "Everything initalized!" << std::endl;
    for(auto it = domain_vec.begin(); it != domain_vec.end(); ++it)
    {
       (*it)->assemble_BR_columns(domain_vec);
    }


    // Initialize the Projector class
     Projector project(domain_vec, neighbour_maxsize);

     //initialize the starting lambda
     project.assemble_lambda_0(domain_vec);
     test_lambda(domain_vec, tol);

     //next construct starting residuum... this will be locally handeled due to B
     //first we will calculate K^-1 f on each domain and save this into a placeholder

     for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
     {
       (*it)->solve((*it)->vec_sol_buffer, (*it)->vec_rhs_local);
     }

     //now we construct the residuum by gathering the data from the buffer on the domain itself and the neighbour domain
     for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
     {
       (*it)->apply_B_matrix((*it)->residuum, domain_vec, -1.);
       //save a copy into another variable, as we need it later
       (*it)->residuum_start_copy();
     }
     test_residuum(domain_vec, tol);



     //next step: generating the starting Fs (BK^-1B^T *lam_0)
     calc_Fs_from_lambda(domain_vec);

     //calculate the residuum: res = res - Fs
     update_residuum(domain_vec, 1.);

     project.project_residuum(domain_vec);



     //apply preconditoning. This again needs communication between neighbours:
     precond(domain_vec);


     //and now project again into w_opt...
     project.project_w(domain_vec);

     //init s_vector
     for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
     {
       (*it)->init_s_vector();
     }
     test_s_vector(domain_vec, tol);

     IndexType max_index = 1000;
     DataType r_0 = 0.;
     DataType s_l = 0.;
     DataType r_1 = gather_residuum_dot(domain_vec);
     DataType nrOopt = norm_res(domain_vec);
     DataType r_alt = norm_res(domain_vec);
     XASSERTM(nrOopt > 0, "norm is zero!");
     DataType gamma = 0.;
     DataType beta = 0.;
     DataType L1_error = 0.;
     DataType CG_tol = 1e-15;

     //we should now have everything ready for the CG-Algorithm

     for(IndexType k(0); k < max_index ; ++k)
     {
       //first step calculate intermediate value Fs = BK^{-1}B^T *s
       calc_Fs(domain_vec);
       test_Fs(domain_vec, tol);

       //now calculate dot(s,Fs)
       s_l = gather_s_dot(domain_vec);
       gamma = r_1/s_l;

       //update lambda : lam = lam + s*gamma
       update_lambda(domain_vec, gamma);
       test_lambda(domain_vec, tol);
       //update residuum : r = r - Fs*gamma
       update_residuum(domain_vec, gamma);
       test_residuum(domain_vec, tol);
       //project residuum
       project.project_residuum(domain_vec);
       //precond
       precond(domain_vec);
       //project preconditioned vektor w
       project.project_w(domain_vec);
       //update r_0 and r_1
       r_0 = r_1;
       r_1 = gather_residuum_dot(domain_vec);
       //calculate beta
       XASSERTM(r_0 != 0., "Division by zero...");
       beta = r_1/r_0;
       //update s_vector : s = w_opt + s*beta
       update_s_vector(domain_vec, beta);

       r_alt = norm_res(domain_vec);

       //calculate relative error:
       L1_error = r_alt/nrOopt;
       std::cout << "CG-Iteration: " << k << " | rel l1_error: " << L1_error << "\n";

       //stopping criteria
       if(L1_error < CG_tol)
       {
         std::cout << "Finished on " << k << "th iteration!\n";
         break;
       }
     }

     //now calculate alpha

     //For this, calculate Fs = BK^{-1}B' * lambda_0
     calc_Fs_from_lambda(domain_vec);

     //then calculate res_copy = Fs - res_copy
     update_residuum_copy(domain_vec, -1.);

     //and now calculate alpha through projector
     LocalVectorType alpha = project.calc_alpha(domain_vec);


     //distribute alpha to the floating domains:
     for(IndexType i(0); i < domain_vec.size() ; ++i)
     {
       domain_vec[i]->alpha = alpha(i);
     }

     //now calclulate the local_solutions
     calc_solution(domain_vec);


    //free Q_fac space...
    for(auto it = domain_vec.begin() ; it != domain_vec.end() ; ++it)
      (*it)->free_umfpack();

    //return local errors
     return print_errors(domain_vec);

  }
}

void error_analysis(FETI::IndexType domain_num, std::vector<FETI::IndexType> lvl, std::vector<FETI::DataType> h0_error,
                     std::vector<FETI::DataType> h1_error, std::vector<FETI::DataType> rate_h0,
                     std::vector<FETI::DataType> rate_h1)
{
  std::cout << "Convergence test for " << domain_num << " sub-domains:\n";
  std::cout << "------------------------------------------\n";
  std::cout << "Level     |      H0 Error      |    rate      |       H1 Error     |      rate \n";
  for(FETI::IndexType i(0); i < lvl.size(); ++i)
  {
    std::cout << lvl[i] << "    |   " << h0_error[i] << "        |    ";
    if(i==0)
      std::cout << "---     |         ";
    else
      std::cout << "  " << rate_h0[i] << "   |    ";
    std::cout << "  " << h1_error[i] << "      |    ";
    if(i==0)
      std::cout << "---             \n";
    else
      std::cout << "  " << rate_h1[i] << "\n ";


  }
}

int main(int argc, char* argv[])
{
  Runtime::initialize(argc, argv);

  std::vector<FETI::IndexType> lvl;
  std::vector<FETI::DataType> h0_error;
  std::vector<FETI::DataType> h1_error;
  std::vector<FETI::DataType> rate_h0;
  std::vector<FETI::DataType> rate_h1;
  FETI::IndexType domain_num = 4*4;
  FETI::IndexType start_level = 3;
  FETI::IndexType end_level = 6;

  for(FETI::IndexType level(start_level); level <= end_level; ++level)
  {
    lvl.push_back(level);
    std::array<FETI::DataType,2> err = FETI::main(level, domain_num);
    if(level > start_level)
    {
      rate_h0.push_back(h0_error.back()/err[0]);
      rate_h1.push_back(h1_error.back()/err[1]);
    }
    else
    {
      rate_h0.push_back(-1.);
      rate_h1.push_back(-1.);
    }
    h0_error.push_back(err[0]);
    h1_error.push_back(err[1]);
  }

  error_analysis(domain_num, lvl, h0_error, h1_error, rate_h0, rate_h1 );

  return Runtime::finalize();
}
