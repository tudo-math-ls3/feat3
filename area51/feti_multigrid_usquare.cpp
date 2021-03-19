// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime
#include <kernel/util/dist.hpp>                            // NEW: for Dist::Comm
#include <kernel/util/simple_arg_parser.hpp>               // NEW: for SimpleArgParser

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
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>              // for L2/H1-error computation
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for LaplaceFunctional
#include <kernel/assembly/mirror_assembler.hpp>            // NEW: for MirrorAssembler
#include <kernel/assembly/grid_transfer.hpp>               // NEW: for GridTransfer

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/vector_mirror.hpp>                  // NEW: for VectorMirror
#include <kernel/lafem/sparse_matrix_factory.hpp>          // for SparseMatrixFactory
#include <kernel/lafem/transfer.hpp>                       // NEW: for Transfer

// not needed, really needed, only for res analysis
// FEAT-Global includes
#include <kernel/global/gate.hpp>                          // NEW: for Global::Gate
#include <kernel/global/filter.hpp>                        // NEW: for Global::Filter
#include <kernel/global/matrix.hpp>                        // NEW: for Global::Matrix
#include <kernel/global/vector.hpp>                        // NEW: for Global::Vector

// FEAT-Solver includes
#include <kernel/solver/umfpack.hpp>                          // for umfpack
#include <kernel/solver/pcg.hpp>                           // for PCG
// #include <kernel/solver/schwarz_precond.hpp>               // NEW: for SchwarzPrecond
// #include <kernel/solver/ilu_precond.hpp>                   // NEW: for ILUPrecond
#include <kernel/solver/richardson.hpp>                    // NEW: for Richardson
#include <kernel/solver/jacobi_precond.hpp>                // NEW: for JacobiPrecond
#include <kernel/solver/multigrid.hpp>                     // NEW: for MultiGrid

#include <kernel/util/stop_watch.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/memory_usage.hpp>

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
  typedef LAFEM::MeanFilter<MemType, DataType, IndexType> MeanFilterType;

  //we create a filter chain with a Diricchlet and a mean Filter
  typedef LAFEM::FilterChain <LocalFilterType, MeanFilterType> FilterChainType;

  // The vector mirror takes the usual memory, data and index types as template parameters:
  typedef LAFEM::VectorMirror<MemType, DataType, IndexType> VectorMirrorType;

  // The next one is new: this is the class that is responsible for the grid transfer.
  // In contrast to the other LAFEM containers that we have typedefed before, the 'Transfer'
  // class template requires a matrix type instead of the usual mem-data-index type triplet.
  typedef LAFEM::Transfer<LocalMatrixType> TransferType;

  // The gate class template is defined in the "Global" namespace and it takes two
  // template arguments: the local vector and the vector mirror types:
  typedef Global::Gate<LocalVectorType, VectorMirrorType> GateType;
  typedef Global::Matrix<LocalMatrixType, VectorMirrorType, VectorMirrorType> GlobalMatrixType;
  typedef Global::Vector<LocalVectorType, VectorMirrorType> GlobalVectorType;
  typedef Global::Filter<LocalFilterType, VectorMirrorType> GlobalFilterType;

  //The two umfpack types we will use:
  typedef Solver::Umfpack RegularUmfpack;
  typedef Solver::UmfpackMean FloatingUmfpack;
  //see synch_vec regarding memory...
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> BufferMain;
  /// the buffer vector type (possibly in device memory) for now the same as buffermain
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> BufferType;

  class Umf
  {
  public:
    virtual ~Umf() = 0;

    virtual void apply(LocalVectorType &l, LocalVectorType const &r) = 0;
  };
  Umf::~Umf(){}

  //wrapper um RegularUmfpack umfpack
  class RegUmf : public Umf
  {
  private:
    RegularUmfpack Q;
  public:
    RegUmf(LocalMatrixType const& loc)
    : Q{loc}
    {
      Q.init();
    }
    virtual ~RegUmf()
    {
      Q.done();
    }
    void apply(LocalVectorType &l, LocalVectorType const &r) override
    {
      Q.apply(l,r);
    }
  }; //class RegUmf

  class FloUmf : public Umf
  {
  private:
    FloatingUmfpack Q;
  public:
    FloUmf(LocalMatrixType const& loc, LocalVectorType const & R_vector)
    : Q{loc, R_vector}
    {
      Q.init();
    }
    FloUmf(LocalMatrixType const& loc, MeanFilterType & R_vector)
    : Q{loc, R_vector}
    {
      Q.init();
    }
    virtual ~FloUmf()
    {
      Q.done();
    }
    void apply(LocalVectorType &l, LocalVectorType const &r) override
    {
      Q.apply(l,r);
    }
  }; //class FloUmf

// This simple Level class will act as a "container" that encapsulates these required objects
  // for each level. We declare all members as "public" and do not implement any member functions
  // (except for the constructor) to keep the program flow as simple as possible.
  class Level
  {
  public:
    // We hold a shared ptr to the root_mesh_node, so we dont create dangling referneces...
    std::shared_ptr<RootMeshNodeType> root_mesh_node;
    // The mesh for this level
    MeshType& mesh;
    // The trafo for this level
    TrafoType trafo;
    // The space for this level
    SpaceType space;

    // The system matrix for this level
    LocalMatrixType matrix;
    // The system filter for this level
    FilterChainType filter_chain;

    // The grid-transfer operators for this level
    TransferType transfer;

    // This constructor will create a mesh, trafo and space based on the mesh factory.
    // Note that the other member variables of this level (matrix, filter, etc.) are
    // not initialized here - this will be performed in the main function below.
    explicit Level(std::shared_ptr<RootMeshNodeType> input_root_mesh_node) :
      root_mesh_node{input_root_mesh_node},
      mesh(*root_mesh_node->get_mesh()),
      trafo(mesh),
      space(trafo)
    {
    }
  }; // class Level

  //a small class which wrapps the local multigrid solver
  class LocalSolver
  {
  public:
//     std::shared_ptr<Umf> test_solver;
    std::shared_ptr<Solver::MultiGridHierarchy<
    LocalMatrixType,          // the system matrix type
    FilterChainType,          // the system filter type
    TransferType         // the transfer operator type
    >> multigrid_hierarchy;

    std::shared_ptr<Solver::MultiGrid<LocalMatrixType, FilterChainType, TransferType>> multigrid;

    std::shared_ptr<Solver::Richardson<LocalMatrixType, FilterChainType>> solver;

    LocalSolver(std::deque<std::shared_ptr<Level>>& levels)
    {
//       //for testing purposes, we will initialise an umfpack solver for the finest level to compare the
//       //"real" solution against he solution of our multigrid
//       Level& test_level = *levels.front();
//       bool _floating = true;
//       if(test_level.filter_chain.at<1>().get_vec_prim().empty())
//         _floating = false;
//       if(_floating)
//         test_solver = std::make_shared<FloUmf>(test_level.matrix, test_level.filter_chain.at<1>());
//       else
//         test_solver = std::make_shared<RegUmf>(test_level.matrix);
      // The solver we want to set up is a "standard" geometric multigrid solver, which
      // 1) uses a simple CG solver as the coarse grid solver
      // 2) uses 4 steps of damped Jacobi as a pre-smoother and post-smoother

      // The first object that we require is a "MultiGridHierarchy" object, which is responsible
      // for managing the level hierarchy used by a multigrid solver. Note that this hierarchy
      // object is not a solver/preconditioner -- we will create that one later on.

      // The MultiGridHierarchy class templates has three template parameters and its only
      // constructor parameter is the size of the level hierarchy, i.e. the total number of
      // levels that we want to create:
      multigrid_hierarchy = std::make_shared<Solver::MultiGridHierarchy<
      LocalMatrixType,          // the system matrix type
      FilterChainType,          // the system filter type
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
        auto jacobi = Solver::new_jacobi_precond(lvl.matrix, lvl.filter_chain);

        // Create a Richardson solver for the smoother
        auto smoother = Solver::new_richardson(lvl.matrix, lvl.filter_chain, DataType(0.8), jacobi);

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
          lvl.filter_chain,     // the system filter for this level
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

//         // We use a simple (unpreconditioned) CG solver as a coarse-grid solver, so let's create one:
//         auto coarse_solver = Solver::new_pcg(lvl.matrix, lvl.filter_chain);
//         // we set an additional stopping criteria for the coarse solver, else the MG has problems...
//         coarse_solver->set_tol_abs_low(1e-12);
        //we will use umfpack as coarse solver:
        //check if we are floating, or not:
        std::shared_ptr<Solver::SolverBase<LAFEM::DenseVector<Mem::Main, double, Index>>> coarse_solver;
        if(lvl.filter_chain.at<1>().get_vec_prim().empty())
          coarse_solver = Solver::new_umfpack(lvl.matrix);
        else
          coarse_solver = Solver::new_umfpack_mean(lvl.matrix, lvl.filter_chain.at<1>());

        // At this point, we could configure the coarse grid solver, i.e. set tolerances and maximum
        // allowed iterations, but we'll just leave it at its default configuration here.

        // Now we need to attach this solver as well as the system matrix and filter to
        // our multigrid hierarchy. This is done by calling the following member function:
        multigrid_hierarchy->push_level(
          lvl.matrix,       // the coarse-level system matrix
          lvl.filter_chain,       // the coarse-level system filter
          coarse_solver     // the coarse-level solver
        );
      }

      // That's it for the coarse level.

      // Next, we need to create a multigrid preconditioner for our hierarchy. This task
      // is quite simple, as we can use one of the convenience functions to obtain a
      // MultiGrid solver object using the hierarchy we have just set up. At this point,
      // we can also choose which multigrid cycle we want to use. We have the choice
      // between V-, F- and W-cycle, but we stick to the simple V-cycle in this example:
      multigrid = Solver::new_multigrid(
        multigrid_hierarchy,          // the multigrid hierarchy object
        Solver::MultiGridCycle::V     // the desired multigrid cycle
      );

      // Alright, now we have a multigrid preconditioner object.

      // However, the 'multigrid' object we have set up so far is just a *single cycle*, i.e.
      // it is merely a preconditioner, so we still need to create a "real solver" around that.
      // As we want to have a "classical multigrid", we just need to use this multigrid object
      // as a preconditioner for a Richardson solver:
      //get ref to finest level:
      Level& lvl_fine = *levels.front();
      solver = Solver::new_richardson(lvl_fine.matrix, lvl_fine.filter_chain, DataType(1), multigrid);

      // Change the solver's plot name (that's what is written to the console) to "Multigrid".
      // This is just for the sake of aesthetics: Without this step, our solver would present
      // itself as "Richardson" ;)
      solver->set_plot_name("Multigrid");

      //last but not least init our multigrid:
      // So, initialize the multigrid hierarchy first:
      multigrid_hierarchy->init();

      // Now we can initialize our solver:
      solver->init();

      // As always, enable the convergence plot:
//       solver->set_plot_mode(Solver::PlotMode::iter);
    }

    ~LocalSolver()
    {
      //Release our solver
      solver->done();

      // Release our multigrid hierarchy; this also needs to be done separately *after*
      // the remaining solver tree was released in the call directly above.
      multigrid_hierarchy->done();
    }

    /**
     * This functions calculates the solution to fine_matrix * output = input
     * Thereby input has to be a filtered vector
     */
    void solve(LocalVectorType& output, const LocalVectorType& input, const LocalMatrixType& fine_matrix, FilterChainType& fine_filter)
    {
      //we will also test if the expected solution is calculated:
//       LocalVectorType comp = output.clone();
//       LocalVectorType output_test = output.clone();
//       test_solver->apply(output_test, input);
      Solver::solve(*solver, output, input, fine_matrix, fine_filter);
      fine_filter.filter_sol(output);
//       comp.axpy(output_test, output, -1.);
//       bool floating = fine_filter.at<1>().get_vec_prim().empty();
//       for(Index i(0); i < comp.size(); ++i)
//       {
//         if(std::abs(comp(i)) > 1e-17)
//         {
//           std::cout << "Floating?" << floating << "\n This was calculated:";
//           std::cout << "Diff vector: " << comp << "\n";
//         }
//       }
    }

  }; //class LocalSolver

  class LocalSystem
  {
  public:
    std::deque<std::shared_ptr<Level>> _levels;
    LocalMatrixType& _matrix_local;
    LocalVectorType _vec_rhs;
    LocalVectorType _vec_sol;
    LocalVectorType _R_vector;
    const bool _floating;
    std::shared_ptr<LocalSolver> _local_solver;

    LocalSystem() = delete;

    template<typename Functional_>
    LocalSystem(std::deque<std::shared_ptr<Level>>&& levels,
                Cubature::DynamicFactory & cubature_factory,
                Functional_& force_funk,
                bool is_floating)
    : _levels{std::move(levels)},
    _matrix_local{_levels.front()->matrix},
    _floating{is_floating}
    {
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // Symbolic linear system assembly
      _vec_sol = _matrix_local.create_vector_r();
      _vec_rhs = _matrix_local.create_vector_l();
      _R_vector = _matrix_local.create_vector_l();


      // Initialise the right-hand-side vector entries to zero.
      _vec_rhs.format();
      _vec_sol.format();
      _R_vector.format(1.);

      // And assemble the local rhs vector:
      Assembly::LinearFunctionalAssembler::assemble_vector(_vec_rhs, force_funk, _levels.front()->space, cubature_factory);

      // Apply the filter onto the system matrix, we of course do this for each level:
      for(auto it = _levels.begin(); it != _levels.end(); ++it)
      {
        (*it)->filter_chain.at<0>().filter_mat((*it)->matrix);
      }
      // Now we impose the boundary conditions into our linear system by applying the local filter
      // onto the local right-hand-side and initial solution vectors
      //as these only exist on our highest level, we only use the first entry of levels:
//       _levels.front()->filter_chain.filter_rhs(_vec_rhs);  //this filter seems to lead to wrong results... i guess this is
      _levels.front()->filter_chain.filter_sol(_vec_sol);    // used somewhere, where we dont want to mean filter
      _levels.front()->filter_chain.at<0>().filter_rhs(_vec_rhs);


      //initialise our local solver
      _local_solver = std::make_shared<LocalSolver>(_levels);

    }

    LocalVectorType create_rhs_buffer() const
    {
      return LocalVectorType(_vec_rhs.size());
    }

    LocalVectorType create_sol_buffer() const
    {
      return LocalVectorType(_vec_sol.size());
    }

    void apply_matrix_local(LocalVectorType & output, const LocalVectorType & input) const
    {
      _matrix_local.apply(output, input);
    }

    void apply_inverse(LocalVectorType & output, const LocalVectorType & input) const
    {
      _local_solver->solve(output, input, _matrix_local, _levels.front()->filter_chain);
    }

    DataType get_rhs_dot_R() const
    {
      if(_floating)
        return _vec_rhs.dot(_R_vector);
      else
        return 0.;
    }

    DataType* get_sol_elements()
    {
      return _vec_sol.elements();
    }

    DataType* get_rhs_elements()
    {
      return _vec_rhs.elements();
    }

    //this function solves for K input = rhs, whereby rhs is the locally assembled member variable rhs
    void solve_starting(LocalVectorType& input) const
    {
      this->apply_inverse(input, _vec_rhs);
    }

    void end_result(LocalVectorType& vec_sol_buffer, LocalVectorType& vec_rhs_buffer, const DataType alpha)
    {
      vec_rhs_buffer.axpy(vec_rhs_buffer, _vec_rhs);
      //calculate solution without kernel
      apply_inverse(vec_sol_buffer, vec_rhs_buffer);
      if(_floating)
      {
        _vec_sol.axpy(_R_vector, vec_sol_buffer, - alpha);
      }
      else
      {
        _vec_sol.axpy(_vec_sol, vec_sol_buffer);
      }
    }

    LocalVectorType& get_sol_ref()
    {
      return _vec_sol;
    }
  }; // LocalSystem

  class FETIGate
  {
  private:
    const Dist::Comm* _comm;
    const IndexType _rank;
    std::vector<int> _gate_ranks;
    std::vector<VectorMirrorType> _gate_mirrors;
    std::vector<IndexType> _mirror_dofs;
    std::vector<int> _gate_signs;
    std::vector<bool> _gate_floating;
    IndexType _gate_ranks_size = 0;
    IndexType _gate_ranks_max_size = 0;
    bool _floating = true;

  public:
    FETIGate(const Dist::Comm& comm) : _comm{&comm}, _rank{IndexType(comm.rank())}
    {
      _mirror_dofs.push_back(0);
    }


    void push(const int neighbor_rank, VectorMirrorType && neighbor_mirror)
    {
      _gate_ranks.push_back(neighbor_rank);
      _gate_mirrors.push_back(std::move(neighbor_mirror));
      _gate_signs.push_back(neighbor_rank > int(_rank) ? 1 : -1);
      _mirror_dofs.push_back(_gate_mirrors.back().num_indices() + _mirror_dofs.back());

    }

    void set_not_floating()
    {
      _floating = false;

    }

    //this function extracts the interface values in input connected to the interface shared with neighbor
    //described through the relativ index neighbor_index

    LocalVectorType extract_edge_values(const LocalVectorType& input, const IndexType neighbor_index) const
    {
      XASSERTM(input.size() == _mirror_dofs.back(), "Interface size does not match!");
      const IndexType start_size = _mirror_dofs.at(neighbor_index);
      const IndexType local_size = _mirror_dofs.at(neighbor_index +1) - start_size;
      LocalVectorType tmp(local_size);
      tmp.format();
      for(IndexType i(0); i < local_size; ++i)
      {
        tmp(i, input(start_size +i));
      }
      return tmp;
    }

    //small helper function
    void extract_all(std::vector<LocalVectorType>& output, const LocalVectorType& input) const
    {
      XASSERTM(input.size() == _mirror_dofs.back(), "Sizes do not match!");
      for(IndexType i(0); i < _gate_ranks_size; ++i)
      {
        output[i] = extract_edge_values(input, i);
      }
    }

    //extract values by performing output = output - extract(input) which is output = P * prev
    void update_by_projection(std::vector<LocalVectorType>& output, const LocalVectorType& input,
                              const std::vector<LocalVectorType>& prev)
    const
    {
      XASSERTM(input.size() == _mirror_dofs.back(), "Sizes do not match!");
      for(IndexType i(0); i < _gate_ranks_size; ++i)
      {
        output[i].axpy(extract_edge_values(input, i), prev[i], -1.);
      }
    }

    //this function scatters the residual values on the egdes regarding the sign in _gate_signs
    void scatter_residual(LocalVectorType & output, const std::vector<LocalVectorType>& res) const
    {
      //format the output
      output.format();
      for(IndexType i(0); i < _gate_ranks_size; ++i)
      {
        _gate_mirrors[i].scatter_axpy(output, res[i], _gate_signs[i]);
      }
    }

    IndexType get_rank() const
    {
      return IndexType(_rank);
    }

    IndexType get_dof_size() const
    {
      return _mirror_dofs.back();
    }

    IndexType get_number_of_domains() const
    {
      return IndexType(_comm->size());
    }

    IndexType get_gate_size() const
    {
      return _gate_ranks_size;
    }

    IndexType get_gate_max_size() const
    {
      return _gate_ranks_max_size;
    }

    bool is_floating() const
    {
      return _floating;
    }

    IndexType mirror_dofs_at(IndexType i) const
    {
      return _mirror_dofs.at(i);
    }

    int signs_at(IndexType i) const
    {
      return _gate_signs.at(i);
    }

    IndexType get_neighbour_rank(IndexType i) const
    {
      return IndexType(_gate_ranks.at(i));
    }

    bool neighbor_floating(IndexType i) const
    {
      return _gate_floating.at(i);
    }

    /**This creates a DenseVector/BufferMain of size n = 2*(neighbour_maxsize+1) representing a row of the global Q matrix, thereby
     * the first n entries represent the Indicies of the non zero entries.
     * Because we will not always have full size vectors these Indicies are encoded the following way:
     * First, always the non diagonal(so Index != domain_num) are saved in ascending order.
     * Then the Diagonal Index is saved, indicating the last vector entry
    */
    BufferMain get_Q_row() const
    {
      IndexType n = (_gate_ranks_max_size +1);
      BufferMain buff(2*n);
      if(!_floating)
      {
        buff(0, DataType(_rank));
        buff(n, 1.);
        return buff;
      }

      IndexType counter = 0;
      for(IndexType i(0); i < _gate_ranks_size; ++i)
      {
        IndexType neighbour = IndexType(_gate_ranks[i]);
        if(!_gate_floating[i])
        { //we skip non floating arrays and dont count them in our buff, so we have to keep a counter of this
          ++counter;
          continue;
        }
        //in case of trivial R the entry coupling with the neighbor is just the number of dofs on the shared
        //interface
        DataType scalar = DataType(_mirror_dofs[i+1] - _mirror_dofs[i]);
        buff(i - counter, DataType(neighbour));
        buff(i - counter + n, -scalar);
      }
      //here the entry is just the number of dofs of the local interface
      buff(_gate_ranks_size - counter, DataType(_rank));
      buff(_gate_ranks_size - counter + n, DataType(_mirror_dofs.back()));
      return buff;

    }

    //this function synchronizes the values in the local vectors input by extracting the values through
    //the mirrors and exchaning those with its neighbors, whereby
    //output[i] = sign*(gate_sign[i]*this->neighbor_input_reduced - gate_sign[i]*neighbor->input_reduced)
    //whereby i denotes a possible edge between two domains and input_reduced is the vector remaining
    //after gathering with mirror
    void apply_B_matrix(std::vector<LocalVectorType>& output, const LocalVectorType& input, int sign = 1) const
    {
      XASSERTM(output.size() == _gate_ranks_size, "Output vector does not match neighbor number!");

      std::size_t n = std::size_t(_gate_ranks_size);
      /// send and receive request vectors
      Dist::RequestVector send_reqs, recv_reqs;
      /// send and receive buffers
      std::vector<BufferMain> send_bufs, recv_bufs;

      // post receives
      recv_reqs.reserve(n);
      recv_bufs.resize(n);

      for(std::size_t i(0); i < n; ++i)
      {
        // create buffer vector in main memory
        recv_bufs.at(i) = BufferMain(_gate_mirrors.at(i).buffer_size(input), LAFEM::Pinning::disabled);
        //buffer size right?

        // post receive
        recv_reqs.push_back(_comm->irecv(recv_bufs.at(i).elements(), recv_bufs.at(i).size(), _gate_ranks.at(i)));
      }

      // post sends
      send_reqs.reserve(n);
      send_bufs.resize(n);

      std::vector<BufferType> buffer_vec;
      for(std::size_t i(0); i < n; ++i)
      {
        // create buffer
        send_bufs.at(i) = BufferMain(_gate_mirrors.at(i).buffer_size(input), LAFEM::Pinning::disabled);
        // gather from mirror
        _gate_mirrors.at(i).gather(send_bufs.at(i), input);
        // post send
        send_reqs.push_back(_comm->isend(send_bufs.at(i).elements(), send_bufs.at(i).size(), _gate_ranks.at(i)));
      }

      // process all pending receives
      for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
      {
        XASSERTM(_gate_signs.at(idx) != 0, "Gate signs are not correctly initialised!");

        //process recv_bufs depending on gate_signs
        if(_gate_signs[idx]*sign > 0)
        {
          output[idx].axpy(recv_bufs.at(idx), send_bufs.at(idx), -1.);
        }
        else
        {
          output[idx].axpy(send_bufs.at(idx), recv_bufs.at(idx), -1.);
        }

      }

      // wait for all sends to finish
      send_reqs.wait_all();
    }

    //This function initialises the starting residual in res by using vec_sol_buffer as a buffer
    void initialise_starting_residual(std::vector<LocalVectorType>& res, LocalVectorType& vec_sol_buffer,
                                      const LocalSystem& system)
    const
    {
      //calculate the local solution
      system.solve_starting(vec_sol_buffer);
      //here negativ sign, because i have done this in test Code...
      apply_B_matrix(res, vec_sol_buffer, -1);
    }

    //test given input vector on whether the values are the same on the same interface across processes
    void test_vector(const std::vector<LocalVectorType>& vec) const
    {
      std::size_t n = std::size_t(_gate_ranks_size);
      /// send and receive request vectors
      Dist::RequestVector send_reqs, recv_reqs;
      /// send and receive buffers
      std::vector<BufferMain> send_bufs, recv_bufs;

      // post receives
      recv_reqs.reserve(n);
      recv_bufs.resize(n);

      for(std::size_t i(0); i < n; ++i)
      {
        // create buffer vector in main memory
        recv_bufs.at(i) = BufferMain(vec.at(i).size(), LAFEM::Pinning::disabled);

        // post receive
        recv_reqs.push_back(_comm->irecv(recv_bufs.at(i).elements(), recv_bufs.at(i).size(), _gate_ranks.at(i)));
      }

      // post sends
      send_reqs.reserve(n);
      send_bufs.resize(n);

      for(std::size_t i(0); i < n; ++i)
      {
        send_bufs.at(i) = BufferMain(vec.at(i).size(), LAFEM::Pinning::disabled);
        send_bufs.at(i).copy(vec.at(i));

        // post send
        send_reqs.push_back(_comm->isend(send_bufs.at(i).elements(), send_bufs.at(i).size(), _gate_ranks.at(i)));
      }
      //the tolerance that should be met... hardcoded because lazy
      const DataType tol_eps = 1e-15;

      // process all pending receives
      for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
      {

        XASSERTM(recv_bufs.at(idx).size() == vec.at(idx).size(), "Sizes do not match!");
        for(IndexType i(0); i < recv_bufs.at(idx).size(); ++i)
        {
          XASSERTM(Math::abs(recv_bufs.at(idx)(i) - vec.at(idx)(i)) < tol_eps , "Values not the same on process"
           + stringify(_comm->rank()) + "Value in recv: " + stringify(recv_bufs.at(idx)(i))
           + " Value in vec: " + stringify(vec.at(idx)(i)));
        }
      }

      // wait for all sends to finish
      send_reqs.wait_all();

    }

    //This function applies B*K-1B^T onto input and saves it into output
    //sol_buffer is used as tempoary vector to comunicate across processes
    //rhs_buffer is used for extracting the edge values
    void calc_sol(std::vector<LocalVectorType>& output, const std::vector<LocalVectorType>& input,
                     const LocalSystem& system, LocalVectorType& sol_buffer, LocalVectorType& rhs_buffer)
    const
    {
      //scatter the edge values into rhs_buffer:
      scatter_residual(rhs_buffer, input);
      //apply local solver
      system.apply_inverse(sol_buffer, rhs_buffer);
      //apply B onto sol_buffer and write it into output
      apply_B_matrix(output, sol_buffer);

    }

    //performs output = BKB^T * input and uses sol_buffer and rhs_buffer as buffers
    void precond(std::vector<LocalVectorType>& output, const std::vector<LocalVectorType>& input,
                     const LocalSystem& system, LocalVectorType& sol_buffer, LocalVectorType& rhs_buffer)
    const
    {
      //scatter the edge values into rhs_buffer:
      scatter_residual(rhs_buffer, input);
      //apply local solver
      system.apply_matrix_local(sol_buffer, rhs_buffer);
      //apply B onto sol_buffer and write it into output
      apply_B_matrix(output, sol_buffer);
    }

    //this function collect the local entry of (BR)^T * input by applying the dot product of input and (1,..,1)^T
    //times the sign in _gate_signs
    DataType gather_R_dot(const std::vector<LocalVectorType>& input) const
    {
      XASSERTM(input.size() == _gate_ranks_size, "Sizes do not match!");
      DataType dot{0};
      if(!_floating)
      {
        return dot; //return 0 if not floating
      }
      //just add up sum(input[i]) * gate_signs[i], as we take R = (1,...,1)^T
      for(IndexType i{0}; i < _gate_ranks_size; ++i)
      {
        DataType tmp{0};
        for(IndexType j{0}; j < input[i].size(); ++j)
        {
          tmp += input[i](j);
        }
        dot += _gate_signs[i]*tmp;
      }
      return dot;
    }

    void end_result(const std::vector<LocalVectorType>& lambda, const DataType alpha, LocalSystem& system,
                    LocalVectorType& vec_sol_buffer, LocalVectorType& vec_rhs_buffer)
    const
    {
      scatter_residual(vec_rhs_buffer, lambda);
      system.end_result(vec_sol_buffer, vec_rhs_buffer, alpha);
    }


    void finalize(/*const bool verbose = false*/)
    {
      // first of all we want to get the maximum number of neighbours over all domains
      // for this we will use iallreduce funtion from dist
      //init buffers, x snedbuf, r recvbuf
      DataType r{0}, x{0};
      _gate_ranks_size = _gate_ranks.size();
      x = DataType(_gate_ranks_size);

      //now reduce the values through the maximum operation
      _comm->allreduce(&x, &r, std::size_t(1), Dist::op_max);
      //now save the value
      _gate_ranks_max_size = IndexType(r);

      XASSERTM(_gate_ranks_size <= _gate_ranks_max_size, "Maxsize is not a maxsize!");

//       //now each prozess should print its value
//       if(verbose)
//       {
//         String msg = "Number of Neigbors of process " + stringify(_comm->rank()) + ":" + " " + stringify(_gate_ranks_size)
//         + " Max Number is " + stringify(_gate_ranks_max_size);
//         _comm->allprint(msg);
//       }

      //now we want to exchange the floating status with the neighbors
      Dist::RequestVector _send_reqs, _recv_reqs;
      std::vector<int> _send_bufs, _recv_bufs;
      const std::size_t n = _gate_ranks_size;

      //post reicives
      _recv_reqs.reserve(n);
      _recv_bufs.resize(n);
      for(std::size_t i(0); i<n; ++i)
      {
        // post receive
        _recv_reqs.push_back(_comm->irecv(&_recv_bufs.at(i), std::size_t(1), _gate_ranks.at(i)));
      }

      // post sends
      _send_reqs.reserve(n);
      _send_bufs.resize(n);
      for(std::size_t i(0); i < n; ++i)
      {
        // gather own floating status
        _send_bufs.at(i) = (_floating ? 1 : 0);
        // post send
        _send_reqs.push_back(_comm->isend(&_send_bufs.at(i), std::size_t(1), _gate_ranks.at(i)));
      }

      //write recv buffs into _gate_floating
      XASSERTM(_gate_floating.size() == 0, "_gate_floating already has values!");
      //resize for right size
      _gate_floating.resize(n);
      for(std::size_t idx(0u); _recv_reqs.wait_any(idx); )
      {
        //we will just write the values directly into target vector
        _gate_floating.at(idx) = (_recv_bufs.at(idx) != 0);
      }

      //wait until all sends finish
      _send_reqs.wait_all();

//       if(verbose){
//         String msg_2 = "Floating status of process " + stringify(_comm->rank()) + "\n";
//         for(IndexType i(0); i < _gate_ranks.size(); ++i)
//         {
//           msg_2 += " Neighbor " + stringify(_gate_ranks[i]) + " is ";
//           msg_2 += _gate_floating[i] ? "floating" : "not floating";
//           msg_2 += " || ";
//         }
//         _comm->allprint(msg_2);
//       }
    }
  }; //class FETIGate

  class Projector
  {
  private:
    const Dist::Comm* _comm;
    LocalMatrixType BR_columns;
    //Q is only non empty if we are on 0 rank
    LocalMatrixType Q_matrix;
    std::shared_ptr<RegUmf> Q_factorised = nullptr;
    mutable LocalVectorType right_side;
    mutable LocalVectorType right_side_buffer;

  public:
    Projector() = default;
    Projector(const FETIGate& gate, const Dist::Comm & comm) : _comm{&comm}
    {
      //assemble BR_columns
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
      //This does not need any information exchange, as we already gathered the information about the floatingness
      //of the neighbour domains
      //for now we just use SparseMatrixFactory to assemble the BR_columns matrix
      LAFEM::SparseMatrixFactory<DataType, IndexType> factory(gate.get_dof_size(), gate.get_number_of_domains());
      if(gate.is_floating())
      {
        for(IndexType i(0); i < gate.get_gate_size(); ++i)
        {
          for(IndexType j(gate.mirror_dofs_at(i)); j < gate.mirror_dofs_at(i+1); ++j)
            factory.add(j, gate.get_rank(), double(gate.signs_at(i)));
        }
      }

      for(IndexType i(0); i < gate.get_gate_size(); ++i)
      {
        if(gate.neighbor_floating(i))
        {
          for(IndexType j(gate.mirror_dofs_at(i)); j < gate.mirror_dofs_at(i+1); ++j)
            factory.add(j, gate.get_neighbour_rank(i), -double(gate.signs_at(i)));
        }
      }
      BR_columns = factory.make_csr();
    }

    void apply_BR_columns(LocalVectorType& output, const LocalVectorType& input) const
    {
      //test if vectors have right size_t
      XASSERTM(input.size() == BR_columns.columns(), "Output vector does not have the right size!");
      XASSERTM(output.size() == BR_columns.rows(), "Output vector does not have the right size!");
      BR_columns.apply(output, input);
    }

    //This function gathers the local scalar, applies Q on the resulting vector and broadcast the result
    //into the local right_side member
    void global_Q_apply(DataType local_scalar) const
    {
      //post the request
      _comm->gather(&local_scalar, std::size_t(1), right_side_buffer.elements(), std::size_t(1), 0);

      //apply Q on right_side on rank 0
      if(_comm->rank() == 0)
      {
        Q_factorised->apply(right_side, right_side_buffer);
      }
      //now broadcast right_side to all processes
      _comm->bcast(right_side.elements(), right_side.size(), 0);
    }

    //this function gathers the Q_matrix onto our rank 0 prozess by locally assembling the rows
    //Q will have 1 entries on the diagonal for non floating domains
    void gather_Q(const FETIGate& gate, const bool verbose = false)
    {
      //assemble the rows of the Q-matrix into our sendbuf
      IndexType n = gate.get_gate_max_size()+1;
      if(verbose)
        _comm->print("Get single Q row...");
      BufferMain send_buf = gate.get_Q_row();
      //our recv_buffer
      BufferMain recv_buf;

      //now we will allocate a buffer of size  2*n * _comm.size() on our core prozess 0
      if(_comm->rank() == 0)
      {
        recv_buf = BufferMain(2*n*(IndexType(_comm->size())));
      }
      //now post the igather request
      _comm->gather(send_buf.elements(), std::size_t(2*n), recv_buf.elements(), std::size_t(2*n), 0);
      if(verbose)
        _comm->print("Construct Q on rank 0...");
      //now construct Q on rank 0, the other prozesses could be send to sleep
      if(_comm->rank() == 0)
      {
        IndexType dom_size = IndexType(_comm->size());
        //setup factory
        LAFEM::SparseMatrixFactory<DataType, IndexType> factory(dom_size, dom_size);
        //iterate over all processes, whereby each process wrote into a 2*n chunk in recv_buffer
        //ordered by its rank
        for(IndexType i(0); i < dom_size ; ++i)
        {
          //thereby, the data is orgenized in two n size chunks, holding index and data off the Q_row entries
          for(IndexType j(0); j < n; ++j)
          {
            factory.add(i, IndexType(recv_buf(i*2*n + j)), recv_buf(i*2*n + n + j));
            //if we read in the diagonal element, we stop, as rest of the data is junk
            if(IndexType(recv_buf(i*2*n + j)) == i)
              break;
          }
        }

        if(verbose)
          _comm->print("Factory call...");
        //after reading in the entries, build the matrix and save it
        Q_matrix = factory.make_csr();
        //initialize some more
        if(verbose)
          _comm->print("Factorize Q_matrix...");
        //factorise the Q_matrix
        Q_factorised = std::make_shared<RegUmf>(Q_matrix);
      }
      //as we possibly need right side on all ranks, init them on all processes
      if(verbose)
        _comm->print("Create right sides...");
      right_side_buffer = LocalVectorType(IndexType(_comm->size()));
      right_side = LocalVectorType(IndexType(_comm->size()));
      right_side_buffer.format();
      right_side.format();
      if(verbose)
        _comm->print("Q_matrix succesfully gathered and initialised!");

    }

    LocalVectorType gather_starting_right_side(const FETIGate& gate, const LocalSystem& system) const
    {
      //negative sign because it makes no difference and this is the way i handleld it before
      DataType send_buf = - system.get_rhs_dot_R();
      //gather, apply Q, broadcast values to right_side
      global_Q_apply(send_buf);

      //apply BR on right side
      LocalVectorType tmp(gate.get_dof_size());
      apply_BR_columns(tmp, right_side);
      return tmp;
    }

    void project(std::vector<LocalVectorType>& output, const std::vector<LocalVectorType>& input,
                 const FETIGate& gate)
    const
    {
      DataType tmp_scalar = gate.gather_R_dot(input);
      //gather, apply Q, broadcast values to right_side
      global_Q_apply(tmp_scalar);
      //apply BR on a tmp vector
      LocalVectorType tmp(gate.get_dof_size());
      apply_BR_columns(tmp, right_side);
      //and now distribute the data back into the output
      gate.update_by_projection(output, tmp, input);
    }

    DataType calc_alpha(const std::vector<LocalVectorType>& input, const FETIGate& gate) const
    {
      DataType tmp_scalar = gate.gather_R_dot(input);

      //gather, apply Q, broadcast values to right_side
      global_Q_apply(tmp_scalar);

      //return the value connected to this rank...
      return right_side(IndexType(_comm->rank()));
    }


  }; //class Projector

  class FETICG
  {
  private:
    const Dist::Comm* _comm;
    std::vector<LocalVectorType> residual;
    std::vector<LocalVectorType> residual_copy;
    std::vector<LocalVectorType> residual_opt;
    std::vector<LocalVectorType> rescond;
    std::vector<LocalVectorType> rescond_opt;
    std::vector<LocalVectorType> Fs_vector;
    std::vector<LocalVectorType> s_vector;
    LocalVectorType vec_sol_buffer;
    LocalVectorType vec_rhs_buffer;
    std::vector<LocalVectorType> lambda;
    const IndexType vec_size;
    DataType gamma = 0.;
    DataType beta = 0.;
    const DataType CG_tol;
    const IndexType maxiter = 1000u;
    DataType alpha = 0.;

  public:
    FETICG(const Dist::Comm & comm, const FETIGate& gate, const LocalSystem& system, DataType tol = 1e-10,
           IndexType niter = 1000u)
    : _comm{&comm},
      vec_size{gate.get_gate_size()},
      CG_tol{tol},
      maxiter{niter}
    {
      IndexType n = gate.get_gate_size();
      for(IndexType i(0); i < n; ++i)
      {
        const IndexType loc_size = gate.mirror_dofs_at(i+1) - gate.mirror_dofs_at(i);
        residual.push_back(LocalVectorType(loc_size, 0.));
        residual_copy.push_back(LocalVectorType(loc_size, 0.));
        residual_opt.push_back(LocalVectorType(loc_size, 0.));
        rescond.push_back(LocalVectorType(loc_size, 0.));
        rescond_opt.push_back(LocalVectorType(loc_size, 0.));
        Fs_vector.push_back(LocalVectorType(loc_size, 0.));
        s_vector.push_back(LocalVectorType(loc_size, 0.));
        lambda.push_back(LocalVectorType(loc_size, 0.));
      }
      vec_sol_buffer = system.create_sol_buffer();
      vec_rhs_buffer = system.create_rhs_buffer();

      vec_sol_buffer.format();
      vec_rhs_buffer.format();
    }

    //this function gathers the dot_product of vector and other across all domains and broadcasts the resulting scalar back
    DataType synch_dot(std::vector<LocalVectorType>& vector,  std::vector<LocalVectorType>& other, bool sqrt = false) const
    {
      XASSERTM(vector.size() == other.size(), "Sizes do not match!");
      DataType x{0}, r{0};
      //gather the dot_product
      for(IndexType i(0); i < vector.size(); ++i)
      {
        XASSERTM(vector[i].size() == other[i].size(), "Sizes do not match!");
        x += vector[i].dot(other[i]);
      }
      //post request-> we want to sum all gathered data
      _comm->allreduce(&x, &r, std::size_t(1), Dist::op_sum);

      //we need to half the value, as each edge was counted twice
      r = r/2.;
      //return value depending wether we want to take square root or not
      return (sqrt ?  Math::sqrt(r) : r);
    }

    //performs target = target + scalar * other
    void update_vector(std::vector<LocalVectorType>& target, std::vector<LocalVectorType>& other, DataType scalar)
    {
      for(IndexType i(0); i < vec_size; ++i)
      {
        target[i].axpy(other[i], target[i], scalar);
      }
    }

    void update_swapped(std::vector<LocalVectorType>& target, std::vector<LocalVectorType>& other, DataType scalar)
    {
      for(IndexType i(0); i < vec_size; ++i)
      {
        target[i].axpy(target[i], other[i], scalar);
      }
    }

    void init(const FETIGate& gate, const LocalSystem& system, const Projector& project)
    {
      _comm->print("Initialising CG-Algorithm...");
      //assemble starting lambda(as this needs projection, Projector handles this)
      const LocalVectorType tmp = project.gather_starting_right_side(gate, system);
      gate.extract_all(lambda, tmp);

      //assemble starting residual(as this only needs exchange between neighbors, Gate will handle this)
      gate.initialise_starting_residual(residual, vec_sol_buffer, system);

      //make a copy of residual for later
      for(IndexType i(0); i < vec_size; ++i)
      {
        residual_copy.at(i).copy(residual.at(i));
      }
//       _comm->print("Residuum initialised!");

      //generate starting Fs = BK^-1B^T * lam
      gate.calc_sol(Fs_vector, lambda, system, vec_sol_buffer, vec_rhs_buffer);


      //udate residual by res = res - Fs
      update_vector(residual, Fs_vector, -1.);


      //project the residual
      project.project(residual_opt, residual, gate);


      //apply precond and save this into rescond
      gate.precond(rescond, residual_opt, system, vec_sol_buffer, vec_rhs_buffer);


      //project the rescond
      project.project(rescond_opt, rescond, gate);


      //init s_vector = rescond_opt
      init_s_vector();

    }

    //Calculate lambda through a CG-Method, we stop if the 2-norm of the relativ projected residual is below our tol
    //after this, lambda will no longer be a valid LocalVectorType as we move it
    bool solve(const FETIGate& gate, const LocalSystem& system, const Projector& project, const bool verbose = false,
                  const bool want_test = false)
    {
      StopWatch watch_proj;
      StopWatch watch_solv;
      //init some needed data
      if(verbose)
        _comm->print("Gathering a bit more data...");
      DataType r_0 = 0.;
      DataType s_l = 0.;
      watch_solv.start();
      DataType r_1 = synch_dot(rescond_opt, residual_opt);
      DataType nrOopt = synch_dot(rescond_opt, rescond_opt, true);
      watch_solv.stop();
      DataType r_alt = nrOopt;
      XASSERTM(r_alt > 0., "Norm is zero or negativ");
      DataType L2_defect = 0.;

      _comm->print("Starting the CG-Algorithm...");

      for(IndexType k(0); k < maxiter ; ++k)
      {
        //first step calculate intermediate value Fs = BK^{-1}B^T *s
        watch_solv.start();
        gate.calc_sol(Fs_vector, s_vector, system, vec_sol_buffer, vec_rhs_buffer);
        //now calculate dot(s,Fs)
        s_l = synch_dot(s_vector, Fs_vector);
        //update gamma
        gamma = r_1/s_l;

        //update lambda by lam = lam + s*gamma
        update_vector(lambda, s_vector, gamma);

        //update residual : r = r - Fs*gamma
        update_vector(residual, Fs_vector, -gamma);
        watch_solv.stop();
        //project residual
        watch_proj.start();
        project.project(residual_opt, residual, gate);
        watch_proj.stop();
        //apply precond and save this into rescond
        watch_solv.start();
        gate.precond(rescond, residual_opt, system, vec_sol_buffer, vec_rhs_buffer);
        watch_solv.stop();
        //project the rescond
        watch_proj.start();
        project.project(rescond_opt, rescond, gate);
        watch_proj.stop();

        //update r_0 and r_1
        watch_solv.start();
        r_0 = r_1;
        r_1 = synch_dot(rescond_opt, residual_opt);
        //calculate beta
        beta = r_1/r_0;
        //update s_vector : s = rescond_opt + s*beta
        update_swapped(s_vector, rescond_opt, beta);

        //calculate standard projected res_error
        r_alt = synch_dot(residual_opt, residual_opt, true);
        //calculate relative error:
        L2_defect = r_alt/nrOopt;
        watch_solv.stop();
        if(verbose)
          _comm->print("CG-Iteration: " + stringify(k) + " | rel l2_defect: " + stringify(L2_defect));

        //stopping criteria
        if(L2_defect < CG_tol)
        {
          if(verbose)
            _comm->print("Finished on " + stringify(k) + "-th iteration!");
          IndexType max_i = 7;
          if(want_test && k >= max_i)
          {
            _comm->print("Test FAILED: Expected less then " + stringify(max_i)
                                 + " Iterations, needed " + stringify(k) + " Iterations.");
            return false;
          }
          break;
        }
      }
      DataType t_proj = watch_proj.elapsed();
      DataType t_solv = watch_solv.elapsed();
      DataType t_max = 0;
      _comm->allreduce(&t_proj, &t_max, std::size_t(1), Dist::op_max);
      _comm->print("CG... Time for projection: " + stringify(t_max) + " s");
      _comm->allreduce(&t_solv, &t_max, std::size_t(1), Dist::op_max);
      _comm->print("CG... Time for solving, includes inter neighbour comm: " + stringify(t_max) + " s");
      return true;
    }

    //This initialises s_0 = h_opt
    void init_s_vector()
    {
      XASSERTM(this->s_vector.size() == this->rescond_opt.size(), "Sizes do not match!");
      for(IndexType i(0); i < this->rescond_opt.size(); ++i)
      {
        this->s_vector[i].copy(this->rescond_opt[i]);
      }
    }

    void calc_alpha(const FETIGate& gate, const LocalSystem& system, const Projector& project)
    {
      //first step calculate intermediate value Fs = BK^{-1}B^T *lambda
        gate.calc_sol(Fs_vector, lambda, system, vec_sol_buffer, vec_rhs_buffer);
      // update residual_copy = Fs - residual_copy
        update_swapped(residual_copy, Fs_vector, -1.);

        //and now calculate local alpha through projector:
        alpha = project.calc_alpha(residual_copy, gate);
    }

    void end_result(const FETIGate& gate, LocalSystem& system)
    {
      _comm->print("Calculating end result...");
      gate.end_result(lambda, alpha, system, vec_sol_buffer, vec_rhs_buffer);
    }

    //testing
    void test_vector(const std::vector<LocalVectorType>& vec, const FETIGate& gate) const
    {
      gate.test_vector(vec);
    }
  }; //class FETICG


  void main(int argc, char* argv[])
  {
    // The very first thing that we need is a so-called "distributed communicator", which
    // represents a communication interface to the set of all processes that have been started
    // for this tutorial. We won't be using the communicator object for any inter-process
    // communication directly in this tutorial, but this communicator will be used "under-the-hood"
    // by all objects that require global communication.
    // The corresponding class is the Dist::Comm class, which is effectively just a wrapper around
    // a MPI_Comm object. We create a new comm object representing the "world" communicator:
    Dist::Comm comm = Dist::Comm::world();

    StopWatch watch;

    //if we only have 1 prozess, we want to abort, as this will generate errors later down the line
    XASSERTM(comm.size() > 1, "Started with only 1 process, please run with more than one process!");

    //first we will parse the input argument
    SimpleArgParser args(argc, argv);

    // So, let's start adding our supported options by calling the parser's 'support'
    // function: The first argument is the name of the option *without* the leading
    // double-hyphen, whereas the second argument is a short description used for the
    // formatting of the argument list description.
    args.support("help", "\nDisplays this help information.\n");
    args.support("level", "<n>\nSets the mesh refinement level.\n");
    args.support("CGplot", "\nDisplay the convergence plot of the CG method.\n");
    args.support("CGtol", "<d>\nSpecifies the tolerance for the stopping criteria of the internal CG-algorithm.\n"
                         "Default value is 1e-8. If set too high could interfere with convergence\n");
    args.support("verbose", "\nDisplay more information about the programm flow.\n");
    args.support("vtk", "<filename>\nSpecifies the filename for the VTK exporter.\n"
      "If this option is not specified, no VTK file will be written\n");
    args.support("CGmaxiter", "<n>\n The maximum iteration step for the Cg-method.\n"
                              "Defaults to 1000");
    args.support("test", "\nStarts in test mode, which disables failing due to insufficient number of process,"
      " while running a predefined problem to compare it to an expected behaviour.\n");

    // We will keep track of whether the caller needs help, which is the case when an
    // invalid option has been supplied or '--help' was given as a command line argument.
    // In this case, we will write out some help information to std::cout later.
    bool need_help = false;

    // Now that we have added all supported options to the parser, we call the 'query_unsupported'
    // function to check whether the user has supplied any unsupported options in the command line:
    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();

    // If the caller did not specify any unsupported (or mistyped) options, the container returned
    // by the query_unsupported function is empty, so we first check for that:
    if(!unsupported.empty())
    {
      // Okay, we have at least one unsupported option.
      // Each entry of the container contains a pair of an int specifying the index of the
      // faulty command line argument as well as the faulty argument itself as a string.
      // We loop over all unsupported arguments and print them out:
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      {
        comm.print(std::cerr, "ERROR: unsupported option #" + stringify((*it).first) + " '--" + stringify((*it).second) + "'");
      }

      // We could abort program execution here, but instead we remember that we have to print
      // the help information containing all supported options later:
      need_help = true;
    }

    // Check whether the caller has specified the '--help' option:
    need_help = need_help || (args.check("help") >= 0);

    // Print help information?
    if(need_help)
    {
      // Okay, let's display some help.
      comm.print("\n");
      comm.print("USAGE: " + stringify(args.get_arg(0)) + " <options> \n");

      // Now we call the 'get_supported_help' function of the argument parser - this
      // will give us a formatted string containing the supported options and their
      // corresponding short descriptions  which we supplied before:
      comm.print("Supported options:");
      comm.print(args.get_supported_help());

    // We abort program execution here:
    return;
    }

    // The first difference in this tutorial is that we do not print console output directly to
    // std::cout, but use the communicator's "print" function instead, which ensures that only
    // one single process (in the communicator) prints the output to the console:
    comm.print("Welcome to FEAT's FETI solver.");

    //check for the simple false true options
    bool CG_plot = (args.check("CGplot") >= 0);
    bool verbose = (args.check("verbose") >= 0);
    bool want_test = (args.check("test") >= 0);
    bool test_passed = true;

    // Our mesh refinement level
    Index level(3);
    // The filename of our VTK file
    String vtk_name("");
    //Our CG_tolerance
    DataType CG_tol(1e-8);
    //Our max iter for the CG algo
    IndexType CG_max_iter(1000u);

    // Now let's start parsing the command line arguments.
    // For this, we call the 'parse' function of the argument parser and supply
    // the name of the desired option as well as the variable(s) to be parsed.
    // This function returns an int which specifies how many parameters have
    // been parsed successfully or which parameter was not parsed, i.e.
    // if the return value is
    // * = 0    , then either the option was not given at all or it was given
    //            but without any parameters.
    // * = n > 0, then the option was given and the first n parameters were
    //            parsed successfully.
    // * = n < 0, if the option was given, but the (-n)-th command line argument
    //            could not be parsed into the variable supplied to the function

    // We first try to parse the option '--level' and then check the return value:
    int iarg_level = args.parse("level", level);
    if(iarg_level < 0)
    {
      // In this case, we have an error, as the corresponding command line
      // argument could not be parsed, so print out an error message:
      comm.print(std::cerr, "ERROR: Failed to parse '" + stringify(args.get_arg(-iarg_level)) + "'");
      comm.print(std::cerr, "as parameter for option '--level");
      comm.print(std::cerr, "Expected: a non-negative integer");

      // and abort our program
      Runtime::abort();
    }

    int iarg_tol = args.parse("CGtol", CG_tol);
    if(iarg_level < 0)
    {
      // In this case, we have an error, as the corresponding command line
      // argument could not be parsed, so print out an error message:
      comm.print(std::cerr, "ERROR: Failed to parse '" + stringify(args.get_arg(-iarg_tol)) + "'");
      comm.print(std::cerr, "as parameter for option '--CGtol");
      comm.print(std::cerr, "Expected: a non-negative double");

      // and abort our program
      Runtime::abort();
    }

    int iarg_iter = args.parse("CGmaxiter", CG_max_iter);
    if(iarg_level < 0)
    {
      // In this case, we have an error, as the corresponding command line
      // argument could not be parsed, so print out an error message:
      comm.print(std::cerr, "ERROR: Failed to parse '" + stringify(args.get_arg(-iarg_iter)) + "'");
      comm.print(std::cerr, "as parameter for option '--CGmaxiter");
      comm.print(std::cerr, "Expected: a non-negative integer");

      // and abort our program
      Runtime::abort();
    }

    // Next, we check for the '--vtk' option, which specifies the filename of the
    // VTK file to be written. If this option is not given or if it has no parameters,
    // then we will not write a VTK file at all, so we just check whether the return
    // value of the parse function is positive:
    bool want_vtk = (args.parse("vtk", vtk_name) > 0);

    //if we want to test, overwrite all given parameters
    if(want_test)
    {
      level = 4;
      want_vtk = false;
      CG_tol = 1e-8;
      CG_max_iter = 1000u;
      CG_plot = false;
      verbose = false;
    }

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
      if(want_test)
      {
        //End the test as successful, small hack, because test functions only with multiple of 4...
        comm.print("Test PASSED: Test started with " + stringify(comm.size()) + " processes, ending test prematurely!");
        return;
      }
      // We pass std::cerr as the first parameter to the 'print' function here:
      comm.print(std::cerr, "ERROR: You must run this tutorial with 1, 4, 16 or 64 processes!");
      Runtime::abort();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Geometry initialization

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
    // but it will also tell us the ranks of all processes that manage our neighbor patches.
    // We will require these ranks to set up the communication "mirrors" that are required for
    // the global simulation. For this, we need to create a std::vector of ints, which will be
    // filled with our neighbor process ranks:
    std::vector<int> neighbor_ranks;

    // Now we can call our generator to obtain our patch mesh-node as well as our neighbor ranks.
    // Moreover, the create function returns an index that corresponds to the refinement level
    // of the global unit-square domain - we will require this for the further refinement below.
    Index lvl = Geometry::UnitCubePatchGenerator<MeshType>::create(
      comm.rank(),          // input:  the rank of this process
      comm.size(),          // input:  the total number of processes
      root_mesh_node,       // output: the root-mesh-node shared pointer
      neighbor_ranks);     // output: the neighbor ranks vector

    // At this point, we have our root mesh node as well as the neighbor ranks.

    // As mentioned before, the generator returned a level index that represents the
    // refinement level of the mesh-node that it has generated. We want to print that out:
    comm.print("\nBase Mesh Level: " + stringify(lvl));

    //As we want to use a multigrid, we will initalise the mesh, space and trafo for each level,
    //for this we will create a deque of the class Level, which later will be moved to the LocalSystem
    std::deque<std::shared_ptr<Level>> levels;

    //As we espacially need the root_mesh_node on our highest level, we use shared_ptr to also use it outside of the Level class


    // Now let's refine the mesh up to the level that was passed as a parameter to this function,
    // assuming that it is greater than the base-mesh level the generator gave us:
    IndexType level_min = lvl;
    IndexType level_max = level_min;
    if(lvl < level)
    {
      comm.print("Refining Mesh to Level " + stringify(level) + "...");
      watch.start();
      for(; lvl < level; ++lvl)
      {
        levels.push_front(std::make_shared<Level>(root_mesh_node));
        root_mesh_node = std::shared_ptr<RootMeshNodeType>(root_mesh_node->refine());
        ++level_max;
      }
      watch.stop();
    }
    //add the current mesh (this is the finest level to the deque)
    levels.push_front(std::make_shared<Level>(root_mesh_node));

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
    // Trafo and Space initialization

    comm.print("Creating Trafo and Space...");

    // The creation of a transformation and a finite element space is identical to the
    // previous tutorials, we just need to get the mesh representing our patch of the domain
    // from the mesh-node, as we created those in Level class, we just use refernces to those:
    MeshType& mesh = levels.front()->mesh;

    // Create the desired finite element space
    SpaceType& space = levels.front()->space;

    // Note that both the transformation and the finite element space are defined merely on
    // the local patch, i.e. the trafo and the space do not know about the fact that we are
    // running a parallel simulation. All global interaction and communication is performed
    // only on the linear algebra level, which are going to set up in the next few steps.

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //We use a pseudo-gate class
    watch.start();
    FETIGate gate(comm);
    watch.stop();

    // For later anaylsis we also init the "standard" GlobalGate
    GateType global_gate(comm);

    for(auto it = neighbor_ranks.begin(); it != neighbor_ranks.end(); ++it)
    {
      // Get the rank of our neighbour process:
      const int neighbor_rank = (*it);

      // As already mentioned before, our mesh-node does not only contain the mesh itself,
      // but also other important stuff, such as the "halos", which describe the overlap of
      // neighboured patches. For the assembly of the mirror, we need to get the halo mesh-part
      // from the mesh-node first:
      const MeshPartType* neighbor_halo = root_mesh_node->get_halo(neighbor_rank);

      // Ensure that we have a halo for this neighbour rank:
      XASSERTM(neighbor_halo != nullptr, "Failed to retrieve neighbour halo!");
      //we create the same vector here and give it to our global_gate.... not needed for FETI itself
      VectorMirrorType mirror_alt;
      // Call the MirrorAssembler to do the dirty work for us:
      Assembly::MirrorAssembler::assemble_mirror(
        mirror_alt,   // the mirror that is to be assembled
        space,              // the FE space for which we want to assemble the mirror
        *neighbor_halo     // the halo mesh-part that the mirror is to be assembled on
      );

      // Once the mirror is assembled, we give it over to our gate.
      global_gate.push(
        neighbor_rank,               // the process rank of the neighbor
        std::move(mirror_alt)   // the mirror for the neighbor
      );

      //Ensure the halo does not contain only one dof, this is for now only for 2D meshes, for 3D meshes edges also have to be checked:
      // Question: I want to check how many elements are in vertices, ie 1 dimensional points. If i want to check if there is only one edge
      // this can be taken out... of course we would then have more communication needed
      if(neighbor_halo->get_num_entities(0) == IndexType(1))
      {
        continue;
      }

      watch.start();
      // Now that we have the halo mesh-part, we can create and assemble the corresponding mirror:
      VectorMirrorType neighbor_mirror;

      // Call the MirrorAssembler to do the dirty work for us:
      Assembly::MirrorAssembler::assemble_mirror(
        neighbor_mirror,   // the mirror that is to be assembled
        space,              // the FE space for which we want to assemble the mirror
        *neighbor_halo     // the halo mesh-part that the mirror is to be assembled on
      );

      // Once the mirror is assembled, we give it over to our gate.
      gate.push(
        neighbor_rank ,                    // the process rank of the neighbor
        std::move(neighbor_mirror)         // the mirror for the neighbor
      );
      watch.stop();

    }

    //internal init for global gate
    global_gate.compile(LocalVectorType(space.get_num_dofs()));

    Cubature::DynamicFactory cubature_factory("auto-degree:5");

    // Again, we use the sine-bubble as a reference solution:
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_function;


    // Create a force functional for our solution:
    Assembly::Common::LaplaceFunctional<decltype(sol_function)> force_functional(sol_function);


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Boundary Condition assembly
    // We create an instance of the global filter for post processing purposes:
    GlobalFilterType filter;

    // And then fetch the internal local filter object:
    LocalFilterType& filter_local = filter.local();

    // Now, we have to assemble the filter.
    // In the case of the unit-filter, which is responsible for Dirichlet boundary conditions,
    // the assembly is quite simple: We just have to assemble the corresponding local filter.
    //As we need a filter on each Level, we loop over all Levels
    bool is_floating = true;
    watch.start();
    for(Index ilevel(level_min); ilevel <= level_max; ++ilevel)
    {
      comm.print("Assembling Matrix and Filter for level " + stringify(ilevel) + "... ");
      // Get a reference to the corresponding level for easier member access
      Level& local_lvl = *levels.at(std::size_t(level_max - ilevel));

      // Assemble the Laplace matrix:
      Assembly::SymbolicAssembler::assemble_matrix_std1(local_lvl.matrix, local_lvl.space);
      Assembly::Common::LaplaceOperator _laplace_operator;
      local_lvl.matrix.format();
      Assembly::BilinearOperatorAssembler::assemble_matrix1(local_lvl.matrix, _laplace_operator, local_lvl.space, cubature_factory);

      // Assemble the unit filter for homogeneous Dirichlet boundary conditions
      Assembly::UnitFilterAssembler<MeshType> unit_asm;

      for(int ibnd(0); ibnd < 2*ShapeType::dimension; ++ibnd)
      {
        // Build the name of the boundary mesh-part and call the mesh-node's "find_mesh_part"
        // function to get a pointer to the mesh-part:
        MeshPartType* bnd_mesh_part = local_lvl.root_mesh_node->find_mesh_part("bnd:" + stringify(ibnd));
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
          is_floating = false;
        }
      }
      //set up Dirichlet Filter
      unit_asm.assemble(local_lvl.filter_chain.at<0>(), local_lvl.space);

      if(is_floating)//only Neumann boundary
      {
        Assembly::MeanFilterAssembler::assemble(local_lvl.filter_chain.at<1>(), local_lvl.space , cubature_factory);
      }
    } // end of level loop
    if(!is_floating)
    {
      gate.set_not_floating();
    }
    watch.stop();

    //create a deep copy of the UnitFilter for the global filter
    filter_local = levels.front()->filter_chain.at<0>().clone();

    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    // Next, we need to assemble the grid transfer (prolongation and restriction) matrices for
    // each pair of consecutive levels.

    // Now loop over all level pairs:
    for(Index ilevel(level_min); ilevel < level_max; ++ilevel)
    {
      comm.print("Assembling Grid Transfer for level " + stringify(ilevel) + "...");

      // Get references to the corresponding coarse and fine levels:
      Level& lvl_fine   = *levels.at(std::size_t(level_max - ilevel - 1));
      Level& lvl_coarse = *levels.at(std::size_t(level_max - ilevel));

      // Now let's get the references to the internal prolongation and restriction matrices
      // of the transfer operator:
      LocalMatrixType& mat_prol = lvl_fine.transfer.get_mat_prol();
      LocalMatrixType& mat_rest = lvl_fine.transfer.get_mat_rest();

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

    //finalize our gate
    if(verbose)
      comm.print("Setting up Gate...");
    TimeStamp ts_1, ts_2;
    double t{0};
    double t_max{0};
    watch.start();
    ts_1.stamp();
    gate.finalize();
    ts_2.stamp();
    watch.stop();
    t = ts_2.elapsed(ts_1);
    comm.allreduce(&t, &t_max, std::size_t(1), Dist::op_max);
    comm.print("Time to init gate: " + stringify(t_max) + " s");


    //now we can initialise the local system by moving our levels deque, after this, we only want to access levels trough local_sys:
    watch.start();
    ts_1.stamp();
    LocalSystem local_sys(std::move(levels), cubature_factory, force_functional, is_floating);
    ts_2.stamp();
    watch.stop();
    t = ts_2.elapsed(ts_1);
    comm.allreduce(&t, &t_max, std::size_t(1), Dist::op_max);
    comm.print("Time to init LocalSystem: " + stringify(t_max) + " s");


    GlobalMatrixType matrix(&global_gate, &global_gate, local_sys._matrix_local.clone(LAFEM::CloneMode::Shallow));

    GlobalVectorType vec_sol = matrix.create_vector_r();
    GlobalVectorType vec_rhs = matrix.create_vector_l();
    GlobalVectorType vec_rhs_buff = matrix.create_vector_l();
    LocalVectorType& vec_sol_local = vec_sol.local();
    LocalVectorType& vec_rhs_local = vec_rhs.local();

    vec_sol_local.copy(local_sys._vec_sol);
    vec_rhs_local.copy(local_sys._vec_rhs);

    vec_rhs.sync_0();


    //now initialise the Projector, this will assemble the local BR_columns matrix
    if(verbose)
      comm.print("Assembling projection...");
    watch.start();
    ts_1.stamp();
    Projector project(gate, comm);
    if(verbose)
      comm.print("Gather Q matrix...");
    project.gather_Q(gate, verbose);
    ts_2.stamp();
    watch.stop();
    t = ts_2.elapsed(ts_1);
    comm.allreduce(&t, &t_max, std::size_t(1), Dist::op_max);
    comm.print("Time to init projector: " + stringify(t_max) + " s");

    if(verbose)
      comm.print("Initialising CG...");
    watch.start();
    ts_1.stamp();
    FETICG CG(comm, gate, local_sys, CG_tol, CG_max_iter);
    CG.init(gate, local_sys, project);
    ts_2.stamp();
    watch.stop();
    t = ts_2.elapsed(ts_1);
    comm.allreduce(&t, &t_max, std::size_t(1), Dist::op_max);
    comm.print("Time to init CG: " + stringify(t_max) + " s");

    watch.start();
    test_passed = CG.solve(gate, local_sys, project, CG_plot, want_test);
    watch.stop();

    //calculate alpha, as CG holds all necessary values, this also will be handled by the CG class
    watch.start();
    ts_1.stamp();
    CG.calc_alpha(gate, local_sys, project);

    //now we will calculate our solution. The solution will be saved in vec_sol in LocalSystem
    CG.end_result(gate, local_sys);
    ts_2.stamp();
    watch.stop();
    t = ts_2.elapsed(ts_1);
    comm.allreduce(&t, &t_max, std::size_t(1), Dist::op_max);
    comm.print("Time to calc final solution: " + stringify(t_max) + " s");

    //final step: post processing
    comm.print("\nComputing errors against reference solution...");

    Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
      local_sys.get_sol_ref(), sol_function, space, cubature_factory);

    // And then we need to synchronize the errors over our communicator to sum up the errors of
    // each patch to obtain the errors over the whole domain:
    errors.synchronize(comm);

    //if we are in test case... test the errors:
    if(want_test)
    {
      DataType max_h0 = 2e-3;
      DataType max_h1 = 2e-1;
      if(errors.norm_h0 >= max_h0)
      {
        comm.print("Test FAILED: H0-Norm " + stringify(errors.norm_h0) + " should be less than "
                                      + stringify(max_h0));
        test_passed = false;
      }
      if(errors.norm_h1 >= max_h1)
      {
        comm.print("Test failed: H1-Norm " + stringify(errors.norm_h1) + " should be less than "
                                      + stringify(max_h1));
        test_passed = false;
      }
    }

    // And let's print the errors to the console; we need to use the "format_string" function here,
    // as the "print" function accepts only String objects as input:
    comm.print(errors.format_string());

    //Calculate the residual of our solution:
    //First copy solution into our vec_sol_local
    vec_sol_local.copy(local_sys.get_sol_ref());
    matrix.apply(vec_rhs_buff, vec_sol);
    //apply filter on buffer
    filter.filter_rhs(vec_rhs_buff);
    vec_rhs_buff.axpy(vec_rhs_buff, vec_rhs, -1.);
    DataType global_res_dot = vec_rhs_buff.norm2();
    DataType max_glob = 4e-9;
    if(want_test && global_res_dot >= max_glob)
    {
      comm.print("Test FAILED: global residual " + stringify(global_res_dot) + " should be less then "
                                       + stringify(max_glob));
      test_passed = false;
    }
    comm.print("The 2 Norm of the residual is: " + stringify(global_res_dot));

    t = watch.elapsed();
    comm.allreduce(&t, &t_max, std::size_t(1), Dist::op_max);
    comm.print("Time for whole programm: " + stringify(t_max) + " s");

    MemoryUsage mem_use;
    std::size_t peak_phys = mem_use.get_peak_physical();
    std::size_t peak_virt = mem_use.get_peak_virtual();

    std::size_t phys_max{0}, phys_sum{0}, virt_max{0}, virt_sum{0};
    comm.allreduce(&peak_phys, &phys_max, std::size_t(1), Dist::op_max);
    comm.allreduce(&peak_phys, &phys_sum, std::size_t(1), Dist::op_sum);
    comm.allreduce(&peak_virt, &virt_max, std::size_t(1), Dist::op_max);
    comm.allreduce(&peak_virt, &virt_sum, std::size_t(1), Dist::op_sum);

    comm.print("Maximum peak pyhsical memory: " + stringify(phys_max) + " B");
    comm.print("Sum peak pyhsical memory: " + stringify(phys_sum) + " B");
    comm.print("Maximum peak virtual memory: " + stringify(virt_max) + " B");
    comm.print("Sum peak virtual memory: " + stringify(virt_sum) + " B");



    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // Finally, we want to export our solution to a (P)VTU file set.

    // In a parallel simulation, each process will write a separate VTU file, which contains the
    // data that is defined on the patch of the corresponding process. Moreover, one process
    // writes a single additional PVTU file, which can be read by ParaView to visualize
    // the whole domain that consists of all patches.

    if(!want_vtk)
    {
      comm.print("Finished!");
      if(want_test && test_passed)
        comm.print("Test PASSED");
      return;
    }
    // Build the VTK filename; we also append the number of processes to the filename:
    vtk_name = vtk_name + stringify(level) + "-n" + stringify(comm.size());

    comm.print("Writing VTK file '" + vtk_name + ".pvtu'...");

    // Create a VTK exporter for our patch mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // Add the vertex-projection of our (local) solution and rhs vectors
    exporter.add_vertex_scalar("sol", local_sys.get_sol_elements());
    exporter.add_vertex_scalar("rhs", local_sys.get_rhs_elements());

    // Finally, write the VTK files by calling the "write" function of the exporter and pass the
    // communicator as a second argument:
    exporter.write(vtk_name, comm);

    // Note: Do not forget the 'comm' argument in the call above as otherwise each process will
    // try to write to the same VTK file, resulting in garbage due to race conditions...

    // That's all, folks.
    comm.print("Finished!");
    if(want_test && test_passed)
        comm.print("Test PASSED");

    return;

  }//main function

}//namespace FETI



int main(int argc, char* argv[])
{
  Runtime::initialize(argc, argv);

  FETI::main(argc, argv);






  return Runtime::finalize();
}