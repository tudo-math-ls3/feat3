#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/solver/legacy_preconditioners.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/precon_wrapper.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/util/dist.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/poisson_mixed.hpp>
#include <control/statistics.hpp>

namespace PoissonMixed2D
{
  using namespace FEAT;

  template
  <
    typename LumpedMatrixA_,
    typename MatrixB_,
    typename MatrixD_,
    typename FilterA_
  >
  class VirtualLumpedSchurMatrix
  {
    public:
      typedef LumpedMatrixA_ LumpedMatrixA;
      typedef MatrixB_ MatrixB;
      typedef MatrixD_ MatrixD;
      typedef FilterA_ FilterA;

      typedef typename LumpedMatrixA_::DataType DataType;

      typedef typename MatrixD::VectorTypeL VectorTypeL;
      typedef typename MatrixB::VectorTypeR VectorTypeR;

      typedef LumpedMatrixA VectorTypeML;
      typedef LumpedMatrixA VectorTypeMR;

      LumpedMatrixA lumped_matrix_a;
      const MatrixB& matrix_b;
      const MatrixD& matrix_d;
      const FilterA& filter_a;

    private:
      mutable VectorTypeML _vec_ml;
      mutable VectorTypeMR _vec_mr;

    public:
      VirtualLumpedSchurMatrix(const LumpedMatrixA& lumped_matrix_a_,
      const MatrixB& matrix_b_,
      const MatrixD& matrix_d_,
      const FilterA& filter_a_) :
        lumped_matrix_a(lumped_matrix_a_.clone(LAFEM::CloneMode::Layout)),
        matrix_b(matrix_b_),
        matrix_d(matrix_d_),
        filter_a(filter_a_),
        _vec_ml(lumped_matrix_a_.clone(LAFEM::CloneMode::Layout)),
        _vec_mr(lumped_matrix_a_.clone(LAFEM::CloneMode::Layout))
        {
          lumped_matrix_a.component_invert(lumped_matrix_a_);
          //lumped_matrix_a.format(2);
        }

      virtual ~VirtualLumpedSchurMatrix()
      {
      }

      VirtualLumpedSchurMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return VirtualSchurMatrix(lumped_matrix_a.clone(mode), matrix_b.clone(mode), matrix_d.clone(mode),
        filter_a.clone(mode));
      }

      VectorTypeL create_vector_l() const
      {
        return matrix_d.create_vector_l();
      }

      VectorTypeR create_vector_r() const
      {
        return matrix_b.create_vector_r();
      }

      /**
       * \brief Gets the total number of columns in this matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of columns
       */
      Index columns() const
      {
        return matrix_b.columns();
      }

      /**
       * \brief Gets the total number of rows in this matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of columns
       */
      Index rows() const
      {
        return matrix_d.columns();
      }

      /**
       * \brief Returns the total number of non-zeros in this matrix.
       *
       * \note This always returns the raw (or POD - Plain Old Data) count, as everything else is ambiguous.
       *
       * \returns The total number of nonzeros in this matrix
       */
      Index used_elements() const
      {
        return lumped_matrix_a.used_elements() + matrix_b.used_elements() + matrix_d.used_elements();
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _vec_ml.bytes() + _vec_mr.bytes();
      }

      void extract_diag(VectorTypeL& diag, bool sync=true) const
      {
        if(diag.get_gate() != nullptr && !diag.get_gate()->_ranks.empty() )
        {
          diag.format();

          typename MatrixB::LocalMatrix matrix_b1(matrix_b.convert_to_1());

          matrix_b1.add_trace_double_mat_mult(diag.local(), matrix_d.local(), lumped_matrix_a.local(), DataType(1));
        }
        else
        {
          matrix_d.local().row_norm2sqr(diag.local(), lumped_matrix_a.local());
        }

        if(sync)
        {
          diag.sync_0();
        }
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        matrix_b.apply(_vec_mr, x);

        filter_a.filter_def(_vec_mr);
        _vec_ml.component_product(_vec_mr, lumped_matrix_a);
        filter_a.filter_cor(_vec_ml);

        matrix_d.apply(r, _vec_ml);
      }

      //auto apply_async(VectorTypeL& r, const VectorTypeR& x) const -> decltype(r.sync_0_async())
      //{
      //}

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const
      {
        // r <- y + alpha*(D A^(-1) B)*x
        matrix_b.apply(_vec_mr, x);

        filter_a.filter_def(_vec_mr);
        _vec_ml.component_product(_vec_mr, lumped_matrix_a);
        filter_a.filter_cor(_vec_ml);

        matrix_d.apply(r, _vec_ml, y, alpha);
      }

      //auto apply_async(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const -> decltype(r.sync_0_async())
      //{
      //  // copy y to r
      //  r.copy(y);

      //  // convert from type-1 to type-0
      //  r.from_1_to_0();

      //  // r <- r + alpha*A*x
      //  _matrix.apply(*r, *x, *r, alpha);

      //  // synchronise r
      //  r.sync_0_async();
      //}

      //LocalMatrix convert_to_1() const
      //{
      //  LocalMatrix locmat = _matrix.clone(LAFEM::CloneMode::Weak);
      //  if((_row_gate != nullptr) && (_col_gate != nullptr))
      //    synch_matrix(locmat, *_row_gate->_comm, _row_gate->_ranks, _row_gate->_mirrors, _col_gate->_mirrors);
      //  return locmat;
      //}

  };

  template
  <
    typename MatrixA_,
    typename MatrixB_,
    typename MatrixD_,
    typename FilterA_
  >
  class VirtualSchurMatrix
  {
    public:
      typedef MatrixA_ MatrixA;
      typedef MatrixB_ MatrixB;
      typedef MatrixD_ MatrixD;
      typedef FilterA_ FilterA;

      typedef typename MatrixA_::DataType DataType;

      typedef typename MatrixD::VectorTypeL VectorTypeL;
      typedef typename MatrixB::VectorTypeR VectorTypeR;

      typedef typename MatrixA::VectorTypeL VectorTypeML;
      typedef typename MatrixA::VectorTypeR VectorTypeMR;

      const MatrixA& matrix_a;
      const MatrixB& matrix_b;
      const MatrixD& matrix_d;
      const FilterA& filter_a;

    private:
      mutable VectorTypeML _vec_ml;
      mutable VectorTypeMR _vec_mr;
      std::shared_ptr<Solver::SolverBase<VectorTypeMR>> _solver_a;

    public:
      VirtualSchurMatrix(const MatrixA& matrix_a_,
      const MatrixB& matrix_b_,
      const MatrixD& matrix_d_,
      const FilterA& filter_a_,
      std::shared_ptr<Solver::SolverBase<VectorTypeMR>> solver_a_) :
        matrix_a(matrix_a_),
        matrix_b(matrix_b_),
        matrix_d(matrix_d_),
        filter_a(filter_a_),
        _vec_ml(matrix_a_.create_vector_l()),
        _vec_mr(matrix_a_.create_vector_r()),
        _solver_a(solver_a_)
        {
        }

      VirtualSchurMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return
          VirtualSchurMatrix(matrix_a.clone(mode), matrix_b.clone(mode),  matrix_d.clone(mode), filter_a.clone(mode));
      }

      VectorTypeL create_vector_l() const
      {
        return matrix_d.create_vector_l();
      }

      VectorTypeR create_vector_r() const
      {
        return matrix_b.create_vector_r();
      }

      /**
       * \brief Gets the total number of columns in this matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of columns
       */
      Index columns() const
      {
        return matrix_b.columns();
      }

      /**
       * \brief Gets the total number of rows in this matrix
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \returns The number of columns
       */
      Index rows() const
      {
        return matrix_d.columns();
      }

      /**
       * \brief Returns the total number of non-zeros in this matrix.
       *
       * \note This always returns the raw (or POD - Plain Old Data) count, as everything else is ambiguous.
       *
       * \returns The total number of nonzeros in this matrix
       */
      Index used_elements() const
      {
        return matrix_a.used_elements() + matrix_b.used_elements() + matrix_d.used_elements();
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _vec_ml.bytes() + _vec_mr.bytes();
      }

      void extract_diag(VectorTypeL& diag, bool sync=true) const
      {

        VectorTypeML lumped_matrix_a(matrix_a.create_vector_l());
        matrix_a.local().lump_rows(lumped_matrix_a.local());
        lumped_matrix_a.sync_0();
        lumped_matrix_a.component_invert(lumped_matrix_a);

        if(diag.get_gate() != nullptr && !diag.get_gate()->_ranks.empty() )
        {
          diag.format();

          typename MatrixB::LocalMatrix matrix_b1(matrix_b.convert_to_1());

          matrix_b1.add_trace_double_mat_mult(diag.local(), matrix_d.local(), lumped_matrix_a.local(), DataType(1));
        }
        else
        {
          matrix_d.local().row_norm2sqr(diag.local(), lumped_matrix_a.local());
        }


        if(sync)
        {
          diag.sync_0();
        }

      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        matrix_b.apply(_vec_mr, x);

        Solver::solve(*_solver_a, _vec_ml, _vec_mr, matrix_a, filter_a);

        matrix_d.apply(r, _vec_ml);
      }

      //auto apply_async(VectorTypeL& r, const VectorTypeR& x) const -> decltype(r.sync_0_async())
      //{
      //}

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const
      {
        // copy y to r
        r.copy(y);

        // r <- y + alpha*(D A^(-1) B)*x
        matrix_b.apply(_vec_mr, x);

        filter_a.filter_def(_vec_mr);
        Solver::solve(*_solver_a, _vec_ml, _vec_mr, matrix_a, filter_a);
        //filter_a.filter_cor(_vec_ml);

        matrix_d.apply(r, _vec_ml, r, alpha);
      }

      //auto apply_async(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const -> decltype(r.sync_0_async())
      //{
      //  // copy y to r
      //  r.copy(y);

      //  // convert from type-1 to type-0
      //  r.from_1_to_0();

      //  // r <- r + alpha*A*x
      //  _matrix.apply(*r, *x, *r, alpha);

      //  // synchronise r
      //  r.sync_0_async();
      //}

      //LocalMatrix convert_to_1() const
      //{
      //  LocalMatrix locmat = _matrix.clone(LAFEM::CloneMode::Weak);
      //  if((_row_gate != nullptr) && (_col_gate != nullptr))
      //    synch_matrix(locmat, *_row_gate->_comm, _row_gate->_ranks, _row_gate->_mirrors, _col_gate->_mirrors);
      //  return locmat;
      //}

      virtual ~VirtualSchurMatrix()
      {
      }
  };

  template<typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our arch types
    typedef Mem::Main MemType;
    typedef double DataType;
    typedef Index IndexType;

    static constexpr int dim = 2;

    // choose our desired analytical solution
    Analytic::Common::ExpBubbleFunction<dim> sol_func;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;

    // define our system level
    typedef Control::PoissonMixedNoneVeloUnitPresSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = domain.size_physical();

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("auto-degree:7");
    //Cubature::DynamicFactory cubature("gauss-lobatto:3");

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    comm.print("Assembling gates...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gates(domain.at(i));
    }

    /* ***************************************************************************************** */

    comm.print("Assembling transfers...");

    for (Index i(0); (i+1) < domain.size_virtual(); ++i)
    {
      system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
      system_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system matrices...");

    for (Index i(0); i < num_levels; ++i)
    {
      // Assemble matrix structure
      system_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
      system_levels.at(i)->assemble_pres_struct(domain.at(i)->space_pres);

      system_levels.at(i)->matrix_a.local().format();
      Assembly::Common::IdentityOperatorBlocked<dim> identity_op;
      Assembly::BilinearOperatorAssembler::assemble_block_matrix1(system_levels.at(i)->matrix_a.local(),
      identity_op, domain.at(i)->space_velo, cubature);

      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
    }

    /* ***************************************************************************************** */

    comm.print("Assembling system filters...");

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local system filters
      //typename SystemLevelType::LocalVeloFilter& fil_loc_v = system_levels.at(i)->filter_velo.local();
      typename SystemLevelType::LocalPresFilter& fil_loc_p = system_levels.at(i)->filter_pres.local();

      // create unit-filter assemblers
      Assembly::UnitFilterAssembler<MeshType> unit_asm_pres;

      // loop over all boundary parts
      std::deque<String> part_names = domain.at(i)->get_mesh_node()->get_mesh_part_names(true);
      for(const auto& name : part_names)
      {
        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);

        // found it?
        XASSERT(mesh_part_node != nullptr);

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        auto* mesh_part = mesh_part_node->get_mesh();
        if(mesh_part == nullptr)
          continue;

        // add to corresponding boundary assembler
        unit_asm_pres.add_mesh_part(*mesh_part);
      }

      // assemble the pressure filter
      unit_asm_pres.assemble(fil_loc_p, domain.at(i)->space_pres, sol_func);
    }

    Statistics::toe_assembly = stamp_ass.elapsed_now();

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalPresVector GlobalSystemVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // create new vector
    GlobalSystemVector vec_sol = the_system_level.matrix_b.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_b.create_vector_r();

    vec_sol.format();
    vec_rhs.format();

    {
      // get the local vector
      typename SystemLevelType::LocalPresVector& vec_f = vec_rhs.local();

      // assemble the force
      Assembly::Common::LaplaceFunctional<decltype(sol_func)> force_func(sol_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_f, force_func, the_domain_level.space_pres, cubature);

      // sync the vector
      vec_rhs.sync_0();
    }

    // and filter it
    the_system_level.filter_pres.filter_sol(vec_sol);
    the_system_level.filter_pres.filter_rhs(vec_rhs);

    typedef VirtualLumpedSchurMatrix
    <
      typename SystemLevelType::GlobalVeloVector,
      typename SystemLevelType::GlobalMatrixBlockB,
      typename SystemLevelType::GlobalMatrixBlockD,
      typename SystemLevelType::GlobalVeloFilter
    > GlobalLumpedSystemMatrix;

    typedef VirtualSchurMatrix
    <
      typename SystemLevelType::GlobalMatrixBlockA,
      typename SystemLevelType::GlobalMatrixBlockB,
      typename SystemLevelType::GlobalMatrixBlockD,
      typename SystemLevelType::GlobalVeloFilter
    > GlobalFullSystemMatrix;

    typedef GlobalLumpedSystemMatrix GlobalSystemMatrix;

    std::deque<GlobalSystemMatrix> matrix_sys;

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("Setting up solver...");

    // PCG ( VCycle ( S: Richardson ( Jacobi )  / C: Richardson ( Jacobi )  )  )
    auto multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy
      <
        GlobalSystemMatrix,
        typename SystemLevelType::GlobalPresFilter,
        typename SystemLevelType::GlobalPresTransfer
      >
    >(domain.size_virtual());

    // push all levels except the coarse most one
    for (Index i(0); (i+1) < num_levels; ++i)
    {
      const SystemLevelType& lvl = *system_levels.at(i);

      // VirtualLumpedSchurMatrix version
      {
        typename SystemLevelType::GlobalVeloVector lumped_matrix_a(lvl.matrix_a.create_vector_l());
        lvl.matrix_a.local().lump_rows(lumped_matrix_a.local());
        lumped_matrix_a.sync_0();

        matrix_sys.emplace_back(lumped_matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo);
      }

      //// VirtualSchurMatrix version
      //{
      //  auto precon_a = Solver::new_jacobi_precond( lvl.matrix_a, lvl.filter_velo, 0.5);

      //  auto inner_solver = Solver::new_pcg(lvl.matrix_a, lvl.filter_velo);
      //  //inner_solver->set_plot_mode(Solver::PlotMode::iter);
      //  //inner_solver->set_plot_name("Inner");

      //  std::shared_ptr<Solver::SolverBase<typename SystemLevelType::GlobalVeloVector>> solver_a = inner_solver;
      //  solver_a->init();

      //  matrix_sys.emplace_back(lvl.matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo, solver_a);
      //}

      auto jac_smoother = Solver::new_jacobi_precond(matrix_sys.at(i), lvl.filter_pres, 0.5);
      auto smoother = Solver::new_richardson(matrix_sys.at(i), lvl.filter_pres, 1.0, jac_smoother);
      smoother->set_min_iter(4);
      smoother->set_max_iter(4);
      multigrid_hierarchy->push_level(matrix_sys.at(i), lvl.filter_pres, lvl.transfer_pres, smoother, smoother, smoother);
    }

    // push the coarse level
    {
      const SystemLevelType& lvl = *system_levels.back();

      // VirtualLumpedSchurMatrix version
      {
        typename SystemLevelType::GlobalVeloVector lumped_matrix_a(lvl.matrix_a.create_vector_l());
        lvl.matrix_a.local().lump_rows(lumped_matrix_a.local());
        lumped_matrix_a.sync_0();

        matrix_sys.emplace_back(lumped_matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo);
      }

      //// VirtualSchurMatrix version
      //{
      //  auto precon_a = Solver::new_jacobi_precond( lvl.matrix_a, lvl.filter_velo, 0.7);

      //  auto inner_solver = Solver::new_pcg(lvl.matrix_a, lvl.filter_velo);
      //  //inner_solver->set_plot_mode(Solver::PlotMode::iter);
      //  //inner_solver->set_plot_name("Inner");

      //  std::shared_ptr<Solver::SolverBase<typename SystemLevelType::GlobalVeloVector>> solver_a = inner_solver;
      //  solver_a->init();

      //  matrix_sys.emplace_back(lvl.matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo, solver_a);
      //}

      auto coarse_precond = Solver::new_jacobi_precond(matrix_sys.back(), lvl.filter_pres, 0.5);
      auto coarse_solver = Solver::new_richardson(matrix_sys.back(), lvl.filter_pres, 1.0, coarse_precond);
      coarse_solver->set_min_iter(4);
      coarse_solver->set_max_iter(4);
      multigrid_hierarchy->push_level(matrix_sys.back(), lvl.filter_pres, coarse_solver);
    }

    //{
    //  const SystemLevelType& lvl = *system_levels.front();
    //  auto precon_a = Solver::new_jacobi_precond( lvl.matrix_a, lvl.filter_velo, 0.7);
    //  auto inner_solver = Solver::new_pcg(lvl.matrix_a, lvl.filter_velo);
    //  std::shared_ptr<Solver::SolverBase<typename SystemLevelType::GlobalVeloVector>> solver_a = inner_solver;
    //  solver_a->init();

    //  GlobalFullSystemMatrix full_matrix_sys(lvl.matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo, solver_a);
    //  typename SystemLevelType::GlobalPresVector diag_a(lvl.matrix_d.create_vector_l());
    //  full_matrix_sys.extract_diag(diag_a);

    //  std::cout << comm.rank() << " full_diag_a " << diag_a.local() << std::endl;

    //  typename SystemLevelType::GlobalVeloVector lumped_matrix_a(lvl.matrix_a.create_vector_r());
    //  lvl.matrix_a.local().lump_rows(lumped_matrix_a.local());
    //  lumped_matrix_a.sync_0();

    //  GlobalLumpedSystemMatrix lumped_matrix_sys(lumped_matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo);
    //  lumped_matrix_sys.extract_diag(diag_a);

    //  std::cout <<  comm.rank() << " lump_diag_a " << diag_a.local() << std::endl;
    //}

    auto mgv = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);
    //auto jac_precond = Solver::new_jacobi_precond(matrix_sys.front(), the_system_level.filter_pres, 0.5);
    //auto solver = Solver::new_pcg(matrix_sys.front(), the_system_level.filter_pres, jac_precond);
    auto solver = Solver::new_richardson(matrix_sys.front(), the_system_level.filter_pres, 1.0, mgv);

    // enable plotting
    if(comm.rank() == 0)
    {
      solver->set_plot_mode(Solver::PlotMode::iter);
    }

    // set tolerance
    solver->set_tol_rel(1E-8);
    solver->set_max_iter(1000);

    // initialise
    multigrid_hierarchy->init();
    solver->init();

    Statistics::reset();

    TimeStamp at;

    // solve
    auto result = Solver::solve(*solver, vec_sol, vec_rhs, matrix_sys, the_system_level.filter_pres);

    if (!Solver::status_success(result))
    {
      comm.print("Solver execution FAILED, with status: " + stringify(result));
    }

    double solver_toe(at.elapsed_now());

    FEAT::Control::Statistics::report(solver_toe, args.check("statistics"), MeshType::ShapeType::dimension,
      system_levels, domain);

    // release solver
    //solver_a->done();
    solver->done();
    multigrid_hierarchy->done();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("no-err") < 0)
    {
      // Compute and print the H0-/H1-errors
      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute
        (vec_sol.local(), sol_func, the_domain_level.space_pres, cubature);

      // synchronise all local errors
      errors.synchronise(comm);

      // print errors
      comm.print("");
      comm.print(errors.format_string());
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./poisson-dirichlet-2d");
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(comm.size());

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

      // project velocity and pressure
      typename SystemLevelType::LocalPresVector vtx_sol, vtx_rhs;
      Assembly::DiscreteVertexProjector::project(vtx_sol, vec_sol.local(), the_domain_level.space_pres);
      Assembly::DiscreteVertexProjector::project(vtx_rhs, vec_rhs.local(), the_domain_level.space_pres);

      // write velocity
      exporter.add_vertex_scalar("sol", vtx_sol.elements());
      exporter.add_vertex_scalar("rhs", vtx_rhs.elements());

      // finally, write the VTK file
      exporter.write(vtk_name, comm);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("test-iter") >= 0)
    {
      int num_iter = (int)solver->get_num_iter();
      int iter_target(0);
      args.parse("test-iter", iter_target);
      if (num_iter < iter_target - 2 || num_iter > iter_target + 2)
      {
        comm.print("FAILED");
        throw InternalError(__func__, __FILE__, __LINE__, "iter count deviation! " + stringify(num_iter) + " vs " + stringify(iter_target));
      }
    }
  }

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

#ifdef FEAT_HAVE_MPI
    comm.print("NUM-PROCS: " + stringify(comm.size()));
#endif

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    args.support("level");
    args.support("no-err");
    args.support("vtk");
    args.support("statistics");
    args.support("min-rank-elems");
    args.support("mesh");
    args.support("test-iter");
    args.support("parti-type");
    args.support("parti-name");
    args.support("parti-rank-elems");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());

      // abort
      FEAT::Runtime::abort();
    }

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
      FEAT::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> AuxSpace;
    //typedef Space::CroRavRanTur::Element<TrafoType> AuxSpace;
    //typedef Space::Lagrange1::Element<TrafoType> SpaceType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1> > SpaceType;

    int lvl_max = 3;
    int lvl_min = 0;
    args.parse("level", lvl_max, lvl_min);

#ifndef DEBUG
    try
#endif
    {
      TimeStamp stamp1;

      // let's create our domain
      comm.print("Preparing domain...");

      // query mesh filename list
      const std::deque<String>& mesh_filenames = args.query("mesh")->second;

      // create our domain control
      typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, AuxSpace, SpaceType> DomainLevelType;
      Control::Domain::PartiDomainControl<DomainLevelType> domain(comm);

      // let the controller parse its arguments
      if(!domain.parse_args(args))
      {
        FEAT::Runtime::abort();
      }

      // read the base-mesh
      domain.read_mesh(mesh_filenames);

      TimeStamp stamp_partition;

      // try to create the partition
      domain.create_partition();

      Statistics::toe_partition = stamp_partition.elapsed_now();

      comm.print("Creating mesh hierarchy...");

      // create the level hierarchy
      domain.create_hierarchy(lvl_max, lvl_min);

      // plot our levels
      comm.print("LVL-MAX: " + stringify(domain.max_level_index()) + " [" + stringify(lvl_max) + "]");
      comm.print("LVL-MIN: " + stringify(domain.min_level_index()) + " [" + stringify(lvl_min) + "]");

      // run our application
      run(args, domain);

      TimeStamp stamp2;

      // get times
      long long time1 = stamp2.elapsed_micros(stamp1);

      // accumulate times over all processes
      long long time2 = time1 * (long long) comm.size();

      // print time
      comm.print("Run-Time: " + stringify(TimeStamp::format_micros(time1, TimeFormat::m_s_m)) + " [" +
        stringify(TimeStamp::format_micros(time2, TimeFormat::m_s_m)) + "]");
    }
#ifndef DEBUG
    catch (const std::exception& exc)
    {
      std::cerr << "ERROR: unhandled exception: " << exc.what() << std::endl;
      FEAT::Runtime::abort();
    }
    catch (...)
    {
      std::cerr << "ERROR: unknown exception" << std::endl;
      FEAT::Runtime::abort();
    }
#endif // DEBUG
  }
} // namespace PoissonMixed2D

int main(int argc, char* argv [])
{
  // initialise
  FEAT::Runtime::initialise(argc, argv);

  PoissonMixed2D::main(argc, argv);
  // okay
  return FEAT::Runtime::finalise();
}
