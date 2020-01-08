// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/global/symmetric_lumped_schur_matrix.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_mixed.hpp>
#include <control/statistics.hpp>

namespace PoissonMixed
{
  using namespace FEAT;

  static void display_help(const Dist::Comm& comm)
  {
    comm.print("poisson-mixed: Solve the Poisson problem in mixed formulation e.g. with Q1~/P0 or Q2/P1disc or Q2/Q1.");
    comm.print("Boundary conditions are either homogeneous Dirichlet or homogeneous Neumann on all outer boundaries.\n");
    comm.print("Usage: poisson-mixed --mesh [PATH_TO_MESH] --level [LVL_MAX LVL_MED LVL_MIN] --bc [BC].\n");
    comm.print("Only LVL_MAX is mandatory.");
    comm.print("--bc needs either 'dirichlet' or 'neumann'");
  }

  template
  <
    typename MatrixA_,
    typename MatrixB_,
    typename MatrixD_,
    typename FilterA_
  >
  class SymmetricSchurMatrix
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
      SymmetricSchurMatrix(const MatrixA& matrix_a_,
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

      virtual ~SymmetricSchurMatrix()
      {
      }

      SymmetricSchurMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return
          SymmetricSchurMatrix(matrix_a.clone(mode), matrix_b.clone(mode),  matrix_d.clone(mode), filter_a.clone(mode));
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

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh and shape types
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainControlType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // choose our desired analytical solution
    // This has homogeneous Dirichlet BCs
    //Analytic::Common::ExpBubbleFunction<dim> sol_func;
    // This has homogeneous Neumann BCs
    Analytic::Common::CosineWaveFunction<dim> sol_func;

    // define our system level
    typedef Control::ScalarMixedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = domain.size_physical();

    String bc_name;
    args.parse("bc", bc_name);

    bool essential_bcs(false);
    if(bc_name == "neumann")
    {
      essential_bcs = true;
      // create system levels
      for (Index i(0); i < num_levels; ++i)
      {
        std::deque<String> neumann_list = domain.at(i)->get_mesh_node()->get_mesh_part_names(true);
        system_levels.push_back(std::make_shared<SystemLevelType>(neumann_list));
      }
    }
    else if(bc_name == "dirichlet")
    {
      // create system levels
      for (Index i(0); i < num_levels; ++i)
      {
        system_levels.push_back(std::make_shared<SystemLevelType>());
      }
    }
    else
    {
        throw InternalError(__func__, __FILE__, __LINE__, "--bc supports only 'dirichlet' or 'neumann', but got "+bc_name);
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

      system_levels.at(i)->matrix_a.local().format();
      Assembly::Common::IdentityOperatorBlocked<dim> identity_op;
      Assembly::BilinearOperatorAssembler::assemble_block_matrix1(system_levels.at(i)->matrix_a.local(),
      identity_op, domain.at(i)->space_velo, cubature);

      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
    }

    /* ***************************************************************************************** */

    if (essential_bcs)
    {
      comm.print("Assembling system filters...");

      for(Index i(0); i < num_levels; ++i)
      {
        // get our local system filters
        typename SystemLevelType::LocalVeloFilter& fil_loc_v = system_levels.at(i)->filter_velo.local();

        for(auto& it: fil_loc_v)
        {
          // Add the MeshPart to the assembler if it is there. There are legimate reasons for it NOT to be
          // there, i.e. we are in parallel and our patch is not adjacent to that MeshPart
          auto* mpp = domain.at(i)->get_mesh_node()->find_mesh_part(it.first);

          // Create empty assembler
          Assembly::SlipFilterAssembler<MeshType> neumann_asm(domain.at(i)->get_mesh());
          if(mpp != nullptr)
          {
            neumann_asm.add_mesh_part(*mpp);
          }
          // Even if the assembler is empty, this call is necessary for the filter to have the right size
          neumann_asm.assemble(it.second, domain.at(i)->space_velo);
        }

        // This shallow-copies the filters to filter_sys
        system_levels.at(i)->compile_system_filter();
        // After assembling to local filters, we need to call this to synchronise e.g. the slip filters
        system_levels.at(i)->assemble_global_filters();
        // Compile the system matrix
        system_levels.at(i)->compile_system_matrix();

      }
    }
    else
    {
      comm.print("Not Assembling system filters because there are no Neumann boundaries");
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

    typedef Global::SymmetricLumpedSchurMatrix
    <
      typename SystemLevelType::GlobalVeloVector,
      typename SystemLevelType::GlobalMatrixBlockB,
      typename SystemLevelType::GlobalMatrixBlockD,
      typename SystemLevelType::GlobalVeloFilter
    > GlobalLumpedSystemMatrix;

    //typedef SymmetricSchurMatrix
    //<
    //  typename SystemLevelType::GlobalMatrixBlockA,
    //  typename SystemLevelType::GlobalMatrixBlockB,
    //  typename SystemLevelType::GlobalMatrixBlockD,
    //  typename SystemLevelType::GlobalVeloFilter
    //> GlobalFullSystemMatrix;

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

      // SymmetricLumpedSchurMatrix version
      {
        typename SystemLevelType::GlobalVeloVector lumped_matrix_a(lvl.matrix_a.create_vector_l());
        lvl.matrix_a.local().lump_rows(lumped_matrix_a.local());
        lumped_matrix_a.sync_0();

        matrix_sys.emplace_back(lumped_matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo);
      }

      //// SymmetricSchurMatrix version
      //{
      //  auto precon_a = Solver::new_jacobi_precond( lvl.matrix_a, lvl.filter_velo, 0.5);

      //  auto inner_solver = Solver::new_pcg(lvl.matrix_a, lvl.filter_velo);
      //  //inner_solver->set_plot_mode(Solver::PlotMode::iter);
      //  //inner_solver->set_plot_name("Inner");

      //  std::shared_ptr<Solver::SolverBase<typename SystemLevelType::GlobalVeloVector>> solver_a = inner_solver;
      //  solver_a->init();

      //  matrix_sys.emplace_back(lvl.matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo, solver_a);
      //}

      auto jac_smoother = Solver::new_jacobi_precond(matrix_sys.at(i), lvl.filter_pres, 0.7);
      auto smoother = Solver::new_richardson(matrix_sys.at(i), lvl.filter_pres, 1.0, jac_smoother);
      smoother->set_min_iter(4);
      smoother->set_max_iter(4);
      multigrid_hierarchy->push_level(matrix_sys.at(i), lvl.filter_pres, lvl.transfer_pres, smoother, smoother, smoother);
    }

    // push the coarse level
    {
      const SystemLevelType& lvl = *system_levels.back();

      // SymmetricLumpedSchurMatrix version
      {
        typename SystemLevelType::GlobalVeloVector lumped_matrix_a(lvl.matrix_a.create_vector_l());
        lvl.matrix_a.local().lump_rows(lumped_matrix_a.local());
        lumped_matrix_a.sync_0();

        matrix_sys.emplace_back(lumped_matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo);
      }

      //// SymmetricSchurMatrix version
      //{
      //  auto precon_a = Solver::new_jacobi_precond( lvl.matrix_a, lvl.filter_velo, 0.7);

      //  auto inner_solver = Solver::new_pcg(lvl.matrix_a, lvl.filter_velo);
      //  //inner_solver->set_plot_mode(Solver::PlotMode::iter);
      //  //inner_solver->set_plot_name("Inner");

      //  std::shared_ptr<Solver::SolverBase<typename SystemLevelType::GlobalVeloVector>> solver_a = inner_solver;
      //  solver_a->init();

      //  matrix_sys.emplace_back(lvl.matrix_a, lvl.matrix_b, lvl.matrix_d, lvl.filter_velo, solver_a);
      //}

      auto coarse_precond = Solver::new_jacobi_precond(matrix_sys.back(), lvl.filter_pres, 0.7);
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
    auto solver = Solver::new_pcg(matrix_sys.front(), the_system_level.filter_pres, mgv);
    //auto solver = Solver::new_richardson(matrix_sys.front(), the_system_level.filter_pres, 1.0, mgv);

    // enable plotting
    solver->set_plot_mode(Solver::PlotMode::iter);

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
      String vtk_name = String("./poisson-mixed");
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
      comm.print("Writing "+vtk_name);
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

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // Add the DomainControl's supported arguments
    Control::Domain::add_supported_pdc_args(args);
    // check command line arguments
    args.support("help", "Displays the help text");
    args.support("bc", " [BC] specifies the homogeneous boundary conditions, either 'dirichlet' or 'neumann'.");
    args.support("level", " [LVL_MAX LVL_MED LVL_MIN] specifies the maximum, intermediate and minimum levels.");
    args.support("mesh", " [PATH_TO_MESH] specifies the mesh to use.");
    args.support("no-err", "Skip the computation of L2/H1 errors vs. the analytical solution.");
    args.support("statistics", "Display execution time and memory statistics.");
    args.support("test-iter", "  [ITER]: Run as a test and FAIL |solver_iterations - ITER| > 2.");
    args.support("vtk", "Write solution and rhs to vtk.");

    if(args.check("help") >= 0)
    {
      PoissonMixed::display_help(comm);
      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());
      // abort
      FEAT::Runtime::abort();
    }

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      PoissonMixed::display_help(comm);

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
      {
        comm.print(std::cerr, "\nERROR: unknown option '--" + (*it).second + "'");
      }

      // abort
      FEAT::Runtime::abort();
    }

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh [PATH_TO_MESH]' is missing!");
      display_help(comm);
      FEAT::Runtime::abort();
    }
    if(args.check("level") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--level [LVL_MAX LVL_MED LVL_MIN]' is missing!");
      display_help(comm);
      FEAT::Runtime::abort();
    }
    if(args.check("bc") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--bc [BC]' is missing!");
      display_help(comm);
      FEAT::Runtime::abort();
    }

    // define our mesh type
    static constexpr int dim = 2;
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> AuxSpace;
    //typedef Space::CroRavRanTur::Element<TrafoType> AuxSpace;
    //typedef Space::Lagrange1::Element<TrafoType> SpaceType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1> > SpaceType;

    // create a time-stamp
    TimeStamp time_stamp;

    // let's create our domain
    comm.print("Preparing domain...");

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, AuxSpace, SpaceType> DomainLevelType;
    bool support_hierarch_partitioning(false);
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, support_hierarch_partitioning);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(args.query("mesh")->second);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("LVL-MAX: " + stringify(domain.max_level_index()) + " [" + stringify(domain.get_desired_level_max()) + "]");
    comm.print("LVL-MIN: " + stringify(domain.min_level_index()) + " [" + stringify(domain.get_desired_level_min()) + "]");

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }


} // namespace PoissonMixed

int main(int argc, char* argv [])
{
  // initialise
  FEAT::Runtime::initialise(argc, argv);

  PoissonMixed::main(argc, argv);
  // okay
  return FEAT::Runtime::finalise();
}
