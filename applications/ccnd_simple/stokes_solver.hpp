// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/// This header/source file pair defines the SteadySolver class, which implements the non-linear Newton-/Picard-solver
/// along with the setup of its linear Multigrid-AmaVanka preconditioner/solver. This class is also responsible for
/// assembling the non-linear Burgers matrices as well as the non-linear defects. This class also contains some helper
/// functions for the assembly of boundary condition filters, the assembly of right-hand-side vectors and the
/// interpolation of initial solution vectors. Finally, this class also offers functionality for the I/O of joined
/// solution vectors.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "base.hpp"
#include "stokes_level.hpp"

#include <kernel/solver/iterative.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/multigrid.hpp>

namespace CCNDSimple
{
  /**
   * \brief Steady-state Navier-Stokes solver class
   *
   * This class provides a simple monolithic steady-state Navier-Stokes solver, which can be used either directly
   * or as a base-class for custom Navier-Stokes solver implementations. This class can also be used to solve each
   * time step update of an unsteady (Navier-)Stokes system.
   *
   * \author Peter Zajac
   */
  class StokesSolver
  {
  public:
    /// our domain controller
    DomainControl& domain;

    /// our communicator
    const Dist::Comm& comm;

    // ----------------
    // input attributes
    // ----------------

    /// solve nonlinear Navier-Stokes or only linear Stokes?
    bool nonlinear_system = true;

    /// use Newton or Picard for non-linear iteration?
    bool newton_solver = true;

    /// deformation tensor or gradient tensor for diffusion?
    bool deform_tensor = true;

    /// use adaptive tolerance for multigrid?
    bool adaptive_tol = true;

    /// plot multigrid iterations during non-linear iteration?
    bool plot_mg_iter = false;

    /// plot linear solver iterations?
    bool plot_linear = true;

    /// plot linear iteration header lines?
    bool plot_linear_header = true;

    /// plot nonlinear iterations?
    bool plot_nonlinear = true;

    /// plot nonlinear iteration header lines?
    bool plot_nonlinear_header = true;

    /// viscosity parameter nu
    DataType nu = DataType(1);

    /// reaction parameter theta (from time stepping scheme)
    DataType theta = DataType(0);

    // streamline diffusion stabilization parameter
    DataType upsam = DataType(0);

    /// min. nonlinear solver iterations
    Index min_nl_iter = Index(1);

    /// max. nonlinear solver iterations
    Index max_nl_iter = Index(25);

    /// min. multigrid iterations
    Index min_mg_iter = Index(1);

    /// max. multigrid iterations
    Index max_mg_iter = Index(25);

    /// number of smoothing steps
    Index smooth_steps = Index(8);

    /// damping parameter for smoother
    DataType smooth_damp = DataType(0.5);

    /// relative tolerance for linear solver
    DataType mg_tol_rel = DataType(1E-3);

    /// absolute tolerance for nonlinear solver
    DataType nl_tol_abs = DataType(1E-8);

    /// non-linear stagnation rate
    DataType nl_stag_rate = DataType(0.95);

    /// use direct solver as coarse grid solver?
    bool direct_coarse_solver = true;

    /// the cubature rule used for assembly
    String cubature = "gauss-legendre:3";

    /// the cubature rule used for post-processing
    String cubature_postproc = "gauss-legendre:4";

    /// Filename to load joined solution vector from
    String load_joined_filename;

    /// Filename to save joined solution vector to
    String save_joined_filename;

    // -----------------
    // output attributes
    // -----------------

    /// vector of all non-linear defect norms (cleared upon each solver call)
    std::vector<DataType> nl_defs;

    /// vector of all multigrid iteration counts (cleared upon each solver call)
    std::vector<Index> mg_iters;

    /// plot line for short plot
    String plot_line;

    /// watch for total nonlinear solver
    StopWatch watch_nonlinear_solve;

    /// watches for matrix/vector/filter assembly
    StopWatch watch_asm_matrix, watch_asm_vector, watch_asm_filter;

    /// watches for symbolic/numeric solver initialization
    StopWatch watch_linsol_create, watch_linsol_sym, watch_linsol_num, watch_linsol_apply;

    /// total time for smoother and coarse grid solver
    double time_mg_smooth = 0.0;
    double time_mg_coarse = 0.0;

    // ----------------
    // state attributes
    // ----------------

    // our actual Stokes levels
    std::deque<std::shared_ptr<StokesLevel>> stokes_levels;

    /// deque of our AmaVanka smoothers
    std::deque<std::shared_ptr<Solver::AmaVanka<
      StokesLevel::LocalSystemMatrix,
      StokesLevel::LocalSystemFilter>>> amavankas;

    /// the multigrid hierarchy
    std::shared_ptr<Solver::MultiGridHierarchy<
      StokesLevel::GlobalSystemMatrix,
      StokesLevel::GlobalSystemFilter,
      StokesLevel::GlobalSystemTransfer>> multigrid_hierarchy;

    /// the multigrid preconditioner object
    std::shared_ptr<Solver::MultiGrid<
      StokesLevel::GlobalSystemMatrix,
      StokesLevel::GlobalSystemFilter,
      StokesLevel::GlobalSystemTransfer>> multigrid_precond;

    /// the iterative solver preconditioned with multigrid
    std::shared_ptr<Solver::IterativeSolver<StokesLevel::GlobalSystemVector>> stokes_solver;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  public:
    /// constructor
    explicit StokesSolver(DomainControl& domain_);

    /// destructor
    virtual ~StokesSolver();

    /// non-copyable
    StokesSolver(const StokesSolver&) = delete;

    /// non-copyable
    StokesSolver& operator=(const StokesSolver&) = delete;

    /// adds all supported arguments to the argument parser
    static void add_supported_args(SimpleArgParser& args);

    /// parses arguments from command line
    virtual bool parse_args(SimpleArgParser& args);

    /// prints configuration to console
    virtual void print_config() const;

    /// creates the actual levels along with gates, muxers and transfers
    virtual void create_levels();

    /**
     * \brief Assembles the in-flow, no-flow, slip and pressure mean boundary conditions
     *
     * \param[in] noflow_parts
     * The name of the no-flow meshparts, separated by whitespaces.
     *
     * \param[in] slip_parts
     * The name of the slip meshparts, separated by whitespaces.
     *
     * \param[in] flow_part
     * The name of the in-flow meshparts, separated by whitespaces
     *
     * \param[in] flow_formula
     * The formula of the inflow profile
     *
     * \param[in] pressure_mean
     * Specifies whether a integral-mean filter for the pressure is to be assembled or not.
     */
    template<typename FlowFunc_>
    void assemble_boundary_conditions(const String& noflow_parts, const String& slip_parts,
      const String& flow_parts, const FlowFunc_& flow_function, const bool pressure_mean)
    {
      watch_asm_filter.start();
      std::deque<String> noflow_deqs = noflow_parts.split_by_charset("|");
      std::deque<String> slip_deqs = slip_parts.split_by_charset("|");
      std::deque<String> flow_deqs = flow_parts.split_by_charset("|");

      String pres_cub = "auto-degree:" + stringify(SpacePresType::local_degree+1);

      for(std::size_t i(0); i < stokes_levels.size(); ++i)
      {
        for(auto it = noflow_deqs.begin(); it != noflow_deqs.end(); ++it)
          stokes_levels.at(i)->assemble_noflow_bc(*it);
        for(auto it = slip_deqs.begin(); it != slip_deqs.end(); ++it)
          stokes_levels.at(i)->assemble_slip_bc(*it);
        for(auto it = flow_deqs.begin(); it != flow_deqs.end(); ++it)
          stokes_levels.at(i)->assemble_flow_bc(*it, flow_function);
        if(pressure_mean)
          stokes_levels.at(i)->assemble_pressure_mean_filter(pres_cub);
        if(!slip_deqs.empty())
          stokes_levels.at(i)->sync_velocity_slip_filters();
        stokes_levels.at(i)->compile_system_filter();
      }
      watch_asm_filter.stop();
    }


    /// compiles the local system matrices for AmaVanka/UMFPACK
    virtual void compile_local_systems();

    /// creates a new global Stokes vector on the finest vector
    virtual GlobalStokesVector create_vector() const;

    /// creates a new global Stokes solution vector on the finest level, formatted to 0 and filtered
    virtual GlobalStokesVector create_sol_vector() const;

    /// creates a new global Stokes right-hand-side vector on the finest level, formatted to 0 and filtered
    virtual GlobalStokesVector create_rhs_vector() const;

    /**
     * \brief Creates a new global Stokes solution vector on the finest level
     *
     * This function interpolates the solution vector from two given analytical functions for the
     * velocity and the pressure.
     *
     * \param[in] function_v
     * The analytic function for the velocity field.
     *
     * \param[in] function_p
     * The analytic function for the pressure function.
     *
     * \returns A new synchronized and filtered Stokes vector that contains the interpolated functions.
     */
    template<typename VeloFunc_, typename PresFunc_>
    GlobalStokesVector create_sol_vector(const VeloFunc_& function_v, const PresFunc_& function_p) const
    {
      // create vector and format it
      GlobalStokesVector vec_sol = create_vector();
      vec_sol.format();

      // interpolate velocity function
      Assembly::Interpolator::project(vec_sol.local().at<0>(), function_v, domain.front()->space_velo);

      // interpolate pressure function
      Assembly::Interpolator::project(vec_sol.local().at<1>(), function_p, domain.front()->space_pres);

      // synchronize and filter
      vec_sol.sync_1();
      stokes_levels.front()->filter_sys.filter_sol(vec_sol);

      return vec_sol;
    }

    /**
     * \brief Loads an initial solution from a joined vector file
     *
     * This function loads the initial solution vector from a file, if the user has specified such
     * a file via the --load-joined-sol <filename> command line option.
     *
     * \param[inout] vec_sol
     * A \transient reference to the initial solution vector.
     *
     * \returns \c true, if the vector was loaded or \c false, if the user did not ask for this.
     */
    virtual bool load_joined_sol_vector(GlobalStokesVector& vec_sol);

    /**
     * \brief Saves the final solution to a joined vector file
     *
     * This function saves the final solution vector to a file, if the user has specified such
     * a file via the --save-joined-sol <filename> command line option.
     *
     * \param[in] vec_sol
     * A \transient reference to the final solution vector.
     */
    virtual void save_joined_sol_vector(const GlobalStokesVector& vec_sol);

    /**
     * \brief Creates a new global Stokes right-hand-side vector on the finest level
     *
     * \param[in] formula_f
     * The analytic function for the right-hand-side of the momentum equation.
     *
     * \param[in] formula_g
     * The analytic function for the right-hand-side of the continuity equation.
     */
    template<typename FuncF_, typename FuncG_>
    GlobalStokesVector create_rhs_vector(const FuncF_& function_f, const FuncG_& function_g) const
    {
      // create vector and format it
      GlobalStokesVector vec_rhs = create_vector();
      vec_rhs.format();

      // assemble force functional for velocity
      Assembly::assemble_force_function_vector(domain.front()->domain_asm, vec_rhs.local().at<0>(),
        function_f, domain.front()->space_velo, cubature);

      // assemble force functional for pressure
      Assembly::assemble_force_function_vector(domain.front()->domain_asm, vec_rhs.local().at<1>(),
        function_g, domain.front()->space_pres, cubature);

      // synchronize and filter
      vec_rhs.sync_0();
      stokes_levels.front()->filter_sys.filter_rhs(vec_rhs);

      return vec_rhs;
    }

    /**
     * \brief Creates a new global Stokes right-hand-side vector on the finest level
     *
     * \param[in] formula_f
     * The analytic function for the right-hand-side of the momentum equation.
     */
    template<typename FuncF_>
    GlobalStokesVector create_rhs_vector(const FuncF_& function_f) const
    {
      // create vector and format it
      GlobalStokesVector vec_rhs = create_vector();
      vec_rhs.format();

      // assemble force functional for velocity
      Assembly::assemble_force_function_vector(domain.front()->domain_asm, vec_rhs.local().at<0>(),
        function_f, domain.front()->space_velo, cubature);

      // synchronize and filter
      vec_rhs.sync_0();
      stokes_levels.front()->filter_sys.filter_rhs(vec_rhs);

      return vec_rhs;
    }

    /**
     * \brief Computes and prints the H^k errors of the solution against a reference solution
     *
     * \param[in] vec_sol
     * A \transient reference to the solution vector whose error is to be computed.
     *
     * \param[in] function_v
     * An analytic function representing the reference velocity field.
     *
     * \param[in] function_p
     * An analytic function representing the reference pressure function.
     */
    template<typename VeloFunc_, typename PresFunc_>
    void compute_errors(const GlobalStokesVector& vec_sol, const VeloFunc_& function_v, const PresFunc_& function_p) const
    {
      // compute velocity errors
      auto err_v = Assembly::integrate_error_function<2>(domain.front()->domain_asm, function_v,
        vec_sol.local().template at<0>(), domain.front()->space_velo, this->cubature_postproc);
      err_v.synchronize(*this->stokes_levels.front()->gate_sys.get_comm());

      // compute pressure errors
      auto err_p = Assembly::integrate_error_function<1>(domain.front()->domain_asm, function_p,
        vec_sol.local().template at<1>(), domain.front()->space_pres, this->cubature_postproc);
      err_p.synchronize(*this->stokes_levels.front()->gate_sys.get_comm());

      Tiny::Vector<DataType, 5> ret;
      ret[0] = Math::sqrt(err_v.norm_h0_sqr);
      ret[1] = Math::sqrt(err_v.norm_h1_sqr);
      ret[2] = Math::sqrt(err_v.norm_h2_sqr);
      ret[3] = Math::sqrt(err_p.norm_h0_sqr);
      ret[4] = Math::sqrt(err_p.norm_h1_sqr);

      String line = "Velocity [H0/H1,H2] / Pressure Errors [H0/H1]:";
      for(int i(0); i < 5; ++i)
        line += stringify_fp_sci(ret[i], 7, 15);
      comm.print(line);
    }

    /**
     * \brief Computes and prints the body forces of the final solution on a specific mesh part
     *
     * \param[in] body_mesh_part
     * The name of the mesh part that represents the body on which the forces are to be computed
     *
     * \param[in] vec_sol
     * A \transient reference to the solution vector whose error is to be computed.
     *
     * \param[in] vec_rhs
     * A \transient reference to the right-hand-side vector.
     *
     * \param[in] scaling_factor
     * A problem-dependent scaling factor for the body forces.
     */
    void compute_body_forces(const String& body_mesh_part, const GlobalStokesVector& vec_sol,
      const GlobalStokesVector& vec_rhs, const DataType scaling_factor = DataType(1));

    /// creates the linear multigrid-AmaVanka solver
    virtual void create_multigrid_solver();

    /**
     * \brief Solves the linear Stokes system
     *
     * \param[inout] vec_sol
     * The solution vector to be updated.
     *
     * \param[in] vec_rhs
     * The right-hand-side vector of the stokes equations
     *
     * \returns
     * true, if the solution was successful, otherwise false
     */
    virtual bool solve_linear(GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs);

    /**
     * \brief Solves the nonlinear Navier-Stokes system
     *
     * \param[inout] vec_sol
     * The solution vector to be updated.
     *
     * \param[in] vec_rhs
     * The right-hand-side vector of the Navier-Stokes equations
     *
     * \returns
     * true, if the solution was successful, otherwise false
     */
    virtual bool solve_nonlinear(GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs);

    /**
     * \brief Assembles the nonlinear defect vector for the Navier-Stokes equations
     *
     * \param[inout] vec_def
     * The defect vector to be assembled. Must be created to correct size, but its contents are irrelevant.
     *
     * \param[in] vec_sol
     * The current solution vector approximation to compute the defect from.
     *
     * \param[in] vec_rhs
     * The right-hand-side vector to compute the defect from.
     *
     * \param[in] filter_def
     * Specifies whether the defect vector is to be filtered or not.
     */
    virtual void assemble_nonlinear_defect(GlobalStokesVector& vec_def, const GlobalStokesVector& vec_sol,
      const GlobalStokesVector& vec_rhs, bool filter_def = true);

    /**
     * \brief Assembles the Burgers matrix blocks A on all levels
     *
     * \param[in] vec_sol
     * The current solution vector approximation to assemble the Burgers matrices form
     */
    virtual void assemble_burgers_matrices(const GlobalStokesVector& vec_sol);

    /**
     * \brief Computes the adaptive multigrid tolerance based on the non-linear defect norms
     */
    virtual void compute_adaptive_mg_tol();

    /**
     * \brief Releases the Stokes solver and frees up (almost) all internal data structures
     */
    virtual void release();

    /// prints runtime information for this object
    virtual void print_runtime(double total_time);

  protected:
    /**
     * \brief Assembles the local Burgers defect vector on a level
     *
     * This method is called by #assemble_nonlinear_defect() to assemble the Burgers part of the
     * defect vector.
     *
     * \param[inout] vec_def
     * The defect vector to be assembled; its contents are already initialized.
     *
     * \param[in] vec_velo
     * The velocity vector for the convective part of the Burgers operator.
     *
     * \param[in] domain_lyr
     * The domain layer that the domain level belongs to.
     *
     * \param[in] domain_lvl
     * The domain level to assemble on; this is typically the finest domain level.
     *
     * \param[in] stokes_lvl
     * The Stokes level to assemble on; this is typically the finest Stokes level.
     */
    virtual void _assemble_local_burgers_defect(LocalVeloVector& vec_def, const LocalVeloVector& vec_velo,
      const DomainLayer& domain_lyr, DomainLevel& domain_lvl, StokesLevel& stokes_lvl);

    /**
     * \brief Assembles the local Burgers matrix on a level
     *
     * This method is called by #assemble_burgers_matrices() to assemble the Burgers part of the
     * matrix on a given level.
     *
     * \param[inout] matrix_a
     * The Burgers matrix to be assembled; its contents are already formatted/initialized.
     *
     * \param[in] vec_velo
     * The velocity vector for the convective part of the Burgers operator.
     *
     * \param[in] domain_lyr
     * The domain layer that the domain level belongs to.
     *
     * \param[in] domain_lvl
     * The domain level to assemble on
     *
     * \param[in] stokes_lvl
     * The Stokes level to assemble on
     */
    virtual void _assemble_local_burgers_matrix(LocalMatrixBlockA& matrix_a, const LocalVeloVector& vec_velo,
      const DomainLayer& domain_lyr, DomainLevel& domain_lvl, StokesLevel& stokes_lvl);

  }; // class StokesSolver
} // namespace CCNDSimple
