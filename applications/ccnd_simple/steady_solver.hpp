// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef FEAT_APPLICATIONS_CCND_SIMPLE_STEADY_SOLVER_HPP
#define FEAT_APPLICATIONS_CCND_SIMPLE_STEADY_SOLVER_HPP 1

#include "stokes_level.hpp"

#include <kernel/assembly/burgers_assembly_job.hpp>

#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/vanka.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/umfpack.hpp>

namespace CCNDSimple
{
  /**
   * \brief Steady-state Navier-Stokes solver class
   *
   * This class provides a simple monolithic steady-state Navier-Stokes solver, which can be used either directly
   * or as a base-class for custom Navier-Stokes solver implementations. This class also acts as a base-class for
   * the unsteady Navier-Stokes solving, supplying the functionality to solve the nonlinear system in each time step.
   *
   * \author Peter Zajac
   */
  class SteadySolver
  {
  public:
    /// our domain controller
    DomainControl& domain;
    /// our communicator
    const Dist::Comm& comm;

    /// solve Navier-Stokes or only Stokes?
    bool navier_stokes = true;
    /// use Newton or Picard for Navier-Stokes?
    bool newton_solver = true;
    /// deformation tensor or gradient tensor for diffusion?
    bool deform_tensor = true;
    /// use adaptive tolerance for multigrid?
    bool adaptive_tol = true;
    /// plot multigrid iterations during non-linear iteration?
    bool plot_mg_iter = false;

    /// plot Stokes solver iterations?
    bool plot_stokes = true;
    /// plot Stokes iterations header lines?
    bool plot_stokes_header = true;
    /// plot Navier-Stokes iterations?
    bool plot_navier = true;
    /// plot Navier-Stokes iterations header lines? (set to false by unsteady derived class)
    bool plot_navier_header = true;

    /// Filename to load joined solution vector from
    String load_joined_filename;
    /// Filename to save joined solution vector to
    String save_joined_filename;

    /// viscosity parameter nu
    DataType nu = DataType(1);
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
    Index smooth_steps = Index(16);
    /// damping parameter for smoother
    DataType smooth_damp = DataType(0.3);
    /// relative tolerance for linear solver
    DataType mg_tol_rel = DataType(1E-3);
    /// absolute tolerance for nonlinear solver
    DataType nl_tol_abs = DataType(1E-8);
    /// non-linear stagnation rate
    DataType nl_stag_rate = DataType(0.95);
    /// use UMPFACK as coarse grid solver?
    bool direct_coarse_solver = true;

    /// the cubature rule used for assembly
    String cubature = "gauss-legendre:3";
    /// the cubature rule used for post-processing
    String cubature_postproc = "gauss-legendre:4";

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

    // vector of all non-linear defect norms
    std::vector<DataType> nl_defs;

    // vector of all multigrid iteration counts
    std::vector<Index> mg_iters;

  public:
    /// constructor
    explicit SteadySolver(DomainControl& domain_);

    /// destructor
    virtual ~SteadySolver();

    /// non-copyable
    SteadySolver(const SteadySolver&) = delete;

    /// non-copyable
    SteadySolver& operator=(const SteadySolver&) = delete;

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
     * \param[in] inflow_part
     * The name of the in-flow meshparts, separated by whitespaces
     *
     * \param[in] inflow_formula
     * The formula of the inflow profile
     *
     * \param[in] pressure_mean
     * Specifies whether a integral-mean filter for the pressure is to be assembled or not.
     */
    virtual void assemble_boundary_conditions(const String& noflow_parts, const String& slip_parts,
      const String& inflow_parts, const String& inflow_formula, const bool pressure_mean);

    /// compiles the local system matrices for the Vanka/UMFPACK
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
     * \param[in] formulae_v
     * The formulae for the components of the velocity field.
     *
     * \param[in] formula_p
     * The formula for the pressure function.
     *
     * \returns A new synchronized and filtered Stokes vector that contains the interpolated functions.
     */
    virtual GlobalStokesVector create_sol_vector(const String& formula_v, const String& formula_p) const;

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
     * The formula for the right-hand-side of the momentum equation. May be an empty string if f = 0.
     *
     * \param[in] formula_g
     * The formula for the right-hand-side of the continuity equation. May be an empty string if g = 0.
     */
    virtual GlobalStokesVector create_rhs_vector(const String& formula_f, const String& formula_g = "") const;

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
     * \brief Computes and prints the body forces of the final solution on a mesh part
     *
     * \note This function prints out the \em raw body forces and not the drag/lift \em coefficients, which are
     * typically used in the DFG95 benchmarks. To obtain these DFG95 drag/lift coefficients,one has to multiply
     * the raw body forces by 2/(rho*U^2*D) in 2D and by 2/(rho*U^2*D*H) in 3D; the factors for the standard benchmarks
     * are given here:
     * - 2D bench1 (steady): 500
     * - 3D bench1 (steady): 50000 / 41
     * - 3D bench7 (steady): 800
     * - 2D bench2 (unsteady): 20
     * - 3D bench2 (unsteady): 2000 / 41
     *
     * \param[in] body_mesh_part
     * The name of the mesh part that represents the body on which the forces are to be computed
     *
     * \param[in] vec_sol
     * A \transient reference to the solution vector whose error is to be computed.
     *
     * \param[in] vec_rhs
     * A \transient reference to the right-hand-side vector.
     */
    void compute_body_forces(const String& body_mesh_part, const GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs);

    /// creates the linear multigrid-Vanka solver
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
    virtual bool solve_stokes(GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs);

    /**
     * \brief Solves the non-linear Navier-Stokes system
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
    virtual bool solve_navier_stokes(GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs);

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
    virtual void assemble_nonlinear_defect(GlobalStokesVector& vec_def, const GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs, bool filter_def = true);

    /**
     * \brief Sets up the Burgers assembly job for the nonlinear defect assembly.
     *
     * This function is called by assemble_nonlinear_defect().
     */
    virtual void setup_defect_burgers_job(Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType>& burgers_def_job);

    /**
     * \brief Assembles the Burgers matrix blocks A on all levels
     *
     * \param[in] vec_sol
     * The current solution vector approximation to assemble the Burgers matrices form
     */
    virtual void assemble_burgers_matrices(const GlobalStokesVector& vec_sol);

    /**
     * \brief Sets up the Burgers assembly job for the nonlinear matrix assembly.
     *
     * This function is called by assemble_burgers_matrices().
     */
    virtual void setup_matrix_burgers_job(Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>& burgers_mat_job);

    /**
     * \brief Computes the adaptive multigrid tolerance based on the non-linear defect norms
     */
    virtual void compute_adaptive_mg_tol();

    /**
     * \brief Releases the Stokes solver and frees up (almost) all internal data structures
     */
    virtual void release();
  }; // class StokesSteadySolver
} // namespace CCNDSimple

#endif // FEAT_APPLICATIONS_CCND_SIMPLE_STEADY_SOLVER_HPP
