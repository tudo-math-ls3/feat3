// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "steady_solver.hpp"

#include <kernel/analytic/parsed_function.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/domain_assembler_helpers.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/solver/direct_stokes_solver.hpp>


namespace CCNDSimple
{

  SteadySolver::SteadySolver(DomainControl& domain_) :
    domain(domain_),
    comm(domain.comm())
  {
  }

  SteadySolver::~SteadySolver()
  {
  }

  void SteadySolver::add_supported_args(SimpleArgParser& args)
  {
    args.support("nu", "<nu>\nSets the viscosity parameter nu.");
    args.support("upsam", "<upsam>\nSets the streamline diffusion stabilization parameter.");
    args.support("max-nl-iter", "<iter>\nSets the maximum number of nonlinear iterations.");
    args.support("min-mg-iter", "<iter>\nSets the minimum number of multigrid iterations.");
    args.support("max-mg-iter", "<iter>\nSets the maximum number of multigrid iterations.");
    args.support("plot-mg-iter", "\nIf given, enabled the output of multigrid iterations to the console during the nonlinear solve.");
    args.support("smooth-steps", "<steps>\nSets the number of Vanka smoothing steps.");
    args.support("smooth-damp", "<damp>\nSets the damping parameter of the Vanka smoother.");
    args.support("mg-tol-rel", "<tol>\nSets the relative tolerance for the multigrid solver.");
    args.support("nl-tol-abs", "<tol>\nSets the absolute tolerance for the nonlinear solver.");
    args.support("nl-stag-rate", "<rate>\nSets the defect stagnation rate for the nonlinear solver.");
    args.support("no-navier", "\nIf given, solve only the linear Stokes equation instead of Navier-Stokes.");
    args.support("no-newton", "\nIf given, use Picard iteration instead of Newton for the nonlinear solve.");
    args.support("no-deform", "\nIf given, use the simple gradient tensor instead of the deformation tensor.");
    args.support("no-umfpack", "\nIf given, use BiCGStab-AmaVanka instead of UMFPACK as a coarse grid solver.");
    args.support("load-joined-sol", "<filename>\nLoads the initial solution vector from a joined vector file.");
    args.support("save-joined-sol", "<filename>\nSaves the final solution vector to a joined vector file.");
  }

  bool SteadySolver::parse_args(SimpleArgParser& args)
  {
    // solve Navier-Stokes or only Stokes?
    navier_stokes = (args.check("no-navier") < 0);
    // use Newton or Picard for Navier-Stokes?
    newton_solver = (args.check("no-newton") < 0);
    // deformation tensor or gradient tensor for diffusion?
    deform_tensor = (args.check("no-deform") < 0);
    // use adaptive tolerance for multigrid?
    adaptive_tol = (args.check("mg-tol-rel") < 0);
    // plot multigrid iterations?
    plot_mg_iter = (args.check("plot-mg-iter") >= 0);

    // viscosity parameter nu
    args.parse("nu", nu);
    // streamline diffusion stabilization
    args.parse("upsam", upsam);
    // max. nonlinear solver iterations
    args.parse("max-nl-iter", max_nl_iter);
    // min. multigrid iterations
    args.parse("min-mg-iter", min_mg_iter);
    // max. multigrid iterations
    args.parse("max-mg-iter", max_mg_iter);
    // number of smoothing steps
    args.parse("smooth-steps", smooth_steps);
    // damping parameter for smoother
    args.parse("smooth-damp", smooth_damp);
    // rel. tolerance for linear solver
    args.parse("mg-tol-rel", mg_tol_rel);
    // absolute tolerance for nonlinear solver
    args.parse("nl-tol-abs", nl_tol_abs);
    // stagnation rate
    args.parse("nl-stag-rate", nl_stag_rate);

    // filename of joined solution vector to read in
    args.parse("load-joined-sol", load_joined_filename);
    // filename of joined solution vector to write out
    args.parse("save-joined-sol", save_joined_filename);

    // use UMPFACK as coarse grid solver?
#ifdef FEAT_HAVE_UMFPACK
    direct_coarse_solver = (domain.back_layer().comm().size() == 1) && (args.check("no-umfpack") < 0);
#else
    direct_coarse_solver = false;
#endif

    return true;
  }

  void SteadySolver::print_config() const
  {
    print_pad(comm, "Nu", stringify(nu));
    print_pad(comm, "upsam", stringify(upsam));
    print_pad(comm, "System", (navier_stokes ? "Navier-Stokes" : "Stokes"));
    print_pad(comm, "Diffusion Tensor", (deform_tensor ? "Deformation" : "Gradient"));
    print_pad(comm, "Nonlinear Solver", (newton_solver ? "Newton" : "Picard"));
    print_pad(comm, "AmaVanka Smoother Steps", stringify(smooth_steps));
    print_pad(comm, "AmaVanka Smoother Damping", stringify(smooth_damp));
    print_pad(comm, "Nonlinear Stagnation Rate", stringify_fp_sci(nl_stag_rate));
    print_pad(comm, "Nonlinear Absolute Tol", stringify_fp_sci(nl_tol_abs));
    print_pad(comm, "Multigrid Relative Tol", (adaptive_tol ? String("adaptive") : stringify_fp_sci(mg_tol_rel)));
    print_pad(comm, "Min Multigrid Iterations", stringify(min_mg_iter));
    print_pad(comm, "Max Multigrid Iterations", stringify(max_mg_iter));
    print_pad(comm, "Max Nonlinear Iterations", stringify(max_nl_iter));
    print_pad(comm, "Coarse Solver", (direct_coarse_solver ? "UMFPACK" : "BiCGStab-AmaVanka"));
    if(load_joined_filename.empty())
      print_pad(comm, "Load Joined Solution", "- N/A -");
    else
      print_pad(comm, "Load Joined Solution", load_joined_filename);
    if(save_joined_filename.empty())
      print_pad(comm, "Save Joined Solution", "- N/A -");
    else
      print_pad(comm, "Save Joined Solution", save_joined_filename);
  }

  void SteadySolver::create_levels()
  {
    bool assemble_splitters = (!load_joined_filename.empty() || !save_joined_filename.empty());
    const Index num_levels = Index(domain.size_physical());

    // create Stokes levels
    for (Index i(0); i < num_levels; ++i)
    {
      stokes_levels.push_back(std::make_shared<StokesLevel>(*domain.at(i)));

      stokes_levels.at(i)->assemble_gates(domain.at(i));

      if(assemble_splitters)
        stokes_levels.at(i)->assemble_base_splitters(domain.at(i));
    }

    // assemble muxers and transfers
    for (Index i(0); (i < domain.size_physical()) && ((i+1) < domain.size_virtual()); ++i)
    {
      stokes_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
      if((i+1) < domain.size_physical())
        stokes_levels.at(i)->assemble_transfers(*stokes_levels.at(i+1), domain.at(i), domain.at(i+1), cubature, true);
      else
        stokes_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature, true);
    }

    // assemble basic stokes matrices on all levels
    for(Index i(0); i < num_levels; ++i)
    {
      stokes_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
      stokes_levels.at(i)->assemble_pres_struct(domain.at(i)->space_pres);
      stokes_levels.at(i)->assemble_velocity_laplace_matrix(cubature, nu, deform_tensor);
      stokes_levels.at(i)->assemble_grad_div_matrices(cubature);
      stokes_levels.at(i)->compile_system_matrix();
    }
  }

  void SteadySolver::assemble_boundary_conditions(const String& noflow_parts, const String& slip_parts,
    const String& inflow_parts, const String& inflow_formula, const bool pressure_mean)
  {
    std::deque<String> noflow_deqs = noflow_parts.split_by_charset("|");
    std::deque<String> slip_deqs = slip_parts.split_by_charset("|");
    std::deque<String> inflow_deqs = inflow_parts.split_by_charset("|");

    String pres_cub = "auto-degree:" + stringify(SpacePresType::local_degree+1);

    for(std::size_t i(0); i < stokes_levels.size(); ++i)
    {
      for(auto it = noflow_deqs.begin(); it != noflow_deqs.end(); ++it)
        stokes_levels.at(i)->assemble_noflow_bc(*it);
      for(auto it = slip_deqs.begin(); it != slip_deqs.end(); ++it)
        stokes_levels.at(i)->assemble_slip_bc(*it);
      for(auto it = inflow_deqs.begin(); it != inflow_deqs.end(); ++it)
        stokes_levels.at(i)->assemble_inflow_bc(*it, inflow_formula);
      if(pressure_mean)
        stokes_levels.at(i)->assemble_pressure_mean_filter(pres_cub);
      if(!slip_deqs.empty())
        stokes_levels.at(i)->sync_velocity_slip_filters();
      stokes_levels.at(i)->compile_system_filter();
    }
  }

  void SteadySolver::compile_local_systems()
  {
    for(std::size_t i(0); i < stokes_levels.size(); ++i)
    {
      //stokes_levels.at(i)->compile_system_filter();
      stokes_levels.at(i)->compile_local_matrix();
    }
  }

  GlobalStokesVector SteadySolver::create_vector() const
  {
    XASSERT(!stokes_levels.empty());
    return stokes_levels.front()->matrix_sys.create_vector_l();
  }

  GlobalStokesVector SteadySolver::create_sol_vector() const
  {
    GlobalStokesVector vec_sol = create_vector();
    vec_sol.format();
    stokes_levels.front()->filter_sys.filter_sol(vec_sol);
    return vec_sol;
  }

  GlobalStokesVector SteadySolver::create_rhs_vector() const
  {
    GlobalStokesVector vec_rhs = create_vector();
    vec_rhs.format();
    stokes_levels.front()->filter_sys.filter_rhs(vec_rhs);
    return vec_rhs;
  }

  GlobalStokesVector SteadySolver::create_sol_vector(const String& formula_v, const String& formula_p) const
  {
    // create vector and format it
    GlobalStokesVector vec_sol = create_vector();
    vec_sol.format();

    // interpolate velocity function
    if(!formula_v.empty())
    {
      Analytic::ParsedVectorFunction<dim, dim> function_v(formula_v);
      Assembly::Interpolator::project(vec_sol.local().at<0>(), function_v, domain.front()->space_velo);
    }

    // interpolate pressure function
    if(!formula_p.empty())
    {
      Analytic::ParsedScalarFunction<dim> function_p(formula_p);
      Assembly::Interpolator::project(vec_sol.local().at<1>(), function_p, domain.front()->space_pres);
    }

    // synchronize and filter
    vec_sol.sync_1();
    stokes_levels.front()->filter_sys.filter_sol(vec_sol);

    return vec_sol;
  }

  GlobalStokesVector SteadySolver::create_rhs_vector(const String& formula_f, const String& formula_g) const
  {
    // create vector and format it
    GlobalStokesVector vec_rhs = create_vector();
    vec_rhs.format();

    // assemble force functional for velocity
    if(!formula_f.empty())
    {
      Analytic::ParsedVectorFunction<dim, dim> function_f(formula_f);
      Assembly::assemble_force_function_vector(domain.front()->domain_asm, vec_rhs.local().at<0>(),
        function_f, domain.front()->space_velo, this->cubature);
    }

    // assemble force functional for pressure
    if(!formula_g.empty())
    {
      Analytic::ParsedScalarFunction<dim> function_g(formula_g);
      Assembly::assemble_force_function_vector(domain.front()->domain_asm, vec_rhs.local().at<1>(),
        function_g, domain.front()->space_pres, this->cubature);
    }

    // synchronize and filter
    vec_rhs.sync_0();
    stokes_levels.front()->filter_sys.filter_rhs(vec_rhs);

    return vec_rhs;
  }

  bool SteadySolver::load_joined_sol_vector(GlobalStokesVector& vec_sol)
  {
    if(load_joined_filename.empty())
      return false;

    comm.print("\nLoading joined solution vector from '" + load_joined_filename + "'...");

    // read and filter
    stokes_levels.front()->base_splitter_sys.split_read_from(vec_sol, load_joined_filename);
    stokes_levels.front()->filter_sys.filter_sol(vec_sol);
    return true;
  }

  void SteadySolver::save_joined_sol_vector(const GlobalStokesVector& vec_sol)
  {
    if(save_joined_filename.empty())
      return;

    comm.print("\nSaving joined solution vector to '" + save_joined_filename + "'...");

    stokes_levels.front()->base_splitter_sys.join_write_out(vec_sol, save_joined_filename);
  }

  void SteadySolver::compute_body_forces(const String& body_mesh_part, const GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs)
  {
    // our local body forces
    Tiny::Vector<DataType, dim> body_forces;
    body_forces.format();

    // fetch the body mesh part
    DomainLevel& the_domain_level = *domain.front();
    const auto* body_part = the_domain_level.get_mesh_node()->find_mesh_part(body_mesh_part);
    if(body_part != nullptr)
    {
      // assemble a unit filter on that mesh-part
      LAFEM::UnitFilter<DataType, IndexType> filter_char;
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      unit_asm.add_mesh_part(*body_part);
      unit_asm.assemble(filter_char, the_domain_level.space_velo);
      LAFEM::SparseVector<DataType, IndexType>& vec_char = filter_char.get_filter_vector();

      // clone the RHS and convert to type 0
      GlobalStokesVector vec_def = vec_rhs.clone(LAFEM::CloneMode::Deep);
      vec_def.from_1_to_0();

      // set up Burgers assembly job for our defect vector
      Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType> burgers_def_job(
        vec_def.local().template at<0>(), vec_sol.local().template at<0>(),
        vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature);

      // call the setup function to set the parameters
      this->setup_defect_burgers_job(burgers_def_job);

      // assemble burgers operator defect
      the_domain_level.domain_asm.assemble(burgers_def_job);

      // compute remainder of velocity defect vector by using matrix B
      stokes_levels.front()->matrix_sys.local().block_b().apply(
        vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -DataType(1));

      // compute dot product of char vector and defect vector
      const auto* velo_def = vec_def.local().first().elements();
      const IndexType* idx = vec_char.indices();
      for(Index i(0); i < vec_char.used_elements(); ++i)
        body_forces += velo_def[idx[i]];
    }

    // synchronize body forces by summing them up
    domain.front().layer().comm().allreduce(body_forces.v, body_forces.v, std::size_t(dim), Dist::op_sum);

    // stringify our body forces
    String out;
    for(int i(0); i < dim; ++i)
      out += stringify_fp_sci(body_forces[i], 16, 25);

    // print body forces
    print_pad(comm, "Body Forces on " + body_mesh_part, out);
  }

  void SteadySolver::create_multigrid_solver()
  {
    // create a multigrid hierarchy
    this->multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<
      StokesLevel::GlobalSystemMatrix,
      StokesLevel::GlobalSystemFilter,
      StokesLevel::GlobalSystemTransfer>>(domain.size_virtual());

    // push levels into multigrid
    for(std::size_t i(0); i < stokes_levels.size(); ++i)
    {
      StokesLevel& lvl = *stokes_levels.at(i);

      if((i+1) < domain.size_virtual())
      {
        auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
        vanka->set_skip_singular(true);
        amavankas.push_back(vanka);
        auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
        schwarz->set_ignore_status(true);
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);
        smoother->set_min_iter(smooth_steps);
        smoother->set_max_iter(smooth_steps);
        smoother->skip_defect_calc(true); // skip defect calculation
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
#ifdef FEAT_HAVE_UMFPACK
      else if(direct_coarse_solver)
      {
        auto cgsolver = Solver::new_direct_stokes_solver(lvl.matrix_sys, lvl.filter_sys);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
#endif //  FEAT_HAVE_UMFPACK
      else
      {
        // create BiCGStab-AmaVanka coarse grid solver
        auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
        vanka->set_skip_singular(true);
        amavankas.push_back(vanka);
        auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
        schwarz->set_ignore_status(true);
        auto cgsolver = Solver::new_bicgstab(lvl.matrix_sys, lvl.filter_sys, schwarz);
        cgsolver->set_max_iter(500);
        cgsolver->set_tol_rel(1e-3);
        //cgsolver->set_plot_mode(Solver::PlotMode::summary);

        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
    }

    // create our multigrid solver
    this->multigrid_precond = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

    // create our solver
    this->stokes_solver = Solver::new_richardson(stokes_levels.front()->matrix_sys,
      stokes_levels.front()->filter_sys, 1.0, this->multigrid_precond);

    stokes_solver->set_plot_name("Multigrid");
    stokes_solver->set_min_iter(min_mg_iter);
    stokes_solver->set_max_iter(max_mg_iter);
    stokes_solver->set_tol_rel(mg_tol_rel);
    stokes_solver->set_min_stag_iter(Index(3));

    // initialize the solver symbolically
    this->multigrid_hierarchy->init_symbolic();
    this->stokes_solver->init_symbolic();
  }

  bool SteadySolver::solve_stokes(GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs)
  {
    // set plot mode to enable output
    stokes_solver->set_plot_mode(plot_stokes ? Solver::PlotMode::iter : Solver::PlotMode::none);

    // initialize solver
    this->multigrid_hierarchy->init_numeric();
    this->stokes_solver->init_numeric();

    // print header
    if(plot_stokes_header)
    {
      //         "Multigrid:  2 : 5.012506e-03 / 4.283371e-04 / 0.028025"
      comm.print(String(this->stokes_solver->get_plot_name().size(), ' ') +
        "  It   Abs. Def.      Rel. Def.      Improve");
      comm.print(String(this->stokes_solver->get_plot_name().size(), '-') +
        "---------------------------------------------");
    }

    // apply solver
    Solver::Status status = stokes_solver->correct(vec_sol, vec_rhs);

    // release solver
    this->multigrid_hierarchy->done_numeric();
    this->stokes_solver->done_numeric();

    // check solver output
    if(!Solver::status_success(status))
    {
      comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\nStatus: " + stringify(status));
      return false;
    }

    return true;
  }

  bool SteadySolver::solve_navier_stokes(GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs)
  {
    // set plot mode to enable output
    stokes_solver->set_plot_mode(plot_navier ? Solver::PlotMode::iter : Solver::PlotMode::none);

    // create two vectors
    GlobalStokesVector vec_cor = create_vector();
    GlobalStokesVector vec_def = create_vector();

    // print header
    if(plot_navier_header)
    {
      //         "Newton:   1: 2.373388e-05 / 4.732559e-03 / 4.733e-03 |   5: 8.7702e-12 / 3.6952e-07 [1.0000e-10]"
      comm.print("         It  Abs. Defect    Rel. Defect    Improve   |  It  Abs. Def.    Rel. Def.   Rel. Tol.");
      comm.print("-----------------------------------------------------+------------------------------------------");
    }

    // nonlinear loop
    for(Index nl_step(0); nl_step <= max_nl_iter; ++nl_step)
    {
      // set up Burgers assembly job for our defect vector
      assemble_nonlinear_defect(vec_def, vec_sol, vec_rhs);

      // compute defect norm
      const DataType def_prev = (nl_defs.empty() ? DataType(1) : nl_defs.back());
      const DataType def_nl = vec_def.norm2();
      const DataType def_improve = (nl_defs.empty() ? DataType(1) : def_nl / def_prev);
      nl_defs.push_back(def_nl);

      // start building nonlinear solver output line
      String line = (newton_solver ? "Newton:" : "Picard:");
      if(plot_navier)
      {
        line += stringify(nl_step).pad_front(4) + ": ";
        line += stringify_fp_sci(def_nl, 6) + " / ";
        line += stringify_fp_sci(def_nl/nl_defs.front()) + " / ";
        line += stringify_fp_sci(def_improve, 3);
      }

      // check for divergence
      if(def_nl > nl_defs.front() * DataType(1E+3))
      {
        comm.print(line + "\n\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
        return false;
      }

      // check for stopping criterions
      if(nl_step >= min_nl_iter)
      {
        if(def_nl < nl_tol_abs)
        {
          if(plot_navier)
            comm.print(line + "\n\nNonlinear solver converged!\n");
          break;
        }
        else if(nl_step >= max_nl_iter)
        {
          if(plot_navier)
            comm.print(line + "\n\nMaximum iterations reached!\n");
          break;
        }
        else if((nl_step >= 3) && (nl_stag_rate*def_prev < def_nl))
        {
          if(plot_navier)
            comm.print(line + "\n\nNonlinear solver stagnated!\n");
          break;
        }
      }

      // reassemble matrix
      assemble_burgers_matrices(vec_sol);

      // compile local systems
      compile_local_systems();

      // specify adaptive tolerance?
      if(adaptive_tol && (nl_step > Index(0)))
        compute_adaptive_mg_tol();
      else
      {
        // disable absolute tolerance
        stokes_solver->set_tol_abs(1E+10);
        stokes_solver->set_tol_rel(mg_tol_rel);
      }

      // initialize linear solver
      multigrid_hierarchy->init_numeric();
      stokes_solver->init_numeric();

      // solve linear system
      Solver::Status status = stokes_solver->apply(vec_cor, vec_def);

      // save MG iterations
      mg_iters.push_back(stokes_solver->get_num_iter());

      // build output line and print it
      if(plot_navier)
      {
        line += String(" | ") + stringify(stokes_solver->get_num_iter()).pad_front(3) + ": "
          + stringify_fp_sci(stokes_solver->get_def_final(), 4) + " / "
          + stringify_fp_sci(stokes_solver->get_def_final() / stokes_solver->get_def_initial(), 4);
        if(adaptive_tol && (nl_step > Index(0)))
          line += String(" [") + stringify_fp_sci(stokes_solver->get_tol_abs(), 4) + "]";
        comm.print(line);
      }

      // release linear solver
      stokes_solver->done_numeric();
      multigrid_hierarchy->done_numeric();

      // check linear solver status
      if(!Solver::status_success(status))
      {
        comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\nStatus: " + stringify(status));
        return false;
      }

      // update solution
      vec_sol.axpy(vec_cor);

      // next non-linear iteration
    }
    // end of Navier-Stokes solve

    return true;
  }

  void SteadySolver::assemble_nonlinear_defect(GlobalStokesVector& vec_def, const GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs, bool filter_def)
  {
    // fetch our finest levels
    DomainLevel& the_domain_level = *domain.front();
    StokesLevel& the_stokes_level = *stokes_levels.front();

    // set up Burgers assembly job for our defect vector
    Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType> burgers_def_job(
      vec_def.local().template at<0>(), vec_sol.local().template at<0>(),
      vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature);

    // call the setup function to set the parameters
    this->setup_defect_burgers_job(burgers_def_job);

    // copy rhs vector and convert it to a 0-type vector
    vec_def.copy(vec_rhs);
    vec_def.from_1_to_0();

    // assemble burgers operator defect
    the_domain_level.domain_asm.assemble(burgers_def_job);

    // compute remainder of defect vector by using matrices B and D
    the_stokes_level.matrix_sys.local().block_b().apply(
      vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
    the_stokes_level.matrix_sys.local().block_d().apply(
      vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);

    // sync and filter
    vec_def.sync_0();
    if(filter_def)
      the_stokes_level.filter_sys.filter_def(vec_def);
  }

  void SteadySolver::setup_defect_burgers_job(Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType>& burgers_def_job)
  {
    burgers_def_job.deformation = deform_tensor;
    burgers_def_job.nu = -nu;
    burgers_def_job.beta = DataType(navier_stokes ? -1 : 0);
  }

  void SteadySolver::assemble_burgers_matrices(const GlobalStokesVector& vec_sol)
  {
    // fetch our finest levels
    StokesLevel& the_stokes_level = *stokes_levels.front();

    // get a clone of the global velocity vector
    // this one will be truncated as we go down the level hierarchy
    GlobalVeloVector vec_conv(&the_stokes_level.gate_velo, vec_sol.local().template at<0>().clone());

    // initialize velocity norm for streamline diffusion (if enabled)
    DataType sd_v_norm = DataType(0);

    // loop over all system levels
    for(std::size_t i(0); i < stokes_levels.size(); ++i)
    {
      // get the local matrix block A on this level
      LocalMatrixBlockA& loc_mat_a = stokes_levels.at(i)->matrix_sys.local().block_a();
      loc_mat_a.format();

      // set up Burgers assembly job
      Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>
        burgers_mat_job(loc_mat_a, vec_conv.local(), domain.at(i)->space_velo, cubature);

      // call the setup function to set the parameters
      setup_matrix_burgers_job(burgers_mat_job);
      if(i == size_t(0))
      {
        burgers_mat_job.set_sd_v_norm(vec_conv);
        sd_v_norm = burgers_mat_job.sd_v_norm;
      }
      else
      {
        // use fine mesh norm
        burgers_mat_job.sd_v_norm = sd_v_norm;
      }

      // assemble our matrix block
      domain.at(i)->domain_asm.assemble(burgers_mat_job);

      // compile the local matrix for the Vanka smoother / UMFPACK solver
      stokes_levels.at(i)->compile_local_matrix();

      // are we there yet?
      if((i+1) >= domain.size_virtual())
        break;

      // does this process have another system level?
      if((i+1) < stokes_levels.size())
      {
        // create a coarse mesh velocity vector
        auto vec_crs = stokes_levels.at(i+1)->matrix_a.create_vector_l();

        // truncate fine mesh velocity vector
        stokes_levels.at(i)->transfer_velo.trunc(vec_conv, vec_crs);

        // the coarse vector is our next convection vector
        vec_conv = std::move(vec_crs);
      }
      else
      {
        // this process is a child, so send truncation to parent
        stokes_levels.at(i)->transfer_velo.trunc_send(vec_conv);
      }
    }
  }

  void SteadySolver::setup_matrix_burgers_job(Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>& burgers_mat_job)
  {
    burgers_mat_job.deformation = deform_tensor;
    burgers_mat_job.nu = nu;
    burgers_mat_job.beta = DataType(navier_stokes ? 1 : 0);
    burgers_mat_job.frechet_beta = DataType(navier_stokes && newton_solver ? 1 : 0);
    burgers_mat_job.sd_delta = upsam;
    burgers_mat_job.sd_nu = nu;
  }


  void SteadySolver::compute_adaptive_mg_tol()
  {
    XASSERT(nl_defs.size() >= std::size_t(2));

    auto it = nl_defs.end();
    const DataType def_nl = *(--it);
    const DataType def_prev = *(--it);
    const DataType def_improve = def_nl / def_prev;

    if(newton_solver)
    {
      // We're using Newton as the nonlinear solver, which optimally should
      // result in quadratic convergence, i.e. let def_{j} and def_{j-1}
      // denote the two previous nonlinear defect norms, then the next
      // defect norm def_{j+1} should fulfill
      //
      //     (def_{j+1} / def_{j}) \approx (def_{j} / def_{j+1})^2
      //
      // If the multiply the above equation by def_{j}, we can therefore
      // estimate the next def_{j+1} based on the two previous defects:
      //
      //      def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})^2
      //
      // We now multiply the approximation by 0.1, which gives us an absolute
      // tolerance for the multigrid solver for this nonlinear iteration.
      // (Note that def_improve := def_{j} / def_{j+1})
      DataType abs_tol = def_nl * def_improve * def_improve * DataType(0.1);
      // We furthermore limit this absolute tolerance to ensure that we do not
      // overshoot the mark by overoptimistic quadratic convergence expectations.
      stokes_solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
      // Also make sure that we gain at least 2 digits.
      stokes_solver->set_tol_rel(1E-2);
    }
    else
    {
      // In the case if Picard iteration, we only expect linear convergence,
      // which (in analogy to Newton) leads us to the following estimate:
      //
      //      def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})
      //
      DataType abs_tol = def_nl * def_improve * DataType(0.1);
      stokes_solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
      stokes_solver->set_tol_rel(1E-2);
    }
  }

  void SteadySolver::release()
  {
    // release solver
    stokes_solver->done_symbolic();
    multigrid_hierarchy->done_symbolic();

    // clear amavankas deque
    amavankas.clear();

    // release solvers
    stokes_solver.reset();
    multigrid_hierarchy.reset();

    // release all levels except for the finest one
    while(stokes_levels.size() > std::size_t(1))
      stokes_levels.pop_back();

    // todo: release matrices, filters, etc on finest level
  }
} // namespace CCNDSimple
