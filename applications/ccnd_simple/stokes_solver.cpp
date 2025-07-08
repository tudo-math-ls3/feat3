// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "stokes_solver.hpp"

#include <kernel/analytic/parsed_function.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
#include <kernel/assembly/burgers_assembly_job.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/solver/direct_stokes_solver.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/bicgstab.hpp>

namespace CCNDSimple
{
  StokesSolver::StokesSolver(DomainControl& domain_) :
    domain(domain_),
    comm(domain.comm())
  {
  }

  StokesSolver::~StokesSolver()
  {
  }

  void StokesSolver::add_supported_args(SimpleArgParser& args)
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
    args.support("stokes", "\nIf given, solve only the linear Stokes equation instead of nonlinear Navier-Stokes.");
    args.support("picard", "\nIf given, use Picard iteration instead of Newton for the nonlinear solve.");
    args.support("deform", "\nIf given, use the deformation tensor instead of the gradient tensor.");
    args.support("no-umfpack", "\nIf given, use BiCGStab-AmaVanka instead of UMFPACK as a coarse grid solver.");
    args.support("load-joined-sol", "<filename>\nLoads the initial solution vector from a joined vector file.");
    args.support("save-joined-sol", "<filename>\nSaves the final solution vector to a joined vector file.");
  }

  bool StokesSolver::parse_args(SimpleArgParser& args)
  {
    // solve Navier-Stokes or only Stokes?
    nonlinear_system = (args.check("stokes") < 0);
    // use Newton or Picard for Navier-Stokes?
    newton_solver = (args.check("picard") < 0);
    // deformation tensor or gradient tensor for diffusion?
    deform_tensor = (args.check("deform") >= 0);
    // use adaptive tolerance for multigrid?
    adaptive_tol = (args.check("mg-tol-rel") < 0);
    // plot multigrid iterations?
    plot_mg_iter = (args.check("plot-mg-iter") >= 0);

    // use UMPFACK as coarse grid solver?
#ifdef FEAT_HAVE_UMFPACK
    direct_coarse_solver = (domain.back_layer().comm().size() == 1) && (args.check("no-umfpack") < 0);
#else
    direct_coarse_solver = false;
#endif

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

    return true;
  }

  void StokesSolver::print_config() const
  {
    print_pad(comm, "Nu", stringify(nu));
    print_pad(comm, "upsam", stringify(upsam));
    print_pad(comm, "System", (nonlinear_system ? "Navier-Stokes" : "Stokes"));
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
    comm.print_flush();
  }

  void StokesSolver::create_levels()
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

    stokes_levels.front()->assemble_velocity_mass_matrix(cubature);
  }

  void StokesSolver::compile_local_systems()
  {
    for(std::size_t i(0); i < stokes_levels.size(); ++i)
    {
      //stokes_levels.at(i)->compile_system_filter();
      stokes_levels.at(i)->compile_local_matrix();
    }
  }

  GlobalStokesVector StokesSolver::create_vector() const
  {
    XASSERT(!stokes_levels.empty());
    return stokes_levels.front()->matrix_sys.create_vector_l();
  }

  GlobalStokesVector StokesSolver::create_sol_vector() const
  {
    GlobalStokesVector vec_sol = create_vector();
    vec_sol.format();
    stokes_levels.front()->filter_sys.filter_sol(vec_sol);
    return vec_sol;
  }

  GlobalStokesVector StokesSolver::create_rhs_vector() const
  {
    GlobalStokesVector vec_rhs = create_vector();
    vec_rhs.format();
    stokes_levels.front()->filter_sys.filter_rhs(vec_rhs);
    return vec_rhs;
  }

  bool StokesSolver::load_joined_sol_vector(GlobalStokesVector& vec_sol)
  {
    if(load_joined_filename.empty())
      return false;

    comm.print("\nLoading joined solution vector from '" + load_joined_filename + "'...");

    // read and filter
    stokes_levels.front()->base_splitter_sys.split_read_from(vec_sol, load_joined_filename);
    stokes_levels.front()->filter_sys.filter_sol(vec_sol);
    return true;
  }

  void StokesSolver::save_joined_sol_vector(const GlobalStokesVector& vec_sol)
  {
    if(save_joined_filename.empty())
      return;

    comm.print("\nSaving joined solution vector to '" + save_joined_filename + "'...");

    stokes_levels.front()->base_splitter_sys.join_write_out(vec_sol, save_joined_filename);
  }

  void StokesSolver::compute_body_forces(const String& body_mesh_part, const GlobalStokesVector& vec_sol,
    const GlobalStokesVector& vec_rhs, const DataType scaling_factor)
  {
    // our local body forces
    Tiny::Vector<DataType, dim> body_forces;
    body_forces.format();

    // clone the RHS and convert to type 0
    GlobalStokesVector vec_def = vec_rhs.clone(LAFEM::CloneMode::Layout);
    vec_def.format();

    // call the setup function to set the parameters
    this->assemble_nonlinear_defect(vec_def, vec_sol, vec_rhs, false);

    // we need a type-0 defect vector here
    vec_def.from_1_to_0();

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
      out += stringify_fp_sci(scaling_factor * body_forces[i], 16, 25);

    // print body forces
    print_pad(comm, "Body Forces on " + body_mesh_part, out);
  }


  void StokesSolver::create_multigrid_solver()
  {
    watch_linsol_create.start();

    // create a multigrid hierarchy
    multigrid_hierarchy = std::make_shared<
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
    multigrid_precond = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

    // create our solver
    stokes_solver = Solver::new_richardson(stokes_levels.front()->matrix_sys,
      stokes_levels.front()->filter_sys, 1.0, multigrid_precond);

    stokes_solver->set_plot_name("Multigrid");
    stokes_solver->set_min_iter(min_mg_iter);
    stokes_solver->set_max_iter(max_mg_iter);
    stokes_solver->set_tol_rel(mg_tol_rel);
    stokes_solver->set_min_stag_iter(Index(3));

    watch_linsol_create.stop();

    // initialize the solver symbolically
    watch_linsol_sym.start();
    multigrid_hierarchy->init_symbolic();
    stokes_solver->init_symbolic();
    watch_linsol_sym.stop();
  }

  bool StokesSolver::solve_linear(GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs)
  {
    // set plot mode to enable output
    stokes_solver->set_plot_mode(plot_linear ? Solver::PlotMode::iter : Solver::PlotMode::none);

    // initialize solver
    watch_linsol_num.start();
    multigrid_hierarchy->init_numeric();
    stokes_solver->init_numeric();
    watch_linsol_num.stop();

    // print header
    if(plot_linear_header)
    {
      //         "Multigrid:  2 : 5.012506e-03 / 4.283371e-04 / 0.028025"
      comm.print(String(stokes_solver->get_plot_name().size(), ' ') +
        "  It   Abs. Def.      Rel. Def.      Improve");
      comm.print(String(stokes_solver->get_plot_name().size(), '-') +
        "---------------------------------------------");
    }

    multigrid_hierarchy->reset_timings();

    // apply solver
    Solver::Status status = stokes_solver->correct(vec_sol, vec_rhs);

    // update solver times
    time_mg_smooth += multigrid_hierarchy->get_time_smooth();
    time_mg_coarse += multigrid_hierarchy->get_time_coarse();

    // release solver
    watch_linsol_num.start();
    multigrid_hierarchy->done_numeric();
    stokes_solver->done_numeric();
    watch_linsol_num.stop();

    // check solver output
    if(!Solver::status_success(status))
    {
      comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\nStatus: " + stringify(status));
      return false;
    }

    return true;
  }

  bool StokesSolver::solve_nonlinear(GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs)
  {
    watch_nonlinear_solve.start();

    // clear statistics
    nl_defs.clear();
    mg_iters.clear();
    plot_line.clear();

    // set plot mode to enable output
    stokes_solver->set_plot_mode(plot_mg_iter ? Solver::PlotMode::iter : Solver::PlotMode::none);

    // create two vectors
    GlobalStokesVector vec_cor = create_vector();
    GlobalStokesVector vec_def = create_vector();

    // print header
    if(plot_nonlinear_header)
    {
      //         "Newton...:   1: 2.373388e-05 / 4.732559e-03 / 4.733e-03 |   5: 8.7702e-12 / 3.6952e-07 [1.0000e-10]"
      comm.print("            It  Abs. Defect    Rel. Defect    Improve   |  It  Abs. Def.    Rel. Def.   Rel. Tol.");
      comm.print("--------------------------------------------------------+------------------------------------------");
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
      String line = (newton_solver ? "Newton...:" : "Picard...:");
      if(plot_nonlinear)
      {
        line += stringify(nl_step).pad_front(3) + " : ";
        line += stringify_fp_sci(def_nl, 6) + " / ";
        line += stringify_fp_sci(def_nl/nl_defs.front()) + " / ";
        line += stringify_fp_sci(def_improve, 3);
      }

      // check for divergence
      if(def_nl > nl_defs.front() * DataType(1E+3))
      {
        comm.print(line + "\n\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
        watch_nonlinear_solve.stop();
        return false;
      }

      // check for stopping criterions
      if(nl_step >= min_nl_iter)
      {
        if(def_nl < nl_tol_abs)
        {
          if(plot_nonlinear)
            comm.print(line + "\n\nNonlinear solver converged!\n");
          break;
        }
        else if(nl_step >= max_nl_iter)
        {
          if(plot_nonlinear)
            comm.print(line + "\n\nMaximum iterations reached!\n");
          break;
        }
        else if((nl_step >= 3) && (nl_stag_rate*def_prev < def_nl))
        {
          if(plot_nonlinear)
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
      watch_linsol_num.start();
      multigrid_hierarchy->init_numeric();
      stokes_solver->init_numeric();
      watch_linsol_num.stop();

      multigrid_hierarchy->reset_timings();

      // solve linear system
      watch_linsol_apply.start();
      Solver::Status status = stokes_solver->apply(vec_cor, vec_def);
      watch_linsol_apply.stop();

      // save MG iterations
      mg_iters.push_back(stokes_solver->get_num_iter());

      // update solver times
      time_mg_smooth += multigrid_hierarchy->get_time_smooth();
      time_mg_coarse += multigrid_hierarchy->get_time_coarse();

      // build output line and print it
      if(plot_nonlinear)
      {
        line += String(" | ") + stringify(stokes_solver->get_num_iter()).pad_front(3) + ": "
          + stringify_fp_sci(stokes_solver->get_def_final(), 4) + " / "
          + stringify_fp_sci(stokes_solver->get_def_final() / stokes_solver->get_def_initial(), 4);
        if(adaptive_tol && (nl_step > Index(0)))
          line += String(" [") + stringify_fp_sci(stokes_solver->get_tol_abs(), 4) + "]";
        comm.print(line);
      }

      // release linear solver
      watch_linsol_num.start();
      stokes_solver->done_numeric();
      multigrid_hierarchy->done_numeric();
      watch_linsol_num.stop();

      // check linear solver status
      if(!Solver::status_success(status))
      {
        comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\nStatus: " + stringify(status));
        watch_nonlinear_solve.stop();
        return false;
      }

      // update solution
      vec_sol.axpy(vec_cor);

      // next non-linear iteration
    }

    // build short plot line
    {
      plot_line += stringify(Math::max(int(nl_defs.size()),1)-1).pad_front(4) + ": ";
      plot_line += stringify_fp_sci(nl_defs.front(), 3) + " > ";
      plot_line += stringify_fp_sci(nl_defs.back(), 3) + " : ";
      Index mgi(0u);
      for(auto& i : mg_iters)
        mgi += i;
      plot_line += stringify(mgi).pad_front(4);
    }

    // end of Navier-Stokes solve
    comm.print_flush();

    watch_nonlinear_solve.stop();

    return true;
  }

  void StokesSolver::assemble_nonlinear_defect(GlobalStokesVector& vec_def, const GlobalStokesVector& vec_sol, const GlobalStokesVector& vec_rhs, bool filter_def)
  {
    watch_asm_vector.start();

    // copy rhs vector and convert it to a 0-type vector
    vec_def.copy(vec_rhs);
    vec_def.from_1_to_0();

    // fetch our finest levels
    StokesLevel& the_stokes_level = *stokes_levels.front();

    // assemble burgers operator defect
    this->_assemble_local_burgers_defect(vec_def.local().template at<0>(), vec_sol.local().template at<0>(),
      domain.front().layer(), *domain.front(), the_stokes_level);

    // compute remainder of defect vector by using matrices B and D
    the_stokes_level.matrix_sys.local().block_b().apply(
      vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
    the_stokes_level.matrix_sys.local().block_d().apply(
      vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);

    // sync and filter
    vec_def.sync_0();
    if(filter_def)
      the_stokes_level.filter_sys.filter_def(vec_def);

    watch_asm_vector.stop();
  }

  void StokesSolver::assemble_burgers_matrices(const GlobalStokesVector& vec_sol)
  {
    watch_asm_matrix.start();

    // fetch our finest levels
    StokesLevel& the_stokes_level = *stokes_levels.front();

    // get a clone of the global velocity vector
    // this one will be truncated as we go down the level hierarchy
    GlobalVeloVector vec_conv(&the_stokes_level.gate_velo, vec_sol.local().template at<0>().clone());

    // loop over all system levels
    for(std::size_t i(0); i < stokes_levels.size(); ++i)
    {
      // get the local matrix block A on this level
      LocalMatrixBlockA& loc_mat_a = stokes_levels.at(i)->matrix_sys.local().block_a();
      loc_mat_a.format();

      // assemble local burgers matrix
      this->_assemble_local_burgers_matrix(loc_mat_a, vec_conv.local(), domain.at(i).layer(), *domain.at(i), *stokes_levels.at(i));

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
    watch_asm_matrix.stop();
  }

  void StokesSolver::compute_adaptive_mg_tol()
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

  void StokesSolver::release()
  {
    // release solver
    watch_linsol_sym.start();
    stokes_solver->done_symbolic();
    multigrid_hierarchy->done_symbolic();
    watch_linsol_sym.stop();

    // clear amavankas deque
    amavankas.clear();

    // release solvers
    stokes_solver.reset();
    multigrid_hierarchy.reset();

    // release all levels except for the finest one
    while(stokes_levels.size() > std::size_t(1))
      stokes_levels.pop_back();
  }

  void StokesSolver::print_runtime(double total_time)
  {
    print_time(comm, "Stokes NonLinear Solver Apply Time", watch_nonlinear_solve.elapsed(), total_time);
    print_time(comm, "Stokes Linear Solver Create Time", watch_linsol_create.elapsed(), total_time);
    print_time(comm, "Stokes Linear Solver Init Symbolic Time", watch_linsol_sym.elapsed(), total_time);
    print_time(comm, "Stokes Linear Solver Init Numeric Time", watch_linsol_num.elapsed(), total_time);
    print_time(comm, "Stokes Linear Solver Apply Time", watch_linsol_apply.elapsed(), total_time);
    print_time(comm, "Stokes Multigrid Smoother Time", time_mg_smooth, total_time);
    print_time(comm, "Stokes Multigrid Coarse Solver Time", time_mg_coarse, total_time);
    print_time(comm, "Stokes Burgers Matrix Assembly Time", watch_asm_matrix.elapsed(), total_time);
    print_time(comm, "Stokes Defect Vector Assembly Time", watch_asm_vector.elapsed(), total_time);
    print_time(comm, "Stokes Filter Assembly Time", watch_asm_filter.elapsed(), total_time);
  }

  void StokesSolver::_assemble_local_burgers_defect(LocalVeloVector& vec_def, const LocalVeloVector& vec_velo,
    const DomainLayer& DOXY(domain_lyr), DomainLevel& domain_lvl, StokesLevel& DOXY(stokes_lvl))
  {
    // set up Burgers assembly job for our defect vector
    Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType> burgers_def_job(
      vec_def, vec_velo, vec_velo, domain_lvl.space_velo, cubature);

    // set the parameters for the Burgers assembly
    burgers_def_job.deformation = deform_tensor;
    burgers_def_job.nu = -nu;
    burgers_def_job.beta = DataType(nonlinear_system ? -1 : 0);
    burgers_def_job.theta = -theta;

    // assemble burgers operator defect
    domain_lvl.domain_asm.assemble(burgers_def_job);
  }

  void StokesSolver::_assemble_local_burgers_matrix(LocalMatrixBlockA& matrix_a, const LocalVeloVector& vec_velo,
    const DomainLayer& DOXY(domain_lyr), DomainLevel& domain_lvl, StokesLevel& stokes_lvl)
  {
    // set up Burgers assembly job
    Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>
      burgers_mat_job(matrix_a, vec_velo, domain_lvl.space_velo, cubature);

    // set the parameters for the Burgers assembly
    burgers_mat_job.deformation = deform_tensor;
    burgers_mat_job.nu = nu;
    burgers_mat_job.beta = DataType(nonlinear_system ? 1 : 0);
    burgers_mat_job.frechet_beta = DataType(nonlinear_system && newton_solver ? 1 : 0);
    burgers_mat_job.sd_delta = upsam;
    burgers_mat_job.sd_nu = nu;
    burgers_mat_job.theta = theta;
    if(burgers_mat_job.sd_delta > DataType(0))
    {
      burgers_mat_job.set_sd_v_norm(vec_velo);
      burgers_mat_job.sd_v_norm = stokes_lvl.gate_sys.sum(burgers_mat_job.sd_v_norm);
    }

    // assemble our matrix block
    domain_lvl.domain_asm.assemble(burgers_mat_job);
  }

} // namespace CCNDSimple
