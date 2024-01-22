// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Steady CCnD solver application base class header
// ------------------------------------------------------------------------------------------------
// This header file defines the CCND::SteadyAppBase class, which can be used as a base class for
// CCnD type applications to outsource commonly used functionality. This comment block describes
// the command line parameters that can be used by all applications derived from this class, unless
// the application overrides the parameters explicitly.
//
// The system is discretized using an isoparametric Q2/P1dc finite element discretization.
// The monolithic nonlinear Oseen systems are solved using an adaptive Newton-Multigrid solver
// with an additive matrix-based Vanka smoother ("AmaVanka") and using UMFPACK (if available)
// as a coarse grid solver. The smoother iteration as well as the outer multigrid iteration is
// usually performed by a Richardson iteration scheme, however, it is also possible to switch to
// a FGMRES(k) solver instead. This application supports recursive partitioning.
//
// ------------------------------------
// Basic Setup and Mandatory Parameters
// ------------------------------------
// This application class defines default values for most of its parameters, however, the following
// parameters are mandatory and always have to be specified explicitly:
//
// --mesh <meshfiles...>
// Specifies the input mesh file(s).
//
// --level <level-max> [levels...] <level-min>
// Specifies the mesh refinement levels in the syntax according to Control::PartiDomainControl.
//
//
// ------------------------
// Optional Mesh Parameters
// ------------------------
// This parameters further define specific parts of the used mesh
//
// --cgal-mesh-chart <chart-name> <filename> [..]
// Specify pairs of chart_name and corresponding filename for cgal based charts that should be explicitly
// added to the mesh atlas before the mesh is parsed.
// This is necessary, if a meshpart references a cgal chart that is not (and can not) be defined in the
// meshfile itself.
//
//
//
// ----------------------------
// System Definition Parameters
// ----------------------------
// Some of the basic parameters of the underlying system can be configured by the following
// parameters.
//
// --nu <nu>
// Specifies the viscosity parameter nu for the diffusion operator. Defaults to 1E-3.
//
// --defo
// If specified, the deformation tensor formulation is used for the diffusion operator,
// otherwise the gradient tensor formulation is used. Defaults to gradient tensor.
//
// --upsam <ups>
// Specifies the stabilization parameter <ups> for the streamline diffusion stabilization.
// Defaults to 0, i.e. unstabilized.
//
// --stokes
// If specified, only the steady-state Stokes equations are solved and the nonlinear iteration
// for solving the Navier-Stokes equations is skipped entirely.
//
//
// -------------------------------
// Solver Configuration Parameters
// -------------------------------
// This section describes the parameters that control the non-linear Newton/Picard solver
// as well as its multigrid preconditioner and its smoother component.
//
// --picard
// If specified, the nonlinear system in each time step will be solved using a simple
// Picard iteration instead of the Newton iteration.
//
// --plot-mg-iter
// If specified, the convergence plot of the multigrid solver in each nonlinear solver iteration
// is printed.
//
// --min-nl-iter <N>
// Specifies the minimum number of nonlinear (Newton/Picard) solver iterations per time step.
// Defaults to 1.
//
// --max-nl-iter <N>
// Specifies the maximum number of nonlinear (Newton/Picard) solver iterations per time step.
// Defaults to 20.
//
// --min-mg-iter <N>
// Specifies the minimum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 1.
//
// --max-mg-iter <N>
// Specifies the maximum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 25.
//
// --smooth-steps <N>
// Specifies the number of pre- and post-smoothing AmaVanka steps. Defaults to 8.
//
// --smooth-damp <omega>
// Specifies the damping parameter for the AmaVanka smoother (if inside a Richardson iteration).
// Defaults to 0.5.
//
// --smooth-gmres
// Use FGMRES-AmaVanka as smoother (instead of Richardson-AmaVanka).
//
// --solve-gmres
// Use FGMRES-Multigrid as solver (instead of Richardson-Multigrid).
//
// --nl-tol-abs <tol>
// Specifies the absolute tolerance for the nonlinear solver. Defaults to 1E-8.
//
// --mg-tol-rel <tol>
// If given, specifies the relative tolerance for the multigrid solver.
// If not given, then the tolerance for the multigrid solver is chosen in an adaptive
// manner depending on the two previous nonlinear solver defects, which is the default case.
// The adaptive tolerance is chosen in each nonlinear iteration by analyzing the nonlinear
// defect improvement in the previous nonlinear (Newton/Picard) solver iteration in the
// following manner: Let def_{j} and def_{j-1} denote the two previous nonlinear defect norms,
// then the next nonlinear defect norm def_{j+1} should approximately fulfill the equation
//
//        (def_{j+1} / def_{j}) \approx (def_{j} / def_{j+1})^C
//
// where C \in {1,2} is the convergence speed of the nonlinear solver, i.e. C=2 for Newton
// and C=1 for Picard iteration. Multiplying the above equation by def_{j} gives us an
// estimate for the next nonlinear defect norm def_{j+1}:
//
//        def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})^C
//
// To obtain an absolute tolerance for the next multigrid solver application, we simply
// multiply the above estimate by 0.1.
//
//
// -----------------------------------------------------
// Initial Solution Read-In and Final Solution Write-Out
// -----------------------------------------------------
//
// --load-initial-sol <filename>
// Specifies that the application should read in the initial (partitioned) solution guess
// from a single binary output file, which was written by a --save-final-sol from a previous run.
// If specified, the solving of the Stokes system to obtain an initial guess is skipped.
//
// --load-initial-sol-joined <filename> [<scale>]
// Specifies that the application should read in the initial joined solution guess from a
// single binary output file, which was written by a --save-final-sol-joined from a previous run.
// The second option argument <scale> is given, then the loaded solution is scaled by that factor,
// which can be used if the loaded solution was computed with a different inflow velocity.
// If specified, the solving of the Stokes system to obtain an initial guess is skipped.
//
// --save-final-sol <filename>
// Specifies that the application should write the final (partitioned) solution (and the
// partitioning) to a single binary output file. The generated output file can be loaded by the
// --load-initial-sol option if the input mesh, the refinement level as well as the partitioning
// (and thus the process count) are identical.
//
// --save-final-sol-joined <filename>
// Specifies that the application should write the final solution into a joined binary output
// file by utilizing the base splitter. The generated output file can be loaded by the
// --load-initial-sol-joined option if the input mesh and refinement level are identical, however,
// the process count and/or partitioning may differ. This feature should only be used with at most
// one parallel domain layer and moderate process counts.
//
// ------------------------
// Miscellaneous Parameters
// ------------------------
// This section describes miscellaneous parameters that do not fit into any other section and
// which do not deserve a custom section of their own.
//
// --vtk [<filename>]
// Specifies that the application should write a VTK visualization output file.
//
// --refine-vtk
// Specifies that the VTK files, which are enabled by the --vtk parameter, should be written on a
// once refined mesh.
//
// --ext-stats
// If given, specifies that the application should output extensive statistics at the end of the
// program run, including detailed MPI timings.
//
// \author Peter Zajac
//
#pragma once
#include "applications/ccnd/ccnd_common.hpp"
#include <iterator>
#ifndef APPLICATIONS_CCND_STEADY_APPBASE_HPP
#define APPLICATIONS_CCND_STEADY_APPBASE_HPP 1

namespace CCND
{
  typedef typename SystemLevel::GlobalSystemMatrix GlobalSystemMatrix;
  typedef typename SystemLevel::GlobalSystemVector GlobalSystemVector;
  typedef typename SystemLevel::GlobalSystemFilter GlobalSystemFilter;
  typedef typename SystemLevel::GlobalSystemTransfer GlobalSystemTransfer;
  typedef typename SystemLevel::GlobalVeloVector GlobalVeloVector;
  typedef typename SystemLevel::GlobalPresVector GlobalPresVector;

  typedef typename SystemLevel::LocalSystemMatrix LocalSystemMatrix;
  typedef typename SystemLevel::LocalSystemVector LocalSystemVector;
  typedef typename SystemLevel::LocalSystemFilter LocalSystemFilter;
  typedef typename SystemLevel::LocalVeloVector LocalVeloVector;
  typedef typename SystemLevel::LocalPresVector LocalPresVector;

  class SteadyAppBase
  {
  public:
    /// our communicator
    const Dist::Comm& comm;
    /// our arg parser
    SimpleArgParser& args;

    /// our domain control
    Control::Domain::PartiDomainControl<DomainLevel> domain;
    /// our system levels
    std::deque<std::unique_ptr<SystemLevel>> system;

    /// our statistics
    BenchmarkStats stats;

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Important Note: All the default values below may be overridden by the constructors of derived classes!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    /// enable fictitious boundary method?
    bool enable_fbm = false;

    /// is this actually an unsteady simulation?
    bool is_unsteady = false;

    /// solve Navier-Stokes rather than just Stokes?
    bool navier = true;
    /// use Newton or Picard iteration?
    bool newton = true;
    /// use deformation tensor instead of gradient tensor?
    bool deformation = false;
    /// use adaptive Multigrid tolerance?
    bool adapt_mg_tol = true;
    /// plot multigrid iterations
    bool plot_mg_iter = false;
    /// enable extended statistics collection?
    bool ext_stats = false;
    /// use FGMRES for smoother
    bool smooth_gmres = false;
    /// use FGMRES as outer solver
    bool solve_gmres = false;
    /// use FGMRES as coarse solver
#ifdef FEAT_HAVE_UMFPACK
    bool coarse_gmres = false;
#else // always use GMRES if UMFPACK is not available
    bool coarse_gmres = true;
#endif
    /// need base splitter for joined load/save?
    bool need_base_splitter = false;
    /// need velocity truncation matrices?
    bool need_velocity_truncation = false;

    /// specifies whether the simulation has a non-zero right-hand-side (always false for unsteady simulations)
    bool homogeneous_rhs = false;

    /// viscosity parameter nu
    DataType nu = DataType(1);
    /// streamline diffusion parameter
    DataType upsam = DataType(0);
    /// minimum non-linear solver iterations
    IndexType min_nl_iter = 1u;
    /// maximum non-linear solver iterations
    IndexType max_nl_iter = 20u;
    /// minimum multigrid iterations per non-linear step
    IndexType min_mg_iter = 1u;
    /// maximum multigrid iterations per non-linear step
    IndexType max_mg_iter = 25u;
    /// number of smoothing steps
    IndexType smooth_steps = 8u;
    /// smoother damping parameter
    DataType smooth_damp = DataType(0.5);
    /// relative tolerance for linear solver
    DataType mg_tol_rel = DataType(1E-2);
    /// absolute tolerance for non-linear solver
    DataType nl_tol_abs = DataType(1E-8);

    /// default filename for all kinds of output files (VTK, checkpoints, etc)
    String default_filename;

    /// cubature rule for transfer matrices
    String cubature_transfer = "gauss-legendre:3";
    String cubature_matrix_a = "gauss-legendre:3";
    String cubature_matrix_b = "gauss-legendre:3";
    String cubature_matrix_m = "gauss-legendre:3";
    String cubature_defect   = "gauss-legendre:3";

    /// name of VTK output file(s)
    String vtk_filename;
    /// write VTK files?
    bool want_vtk = false;
    /// write VTK on refined mesh?
    bool refine_vtk = false;
    /// refined mesh for VTK output
    std::unique_ptr<MeshType> refined_mesh_vtk;

    /// the list of all mesh files that have been read
    std::deque<String> mesh_file_names;
    String mesh_file_type;

    //create a new cgal_chart
    bool create_cgal_chart = false;
    std::deque<String> cgal_chart_pairs;

    /// a handful of stopwatches
    StopWatch watch_create_domain, watch_create_system, watch_create_solver;
    StopWatch watch_stokes_solve, watch_nonlin_loop, watch_nonlin_def_asm, watch_nonlin_mat_asm, watch_nonlin_solver_init, watch_nonlin_solver_apply;
    StopWatch watch_sol_init_read, watch_sol_final_write, watch_sol_analysis, watch_total_run, watch_vtk;

    /// our global vectors
    GlobalSystemVector vec_sol, vec_rhs, vec_def, vec_cor, vec_def_unsynced;

    /// our multigrid hierarchy
    std::shared_ptr<Solver::MultiGridHierarchy<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>> multigrid_hierarchy;

    /// our AmaVanka smoothers
    std::deque<std::shared_ptr<Solver::AmaVanka<LocalSystemMatrix, LocalSystemFilter>>> ama_vankas;

    /// our multigrid solver
    std::shared_ptr<Solver::MultiGrid<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>> multigrid;

    /// our outer iterative solver
    std::shared_ptr<Solver::IterativeSolver<GlobalSystemVector>> solver_iterative;

    /// our solver
    std::shared_ptr<Solver::SolverBase<GlobalSystemVector>> solver;

    /// vector of all non-linear defects
    std::vector<DataType> nonlinear_defects;

  public:
    explicit SteadyAppBase(const Dist::Comm& comm_, SimpleArgParser& args_) :
      comm(comm_),
      args(args_),
      domain(comm, true)
    {
      watch_total_run.start();

      Control::Domain::add_supported_pdc_args(args);
      args.support("level");
      args.support("vtk");
      args.support("refine-vtk");
      args.support("mesh");
      args.support("cgal-mesh-chart");
      args.support("nu");
      args.support("upsam");
      args.support("min-nl-iter");
      args.support("max-nl-iter");
      args.support("min-mg-iter");
      args.support("max-mg-iter");
      args.support("plot-mg-iter");
      args.support("smooth-steps");
      args.support("smooth-damp");
      args.support("smooth-gmres");
      args.support("solve-gmres");
      args.support("coarse-gmres");
      args.support("mg-tol-rel");
      args.support("nl-tol-abs");
      args.support("picard");
      args.support("defo");
      args.support("stokes");
      args.support("ext-stats");
      args.support("load-initial-sol");
      args.support("load-initial-sol-joined");
      args.support("save-final-sol");
      args.support("save-final-sol-joined");
    }

    virtual ~SteadyAppBase()
    {
    }

    // no move, no copy, no problems
    SteadyAppBase(SteadyAppBase&&) = delete;
    SteadyAppBase(const SteadyAppBase&) = delete;
    SteadyAppBase& operator=(SteadyAppBase&&) = delete;
    SteadyAppBase& operator=(const SteadyAppBase&) = delete;

    virtual void check_args()
    {
      // check for unsupported options
      auto unsupported = args.query_unsupported();
      if (!unsupported.empty())
      {
        // print all unsupported options to cerr
        for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
          comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

        // abort
        FEAT::Runtime::abort();
      }

      if(args.check("mesh") < 1)
      {
        comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
        FEAT::Runtime::abort();
      }
      if(args.check("level") < 1)
      {
        comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
        FEAT::Runtime::abort();
      }
    }

    virtual void parse_args()
    {
      // parse simple bools
      navier = (args.check("stokes") < 0);
      newton = (args.check("picard") < 0);
      deformation = (args.check("defo") >= 0);
      adapt_mg_tol = (args.check("mg-tol-rel") < 0);
      ext_stats = (args.check("ext-stats") >= 0);
      plot_mg_iter = (args.check("plot-mg-iter") >= 0);
      coarse_gmres = (args.check("coarse-gmres") >= 0);
      smooth_gmres = (args.check("smooth-gmres") >= 0);
      solve_gmres = (args.check("solve-gmres") >= 0);

      // parse other parameters
      args.parse("nu", nu);
      args.parse("upsam", upsam);
      args.parse("min-nl-iter", min_nl_iter);
      args.parse("max-nl-iter", max_nl_iter);
      args.parse("min-mg-iter", min_mg_iter);
      args.parse("max-mg-iter", max_mg_iter);
      args.parse("mg-tol-rel", mg_tol_rel);
      args.parse("nl-tol-abs", nl_tol_abs);
      args.parse("smooth-steps", smooth_steps);
      args.parse("smooth-damp", smooth_damp);

      // want to write VTK?
      refine_vtk = (args.check("refine-vtk") >= 0);
      args.parse("vtk", vtk_filename);
      want_vtk = (args.check("vtk") >= 0) || refine_vtk;

      // need base splitter for joined save/load
      need_base_splitter = need_base_splitter || (args.check("save-final-sol-joined") >= 0) || (args.check("load-initial-sol-joined") >= 0);

      // need velocity truncation for Navier-Stokes
      need_velocity_truncation = need_velocity_truncation || navier;

      // enable extended statistics if desired
      if(ext_stats)
        Statistics::enable_solver_expressions = true;

      // let the domain parse the arguments
      domain.parse_args(args);

      // set desired levels
      domain.set_desired_levels(args.query("level")->second);

      // get mesh files to read in
      mesh_file_names = args.query("mesh")->second;

      create_cgal_chart = (args.check("cgal-mesh-chart") >= 0);
      if(create_cgal_chart)
      {
        cgal_chart_pairs = args.query("cgal-mesh-part")->second;
        if(cgal_chart_pairs.size() == 0 || cgal_chart_pairs.size()%2 != 0)
        {
          XABORTM("Either chart pair names are empty or a part of a pair is missing");
        }
      }

      // keep base levels if we need to save/load joined solution vectors
      if(need_base_splitter)
        domain.keep_base_levels();
    }

    virtual void print_problem()
    {
      comm.print("\nProblem Parameters:");
      comm.print(String("Dimension").pad_back(pad_len, pad_char) + ": " + stringify(dim));
      comm.print(String("Mesh Files").pad_back(pad_len, pad_char) + ": '" + stringify_join(mesh_file_names, "' , '") + "'");
      comm.print(String("Mesh Type").pad_back(pad_len, pad_char) + ": " + mesh_file_type);
      if(create_cgal_chart)
        comm.print(String("CgalMesh Files/Charts").pad_back(pad_len, pad_char) + ": '" + stringify_join(cgal_chart_pairs, "' , '") + "'");
      comm.print(String("Desired Levels").pad_back(pad_len, pad_char) + ": " + domain.format_desired_levels());
      comm.print(String("Chosen Levels").pad_back(pad_len, pad_char) + ": " + domain.format_chosen_levels());
#ifdef FEAT_CCND_APP_ISOPARAM
      comm.print(String("Transformation").pad_back(pad_len, pad_char) + ": Isoparametric:2");
#else
      comm.print(String("Transformation").pad_back(pad_len, pad_char) + ": Standard");
#endif

      comm.print(String("System").pad_back(pad_len, pad_char) + ": steady " + (navier ? "Navier-Stokes" : "Stokes"));
      comm.print(String("Diffusion Tensor").pad_back(pad_len, pad_char) + ": " + (deformation ? "Deformation" : "Gradient"));
      comm.print(String("Fictitious Boundary").pad_back(pad_len, pad_char) + ": " + (enable_fbm ? "yes" : "no"));
      comm.print(String("Nu").pad_back(pad_len, pad_char) + ": " + stringify(nu));
      comm.print(String("Upsam").pad_back(pad_len, pad_char) + ": " + stringify(upsam));
      comm.print(String("Nonlinear Solver").pad_back(pad_len, pad_char) + ": " + (newton ? "Newton" : "Picard"));
      comm.print(String("AmaVanka Smoother Steps").pad_back(pad_len, pad_char) + ": " + stringify(smooth_steps));
      comm.print(String("AmaVanka Smoother Damping").pad_back(pad_len, pad_char) + ": " + stringify(smooth_damp));
      comm.print(String("Smoother Iteration").pad_back(pad_len, pad_char) + ": " + (smooth_gmres ? "FGMRES" : "Richardson"));
      comm.print(String("Nonlinear Absolute Tol").pad_back(pad_len, pad_char) + ": " + stringify_fp_sci(nl_tol_abs));
      comm.print(String("Multigrid Relative Tol").pad_back(pad_len, pad_char) + ": " + (adapt_mg_tol ? String("adaptive") : stringify_fp_sci(mg_tol_rel)));
      comm.print(String("Min Multigrid Iterations").pad_back(pad_len, pad_char) + ": " + stringify(min_mg_iter));
      comm.print(String("Max Multigrid Iterations").pad_back(pad_len, pad_char) + ": " + stringify(max_mg_iter));
      comm.print(String("Max Nonlinear Iterations").pad_back(pad_len, pad_char) + ": " + stringify(max_nl_iter));
      comm.print(String("Multigrid Iteration").pad_back(pad_len, pad_char) + ": " + (solve_gmres ? "FGMRES" : "Richardson"));
      if(coarse_gmres)
        comm.print(String("Multigrid Coarse Solver").pad_back(pad_len, pad_char) + ": FGMRES-AmaVanka");
      else
        comm.print(String("Multigrid Coarse Solver").pad_back(pad_len, pad_char) + ": Direct Solver (UMFPACK)");
    }

    virtual void create_domain()
    {
      watch_create_domain.start();

      // Our mesh file reader
      Geometry::MeshFileReader mesh_reader;

      // start reading mesh files
      mesh_reader.add_mesh_files(comm, mesh_file_names);

      // read the mesh file root markups and write mesh type
      mesh_reader.read_root_markup();
      mesh_file_type = mesh_reader.get_meshtype_string();

      //before we create the mesh in our domain, we need to add our charts if required
      if(create_cgal_chart)
      {
        //requires create_cgal_chart to be divisble by 2
        if(cgal_chart_pairs.size()%2 != 0)
        {
          XABORTM("CGAL chart pair names are not pairs.");
        }
#if defined(FEAT_HAVE_CGAL) && (FEAT_CCND_APP_DIM == 3)
        for(auto it = cgal_chart_pairs.begin(); it != cgal_chart_pairs.end(); std::advance(it, 2))
        {
          //first the chart name, then the filename
          if(!domain.get_atlas().add_mesh_chart(*it, Geometry::Atlas::CGALSurfaceMesh<MeshType>::create_cgal_surface_mesh(*(std::next(it)))))
          {
            XABORTM("Could not add cgal surface mesh " + *it + " " + *(std::next(it)));
          }
        }
#elif defined(FEAT_HAVE_CGAL)
        XABORTM("Trying to create CGALSurfaceMesh with dimension < 3");
#else
        XABORTM("Trying to create CGALSurfaceMesh without CGAL");
#endif
      }

      // try to create the domain
      domain.create(mesh_reader);

      // print partitioning info
      comm.print(domain.get_chosen_parti_info());

      // always use FGMRES as coarse solver if coarse layer has more than 1 processes
      coarse_gmres |= (domain.back_layer().comm().size() > 1);

      // do we need a refined mesh for VTK output?
      if(want_vtk && refine_vtk)
      {
        Geometry::StandardRefinery<MeshType> refinery(domain.front()->get_mesh());
        refined_mesh_vtk = refinery.make_unique();
      }

      watch_create_domain.stop();
    }

    virtual void add_all_isoparam_parts()
    {
#ifdef FEAT_CCND_APP_ISOPARAM
      watch_create_domain.start();
      for(std::size_t i(0); i < domain.size_physical(); ++i)
      {
        DomainLevel& dom_lvl = *domain.at(i);

        // get all mesh-parts in our mesh node
        const MeshNodeType& mesh_node = *dom_lvl.get_mesh_node();
        std::deque<String> part_names = mesh_node.get_mesh_part_names(true);

        // loop over all mesh parts
        for(const auto& p_n : part_names)
        {
          // get the mesh part and its charts
          const MeshPartType* part = mesh_node.find_mesh_part(p_n);
          const ChartBaseType* chart = mesh_node.find_mesh_part_chart(p_n);

          // note: part may be a nullptr if the mesh part does not intersect with the patch of this process, which
          // is not an error case, since the mesh part exists on at least one other process -- just not this one

          // if both exist, then we're adding them to the trafo
          if((part != nullptr) && (chart != nullptr))
            dom_lvl.trafo.add_meshpart_chart(*part, *chart);
        }
      }
      watch_create_domain.stop();
#endif // FEAT_CCND_SIMPLE_ISOPARAM
    }

    virtual void create_system(bool need_velo_mass = false)
    {
      watch_create_system.start();
      const Index num_levels(domain.size_physical());

      // allocate system levels
      system.resize(num_levels);
      for (Index i(0); i < num_levels; ++i)
        system.at(i).reset(new SystemLevel());

      // create system levels
      for (Index i(0); i < num_levels; ++i)
      {
        // get domain and system levels
        DomainLevel& dom_lvl = *domain.at(i);
        SystemLevel& sys_lvl = *system.at(i);

        // compile domain assembler
        dom_lvl.domain_asm.compile_all_elements();

        // assemble gates
        sys_lvl.assemble_gates(domain.at(i));

        // assemble trafo and muxers
        if((i+1) < domain.size_virtual())
        {
          sys_lvl.assemble_coarse_muxers(domain.at(i+1));
          sys_lvl.assemble_transfers(domain.at(i), domain.at(i+1), cubature_transfer);
        }

        // assemble matrix structures
        sys_lvl.assemble_velo_struct(dom_lvl.space_velo);
        sys_lvl.assemble_pres_struct(dom_lvl.space_pres);

        // assemble matrices
        sys_lvl.assemble_velocity_laplace_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, nu, deformation, cubature_matrix_a);
        if(enable_fbm || is_unsteady || need_velo_mass)
          sys_lvl.assemble_velocity_mass_matrix(dom_lvl.domain_asm, dom_lvl.space_velo, cubature_matrix_m);
        sys_lvl.assemble_grad_div_matrices(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b);

        // compile the system matrix
        sys_lvl.compile_system_matrix();
      }

      // need velocity truncation matrices?
      if(need_velocity_truncation)
      {
        for (Index i(0); i < num_levels; ++i)
        {
          if(i+1 < num_levels)
            system.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature_transfer, system.at(i+1).get());
          else if(i+1 < domain.size_virtual())
            system.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature_transfer);
        }
      }

      // assemble base splitter on finest level if required
      if(need_base_splitter)
        system.front()->assemble_base_splitters(domain.front());
      watch_create_system.stop();
    }

    virtual void collect_system_statistics()
    {
      // collect some finest-level statistics
      {
        auto tv = system.front()->gate_velo._freqs.clone(LAFEM::CloneMode::Deep);
        tv.format(1.0);
        Index velo_dofs = Index(system.front()->gate_velo.dot(tv, tv));
        Index locp_dofs = system.front()->gate_pres._freqs.size();
        Index pres_dofs = Index(system.front()->gate_pres.sum(DataType(locp_dofs)));
        stats.counts[Counts::velo_dofs] = velo_dofs/dim;
        stats.counts[Counts::pres_dofs] = pres_dofs;
        stats.counts[Counts::total_dofs] = velo_dofs+pres_dofs;
        stats.counts[Counts::elements] = domain.front()->get_mesh().get_num_elements();
      }

      // accumulate sizes
      for(std::size_t i(0); i < domain.size_physical(); ++i)
      {
        stats.bytes[Bytes::mesh] += domain.at(i)->get_mesh_node()->bytes();
        stats.bytes[Bytes::gate] += system.at(i)->gate_sys.bytes();
        stats.bytes[Bytes::muxer] += system.at(i)->coarse_muxer_sys.bytes();
        stats.bytes[Bytes::matrix] += system.at(i)->matrix_sys.local().bytes();
        stats.bytes[Bytes::transfer] += system.at(i)->transfer_sys.bytes();
        const auto& loc_a = system.at(i)->matrix_sys.local().block_a();
        const auto& loc_b = system.at(i)->matrix_sys.local().block_b();
        const auto& loc_d = system.at(i)->matrix_sys.local().block_d();
        stats.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_a.used_elements() + loc_a.rows() + Index(1));
        stats.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_b.used_elements() + loc_b.rows() + Index(1));
        stats.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_d.used_elements() + loc_d.rows() + Index(1));
        stats.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_a.template used_elements<LAFEM::Perspective::pod>());
        stats.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_b.template used_elements<LAFEM::Perspective::pod>());
        stats.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_d.template used_elements<LAFEM::Perspective::pod>());
      }

      {
        // count non-zeros in a and b
        stats.counts[Counts::nnze_a] = system.front()->matrix_sys.local().block_a().used_elements();
        stats.counts[Counts::nnze_b] = system.front()->matrix_sys.local().block_b().used_elements();
        stats.counts[Counts::nnze_total] = system.front()->matrix_sys.local().template used_elements<LAFEM::Perspective::pod>();
      }
    }

    virtual void compile_system_matrices()
    {
      watch_create_system.start();
      for(std::size_t i(0); i < domain.size_physical(); ++i)
      {
        if(enable_fbm)
        {
          system.at(i)->filter_interface_fbm.filter_weak_matrix_rows(system.at(i)->matrix_a.local(), system.at(i)->velo_mass_matrix.local());
        }
        system.at(i)->compile_system_matrix();
        system.at(i)->compile_local_matrix();
      }
      watch_create_system.stop();
    }

    virtual void create_vectors()
    {
      watch_create_system.start();
      // get our global solve matrix and filter
      GlobalSystemMatrix& matrix = system.front()->matrix_sys;

      // create new vectors
      vec_sol = matrix.create_vector_r();
      vec_rhs = matrix.create_vector_r();
      vec_def = matrix.create_vector_r();
      vec_cor = matrix.create_vector_r();
      vec_def_unsynced = matrix.create_vector_r();

      // format the vectors
      vec_sol.format();
      vec_rhs.format();
      watch_create_system.stop();
    }

    virtual bool load_initial_solution()
    {
      // load distributed solution?
      if(args.check("load-initial-sol") > 0)
      {
        watch_sol_init_read.start();

        String save_name;
        args.parse("load-initial-sol", save_name);
        comm.print("\nReading (partitioned) initial solution from '" + save_name + "'...");

        // read file into streams
        BinaryStream bs_com, bs_sol;
        DistFileIO::read_combined(bs_com, bs_sol, save_name, comm);

        // parse vector from stream
        vec_sol.local().read_from(LAFEM::FileMode::fm_binary, bs_sol);

        // apply our solution filter in case the BCs have changed
        system.front()->filter_sys.filter_sol(vec_sol);

        watch_sol_init_read.stop();
        return true;
      }

      // load joined solution?
      if(args.check("load-initial-sol-joined") > 0)
      {
        watch_sol_init_read.start();

        String save_name;
        DataType load_scale = DataType(1);
        args.parse("load-initial-sol-joined", save_name, load_scale);
        comm.print("\nReading joined initial solution from '" + save_name + "'...");

        // parse vector from file
        system.front()->base_splitter_sys.split_read_from(vec_sol, save_name);

        // scale vector if desired
        if(Math::abs(load_scale - DataType(1)) > DataType(1E-8))
        {
          comm.print("\nScaling joined initial solution by " + stringify(load_scale));
          vec_sol.scale(vec_sol, load_scale);
        }

        // apply our solution filter in case the BCs have changed
        system.front()->filter_sys.filter_sol(vec_sol);

        watch_sol_init_read.stop();
        return true;
      }

      // nope
      return false;
    }

    virtual void save_final_solution()
    {
      if(args.check("save-final-sol") >= 0)
      {
        watch_sol_final_write.start();

        String save_name;
        if(args.parse("save-final-sol", save_name) < 1)
        {
          save_name = default_filename;
          save_name += "-lvl" + stringify(domain.front()->get_level_index());
          save_name += "-n" + stringify(comm.size()) + ".sol";
        }

        comm.print("\nWriting (partitioned) solution to '" + save_name + "'");

        // serialize the partitioning
        std::vector<char> buf_pdc = domain.serialize_partitioning();

        // serialize solution vector
        BinaryStream bs_sol;
        vec_sol.local().write_out(LAFEM::FileMode::fm_binary, bs_sol);

        // write to combined output file
        DistFileIO::write_combined(buf_pdc, bs_sol.container(), save_name, comm);

        watch_sol_final_write.stop();
      }

      /* ***************************************************************************************** */
      /* ***************************************************************************************** */
      /* ***************************************************************************************** */

      if(args.check("save-final-sol-joined") >= 0)
      {
        watch_sol_final_write.start();

        String save_name;
        if(args.parse("save-final-sol-joined", save_name) < 1)
        {
          save_name = default_filename;
          save_name += "-joined-lvl" + stringify(domain.front()->get_level_index()) + ".bin";
          // don't include number of processes, because its invariant
        }

        comm.print("\nWriting joined solution to '" + save_name + "'");

        // write out in binary format
        system.front()->base_splitter_sys.join_write_out(vec_sol, save_name);

        watch_sol_final_write.stop();
      }
    }

    virtual void create_solver()
    {
      watch_create_solver.start();

      const GlobalSystemMatrix& matrix_sys = system.front()->matrix_sys;
      const GlobalSystemFilter& filter_sys = system.front()->filter_sys;

      // do we have just a single level?
      if(domain.size_virtual() == std::size_t(1))
      {
        // create direct solver or FMGRES?
#ifdef FEAT_HAVE_UMFPACK
        if(!coarse_gmres)
        {
          comm.print("\nINFO: only 1 level chosen; creating single grid solver: DirectStokesSolver");
          solver = Solver::new_direct_stokes_solver(matrix_sys, filter_sys);
        }
        else
#endif //  FEAT_HAVE_UMFPACK
        {
          comm.print("\nINFO: only 1 level chosen; creating single grid solver: FGMRES-AmaVanka");
          auto vanka = Solver::new_amavanka(system.front()->local_matrix_sys, filter_sys.local());
          vanka->set_skip_singular(true);
          ama_vankas.push_back(vanka);
          auto schwarz = Solver::new_schwarz_precond(vanka, filter_sys);
          solver = solver_iterative = Solver::new_fgmres(matrix_sys, filter_sys, 16, 0.0, schwarz);

          // configure solver
          solver_iterative->set_plot_name("FGMRES-AmaVanka");
          solver_iterative->set_min_iter(min_mg_iter);
          solver_iterative->set_max_iter(max_mg_iter);
          solver_iterative->set_tol_rel(mg_tol_rel);
          solver_iterative->set_min_stag_iter(Index(3));
        }

        watch_create_solver.stop();
        return;
      }

      // create new hierarchy
      multigrid_hierarchy =
        std::make_shared<Solver::MultiGridHierarchy<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>>(domain.size_virtual());

      // push levels into multigrid
      for(std::size_t i(0); i < system.size(); ++i)
      {
        SystemLevel& lvl = *system.at(i);

        if((i+1) < domain.size_virtual())
        {
          // create Schwarz-AmaVanka
          auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
          vanka->set_skip_singular(true);
          ama_vankas.push_back(vanka);
          auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);

          // create smoother: either GMRES or Richardson
          std::shared_ptr<Solver::IterativeSolver<GlobalSystemVector>> smoother;
          if(smooth_gmres)
            smoother = Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, smooth_steps, 0.0, schwarz);
          else
            smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);
          smoother->set_min_iter(smooth_steps);
          smoother->set_max_iter(smooth_steps);
          smoother->skip_defect_calc(true); // skip defect calculation if possible
          multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
        }
#ifdef FEAT_HAVE_UMFPACK
        else if(!coarse_gmres)
        {
          // create UMFPACK coarse grid solver
          auto coarse_solver = Solver::new_direct_stokes_solver(lvl.matrix_sys, lvl.filter_sys);
          multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, coarse_solver);
        }
#endif //  FEAT_HAVE_UMFPACK
        else
        {
          // create FGMRES-AmaVanka coarse grid solver
          auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
          vanka->set_skip_singular(true);
          ama_vankas.push_back(vanka);
          auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
          //auto coarse_solver = Solver::new_bicgstab(lvl.matrix_sys, lvl.filter_sys, schwarz);
          auto coarse_solver = Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, 16, 0.0, schwarz);
          coarse_solver->set_max_iter(500);
          coarse_solver->set_tol_rel(1e-3);
          //coarse_solver->set_plot_mode(Solver::PlotMode::summary);

          multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, coarse_solver);
        }
      }

      // create our multigrid solver
      multigrid = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

      // create our solver
      if(solve_gmres)
        solver = solver_iterative = Solver::new_fgmres(matrix_sys, filter_sys, 16, 0.0, multigrid);
      else
        solver = solver_iterative = Solver::new_richardson(matrix_sys, filter_sys, 1.0, multigrid);

      // configure iterative solver
      if(solver_iterative)
      {
        // set solver name
        if(solve_gmres)
          solver_iterative->set_plot_name("FGMRES-MG");
        else
          solver_iterative->set_plot_name("Multigrid");

        // configure solver
        solver_iterative->set_min_iter(min_mg_iter);
        solver_iterative->set_max_iter(max_mg_iter);
        solver_iterative->set_tol_rel(mg_tol_rel);
        solver_iterative->set_min_stag_iter(Index(3));
      }
      watch_create_solver.stop();
    }

    virtual void init_solver_symbolic()
    {
      watch_nonlin_solver_init.start();
      if(multigrid_hierarchy)
        multigrid_hierarchy->init_symbolic();
      solver->init_symbolic();
      watch_nonlin_solver_init.stop();
    }

    virtual void init_solver_numeric()
    {
      watch_nonlin_solver_init.start();
      if(multigrid_hierarchy)
        multigrid_hierarchy->init_numeric();
      solver->init_numeric();
      watch_nonlin_solver_init.stop();
    }

    virtual void done_solver_numeric()
    {
      solver->done_numeric();
      if(multigrid_hierarchy)
        multigrid_hierarchy->done_numeric();

      FEAT::Statistics::compress_solver_expressions();
    }

    virtual void done_solver_symbolic()
    {
      solver->done_symbolic();
      if(multigrid_hierarchy)
        multigrid_hierarchy->done_symbolic();
    }

    virtual bool solve_stokes()
    {
      comm.print("\nSolving Stokes system...");

      init_solver_numeric();

      Solver::Status stokes_status = Solver::Status::undefined;

      // solve Stokes system
      if(solver_iterative)
      {
        solver_iterative->set_plot_mode(Solver::PlotMode::iter);

        watch_stokes_solve.start();
        stokes_status = solver_iterative->correct(vec_sol, vec_rhs);
        watch_stokes_solve.stop();

        stats.counts[Counts::linsol_iter] = solver_iterative->get_num_iter();
      }
      else
      {
        watch_stokes_solve.start();
        // solve
        stokes_status = Solver::solve(*solver, vec_sol, vec_rhs, system.front()->matrix_sys, system.front()->filter_sys);

        // compute defect after solve manually
        system.front()->matrix_sys.apply(vec_def, vec_sol, vec_rhs, -DataType(1));
        system.front()->filter_sys.filter_def(vec_def);
        //if(enable_fbm)
          //system.front()->filter_interface_fbm.filter_def(vec_def.local().first());

        comm.print("Defect after direct solve: " + stringify_fp_sci(vec_def.norm2()));

        watch_stokes_solve.stop();
      }

      // release solvers
      done_solver_numeric();

      if(!Solver::status_success(stokes_status))
      {
        comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
        //if(testmode)
          //comm.print("Test-Mode: FAILED");
        return false;
      }

      // success
      return true;
    }

    virtual void compute_unsynced_defect()
    {
      // In this case, let's compute the final unsynchronized defect for volumetric body forces computation
      system.front()->matrix_sys.apply(vec_def_unsynced, vec_sol, vec_rhs, -DataType(1));
      vec_def_unsynced.from_1_to_0();
    }

    virtual void setup_burgers_defect_job(Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType>& burgers_def_job)
    {
      burgers_def_job.deformation = deformation;
      burgers_def_job.nu = -nu;
      burgers_def_job.beta = -DataType(1);
    }

    virtual void initialize_nonlinear_defect()
    {
      if(!homogeneous_rhs)
      {
        vec_def.copy(vec_rhs);
        vec_def.from_1_to_0();
      }
      else
        vec_def.format();
    }

    virtual void assemble_nonlinear_defect()
    {
      // assemble nonlinear defect vector
      watch_nonlin_def_asm.start();

      initialize_nonlinear_defect();

      // if FBM is enabled, then we temporarily assemble the Burgers defect into vec_def_unsynced
      if(enable_fbm)
        vec_def_unsynced.format();

      // set up Burgers assembly job for our defect vector
      Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType> burgers_def_job(
        enable_fbm ? vec_def_unsynced.local().template at<0>() : vec_def.local().template at<0>(),
        vec_sol.local().template at<0>(),
        vec_sol.local().template at<0>(), domain.front()->space_velo, cubature_defect);
      setup_burgers_defect_job(burgers_def_job);

      // assemble burgers operator defect
      domain.front()->domain_asm.assemble(burgers_def_job);

      // apply FBM filter if required
      if(enable_fbm)
      {
        system.front()->apply_fbm_filter_to_def(vec_def_unsynced, vec_sol, -DataType(1));
        vec_def.axpy(vec_def, vec_def_unsynced);
      }

      // compute remainder of defect vector
      system.front()->matrix_sys.local().block_b().apply(
        vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
      system.front()->matrix_sys.local().block_d().apply(
        vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);
      // store the unsynced and unfiltered defect for later body forces computation
      vec_def_unsynced.copy(vec_def);
      // synchronize and filter the defect
      vec_def.sync_0();
      //if(enable_fbm)
        //system.front()->filter_interface_fbm.filter_def(vec_def.local().first());
      system.front()->filter_sys.filter_def(vec_def);
      watch_nonlin_def_asm.stop();
    }

    virtual void setup_burgers_matrix_job(Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>& burgers_mat_job)
    {
      burgers_mat_job.deformation = deformation;
      burgers_mat_job.nu = nu;
      burgers_mat_job.beta = DataType(1);
      burgers_mat_job.frechet_beta = DataType(newton ? 1 : 0);
      burgers_mat_job.sd_delta = upsam;
      burgers_mat_job.sd_nu = nu;
    }

    virtual void assemble_nonlinear_matrices()
    {
      // assemble burgers matrices on all levels
      watch_nonlin_mat_asm.start();
      {
        // get a clone of the global velocity vector
        GlobalVeloVector vec_conv(&system.front()->gate_velo, vec_sol.local().template at<0>().clone());

        // initialize velocity norm for streamline diffusion (if enabled)
        DataType sd_v_norm = DataType(0);

        // loop over all system levels
        for(std::size_t i(0); i < system.size(); ++i)
        {
          auto& loc_mat_a = system.at(i)->matrix_sys.local().block_a();
          loc_mat_a.format();

          // set up Burgers assembly job
          Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>
            burgers_mat_job(loc_mat_a, vec_conv.local(), domain.at(i)->space_velo, cubature_matrix_a);
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
          setup_burgers_matrix_job(burgers_mat_job);

          // assemble our system matrix
          domain.at(i)->domain_asm.assemble(burgers_mat_job);

          // restrict our convection vector
          if((i+1) >= domain.size_virtual())
            break;

          // does this process have another system level?
          if((i+1) < system.size())
          {
            // create a coarse mesh velocity vector
            auto vec_crs = system.at(i+1)->matrix_a.create_vector_l();

            // truncate fine mesh velocity vector
            system.at(i)->transfer_velo.trunc(vec_conv, vec_crs);

            // the coarse vector is our next convection vector
            vec_conv = std::move(vec_crs);
          }
          else
          {
            // this process is a child, so send truncation to parent
            system.at(i)->transfer_velo.trunc_send(vec_conv);
          }
        }
      }
      watch_nonlin_mat_asm.stop();
    }

    virtual void choose_nonlinear_tolerances(bool first_iteration, DataType def_nl, DataType def_improve)
    {
      // specify adaptive tolerance?
      if(adapt_mg_tol && !first_iteration)
      {
        if(newton)
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
          solver_iterative->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
          // Also make sure that we gain at least 2 digits.
          solver_iterative->set_tol_rel(1E-2);
        }
        else
        {
          // In the case if Picard iteration, we only expect linear convergence,
          // which (in analogy to Newton) leads us to the following estimate:
          //
          //      def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})
          //
          DataType abs_tol = def_nl * def_improve * DataType(0.1);
          solver_iterative->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
          solver_iterative->set_tol_rel(1E-2);
        }
      }
      else
      {
        // disable absolute tolerance
        solver_iterative->set_tol_abs(1E+10);
        solver_iterative->set_tol_rel(mg_tol_rel);
      }
    }

    virtual bool solve_navier_stokes(String line_prefix = "")
    {
      bool single_line_plot = !line_prefix.empty();

      if(!single_line_plot)
        comm.print("\nSolving Navier-Stokes system...");

      if(!plot_mg_iter && (solver_iterative))
        solver_iterative->set_plot_mode(Solver::PlotMode::none);

      // vector of all non-linear defect norms

      watch_nonlin_loop.start();

      // nonlinear loop
      for(Index nl_step(0); nl_step <= max_nl_iter; ++nl_step)
      {
        stats.counts[Counts::nonlin_iter] = nl_step;

        // compute defect
        assemble_nonlinear_defect();

        // compute defect norm
        const DataType def_prev = (nonlinear_defects.empty() ? DataType(1) : nonlinear_defects.back());
        const DataType def_nl = vec_def.norm2();
        const DataType def_improve = def_nl / def_prev;
        nonlinear_defects.push_back(def_nl);

        String line = line_prefix;
        if(!single_line_plot)
          line = (newton ? "Newton: " : "Picard: ");
        line += stringify(nl_step).pad_front(2) + ": ";
        line += stringify_fp_sci(def_nl, 6) + " / ";
        line += stringify_fp_sci(def_nl/nonlinear_defects.front()) + " / ";
        line += stringify_fp_sci(def_nl/def_prev, 3);

        if(def_nl > nonlinear_defects.front() * DataType(1E+3))
        {
          comm.print(line + "\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
          //if(testmode)
            //comm.print("Test-Mode: FAILED");
          return false;
        }
        else if(nl_step < min_nl_iter)
        {
          // nothing to do here; this else-case exists just to ensure
          // that none of the following cases fires and breaks the loop
        }
        else if(def_nl < nl_tol_abs)
        {
          if(single_line_plot)
            comm.print(line + " > converged!");
          else
            comm.print(line + "\nNonlinear solver converged!\n");
          break;
        }
        else if(nl_step >= max_nl_iter)
        {
          if(single_line_plot)
            comm.print(line + " > maximum iterations reached!");
          else
            comm.print(line + "\nMaximum iterations reached!\n");
          break;
        }
        else if((nl_step >= 3) && (DataType(0.97)*def_prev < def_nl))
        {
          if(single_line_plot)
            comm.print(line + " > stagnated!");
          else
            comm.print(line + "\nNonlinear solver stagnated!\n");
          break;
        }

        // assemble nonlinear matrices
        assemble_nonlinear_matrices();

        // compile local matrices
        compile_system_matrices();

        // initialize linear solver
        init_solver_numeric();

        // choose nonlinear tolerances if using an iterative solver
        if(solver_iterative)
          choose_nonlinear_tolerances(nl_step == Index(0), def_nl, def_improve);

        // solve linear system
        watch_nonlin_solver_apply.start();
        Solver::Status status = solver->apply(vec_cor, vec_def);
        watch_nonlin_solver_apply.stop();

        // get solver statistics
        if(solver_iterative)
        {
          stats.counts[Counts::linsol_iter] += solver_iterative->get_num_iter();
          line += String(" | ") + stringify(solver_iterative->get_num_iter()).pad_front(3) + ": "
            + stringify_fp_sci(solver_iterative->get_def_final(), 4) + " / "
            + stringify_fp_sci(solver_iterative->get_def_final() / solver_iterative->get_def_initial(), 4);
          if(adapt_mg_tol && (nl_step > Index(0)))
            line += String(" [") + stringify_fp_sci(solver_iterative->get_tol_abs(), 4) + "]";
        }
        comm.print(line);

        // release linear solver
        done_solver_numeric();

        if(!Solver::status_success(status))
        {
          comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
          //if(testmode)
            //comm.print("Test-Mode: FAILED");
          return false;
        }

        // update solution
        vec_sol.axpy(vec_cor, vec_sol, DataType(1));

        // next non-linear iteration
      }

      watch_nonlin_loop.stop();

      return true;
    }

    virtual void collect_solver_statistics()
    {
      // save timings
      stats.times[Times::stokes_solve] = watch_stokes_solve.elapsed();
      stats.times[Times::nonlin_total] = watch_nonlin_loop.elapsed();
      stats.times[Times::nonlin_asm_def] = watch_nonlin_def_asm.elapsed();
      stats.times[Times::nonlin_asm_mat] = watch_nonlin_mat_asm.elapsed();
      stats.times[Times::linsol_init] = watch_nonlin_solver_init.elapsed();
      stats.times[Times::linsol_apply] = watch_nonlin_solver_apply.elapsed();

      // get multigrid timings
      if(multigrid_hierarchy)
      {
        stats.times[Times::mg_defect] = multigrid_hierarchy->get_time_defect();
        stats.times[Times::mg_smooth] = multigrid_hierarchy->get_time_smooth();
        stats.times[Times::mg_coarse] = multigrid_hierarchy->get_time_coarse();
        stats.times[Times::mg_transfer] = multigrid_hierarchy->get_time_transfer();
      }

      // accumulate vanka sizes
      if(!ama_vankas.empty())
        stats.counts[Counts::vanka_data] = Index(ama_vankas.front()->data_size());
      stats.bytes[Bytes::vanka] = 0ull;

      // accumulate vanka timings
      for(auto& v : ama_vankas)
      {
        stats.times[Times::vanka_init_sym] += v->time_init_symbolic();
        stats.times[Times::vanka_init_num] += v->time_init_numeric();
        stats.times[Times::vanka_apply] += v->time_apply();
        stats.bytes[Bytes::vanka] += v->bytes();
      }
    }

    virtual void write_vtk()
    {
      if(!want_vtk)
        return;

      watch_vtk.start();

      // build VTK name
      String vtk_name = vtk_filename;
      if(vtk_filename.empty())
      {
        vtk_name = default_filename;
        vtk_name += "-lvl" + stringify(domain.front()->get_level_index());
        vtk_name += "-n" + stringify(comm.size());
      }

      comm.print("Writing VTK output to '" + vtk_name + ".pvtu'\n");

      // get the mesh for the VTK output
      const MeshType* mesh = &domain.front()->get_mesh();
      if(refined_mesh_vtk)
        mesh = refined_mesh_vtk.get();

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(*mesh);

      // write velocity
      exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());
      exporter.add_vertex_vector("rhs_v", vec_rhs.local().template at<0>());

      // project pressure
      LAFEM::DenseVector<DataType, Index> vtx_p, vtx_q;
      if(refined_mesh_vtk)
      {
        Assembly::DiscreteCellProjector::project_refined(vtx_p, vec_sol.local().template at<1>(), domain.front()->space_pres);
        Assembly::DiscreteCellProjector::project_refined(vtx_q, vec_rhs.local().template at<1>(), domain.front()->space_pres);
      }
      else
      {
        Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), domain.front()->space_pres, "gauss-legendre:2");
        Assembly::DiscreteCellProjector::project(vtx_q, vec_rhs.local().template at<1>(), domain.front()->space_pres, "gauss-legendre:2");
      }

      // write pressure
      exporter.add_cell_scalar("pressure", vtx_p.elements());
      exporter.add_cell_scalar("rhs_p", vtx_q.elements());

      // write FBM masks if
      if(enable_fbm)
      {
        if(refined_mesh_vtk)
        {
          exporter.add_vertex_scalar("fbm_mask_v", system.front()->fbm_mask_velo.data());

          const std::vector<int>& mask_p = domain.front()->fbm_assembler->get_fbm_mask_vector(dim);
          const int nc = 1 << dim;
          std::vector<int> mask_p_ref;
          mask_p_ref.reserve(mask_p.size() * std::size_t(nc));
          for(int i : mask_p)
            for(int k(0); k < nc; ++k)
              mask_p_ref.push_back(i);
          exporter.add_cell_scalar("fbm_mask_p", mask_p_ref.data());
        }
        else
        {
          exporter.add_vertex_scalar("fbm_mask_v", domain.front()->fbm_assembler->get_fbm_mask_vector(0).data());
          exporter.add_cell_scalar("fbm_mask_p", domain.front()->fbm_assembler->get_fbm_mask_vector(dim).data());
        }
      }

      // finally, write the VTK file
      exporter.write(vtk_name, comm);

      watch_vtk.stop();
    }

    virtual void collect_final_statistics()
    {
      stats.times[Times::create_domain] = watch_create_domain.elapsed();
      stats.times[Times::create_system] = watch_create_system.elapsed();
      stats.times[Times::create_solver] = watch_create_solver.elapsed();
      stats.times[Times::sol_analysis] = watch_sol_analysis.elapsed();
      stats.times[Times::vtk_write] = watch_vtk.elapsed();
      stats.times[Times::sol_init_read] = watch_sol_init_read.elapsed();
      stats.times[Times::sol_final_write] = watch_sol_final_write.elapsed();

      watch_total_run.stop();
      stats.times[Times::total_run] = watch_total_run.elapsed();

      {
        MemoryUsage mi;
        stats.bytes[Bytes::peak_p] = mi.get_peak_physical();
        stats.bytes[Bytes::peak_v] = mi.get_peak_virtual();
      }

      stats.sync(comm);
    }

    virtual void print_statistics()
    {
      // print basic statistics
      comm.print(stats.format());

      // print multigrid timings
      if((comm.rank() == 0) && (multigrid_hierarchy))
      {
        comm.print("Multigrid Timings:");
        comm.print("              Defect /   Smoother /   Transfer /     Coarse");
        comm.print("Overall : " +
          stringify_fp_fix(multigrid_hierarchy->get_time_defect(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_smooth(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_transfer(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_coarse(), 3, 10));
        for(int i(0); i < int(multigrid_hierarchy->size_physical()); ++i)
        {
          comm.print("Level " + stringify(i).pad_front(2) + ": " +
            stringify_fp_fix(multigrid_hierarchy->get_time_defect(i), 3, 10) + " / " +
            stringify_fp_fix(multigrid_hierarchy->get_time_smooth(i), 3, 10) + " / " +
            stringify_fp_fix(multigrid_hierarchy->get_time_transfer(i), 3, 10) + " / " +
            stringify_fp_fix(multigrid_hierarchy->get_time_coarse(i), 3, 10));
        }
      }

      // print extended statistics if desired
      if(ext_stats)
      {
        comm.print("\n");
        comm.print(FEAT::Statistics::get_formatted_flops(stats.times[Times::stokes_solve] + stats.times[Times::linsol_apply], comm.size()));
        comm.print(FEAT::Statistics::get_formatted_times(stats.times[Times::stokes_solve] + stats.times[Times::linsol_apply]));
        comm.print(FEAT::Statistics::get_formatted_solver_internals("default"));
        comm.print("\n");
        comm.print(FEAT::Statistics::get_formatted_solver_tree("default").trim());
      }

      comm.print("\nTotal Runtime: " + watch_total_run.elapsed_string(TimeFormat::h_m_s_m)
        + " [" + watch_total_run.elapsed_string(TimeFormat::s_m) + " seconds]");

    }
  }; // class SteadyAppBase
} // namespace CCND

#endif // APPLICATIONS_CCND_STEADY_APPBASE_HPP
