// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Steady CCnD solver for DFG95 Flow-Around-A-Cylinder Benchmarks
// ------------------------------------------------------------------------------------------------
// This application implements a parallel steady CCND solver, which is pre-configured to solve
// the infamous unsteady "flow-around-a-cylinder" benchmark problem, which is defined in
//
//     M. Schaefer and S. Turek: Benchmark Computations of Laminar Flow Around a Cylinder
//
// The system is discretized using an isoparametric Q2/P1dc finite element discretization.
// The monolithic nonlinear Oseen systems are solved using an adaptive Newton-Multigrid solver
// with an additive matrix-based Vanka smoother ("AmaVanka") and using UMFPACK (if available)
// as a coarse grid solver. This application supports recursive partitioning.
//
//
// ------------------------------------
// Basic Setup and Mandatory Parameters
// ------------------------------------
// This application defines default values for most of its parameters, however, three parameters
// are mandatory and always have to be specified explicitly:
//
// --mesh <meshfiles...>
// Specifies the input mesh file(s).
//
// --level <level-max> [levels...] <level-min>
// Specifies the mesh refinement levels in the syntax according to Control::PartiDomainControl.
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
// --v-max <vmax>
// Specifies the maximum velocity of the inflow boundary condition.
// Defaults to 1.5 in 2D and 2.25 in 3D.
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
// Defaults to 10.
//
// --min-mg-iter <N>
// Specifies the minimum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 1.
//
// --max-mg-iter <N>
// Specifies the maximum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 5.
//
// --smooth-steps <N>
// Specifies the number of pre- and post-smoothing AmaVanka steps. Defaults to 8.
//
// --smooth-damp <omega>
// Specifies the damping parameter for the AmaVanka smoother. Defaults to 0.7.
//
// --no-umfpack
// If specified, the multigrid solver will use a BiCGStab-AmaVanka solver as the coarse grid solver
// instead of the UMFPACK direct solver. Note that UMFPACK is only used if it is included in the
// build id and if the coarse system is solved on a single process.
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
// --save-sol <filename>
// Specifies that the application should write the final (partitioned) solution (and the
// partitioning) to a single binary output file. The output file can be loaded by the --load-sol
// option if the input mesh, the refinement level as well as the partitioning (and thus the
// process count) are identical.
//
// --save-joined-sol <filename>
// Specifies that the application should write the final solution into a joined binary output
// file by utilizing the base splitter. The output file can be loaded by the --load-joined-sol
// option if the input mesh and refinement level are identical, however, the process count
// and/or partitioning may differ. This feature should only be used with at most one parallel
// domain layer and moderate process counts.
//
// --load-sol <filename>
// Specifies that the application should read in the initial (partitioned) solution guess
// from a single binary output file, which was written by a --save-sol from a previous run.
// If specified, the solving of the Stokes system to obtain an initial guess is skipped.
//
// --load-joined-sol <filename>
// Specifies that the application should read in the initial joined solution guess from a
// single binary output file, which was written by a --save-joined-sol from a previous run.
// If specified, the solving of the Stokes system to obtain an initial guess is skipped.
//
//
// ------------------------
// Miscellaneous Parameters
// ------------------------
// This section describes miscellaneous parameters that do not fit into any other section and
// which do not deserve a custom section of their own.
//
// --vtk <filename> [<refined-filename>]
// Specifies that the application should write a VTK visualization output file. The second
// optional parameter specifies the filename for a VTK output file on a once refined mesh,
// if given. Note: it is possible to output only the refined VTKs by passing a whitespace
// string as the first filename, e.g.: --vtk " " myvtk
//
// --ext-stats
// If given, specifies that the application should output extensive statistics at the end of the
// program run, including detailed MPI timings.
//
// --test-mode
// If given, specifies that the application should run in test-mode rather than its normal mode.
// In test mode, the application only perform 3 time steps and writes some additional output which
// is parsed by the test system.
//
//
// \author Peter Zajac
//

// include common CCND header
#include "ccnd_common.hpp"

namespace DFG95
{
  template<typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::PartiDomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainLevelType::SpaceVeloType SpaceVeloType;
    typedef typename DomainLevelType::SpacePresType SpacePresType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // define our system level
    typedef NavierStokesBlockedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    BenchmarkSummary<DataType, dim> summary;
    BenchmarkStats statistics;
    StopWatch watch_total_run;
    watch_total_run.start();

    /* ****************************************************************************************** */

    const bool navier = (args.check("stokes") < 0);
    const bool newton = (args.check("picard") < 0);
    const bool defo = (args.check("defo") >= 0);
    const bool adapt_tol = (args.check("mg-tol-rel") < 0);
    const bool testmode = (args.check("test-mode") >= 0);
    const bool ext_stats = (args.check("ext-stats") >= 0);
    const bool plot_mg_iter = (args.check("plot-mg-iter") >= 0);

    // viscosity parameter
    const DataType nu = parse(args, "nu", DataType(1e-3));
    // maximum velocity, default: 2D: 0.3, 3D: 0.45
    const DataType v_max = parse(args, "v-max", DataType(dim) * DataType(0.15));
    // streamline diffusion parameter
    const DataType upsam = parse(args, "upsam", DataType(0));
    // max. nonlinear solver iterations
    const Index max_nl_iter = parse(args, "max-nl-iter", Index(20));
    // min. multigrid iterations
    const Index min_mg_iter = parse(args, "min-mg-iter", Index(1));
    // max. multigrid iterations
    const Index max_mg_iter = parse(args, "max-mg-iter", Index(25));
    // number of smoothing steps
    const Index smooth_steps = parse(args, "smooth-steps", Index(8));
    // damping parameter for smoother
    const DataType smooth_damp = parse(args, "smooth-damp", DataType(0.5));
    // rel. tolerance for linear solver
    const DataType mg_tol_rel = parse(args, "mg-tol-rel", DataType(1E-2));
    // absolute tolerance for nonlinear solver
    const DataType nl_tol_abs = parse(args, "nl-tol-abs", DataType(1E-8));


#ifdef FEAT_HAVE_UMFPACK
    const bool umf_cgs = (domain.back_layer().comm().size() == 1) && (args.check("no-umfpack") < 0);
#else
    const bool umf_cgs = false;
#endif

    {
      static constexpr std::size_t pl = 30u;
      static constexpr char pc = '.';
      comm.print("\nProblem Parameters:");
      comm.print(String("Nu").pad_back(pl, pc) + ": " + stringify(nu));
      comm.print(String("V-Max").pad_back(pl, pc) + ": " + stringify(v_max));
      comm.print(String("Upsam").pad_back(pl, pc) + ": " + stringify(upsam));
      comm.print(String("System").pad_back(pl, pc) + ": " + (navier ? "Navier-Stokes" : "Stokes"));
      comm.print(String("Tensor").pad_back(pl, pc) + ": " + (defo ? "Deformation" : "Gradient"));
      comm.print(String("Nonlinear Solver").pad_back(pl, pc) + ": " + (newton ? "Newton" : "Picard"));
      comm.print(String("AmaVanka Smoother Steps").pad_back(pl, pc) + ": " + stringify(smooth_steps));
      comm.print(String("AmaVanka Smoother Damping").pad_back(pl, pc) + ": " + stringify(smooth_damp));
      comm.print(String("Nonlinear Absolute Tol").pad_back(pl, pc) + ": " + stringify_fp_sci(nl_tol_abs));
      comm.print(String("Multigrid Relative Tol").pad_back(pl, pc) + ": " + (adapt_tol ? String("adaptive") : stringify_fp_sci(mg_tol_rel)));
      comm.print(String("Min Multigrid Iterations").pad_back(pl, pc) + ": " + stringify(min_mg_iter));
      comm.print(String("Max Multigrid Iterations").pad_back(pl, pc) + ": " + stringify(max_mg_iter));
      comm.print(String("Max Nonlinear Iterations").pad_back(pl, pc) + ": " + stringify(max_nl_iter));
      if(umf_cgs)
        comm.print(String("Coarse Solver").pad_back(pl, pc) + ": UMFPACK");
      else
        comm.print(String("Coarse Solver").pad_back(pl, pc) + ": BiCGStab-AmaVanka");
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = Index(domain.size_physical());

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    // cubature for assembly
    const String cubature("gauss-legendre:3");

    // cubature for post-processing
    Cubature::DynamicFactory cubature_postproc("gauss-legendre:5");

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    // assemble gates, muxers and transfers
    for (Index i(0); i < num_levels; ++i)
    {
      domain.at(i)->domain_asm.compile_all_elements();
      system_levels.at(i)->assemble_gates(domain.at(i));

      if((i+1) < domain.size_virtual())
      {
        system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
        system_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature);
      }
    }

    // assemble base splitter on finest level if required
    if((args.check("save-joined-sol") >= 0) || (args.check("load-joined-sol") >= 0))
      system_levels.front()->assemble_base_splitters(domain.front());

    if(navier)
    {
      // assemble velocity truncation operators -- we need those for the assembly of the
      // non-linear burgers operators on the coarser levels
      for (Index i(0); i < num_levels; ++i)
      {
        if(i+1 < num_levels)
          system_levels.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature, system_levels.at(i+1).get());
        else if(i+1 < domain.size_virtual())
          system_levels.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature);
      }
    }

    // collect some finest-level statistics
    {
      auto tv = system_levels.front()->gate_velo._freqs.clone(LAFEM::CloneMode::Deep);
      tv.format(1.0);
      Index velo_dofs = Index(system_levels.front()->gate_velo.dot(tv, tv));
      Index locp_dofs = system_levels.front()->gate_pres._freqs.size();
      Index pres_dofs = Index(system_levels.front()->gate_pres.sum(DataType(locp_dofs)));
      statistics.counts[Counts::velo_dofs] = velo_dofs/dim;
      statistics.counts[Counts::pres_dofs] = pres_dofs;
      statistics.counts[Counts::total_dofs] = velo_dofs+pres_dofs;
      statistics.counts[Counts::elements] = domain.front()->get_mesh().get_num_elements();
    }

    /* ***************************************************************************************** */

    // assemble basic matrices
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
      system_levels.at(i)->assemble_pres_struct(domain.at(i)->space_pres);
      system_levels.at(i)->assemble_velocity_laplace_matrix(domain.at(i)->domain_asm, domain.at(i)->space_velo, cubature, nu, defo);
      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->domain_asm, domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
      system_levels.at(i)->compile_system_matrix();
    }

    /* ***************************************************************************************** */

    // our inflow BC function
    SteadyInflowFunction<dim> inflow_func(v_max);

    // the names of the mesh parts on which to assemble
    std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local velocity filter
      auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow, unit_asm_noflow;

      // loop over all boundary parts except for the right one, which is outflow
      for(const auto& name : part_names)
      {
        // skip non-boundary mesh-parts
        if(!name.starts_with("bnd:"))
          continue;

        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
        XASSERT(mesh_part_node != nullptr);

        auto* mesh_part = mesh_part_node->get_mesh();
        if (mesh_part != nullptr)
        {
          if(name == "bnd:l")
          {
            // inflow
            unit_asm_inflow.add_mesh_part(*mesh_part);
          }
          else if(name != "bnd:r")
          {
            // outflow
            unit_asm_noflow.add_mesh_part(*mesh_part);
          }
        }
      }

      // assemble the filters
      unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, inflow_func);
      unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo);

      // finally, compile the system filter
      system_levels.at(i)->compile_system_filter();
    }

    // finally, compile the local type-1 matrices
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->compile_local_matrix();
    }

    Statistics::toe_assembly = stamp_ass.elapsed_now();

    // accumulate sizes
    for(Index i(0); i < num_levels; ++i)
    {
      statistics.bytes[Bytes::mesh] += domain.at(i)->get_mesh_node()->bytes();
      statistics.bytes[Bytes::gate] += system_levels.at(i)->gate_sys.bytes();
      statistics.bytes[Bytes::muxer] += system_levels.at(i)->coarse_muxer_sys.bytes();
      statistics.bytes[Bytes::matrix] += system_levels.at(i)->matrix_sys.local().bytes();
      statistics.bytes[Bytes::transfer] += system_levels.at(i)->transfer_sys.bytes();
      const auto& loc_a = system_levels.at(i)->matrix_sys.local().block_a();
      const auto& loc_b = system_levels.at(i)->matrix_sys.local().block_b();
      const auto& loc_d = system_levels.at(i)->matrix_sys.local().block_d();
      statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_a.used_elements() + loc_a.rows() + Index(1));
      statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_b.used_elements() + loc_b.rows() + Index(1));
      statistics.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_d.used_elements() + loc_d.rows() + Index(1));
      statistics.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_a.template used_elements<LAFEM::Perspective::pod>());
      statistics.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_b.template used_elements<LAFEM::Perspective::pod>());
      statistics.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_d.template used_elements<LAFEM::Perspective::pod>());
    }

    /* ***************************************************************************************** */

    // get our global solver system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalSystemTransfer GlobalSystemTransfer;

    // get our local system types
    typedef typename SystemLevelType::LocalMatrixBlockA LocalMatrixBlockA;
    typedef typename SystemLevelType::LocalVeloVector LocalVeloVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // get our global solve matrix and filter
    GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
    GlobalSystemFilter& filter = the_system_level.filter_sys;

    // create new vectors
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_def = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_cor = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_def_unsynced = the_system_level.matrix_sys.create_vector_r();

    // format the vectors
    vec_sol.format();
    vec_def.format();

    // and filter them
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_def);

    {
      // count non-zeros in a and b
      statistics.counts[Counts::nnze_a] = the_system_level.matrix_sys.local().block_a().used_elements();
      statistics.counts[Counts::nnze_b] = the_system_level.matrix_sys.local().block_b().used_elements();
      statistics.counts[Counts::nnze_total] = the_system_level.matrix_sys.local().template used_elements<LAFEM::Perspective::pod>();
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    auto* mesh_part_bnd_c = the_domain_level.get_mesh_node()->find_mesh_part("bnd:c");
    auto* mesh_part_inner_u = the_domain_level.get_mesh_node()->find_mesh_part("inner:u");
    auto* mesh_part_inner_l = the_domain_level.get_mesh_node()->find_mesh_part("inner:l");

    // reference values for benchmark values
    const DataType ref_drag = DataType(dim == 2 ? 5.57953523384 : 6.18533);
    const DataType ref_lift = DataType(dim == 2 ? 0.010618948146 : 0.009401);
    const DataType ref_d_p = DataType(dim == 2 ? 0.11752016697 : 0.170826996);

    // create trace assembler for body force assembly
    Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> body_force_asm(the_domain_level.trafo);
    if(mesh_part_bnd_c != nullptr)
      body_force_asm.add_mesh_part(*mesh_part_bnd_c);
    body_force_asm.compile();

    // create trace assembler for upper flux
    Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> flux_asm_u(the_domain_level.trafo);
    if(mesh_part_inner_u != nullptr)
      flux_asm_u.add_mesh_part(*mesh_part_inner_u);
    flux_asm_u.compile();

    // create trace assembler for lower flux
    Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> flux_asm_l(the_domain_level.trafo);
    if(mesh_part_inner_l != nullptr)
      flux_asm_l.add_mesh_part(*mesh_part_inner_l);
    flux_asm_l.compile();

    // unmap pressure evaluation points p_a and p_e
    Trafo::InverseMappingData<DataType, dim> point_iv_a, point_iv_e;
    {
      typedef Trafo::InverseMapping<typename SpacePresType::TrafoType, DataType> InvMappingType;
      InvMappingType inv_mapping(the_domain_level.trafo);

      // reference pressure points
      typename InvMappingType::ImagePointType v_a, v_e;
      if(dim == 2)
      {
        v_a[0] = DataType(0.15);
        v_e[0] = DataType(0.25);
        v_a[1] = v_e[1] = DataType(0.2);
      }
      else
      {
        v_a[0] = DataType(0.45);
        v_e[0] = DataType(0.55);
        v_a[1] = v_e[1] = DataType(0.2);
        v_a[2] = v_e[2] = DataType(0.205);
      }

      // unmap points
      point_iv_a = inv_mapping.unmap_point(v_a, true);
      point_iv_e = inv_mapping.unmap_point(v_e, true);
    }

    // create characteristic function vector for circle boundary
    // this is needed for the volumetric drag/lift computation
    LAFEM::DenseVector<Mem::Main, DataType, IndexType> vec_char(vec_sol.local().template at<0>().size());
    vec_char.format(DataType(0));
    if(mesh_part_bnd_c != nullptr)
    {
      LAFEM::UnitFilter<Mem::Main, DataType, IndexType> filter_char;
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      unit_asm.add_mesh_part(*mesh_part_bnd_c);
      unit_asm.assemble(filter_char, the_domain_level.space_velo);
      filter_char.get_filter_vector().format(DataType(1));
      filter_char.filter_sol(vec_char);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a multigrid solver
    auto multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>>(domain.size_virtual());

    // array of Vanka pointers - this is only required to collect the memory usage
    // statistics of the Vankas, as we need the memory usage after factorization
    std::deque<
      std::shared_ptr<
        Solver::AmaVanka<
          typename SystemLevelType::LocalSystemMatrix,
          typename SystemLevelType::LocalSystemFilter>>> ama_vankas;

    // push levels into multigrid
    for(std::size_t i(0); i < system_levels.size(); ++i)
    {
      SystemLevelType& lvl = *system_levels.at(i);

      if((i+1) < domain.size_virtual())
      {
        auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
        ama_vankas.push_back(vanka);
        auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);
        smoother->set_min_iter(smooth_steps);
        smoother->set_max_iter(smooth_steps);
        smoother->skip_defect_calc(true); // skip defect calculation
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
#ifdef FEAT_HAVE_UMFPACK
      else if(umf_cgs)
      {
        // create UMFPACK coarse grid solver
        auto umfpack = Solver::new_generic_umfpack(lvl.local_matrix_sys);
        auto cgsolver = Solver::new_schwarz_precond(umfpack, lvl.filter_sys);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
#endif //  FEAT_HAVE_UMFPACK
      else
      {
        // create BiCGStab-AmaVanka coarse grid solver
        auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
        ama_vankas.push_back(vanka);
        auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
        auto cgsolver = Solver::new_bicgstab(lvl.matrix_sys, lvl.filter_sys, schwarz);
        cgsolver->set_max_iter(500);
        cgsolver->set_tol_rel(1e-3);
        //cgsolver->set_plot_mode(Solver::PlotMode::summary);

        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
    }

    // create our multigrid solver
    auto multigrid = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

    // create our solver
    auto solver = Solver::new_richardson(matrix, filter, 1.0, multigrid);

    solver->set_plot_name("Multigrid");

    solver->set_min_iter(min_mg_iter);
    solver->set_max_iter(max_mg_iter);
    solver->set_tol_rel(mg_tol_rel);
    solver->set_min_stag_iter(Index(3));

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a few watches
    StopWatch watch_stokes_solve;
    StopWatch watch_nonlin_loop;
    StopWatch watch_nonlin_def_asm;
    StopWatch watch_nonlin_mat_asm;
    StopWatch watch_nonlin_solver_init;
    StopWatch watch_nonlin_solver_apply;

    // initialize solver
    watch_nonlin_solver_init.start();
    multigrid_hierarchy->init_symbolic();
    solver->init_symbolic();
    watch_nonlin_solver_init.stop();

    // accumulate vanka sizes
    statistics.counts[Counts::vanka_data] = Index(ama_vankas.front()->data_size());
    statistics.bytes[Bytes::vanka] = 0ull;
    for(auto& v : ama_vankas)
      statistics.bytes[Bytes::vanka] += v->bytes();

    // load distributed solution?
    if(args.check("load-sol") > 0)
    {
      StopWatch watch_checkpoint;
      watch_checkpoint.start();

      String save_name;
      args.parse("load-sol", save_name);
      comm.print("\nReading (partitioned) initial solution from '" + save_name + "'...");

      // read file into streams
      BinaryStream bs_com, bs_sol;
      DistFileIO::read_combined(bs_com, bs_sol, save_name, comm);

      // parse vector from stream
      vec_sol.local().read_from(LAFEM::FileMode::fm_binary, bs_sol);

      // apply our solution filter in case the BCs have changed
      the_system_level.filter_sys.filter_sol(vec_sol);

      watch_checkpoint.stop();
      statistics.times[Times::checkpoint] += watch_checkpoint.elapsed();
    }
    else if(args.check("load-joined-sol") > 0)
    {
      StopWatch watch_checkpoint;
      watch_checkpoint.start();

      String save_name;
      args.parse("load-joined-sol", save_name);
      comm.print("\nReading joined initial solution from '" + save_name + "'...");

      // parse vector from file
      the_system_level.base_splitter_sys.split_read_from(vec_sol, save_name);

      // apply our solution filter in case the BCs have changed
      the_system_level.filter_sys.filter_sol(vec_sol);

      watch_checkpoint.stop();
      statistics.times[Times::checkpoint] += watch_checkpoint.elapsed();
    }
    else // solve Stokes to obtain initial solution
    {
      comm.print("\nSolving Stokes system...");

      watch_nonlin_solver_init.start();
      multigrid_hierarchy->init_numeric();
      solver->init_numeric();
      watch_nonlin_solver_init.stop();

      // solve Stokes system
      solver->set_plot_mode(Solver::PlotMode::iter);
      watch_stokes_solve.start();
      Solver::Status stokes_status = Solver::solve(*solver, vec_sol, vec_def, matrix, filter);
      watch_stokes_solve.stop();

      statistics.counts[Counts::linsol_iter] = solver->get_num_iter();

      // release solvers
      solver->done_numeric();
      multigrid_hierarchy->done_numeric();

      FEAT::Statistics::compress_solver_expressions();

      if(!Solver::status_success(stokes_status))
      {
        comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
        if(testmode)
          comm.print("Test-Mode: FAILED");
        return;
      }
    }

    // Don't we solve Navier-Stokes?
    if(!navier)
    {
      // In this case, let's compute the final unsynchronized defect for volumetric body forces computation
      matrix.apply(vec_def_unsynced, vec_sol, vec_def, -DataType(1));
      vec_def_unsynced.from_1_to_0();
    }

    // --------------------------------------------------------------------------------------------
    // SOLVE NAVIER-STOKES (IF DESIRED)
    // --------------------------------------------------------------------------------------------
    if(navier)
    {
      comm.print("\nSolving Navier-Stokes system...");

      if(!plot_mg_iter)
        solver->set_plot_mode(Solver::PlotMode::none);

      // vector of all non-linear defect norms
      std::vector<DataType> nl_defs;

      watch_nonlin_loop.start();

      // nonlinear loop
      for(Index nl_step(0); nl_step <= max_nl_iter; ++nl_step)
      {
        statistics.counts[Counts::nonlin_iter] = nl_step;

        // set up Burgers assembly job for our defect vector
        Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType> burgers_def_job(
          vec_def.local().template at<0>(), vec_sol.local().template at<0>(),
          vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature);
        burgers_def_job.deformation = defo;
        burgers_def_job.nu = -nu;
        burgers_def_job.beta = -DataType(1);

        // assemble nonlinear defect vector
        watch_nonlin_def_asm.start();
        vec_def.format();
        // assemble burgers operator defect
        the_domain_level.domain_asm.assemble(burgers_def_job);
        // compute remainder of defect vector
        the_system_level.matrix_sys.local().block_b().apply(
          vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
        the_system_level.matrix_sys.local().block_d().apply(
          vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);
        // store the unsynced and unfiltered defect for later body forces computation
        vec_def_unsynced.copy(vec_def);
        // synchronize and filter the defect
        vec_def.sync_0();
        filter.filter_def(vec_def);
        watch_nonlin_def_asm.stop();

        // compute defect norm
        const DataType def_prev = (nl_defs.empty() ? DataType(1) : nl_defs.back());
        const DataType def_nl = vec_def.norm2();
        const DataType def_improve = def_nl / def_prev;
        nl_defs.push_back(def_nl);


        String line = (newton ? "Newton: " : "Picard: ");
        line += stringify(nl_step).pad_front(2) + ": ";
        line += stringify_fp_sci(def_nl, 6) + " / ";
        line += stringify_fp_sci(def_nl/nl_defs.front()) + " / ";
        line += stringify_fp_sci(def_nl/def_prev, 3);

        if(def_nl > nl_defs.front() * DataType(1E+3))
        {
          comm.print(line);
          comm.print("\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
          if(testmode)
            comm.print("Test-Mode: FAILED");
          return;
        }
        else if(def_nl < nl_tol_abs)
        {
          comm.print(line);
          comm.print("\nNonlinear solver converged!");
          break;
        }
        else if(nl_step >= max_nl_iter)
        {
          comm.print(line);
          comm.print("\nMaximum iterations reached!");
          break;
        }
        else if((nl_step >= 3) && (DataType(0.95)*def_prev < def_nl))
        {
          comm.print(line);
          comm.print("\nNonlinear solver stagnated!");
          break;
        }

        // assemble burgers matrices on all levels
        watch_nonlin_mat_asm.start();
        {
          // get a clone of the global velocity vector
          typename SystemLevelType::GlobalVeloVector vec_conv(
            &the_system_level.gate_velo, vec_sol.local().template at<0>().clone());

          // initialize velocity norm for streamline diffusion (if enabled)
          DataType sd_v_norm = DataType(0);

          // loop over all system levels
          for(std::size_t i(0); i < system_levels.size(); ++i)
          {
            auto& loc_mat_a = system_levels.at(i)->matrix_sys.local().block_a();
            loc_mat_a.format();

            // set up Burgers assembly job
            Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType, LocalVeloVector>
              burgers_mat_job(loc_mat_a, vec_conv.local(), domain.at(i)->space_velo, cubature);

            burgers_mat_job.deformation = defo;
            burgers_mat_job.nu = nu;
            burgers_mat_job.beta = DataType(1);
            burgers_mat_job.frechet_beta = DataType(newton ? 1 : 0);
            burgers_mat_job.sd_delta = upsam;
            burgers_mat_job.sd_nu = nu;
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

            // assemble our system matrix
            domain.at(i)->domain_asm.assemble(burgers_mat_job);
            system_levels.at(i)->compile_local_matrix();

            // restrict our convection vector
            if((i+1) >= domain.size_virtual())
              break;

            // does this process have another system level?
            if((i+1) < system_levels.size())
            {
              // create a coarse mesh velocity vector
              auto vec_crs = system_levels.at(i+1)->matrix_a.create_vector_l();

              // truncate fine mesh velocity vector
              system_levels.at(i)->transfer_velo.trunc(vec_conv, vec_crs);

              // the coarse vector is our next convection vector
              vec_conv = std::move(vec_crs);
            }
            else
            {
              // this process is a child, so send truncation to parent
              system_levels.at(i)->transfer_velo.trunc_send(vec_conv);
            }
          }
        }
        watch_nonlin_mat_asm.stop();

        // initialize linear solver
        watch_nonlin_solver_init.start();
        multigrid_hierarchy->init_numeric();
        solver->init_numeric();
        watch_nonlin_solver_init.stop();

        // specify adaptive tolerance?
        if(adapt_tol && (nl_step > Index(0)))
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
            solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
            // Also make sure that we gain at least 2 digits.
            solver->set_tol_rel(1E-2);
          }
          else
          {
            // In the case if Picard iteration, we only expect linear convergence,
            // which (in analogy to Newton) leads us to the following estimate:
            //
            //      def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})
            //
            DataType abs_tol = def_nl * def_improve * DataType(0.1);
            solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
            solver->set_tol_rel(1E-2);
          }
        }
        else
        {
          // disable absolute tolerance
          solver->set_tol_abs(1E+10);
          solver->set_tol_rel(mg_tol_rel);
        }

        // solve linear system
        watch_nonlin_solver_apply.start();
        Solver::Status status = solver->apply(vec_cor, vec_def);
        watch_nonlin_solver_apply.stop();

        statistics.counts[Counts::linsol_iter] += solver->get_num_iter();

        line += String(" | ") + stringify(solver->get_num_iter()).pad_front(3) + ": "
          + stringify_fp_sci(solver->get_def_final(), 4) + " / "
          + stringify_fp_sci(solver->get_def_final() / solver->get_def_initial(), 4);
        if(adapt_tol && (nl_step > Index(0)))
          line += String(" [") + stringify_fp_sci(solver->get_tol_abs(), 4) + "]";
        comm.print(line);

        // release linear solver
        solver->done_numeric();
        multigrid_hierarchy->done_numeric();

        if(!Solver::status_success(status))
        {
          comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
          if(testmode)
            comm.print("Test-Mode: FAILED");
          return;
        }

        // update solution
        vec_sol.axpy(vec_cor, vec_sol, DataType(1));

        FEAT::Statistics::compress_solver_expressions();
        // next non-linear iteration
      }

      watch_nonlin_loop.stop();

      // end of Navier-Stokes solve
    }

    // release solver
    solver->done_symbolic();
    multigrid_hierarchy->done_symbolic();

    // save timings
    statistics.times[Times::stokes_solve] = watch_stokes_solve.elapsed();
    statistics.times[Times::nonlin_total] = watch_nonlin_loop.elapsed();
    statistics.times[Times::nonlin_asm_def] = watch_nonlin_def_asm.elapsed();
    statistics.times[Times::nonlin_asm_mat] = watch_nonlin_mat_asm.elapsed();
    statistics.times[Times::linsol_init] = watch_nonlin_solver_init.elapsed();
    statistics.times[Times::linsol_apply] = watch_nonlin_solver_apply.elapsed();

    // get multigrid timings
    statistics.times[Times::mg_defect] = multigrid_hierarchy->get_time_defect();
    statistics.times[Times::mg_smooth] = multigrid_hierarchy->get_time_smooth();
    statistics.times[Times::mg_coarse] = multigrid_hierarchy->get_time_coarse();
    statistics.times[Times::mg_transfer] = multigrid_hierarchy->get_time_transfer();

    // accumulate vanka timings
    for(auto& v : ama_vankas)
    {
      statistics.times[Times::vanka_init_sym] += v->time_init_symbolic();
      statistics.times[Times::vanka_init_num] += v->time_init_numeric();
      statistics.times[Times::vanka_apply] += v->time_apply();
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    StopWatch watch_sol_analysis;
    watch_sol_analysis.start();

    // compute drag & lift coefficients by line integration using the trace assembler
    {
      BenchBodyForceAccumulator<DataType> body_force_accum(defo, nu, v_max);
      body_force_asm.assemble_flow_accum(
        body_force_accum,
        vec_sol.local().template at<0>(),
        vec_sol.local().template at<1>(),
        the_domain_level.space_velo,
        the_domain_level.space_pres,
        cubature_postproc);
      body_force_accum.sync(comm);

      summary.drag_coeff_line = body_force_accum.drag;
      summary.lift_coeff_line = body_force_accum.lift;
      summary.drag_err_line = Math::abs((body_force_accum.drag - ref_drag) / ref_drag);
      summary.lift_err_line = Math::abs((body_force_accum.lift - ref_lift) / ref_lift);
    }

    // compute drag & lift coefficients via volume integration from unsynchronized final defect
    {
      Tiny::Vector<DataType, dim> bdf;
      assemble_bdforces_vol<DataType, dim>(bdf, vec_def_unsynced.local().first(), vec_char);

      const DataType dpf2 = DataType(2) / (dim == 2 ?
        DataType(0.100)*Math::sqr(v_max*(DataType(2)/DataType(3))) : // = 2 / (rho * U^2 * D)
        DataType(0.041)*Math::sqr(v_max*(DataType(4)/DataType(9)))); // = 2 / (rho * U^2 * D * H)

      summary.drag_coeff_vol = the_system_level.gate_sys.sum(dpf2 * bdf[0]);
      summary.lift_coeff_vol = the_system_level.gate_sys.sum(dpf2 * bdf[1]);

      summary.drag_err_vol = Math::abs((summary.drag_coeff_vol - ref_drag) / ref_drag);
      summary.lift_err_vol = Math::abs((summary.lift_coeff_vol - ref_lift) / ref_lift);
    }

    // compute pressure difference
    {
      // evaluate pressure
      auto pval_a = Assembly::DiscreteEvaluator::eval_fe_function(
        point_iv_a, vec_sol.local().template at<1>(), the_domain_level.space_pres);
      auto pval_e = Assembly::DiscreteEvaluator::eval_fe_function(
        point_iv_e, vec_sol.local().template at<1>(), the_domain_level.space_pres);

      // compute pressure mean
      const DataType p_a = pval_a.mean_value_dist(comm);
      const DataType p_e = pval_e.mean_value_dist(comm);
      const DataType d_p = p_a - p_e;

      // compute error to reference values
      summary.pres_diff = d_p;
      summary.pres_err = Math::abs((d_p - ref_d_p) / ref_d_p);
    }

    // compute flow through upper region
    {
      XFluxAccumulator<DataType> flux_accum_u;
      flux_asm_u.assemble_flow_accum(
        flux_accum_u,
        vec_sol.local().template at<0>(),
        vec_sol.local().template at<1>(),
        the_domain_level.space_velo,
        the_domain_level.space_pres,
        cubature_postproc);
      flux_accum_u.sync(comm);

      summary.flux_upper = flux_accum_u.flux / DataType(2);
    }

    // compute flow through lower region
    {
      XFluxAccumulator<DataType> flux_accum_l;
      flux_asm_l.assemble_flow_accum(
        flux_accum_l,
        vec_sol.local().template at<0>(),
        vec_sol.local().template at<1>(),
        the_domain_level.space_velo,
        the_domain_level.space_pres,
        cubature_postproc);
      flux_accum_l.sync(comm);

      summary.flux_lower = flux_accum_l.flux / DataType(2);
    }

    // perform analysis of velocity field
    {
      Assembly::VelocityInfo<DataType, dim> vi = Assembly::VelocityAnalyser::compute(
        vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature_postproc);
      vi.synchronize(comm);

      summary.velo_info = vi;
    }

    watch_sol_analysis.stop();
    statistics.times[Times::sol_analysis] = watch_sol_analysis.elapsed();

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if(args.check("save-sol") >= 0)
    {
      StopWatch watch_checkpoint;
      watch_checkpoint.start();

      String save_name;
      if(args.parse("save-sol", save_name) < 1)
      {
        save_name = String("dfg95-cc") + stringify(dim) + "d-bench1";
        save_name += "-lvl" + stringify(the_domain_level.get_level_index());
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

      watch_checkpoint.stop();
      statistics.times[Times::checkpoint] += watch_checkpoint.elapsed();
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if(args.check("save-joined-sol") >= 0)
    {
      StopWatch watch_checkpoint;
      watch_checkpoint.start();

      String save_name;
      if(args.parse("save-joined-sol", save_name) < 1)
      {
        save_name = String("dfg95-cc") + stringify(dim) + "d-bench1";
        save_name += "-joined-lvl" + stringify(the_domain_level.get_level_index()) + ".bin";
        // don't include number of processes, because its invariant
      }

      comm.print("\nWriting joined solution to '" + save_name + "'");

      // write out in binary format
      the_system_level.base_splitter_sys.join_write_out(vec_sol, save_name);

      watch_checkpoint.stop();
      statistics.times[Times::checkpoint] += watch_checkpoint.elapsed();
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    //*
    if(args.check("vtk") >= 0)
    {
      StopWatch watch_vtk;
      watch_vtk.start();

      // build VTK name
      String vtk_name, vtk_name2;
      int npars = args.parse("vtk", vtk_name, vtk_name2);
      if(npars < 1)
      {
        vtk_name = String("dfg95-cc") + stringify(dim) + "d-bench1";
        vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
        vtk_name += "-n" + stringify(comm.size());
      }

      comm.print("\nWriting VTK output to '" + vtk_name + ".pvtu'");

      if(!vtk_name.empty())
      {
        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

        // project velocity and pressure
        exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());

        // project pressure
        Cubature::DynamicFactory cub("gauss-legendre:2");
        LAFEM::DenseVector<Mem::Main, DataType, Index> vtx_p;
        Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), the_domain_level.space_pres, cub);

        // write pressure
        exporter.add_cell_scalar("pressure", vtx_p.elements());

        // finally, write the VTK file
        exporter.write(vtk_name, comm);
      }

      // write refined VTK?
      if(!vtk_name2.empty())
      {
        comm.print("Writing refined VTK output to '" + vtk_name2 + ".pvtu'");

        // refine mesh
        Geometry::StandardRefinery<MeshType> refinery(the_domain_level.get_mesh());
        MeshType mesh(refinery);

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(mesh);

        // project velocity and pressure
        exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());

        // finally, write the VTK file
        exporter.write(vtk_name2, comm);
      }

      watch_vtk.stop();
      statistics.times[Times::vtk_write] = watch_vtk.elapsed();
    }

    watch_total_run.stop();
    statistics.times[Times::total_run] = watch_total_run.elapsed();

    {
      MemoryUsage mi;
      statistics.bytes[Bytes::peak_p] = mi.get_peak_physical();
      statistics.bytes[Bytes::peak_v] = mi.get_peak_virtual();
    }

    statistics.sync(comm);

    comm.print(String("\n") + String(80, '=') + "\n");
    comm.print(summary.format());
    comm.print(statistics.format());

    // print multigrid timings
    if(comm.rank() == 0)
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
      comm.print(FEAT::Statistics::get_formatted_flops(statistics.times[Times::stokes_solve] + statistics.times[Times::linsol_apply], comm.size()));
      comm.print(FEAT::Statistics::get_formatted_times(statistics.times[Times::stokes_solve] + statistics.times[Times::linsol_apply]));
      comm.print(FEAT::Statistics::get_formatted_solver_internals("default"));
      comm.print("\n");
      comm.print(FEAT::Statistics::get_formatted_solver_tree("default").trim());
    }

    if(testmode)
      comm.print("\nTest-Mode: PASSED");
  }

  /*template<int dim_>
  void run_dim(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
    // define our mesh type
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // create a time-stamp
    TimeStamp time_stamp;

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(mesh_reader);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());
    comm.print("Transformation: Standard");

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }*/

  template<int dim_>
  void run_dim_iso(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
    // define our mesh type
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, DataType> MeshType;
    typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // create a time-stamp
    TimeStamp time_stamp;

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);

    // if we want to write/read a joined solution, we also need to tell the domain control to keep the base levels
    if((args.check("save-joined-sol") >= 0) || (args.check("load-joined-sol") >= 0))
      domain.keep_base_levels();

    domain.create(mesh_reader);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());
    comm.print("Transformation: Isoparametric:2");

    // get circle chart
    auto chart = domain.get_atlas().find_mesh_chart("circle");

    // loop over all physical levels and add meshpart charts
    for(Index ilvl(0); ilvl < domain.size_physical(); ++ilvl)
    {
      auto mesh_node = domain.at(ilvl)->get_mesh_node();
      const auto mpart = mesh_node->find_mesh_part("bnd:c");
      if((chart != nullptr) && (mpart != nullptr))
        domain.at(ilvl)->trafo.add_meshpart_chart(*mpart, *chart);
    }

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("\nRun-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }

  void main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));
    comm.print("Floating Point Type: " + String(fp_typename) + " precision");

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("level");
    args.support("vtk");
    args.support("mesh");
    args.support("nu");
    args.support("upsam");
    args.support("max-nl-iter");
    args.support("min-mg-iter");
    args.support("max-mg-iter");
    args.support("plot-mg-iter");
    args.support("smooth-steps");
    args.support("smooth-damp");
    args.support("mg-tol-rel");
    args.support("nl-tol-abs");
    //args.support("no-adapt");
    args.support("picard");
    args.support("v-max");
    args.support("defo");
    args.support("stokes");
    args.support("test-mode");
    args.support("ext-stats");
    args.support("no-umfpack");
    args.support("load-sol");
    args.support("load-joined-sol");
    args.support("save-sol");
    args.support("save-joined-sol");
    //args.support("isoparam");

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

    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run 2D or 3D ?
    //if(args.check("isoparam") >= 0)
    {
      if(mesh_type == "conformal:hypercube:2:2")
        run_dim_iso<2>(args, comm, mesh_reader);
      else if(mesh_type == "conformal:hypercube:3:3")
        run_dim_iso<3>(args, comm, mesh_reader);
      else
      {
        comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
        FEAT::Runtime::abort();
      }
    }
    /*else
    {
      if(mesh_type == "conformal:hypercube:2:2")
        run_dim<2>(args, comm, mesh_reader);
      else if(mesh_type == "conformal:hypercube:3:3")
        run_dim<3>(args, comm, mesh_reader);
      else
      {
        comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
        FEAT::Runtime::abort();
      }
    }*/
  }
} // namespace DFG95

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialize(argc, argv);
  try
  {
    DFG95::main(argc, argv);
  }
  catch(std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return FEAT::Runtime::finalize();
}
