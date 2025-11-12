// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// PURPOSE
// -------
// This application demonstrates how to assemble a 2D/3D Pressure Poisson equation in mixed formulation by
// using an inf-sup stable finite element pair such as Q2/P1dc on a mesh that is read in from a mesh file and
// afterwards solve it using Trilinos FROSch solver. To keep things simple, this application does *NOT* use
// FEAT's internal solver (wrapper) classes, but instead exposes all calls required to use HYPRE's BoomerAMG
// solver directly in the main function. This application makes use of the so-called "pre-multiplied discontinuous
// diagonal Schur complement matrix" to discretize the Pressure Poisson operator based the lumped velocity mass
// matrix and the B/D matrices from the Stokes system. This matrix class also offers the conversion of the
// system into an algebraic DOF partitioning form, which effectively describes a row-wise partitioned distributed
// CSR matrix system that can be for example directly converted into HYPRE's internal ParCSR IJ-Matrix structures.
//
// USAGE
// -----
// run this application with MPI with an arbitrary number of MPI ranks
//
// command line options:
//
// --mesh <path-to-meshfile>
// MANDATORY: This command line option specifies the path to the mesh file to be read in.
// The most commonly used mesh file will be "data/meshes/unit-square-quad.xml", but in principle
// any other 2D quadrilateral mesh file could be used; these mesh files all have "quad" in their names.
//
// --level <level-max>
// MANDATORY: This command line option specifies the maximum refinement level <level-max>, which
// defines how many times the input mesh has to be refined regularly to obtain the final mesh that
// the PDE is to be solved on.
//
// --neumann <meshparts...>
// Optional: Specifies the names of the mesh-parts which make up the Neumann boundary region(s).
// If specified without any meshparts names, then ALL boundary meshparts are treated as Neumann boundaries.
//
// --vtk [<filename>]
// Optional: Specifies that the application should write out PVTU/VTU files at the end that can be
// used to visualize the solutions.
//
// --dump
// Optional: Enables additional debug output about the partitioning.
//

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/global/alg_dof_parti.hpp>
#include <kernel/global/pmdcdsc_matrix.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/fgmres.hpp>

#include <control/domain/unit_cube_domain_control.hpp>
#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_mixed.hpp>
#include <control/statistics.hpp>

#include <kernel/solver/frosch.hpp>
#include <vector>

#ifndef FEAT_HAVE_MPI
#error This application must be configured and compiled with MPI support enabled.
#endif

#ifndef SCALEXA_DIM
#error You have to define the dimension of your app
#define SCALEXA_DIM 2
#endif

using DataType = double;
using IndexType = FEAT::Index;

namespace ScalexaPoisson
{
  using namespace FEAT;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Use the following typedefs to define the shape type of the mesh as well as the finite element space

  // select the element shape to be used; must be one of:
  //typedef Shape::Triangle ShapeType;    // 2D Triangle Elements
  #if SCALEXA_DIM == 2
  typedef Shape::Quadrilateral ShapeType; // 2D Quadrilateral Elements
  #else
  typedef Shape::Hexahedron ShapeType;  // 3D Hexahedral Elements
  #endif
  //typedef Shape::Tetrahedron ShapeType; // 3D Tetrahedral Elements

  // define the mesh type; this is always
  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  // define the trafo type; this is usually
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  // select an inf-sup stable finite element space pair to be used:
  typedef Space::Lagrange2::Element<TrafoType>       VeloSpaceType; // P2 or Q2 finite elements
  typedef Space::Discontinuous::ElementP1<TrafoType> PresSpaceType; // P1dc finite elements

  // define the domain level type; for this application it is typically
  typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, VeloSpaceType, PresSpaceType> DomainLevelType;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // this main function is called by the actual main function at the end of this source file
  void main(int argc, char* argv [])
  {
    // stop watch for time measurement
    StopWatch watch_total;
    watch_total.start();

    // create world communicator (this is our MPI_Comm wrapper object)
    Dist::Comm comm(Dist::Comm::world());
    comm.print(String(100u, '*'));

    // dump system call
    {
      String s("Arguments: ");
      s.append(argv[0]);
      for(int i(1); i < argc; ++i)
        s.append(" ").append(argv[i]);
      comm.print(s);
    }

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("mesh");
    args.support("level");
    args.support("dump");
    args.support("vtk");
    args.support("neumann");

    // FROSch
    Solver::Trilinos::add_supported_fpl_args(args);

    args.support("krylov_dim");
    args.support("krylov_tol_rel");
    args.support("krylov_tol_abs");
    args.support("krylov_inner_res_scale");
    args.support("krylov_maxit");
    args.support("krylov_minit");
    args.support("krylov_info");

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
    if(args.check("level") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
      FEAT::Runtime::abort();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    StopWatch watch_domain_setup;
    watch_domain_setup.start();

    // create our domain control
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, false);

    // parse command line arguments
    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);

    // create domain partitioning; this is where the magic happens
    domain.create(args.query("mesh")->second);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());

    // dump domain info if desired
    if(args.check("dump") >= 0)
    {
      comm.print("\nDomain Layers:");
      comm.allprint(domain.dump_layers());
      comm.print("\nDomain Layer Levels:");
      comm.allprint(domain.dump_layer_levels());
      comm.print("\nDomain Virtual Levels:");
      comm.allprint(domain.dump_virt_levels());
    }

    watch_domain_setup.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    StopWatch watch_system_asm, watch_pmdcdsc_init;
    watch_system_asm.start();

    // define our arch types
    typedef double DataType; // usually double
    typedef Index IndexType; // usually 64 bit unsigned int

    // get the dimension
    static constexpr int dim = ShapeType::dimension;

    // define our system level
    typedef Control::ScalarMixedSystemLevel<dim, DataType, IndexType> SystemLevelType;

    // get the finest domain level
    DomainLevelType& the_domain_level = *domain.front();

    // get the list of all boundary mesh-parts in the mesh
    std::set<String> all_meshpart_names, bnd_meshpart_names;
    for(const String& name : the_domain_level.get_mesh_node()->get_mesh_part_names(true))
    {
      all_meshpart_names.insert(name);
      if(name.starts_with("bnd:"))
        bnd_meshpart_names.insert(name);
      /* Gendie */
      else if(name.starts_with("inflow") || name.starts_with("newwall") || name.starts_with("z+"))
        bnd_meshpart_names.insert(name);
    }

    // determine which boundary mesh-parts are Dirichlet and which are Neumann
    std::set<String> boundary_dirichlet, boundary_neumann;
    if(args.check("neumann") < 0)
    {
      // all boundary mesh-parts belong to Dirichlet boundary
      boundary_dirichlet = bnd_meshpart_names;
    }
    else if(args.check("neumann") == 0)
    {
      // all boundary mesh-part belong to Neumann boundary
      boundary_neumann = bnd_meshpart_names;
    }
    else
    {
      // mixed Dirichlet-Neumann boundary conditions
      for(const String& name : args.query("neumann")->second)
      {
        if(all_meshpart_names.find(name) == all_meshpart_names.end())
        {
          comm.print(String("ERROR: Neumann meshpart '") + name + String("' not found!"));
          comm.print(String("Meshpart names for the given mesh are: '") + stringify_join(all_meshpart_names, "', '") + String("'"));
          Runtime::abort();
        }
        boundary_neumann.insert(name);
        auto it = bnd_meshpart_names.find(name);
        if(it != bnd_meshpart_names.end())
          bnd_meshpart_names.erase(it);
      }
      boundary_dirichlet = bnd_meshpart_names;
    }

    if(boundary_dirichlet.empty())
      comm.print("Dirichlet Boundaries: - none -");
    else
      comm.print(String("Dirichlet Boundaries: ") + stringify_join(boundary_dirichlet, " "));
    if(boundary_neumann.empty())
      comm.print("Neumann   Boundaries: - none -");
    else
      comm.print(String("Neumann   Boundaries: ") + stringify_join(boundary_neumann, " "));

    // create a system level for the finest domain level
    SystemLevelType the_system_level;

    // define our cubature rule
    const String cubature("gauss-legendre:" + stringify(VeloSpaceType::local_degree+1));

    // prepare domain assembler
    the_domain_level.domain_asm.compile_all_elements();

    // assemble gate; this is the main communication structure
    the_system_level.assemble_gates(domain.front());

    // assemble velocity mass matrix
    the_system_level.assemble_velo_struct(the_domain_level.space_velo);
    {
      Assembly::Common::IdentityOperatorBlocked<dim> identity_op;
      Assembly::assemble_bilinear_operator_matrix_1(the_domain_level.domain_asm, the_system_level.matrix_a.local(),
        identity_op, the_domain_level.space_velo, cubature);
    }

    // assemble pressure gradient/velocity divergence matrices
    the_system_level.assemble_grad_div_matrices(the_domain_level.space_velo, the_domain_level.space_pres, cubature);

    // create a velocity unit filter
    LAFEM::UnitFilterBlocked<DataType, IndexType, dim> local_velocity_filter(the_domain_level.space_velo.get_num_dofs());

    // loop over all Neumann boundary parts
    {
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      for(const String& name : boundary_neumann)
      {
        // assemble velocity slip filter for pressure Neumann boundary
        const auto* mpart = the_domain_level.get_mesh_node()->find_mesh_part(name);
        if(mpart)
          unit_asm.add_mesh_part(*mpart);
      }
      unit_asm.assemble(local_velocity_filter, the_domain_level.space_velo);
    }

    // compile filters
    the_system_level.compile_system_filter();
    the_system_level.assemble_global_filters();

    // compile system matrix
    the_system_level.compile_system_matrix();

    // get the global matrix types for B and D (=B^T)
    typedef typename SystemLevelType::GlobalMatrixBlockB GlobalMatrixBlockB;
    typedef typename SystemLevelType::GlobalMatrixBlockD GlobalMatrixBlockD;
    typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;
    typedef typename SystemLevelType::LocalPresVector LocalPresVector;
    //typedef typename SystemLevelType::LocalScalarMatrix LocalScalarMatrix;

    // invert lumped velocity mass matrix
    the_system_level.lumped_matrix_a.component_invert(the_system_level.lumped_matrix_a, DataType(1));

    // apply velocity slip filters onto the lumped matrix
    local_velocity_filter.filter_def(the_system_level.lumped_matrix_a.local());
    //local_velocity_filter.filter_offdiag_row_mat(the_system_level.matrix_b.local());

    // define a global Pre-Multiplied Discontinuous Diagonal Schur-Complement Matrix
    Global::PMDCDSCMatrix<GlobalMatrixBlockB, GlobalMatrixBlockD> system_matrix(
      the_system_level.lumped_matrix_a, the_system_level.matrix_b, the_system_level.matrix_d);

    watch_system_asm.stop();
    watch_pmdcdsc_init.start();

    // initialize the PMDCDSC matrix
    system_matrix.init();

    watch_pmdcdsc_init.stop();
    watch_system_asm.start();

    // create four new vectors
    GlobalPresVector vec_sol = system_matrix.create_vector_r();
    GlobalPresVector vec_rhs = system_matrix.create_vector_r();
    GlobalPresVector vec_def = system_matrix.create_vector_r();

    // format vectors to zero
    vec_sol.format();
    vec_rhs.format();

    // assemble the constant force functional
    Analytic::Common::ConstantFunction<ShapeType::dimension, DataType> one_func(DataType(1));
    Assembly::assemble_force_function_vector(the_domain_level.domain_asm, vec_rhs.local(),
      one_func, the_domain_level.space_pres, cubature);

    // synchronize the rhs vector
    vec_rhs.sync_0();

    watch_system_asm.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //
    // >> S C A L E X A  -  S C A L E X A  -  S C A L E X A  -  S C A L E X A  -  S C A L E X A  -  S C A L E X A >> //
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    comm.print(String(100u, '-'));

    // two stop watches for time measurement
    StopWatch watch_adp_convert;

#if 0

    // Note: by default, IndexType is a typedef for unsigned long and DataType is a typedef for double

    watch_adp_convert.start();

    Index num_total_dofs(0u), first_owned_dof(0u);

    // initialize ADP for this type of system matrix and obtain our local matrix slice
    LocalScalarMatrix local_matrix_slice = system_matrix.asm_adp_symbolic(first_owned_dof, num_total_dofs);

    // the number of global DOFS that this process owns:
    const IndexType num_owned_dofs = local_matrix_slice.rows();

    // perform numeric initialization
    system_matrix.asm_adp_numeric(local_matrix_slice);

    watch_adp_convert.stop();

    // print basic information
    comm.print("Basic Information of Algebraic DOF Partitioning on all ranks:");
    comm.print("Total Number of Global DOFs: " + stringify(num_total_dofs));
    if(comm.size() < 1000) // don't print on more than 1000 MPI ranks
      comm.allprint("First Owned DOF:" + stringify(first_owned_dof).pad_front(9) +
        "   Number of Owned DOFs:" + stringify(num_owned_dofs).pad_front(9));

    // get the basic information about the partitioned CSR matrix for this process
    // note: in our format, all entries in a row are sorted in ascending column order!

    // the number of rows for our CSR partition; this should coincide with num_owned_dofs
    const IndexType pcsr_num_rows = local_matrix_slice.rows();

    // the number of columns for our CSR partition; this should coincide with num_total_dofs
    const IndexType pcsr_num_cols = local_matrix_slice.columns();

    // the number of non-zeros for out CSR partition
    const IndexType pcsr_num_nze = local_matrix_slice.used_elements();

    // print the basic PCSR information
    if(comm.size() < 1000) // don't print on more than 1000 MPI ranks
      comm.allprint("#Rows: " + stringify(pcsr_num_rows).pad_front(9) + "  #Cols:" +
        stringify(pcsr_num_cols).pad_front(9) + "  #NZEs:" + stringify(pcsr_num_nze).pad_front(9));
    comm.print("");
#endif
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //
    //  F R O S C H  -  F R O S C H  -  F R O S C H  -  F R O S C H  -  F R O S C H  -  F R O S C H  -  F R O S C H  //
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    StopWatch watch_frosch_setup, watch_frosch_solve;
    watch_frosch_setup.start();

    comm.print(String(100u, '-'));

    std::shared_ptr<Solver::Trilinos::FROSchParameterList> params = std::make_shared<Solver::Trilinos::FROSchParameterList>(comm, dim, Solver::Trilinos::FROSchParameterList::PRESSUREPOISSON);

    if(!params->parse_args(args))
      XASSERTM(false, "Parse went wrong");
    params->create_core();

    Index krylov_dim = 10;
    DataType tol_rel = 1E-8, tol_abs = 1E6;
    DataType inner_res_scale = 0;
    Index krylov_maxit = 1000;
    Index krylov_minit = 10;
    bool krylov_info = false;
    args.parse("krylov_dim", krylov_dim);
    args.parse("krylov_tol_rel", tol_rel);
    args.parse("krylov_tol_abs", tol_abs);
    args.parse("krylov_inner_res_scale", inner_res_scale);
    args.parse("krylov_maxit", krylov_maxit);
    args.parse("krylov_minit", krylov_minit);
    args.parse("krylov_info", krylov_info);

    auto frosch_precond = Solver::new_frosch(system_matrix, the_system_level.filter_pres, *params);
#if 1
    // auto frosch_solver  = Solver::new_gmres(system_matrix,
    //                                         the_system_level.filter_pres,
    //                                         krylov_dim,
    //                                         inner_res_scale,
    //                                         frosch_precond);
    // frosch_solver->set_plot_name("GMRES with FROSch");
    auto frosch_solver  = Solver::new_fgmres(system_matrix,
                                             the_system_level.filter_pres,
                                             krylov_dim,
                                             inner_res_scale,
                                             frosch_precond);
    frosch_solver->set_plot_name("FGMRES with FROSch");
#else
    auto frosch_solver  = Solver::new_pcg(system_matrix,
                                          the_system_level.filter_pres,
                                          frosch_precond);
    frosch_solver->set_plot_name("PCG with FROSch");
#endif

    // enable plotting
    if(krylov_info)
      frosch_solver->set_plot_mode(Solver::PlotMode::iter);

    // set tolerance
    frosch_solver->set_tol_rel(tol_rel);
    frosch_solver->set_tol_abs(tol_abs);
    frosch_solver->set_max_iter(krylov_maxit);

    frosch_solver->init();

    Statistics::reset();

    comm.print("");

    watch_frosch_setup.stop();
    watch_frosch_solve.start();

    /* __TODO__: Ist hier the_system_level.matrix_sys, the_system_level.filter_sys richtig? */
//     result = Solver::solve(*frosch_solver, vec_sol, vec_rhs,
//                       the_system_level.matrix_sys, the_system_level.filter_sys);
    auto result = Solver::solve(*frosch_solver, vec_sol, vec_rhs,
                                system_matrix, the_system_level.filter_pres);
    watch_frosch_solve.stop();

    if (!Solver::status_success(result)) {
      comm.print("FROSch Solver execution FAILED, with status: " + stringify(result));
    }

    if (!Solver::status_success(result)) {
      comm.print("FROSch Solver execution FAILED, with status: " + stringify(result));
    }

    // release solver
    frosch_solver->done();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
    //  F R O S C H  -  F R O S C H  -  F R O S C H  -  F R O S C H  -  F R O S C H  -  F R O S C H  -  F R O S C H  //
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // compute and print final defect norm
    system_matrix.apply(vec_def, vec_sol, vec_rhs, -1.0);
    comm.print("\nFinal Defect Norm: " + stringify_fp_sci(vec_def.norm2()));
    comm.print(String(100u, '-'));

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
    // << S C A L E X A  -  S C A L E X A  -  S C A L E X A  -  S C A L E X A  -  S C A L E X A  -  S C A L E X A << //
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // write VTK output if desired
    if (args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./scalexa-poisson-08");
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(comm.size());
      args.parse("vtk", vtk_name);

      comm.print("\nWriting VTU output to '" + vtk_name + ".pvtu'...");

      // Refine the mesh one time
      auto refined_node = the_domain_level.get_mesh_node()->refine_unique();

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(*refined_node->get_mesh());

      // refine pressure solution
      Cubature::DynamicFactory cubature_factory("refine:midpoint");
      Cubature::Rule<ShapeType, DataType, DataType> cubature_rule(Cubature::ctor_factory, cubature_factory);
      const int num_pts = cubature_rule.get_num_points();

      const Index num_elems = the_domain_level.space_pres.get_mesh().get_num_elements();
      std::vector<double> vsol(num_elems*Index(num_pts));

      // assembly traits
      typedef Assembly::AsmTraits1<DataType, PresSpaceType, TrafoTags::none, SpaceTags::value> AsmTraits;

      // fetch the trafo
      const typename AsmTraits::TrafoType& trafo = the_domain_level.space_pres.get_trafo();

      // create a trafo evaluator
      typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

      // create a space evaluator and evaluation data
      typename AsmTraits::SpaceEvaluator space_eval(the_domain_level.space_pres);

      // create a dof-mapping
      typename AsmTraits::DofMapping dof_mapping(the_domain_level.space_pres);

      // create trafo evaluation data
      typename AsmTraits::TrafoEvalData trafo_data;

      // create space evaluation data
      typename AsmTraits::SpaceEvalData space_data;

      // create local vector data
      typename AsmTraits::template TLocalVector<DataType> local_vector;

      // create matrix scatter-axpy
      typename LocalPresVector::GatherAxpy gather_axpy(vec_sol.local());

      // loop over all elements
      for(Index ielem(0), i(0); ielem < num_elems; ++ielem)
      {
        // format local vector
        local_vector.format();

        // initialize dof-mapping
        dof_mapping.prepare(ielem);

        // gather local vector data
        gather_axpy(local_vector, dof_mapping);

        // finish dof-mapping
        dof_mapping.finish();

        // prepare trafo evaluator
        trafo_eval.prepare(ielem);

        // prepare space evaluator
        space_eval.prepare(trafo_eval);

        // fetch number of local dofs
        const int num_loc_dofs = space_eval.get_num_local_dofs();

        for(int k(0); k < num_pts; ++k, ++i)
        {
          // compute trafo data
          trafo_eval(trafo_data, cubature_rule.get_point(k));

          // compute basis function data
          space_eval(space_data, trafo_data);

          // compute function value
          DataType value = DataType(0);
          for(int j(0); j < num_loc_dofs; ++j)
            value += local_vector[j] * space_data.phi[j].value;

          // save value
          vsol[i] = value;
        }

        // finish evaluators
        space_eval.finish();
        trafo_eval.finish();
      }

      // export pressure
      exporter.add_cell_scalar("solution", vsol.data());

      // finally, write the VTK file
      exporter.write(vtk_name, comm);
    }

    // print elapsed runtime
    watch_total.stop();
    const double time_total = watch_total.elapsed();
    comm.print("Time for Domain Setup.........: " + watch_domain_setup.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_domain_setup.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for System Assembly......: " + watch_system_asm.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_system_asm.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for PMDCDSC Matrix setup.: " + watch_pmdcdsc_init.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_pmdcdsc_init.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for ADP system conversion: " + watch_adp_convert.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_adp_convert.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for FROsch solver setup..: " + watch_frosch_setup.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_frosch_setup.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for FROsch solver apply..: " + watch_frosch_solve.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_frosch_solve.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Total Runtime.................: " + watch_total.elapsed_string().pad_front(10) + " seconds");
  }
} // namespace ScalexaPoisson

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialize(argc, argv);
  try
  {
    ScalexaPoisson::main(argc, argv);
  }
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
  return FEAT::Runtime::finalize();
}
