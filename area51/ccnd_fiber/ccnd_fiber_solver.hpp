#pragma once
#ifndef CCND_FIBER_STEADY_SOLVER_HPP
#define CCND_FIBER_STEADY_SOLVER_HPP 1

#include <area51/ccnd_fiber/ccnd_fiber_common.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <area51/ccnd_fiber/ccnd_steady_function.hpp> //for steady boundary conditions and test functions
#include <kernel/assembly/common_functionals.hpp>

namespace CCND_FIBER
{
  template<typename DomainLevel_>
  class CCNDFiberSolver;

  /**
     * \brief Interface for a paralizable local MPSC solver for the extended, time dependent, incompressible Navier-Stokes equations.
     *
     *
     * \tparam DomainLevel_ The DomainLevel on which the partioning of the mesh is based, refer to \ref FEAT::Domain::Control.
     *
     * This class represents an interface to solve \f[(u,v)\f] for the equations
     *
     * \f[\delta_t v + \rho(v\cdot \nabla) v - \nabla \cdot \sigma(v) + \nabla p = f\f] in \f[\Omega\f]
     * \f[\nabla v = 0 \f] in \f[\Omega\f]
     * \f[v = v_0\f] on \f[\Gamma_I\f]
     * \f[v = v_{noslip}\f] on \f[\Gamma_{noslip}\f]
     * \f[\sigma(v) \cdot n - pn = h\f] on \f[\partial\Omega\setminus (\Gamma_I\cup\Gamma_{noslip})\f]
     *
     * where \f[sigma(v) = 2\mu(D(v) + N_s(A\cdot D(v) + D(v) \cdot A) + N_p \underline{A} : D(v))\f].
     *
     * The solver itself can be adjusted through parameters. This can either be set through the connected methods or initialized through an ArgsParser, see \ref FEAT::SimpleArgParser.
     * The most important ones are listed below as their ArgsParser argument:
     *
     *
     * ------------------------------------
     * Basic Setup and Mandatory Parameters
     * ------------------------------------
     * While not directly a part of our solver, this are the two mandatory parameters you have to give your MeshReade resp. PartiDomainControl,
     * so we list these for the sake of completness.
     *
     * --mesh <meshfiles...>
     * Specifies the input mesh file(s).
     *
     * --level <level-max> [levels...] <level-min>
     * Specifies the mesh refinement levels in the syntax according to Control::PartiDomainControl.
     *
     *
     * ----------------------------
     * System Definition Parameters
     * ----------------------------
     * Some of the basic parameters of the underlying system can be configured by the following
     * parameters.
     *
     * --mu <mu>
     * Specifies the viscosity parameter mu for the diffusion operator. Defaults to 1E-1.
     *
     * --v-max <vmax>
     * Specifies the maximum velocity of the inflow boundary condition.
     * Defaults to 1.5 in 2D and 2.25 in 3D.
     *
     * --n-s <ns>
     * Specifies the shear rate number N_s. Defaults to 0.
     *
     * --n-p <np>
     * Specifies the particle number N_p. Defaults to 1.
     *
     *
     * ------------------------
     * Time Stepping Parameters
     * ------------------------
     * This section describes the parameters that control the time stepping scheme.
     *
     * --delta-t <dt>
     * Specifies the time step size. Defaults to 1E-2.
     *
     * --t-expo <0|1|2>
     * Specifies the time step extrapolation order, which is used to define an initial guess u_k for
     * the non-linear solver in each time step k, where:
     *   --t-expo 0 use constant extrapolation from previous time step, i.e.
     *              u_k := u_{k-1}
     *   --t-expo 1 use linear extrapolation from two previous time steps, i.e.
     *              u_k := 2*u_{k-1} - u_{k-2}
     *   --t-expo 2 use quadratic extrapolation from three previous time steps, i.e.
     *              u_k := 3*u_{k-1} - 3*u_{k-2}) + u_{k-3}
     * Defaults to 2.
     *
     * --euler
     * Specifies if an implicit Euler instead of BDF2 shall be used.
     *
     * -------------------------------
     * Solver Configuration Parameters
     * -------------------------------
     * This section describes the parameters that control the non-linear Newton/Picard solver
     * as well as its multigrid preconditioner and its smoother component.
     *
     * --verbose
     * Print more information to terminal-
     *
     * --picard
     * If specified, the nonlinear system in each time step will be solved using a simple
     * Picard iteration instead of the Newton iteration.
     *
     * --min-nl-iter <N>
     * Specifies the minimum number of nonlinear (Newton/Picard) solver iterations per time step.
     * Defaults to 1.
     *
     * --max-nl-iter <N>
     * Specifies the maximum number of nonlinear (Newton/Picard) solver iterations per time step.
     * Defaults to 20.
     *
     * --min-mg-iter <N>
     * Specifies the minimum number of multigrid iterations per nonlinear solver iteration.
     * Defaults to 1.
     *
     * --max-mg-iter <N>
     * Specifies the maximum number of multigrid iterations per nonlinear solver iteration.
     * Defaults to 25.
     *
     * --no-umfpack
     * If specified, the multigrid solver will use a BiCGStab-AmaVanka solver as the coarse grid solver
     * instead of the UMFPACK direct solver. Note that UMFPACK is only used if it is included in the
     * build id and if the coarse system is solved on a single process.
     *
     * --nl-tol-abs <tol>
     * Specifies the absolute tolerance for the nonlinear solver. Defaults to 1E-8.
     *
     * --mg-tol-rel <tol>
     * If given, specifies the relative tolerance for the multigrid solver.
     * If not given, then the tolerance for the multigrid solver is chosen in an adaptive
     * manner depending on the two previous nonlinear solver defects, which is the default case.
     * The adaptive tolerance is chosen in each nonlinear iteration by analyzing the nonlinear
     * defect improvement in the previous nonlinear (Newton/Picard) solver iteration in the
     * following manner: Let def_{j} and def_{j-1} denote the two previous nonlinear defect norms,
     * then the next nonlinear defect norm def_{j+1} should approximately fulfill the equation
     *
     *        (def_{j+1} / def_{j}) \approx (def_{j} / def_{j+1})^C
     *
     * where C \in {1,2} is the convergence speed of the nonlinear solver, i.e. C=2 for Newton
     * and C=1 for Picard iteration. Multiplying the above equation by def_{j} gives us an
     * estimate for the next nonlinear defect norm def_{j+1}:
     *
     *        def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})^C
     *
     * To obtain an absolute tolerance for the next multigrid solver application, we simply
     * multiply the above estimate by 0.1.
     *
     * \author Maximilian Esser
     */
  template<typename DomainLevel_>
  class CCNDUnsteadyInterface
  {
  public:
    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainLevelType::SpaceVeloType SpaceVeloType;
    typedef typename DomainLevelType::SpacePresType SpacePresType;
    //     typedef Space::Lagrange1::Element<typename DomainLevelType::TrafoType> SpaceOrientationType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim*(dim+1)/2> DenseVectorBlocked2ndMoment;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim*(dim+1)*(dim+2)*(dim+3)/24> DenseVectorBlocked4thMoment;

    // define our system level
    typedef ModNavierStokesBlockedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    // get our global solver system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalSystemTransfer GlobalSystemTransfer;
  public:
    TimeStamp stamp_init_time;
    double init_time;
    Cubature::DynamicFactory cubature;
    const String InflowFacetName;
    const String OutflowFacetName;
    const Dist::Comm& _comm;
    Index time_step;
    //statistics for Solver data
    BenchmarkStats statistics;
    CCNDFiberSolver<DomainLevel_> _inner_solver;
    GlobalSystemVector _vec_sol;
    GlobalSystemVector _vec_rhs;
    GlobalSystemVector _prev_sol1;
    GlobalSystemVector _prev_sol2;
    DataType delta_t;
    //DataType steady_tol;
    Index t_expo;
    StopWatch watch_total_run;
    StopWatch watch_vtk_write;
    bool euler;
    bool verbose;

    public:
    template<typename DirichletFunctionType, typename RHSFunctionType, typename NeumannFunctionType, typename NoSlipFunctionType>
    CCNDUnsteadyInterface(const Dist::Comm& comm, SimpleArgParser& args, std::unique_ptr<Control::Domain::PartiDomainControl<DomainLevel_>>&& domain_ptr,
                        DirichletFunctionType& diri_func, RHSFunctionType& rhs_func, NeumannFunctionType& neumann_func, NoSlipFunctionType& no_slip_func, String cubature_name = "auto-degree:6",
                        String inflow_name = "bnd:io:big", String outflow_name = "bnd:io:small")
    : cubature(cubature_name),
    InflowFacetName(inflow_name),
    OutflowFacetName(outflow_name),
    _comm(comm),
    _inner_solver(comm, args, std::move(domain_ptr), diri_func, no_slip_func, InflowFacetName, OutflowFacetName, cubature_name, statistics)
    {
      init_time = stamp_init_time.elapsed_now();
      watch_total_run.start();
      // fetch our finest levels
      DomainLevelType& the_domain_level = *(_inner_solver._domain_ptr->front());
      SystemLevelType& the_system_level = *(_inner_solver.system_levels.front());
      SpaceVeloType& velo_space = the_domain_level.space_velo;
//       SpacePresType& pres_space = the_domain_level.space_pres;

      time_step = Index(1);
      //parse time-stepping parameters
      delta_t = parse(args, "delta-t", DataType(0.01));
      //steady_tol = parse(args, "steady-tol", DataType(1E-3));
      t_expo = parse(args, "t-expo", Index(2));
      euler = (args.check("euler") >= 0);
      verbose = (args.check("verbose") >= 0);
//       euler = true;
      //init _vec_sol and _vec_rhs
      _vec_sol = the_system_level.matrix_sys.create_vector_r();
      _vec_rhs = the_system_level.matrix_sys.create_vector_r();
      auto vec_temp = the_system_level.matrix_sys.create_vector_r();

      _vec_sol.format();
      _vec_rhs.format();
      vec_temp.format();

      //create a linear functional assembler
      Assembly::Common::ForceFunctional<RHSFunctionType> functional(rhs_func);
      //write this into vec_rhs for now...
      Assembly::LinearFunctionalAssembler::assemble_vector(_vec_rhs.local().template at<0>(), functional, velo_space, cubature);

      //create vector for neumann boundary condition on right side
      //first we get the facet
      auto* mesh_part_bnd_r = the_domain_level.get_mesh_node()->find_mesh_part(OutflowFacetName);

      //and now use the trace assembler to assemble the additional part for our rhs
      Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> neumann_boundary_asm(velo_space.get_trafo());
      if(mesh_part_bnd_r != nullptr)
        neumann_boundary_asm.add_mesh_part(*mesh_part_bnd_r);
      neumann_boundary_asm.compile();


      //create vector of the right size
      //we will just use vec_def for this

      //create the associated lin func assembler
      Assembly::Common::ForceFunctional<NeumannFunctionType> neumann_functional(neumann_func);
      //and now we should be able to compile our vector
      neumann_boundary_asm.assemble_functional_vector(_vec_sol.local().template at<0>(), neumann_functional, velo_space, cubature);

      //and now add this vector onto cor
      vec_temp.axpy(_vec_sol, _vec_rhs, DataType(1.));
      //multiply with rho, save this into rhs
      _vec_rhs.scale(vec_temp, _inner_solver.rho);
      //set help vectors to zero
      _vec_sol.format();


      //sync the rhs vector! no? should make no difference because we format in new time step...
      _vec_rhs.sync_0();



      //filter
      the_system_level.filter_sys.filter_sol(_vec_sol);
      the_system_level.filter_sys.filter_rhs(_vec_rhs);

      //clone our previous solution vectors
      _prev_sol1.clone(_vec_sol);
      _prev_sol2.clone(_vec_sol);

      watch_total_run.stop();

    }

    template<typename DirichletFunctionType, typename NoSlipFunctionType >
    void update_filters(DirichletFunctionType& diri_func, NoSlipFunctionType& noslip_func)
    {
      watch_total_run.start();
      _inner_solver.update_filters(diri_func, noslip_func, InflowFacetName, OutflowFacetName);
      //update our solution
//       _inner_solver.system_levels.front()->filter_sys.filter_sol(_vec_sol);  //this should not be done before copying into previous timesteps...
      watch_total_run.stop();
    }

    template<typename RHSFunctionType, typename NeumannFunctionType>
    void update_rhs(RHSFunctionType& rhs_func, NeumannFunctionType& neumann_func)
    {
      watch_total_run.start();

      DomainLevelType& the_domain_level = *(_inner_solver._domain_ptr->front());
      SystemLevelType& the_system_level = *(_inner_solver.system_levels.front());
      SpaceVeloType& velo_space = the_domain_level.space_velo;

      auto vec_temp = the_system_level.matrix_sys.create_vector_r();
      vec_temp.format();
      _vec_rhs.format();

      //create a linear functional assembler
      Assembly::Common::ForceFunctional<RHSFunctionType> functional(rhs_func);
      //write this into vec_rhs for now...
      Assembly::LinearFunctionalAssembler::assemble_vector(_vec_rhs.local().template at<0>(), functional, velo_space, cubature);

      //create vector for neumann boundary condition on right side
      //first we get the facet
      auto* mesh_part_bnd_r = the_domain_level.get_mesh_node()->find_mesh_part(OutflowFacetName);

      //and now use the trace assembler to assemble the additional part for our rhs
      Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> neumann_boundary_asm(velo_space.get_trafo());
      if(mesh_part_bnd_r != nullptr)
        neumann_boundary_asm.add_mesh_part(*mesh_part_bnd_r);
      neumann_boundary_asm.compile();

      //create vector of the right size
      //we will just use vec_def for this

      //create the associated lin func assembler
      Assembly::Common::ForceFunctional<NeumannFunctionType> neumann_functional(neumann_func);
      //and now we should be able to compile our vector
      neumann_boundary_asm.assemble_functional_vector(vec_temp.local().template at<0>(), neumann_functional, velo_space, cubature);

      //and now add this vector onto cor
      vec_temp.axpy(vec_temp, _vec_rhs, DataType(1.));
      //multiply with rho, save this into rhs
      _vec_rhs.scale(vec_temp, _inner_solver.rho);

      //sync the rhs vector!
//       _vec_rhs.sync_0();

      //filter
//       the_system_level.filter_sys.filter_rhs(_vec_rhs);
      watch_total_run.stop();
    }



    /**
     * \brief This solves the standard steady incompressible Navier-Stokes equations, i.e. N_s = N_p = 0, and writes the solution as projection into the OrientSpace.
     *
     * \warning This only works as inteded if called directly after initialization, so without solve_time_step calls before.
     *
     *
     * \tparam SpaceOrientationType The type of the space the solution should be projected into.
     *
     * \param[out] vec_sol_orient A velocity vector matching orient_space in which the velocity solution is to be projected.
     *
     * \param[in] orient_space The space in which the velocity solution should be represented.
     *
     *
     * \author Maximilian Esser
     */

    template<typename SpaceOrientationType>
    Solver::Status solve_basic_navier(LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim>& vec_sol_orient, SpaceOrientationType& orient_space)
    {
      XASSERTM(vec_sol_orient.size() == orient_space.get_num_dofs(), "ERROR: Space Orientation type does not fit orient vector");
      watch_total_run.start();

      if(verbose)
        _comm.print("Solving Basic Navier-Stokes system");
      //now simply call the inner solver function
      Solver::Status status = _inner_solver.solve_basic_navier(_vec_sol, _vec_rhs);

      vec_sol_orient.format();
      //and now project the velocity solution into vec_sol_orient
      Assembly::FEInterpolator<SpaceOrientationType, SpaceVeloType>::interpolate(vec_sol_orient, _vec_sol.local().template at<0>(), orient_space, _inner_solver._domain_ptr->front()->space_velo);
      watch_total_run.stop();

      return status;
    }

    /**
     * \brief This solves one timestep of the extended incompressible Navier-Stokes equations and writes the solution as projection into the OrientSpace.
     *
     * \warning Do not change the inner _vec_sol vectors and the time_step variable between different calls of solve_time_step.
     *
     * \warning If RHS or boundary conditions are time dependent, you have to call the appropriate update function before calling this function.
     *
     *
     * \tparam SpaceOrientationType The type of the space the solution should be projected into.
     *
     * \param[out] vec_sol_orient A velocity vector matching orient_space in which the velocity solution is to be projected.
     *
     * \param[in] orient_space The space in which the velocity solution should be represented.
     *
     * \param[in] second_moment A blocked vector representing the second moment orientation tensor as FE function in orient_space.
     *
     * \param[in] fourth_moment A blocked vector representing the fourth moment orientation tensor as FE function in orient_space.
     *
     * \note Due to higher symmetry in the orientation tensors, the entries are, depending on the dimension, encoded in a special way desribed below:
     *
     *            2D                                                                         3D
     *                                                second_moment Aij -> w(k)
     *
     *        w(0) = A11                                        |                       w(0) = A11
     *        w(1) = A22                                        |                       w(1) = A22
     *        w(2) = A12                                        |                       w(2) = A33
     *                                                          |                       w(3) = A12
     *                                                          |                       w(4) = A13
     *                                                          |                       w(5) = A23
     *                                                          |
     *                                                          |
     *
     *                                                fourth_moment Aijmn -> t(k)
     *                                                          |
     *                                                          |
     *        t(0) = A1111                                      |                       t(0)  = A1111
     *        t(1) = A2222                                      |                       t(1)  = A2222
     *        t(2) = A1112                                      |                       t(2)  = A3333
     *        t(3) = A2221                                      |                       t(3)  = A1112
     *        t(4) = A1122                                      |                       t(4)  = A1113
     *                                                          |                       t(5)  = A2221
     *                                                          |                       t(6)  = A2223
     *                                                          |                       t(7)  = A3331
     *                                                          |                       t(8)  = A3332
     *                                                          |                       t(9)  = A1122
     *                                                          |                       t(10) = A1133
     *                                                          |                       t(11) = A2233
     *                                                          |                       t(12) = A1123
     *                                                          |                       t(13) = A2213
     *                                                          |                       t(14) = A3312
     *
     *
     *
     *
     * \author Maximilian Esser
     */
    template<typename SpaceOrientationType>
    Solver::Status solve_time_step(LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim>& vec_sol_orient, const SpaceOrientationType& orient_space, const DenseVectorBlocked2ndMoment& second_moment, const DenseVectorBlocked4thMoment& fourth_moment)
    {
      XASSERTM(vec_sol_orient.size() == orient_space.get_num_dofs(), "ERROR: Space Orientation type does not fit orient vector");
      watch_total_run.start();

      if(verbose)
        _comm.print("Solving time step " + stringify(time_step) + " of the extended Navier-Stokes equations");

      //get ref to the system level
      SystemLevelType& the_system_level = *(_inner_solver.system_levels.front());

      // extrapolate initial u[k] from u[k-1] and u[k-2] for this time-step
      //we have to prepare a new rhs vector according to the chosen time-stepping scheme, which will be done later, but we can use this vector for extrapolation
      auto vec_rhs_time_stepping = _vec_rhs.clone(LAFEM::CloneMode::Deep);

      if(time_step > Index(1))
      {
        auto& _prev_sol3 = vec_rhs_time_stepping; // abuse RHS vector, which is formatted below anyways

        // shift solution vectors
        _prev_sol3.copy(_prev_sol2); // u_{k-3}
        _prev_sol2.copy(_prev_sol1); // u_{k-2}
        _prev_sol1.copy(_vec_sol);   // u_{k-1}

        // extrapolate to current time step
        if((t_expo >= Index(2)) && (time_step > Index(2)))
        {
          // perform quadratic extrapolation of solution to current timestep
          // u_{k} := 3*u_{k-1} - 3*u_{k-2}) + u_{k-3}
          _vec_sol.axpy(_prev_sol1, _prev_sol3, DataType(3));
          _vec_sol.axpy(_prev_sol2, _vec_sol, -DataType(3));
        }
        else if((t_expo >= Index(1)) && (time_step > Index(1)))
        {
          // perform linear extrapolation of solution to current timestep
          // u_{k} := 2*u_{k-1} - u_{k-2}
          _vec_sol.scale(_prev_sol1, DataType(2));
          _vec_sol.axpy(_prev_sol2, _vec_sol, -DataType(1));
        }
        // else-case is constant extrapolation, which is done already by copy
      }

      vec_rhs_time_stepping.format();
      //apply solution filter, before calling this solving method, the filters should have been updated...
      the_system_level.filter_sys.filter_sol(_vec_sol);




      //first of all, is this the first time_step?
      if(euler || time_step == Index(1))
      {
        //we use an implicit Euler in any case
        _inner_solver.alpha = DataType(1) / delta_t;
        // f_k := 1/dt * M * u_{k-1}
        the_system_level.local_velo_mass_matrix.apply(
          vec_rhs_time_stepping.local().template at<0>(),
          _prev_sol1.local().template at<0>());
        vec_rhs_time_stepping.axpy(vec_rhs_time_stepping, _vec_rhs, DataType(1) / delta_t);
      }
      else
      {
        // we're beyond the first time step ==> BDF(2)
        // First, adjust the mass matrix parameter in the burgers assemblers
        _inner_solver.alpha = DataType(1.5) / delta_t;

        auto& vec_temp = _inner_solver.vec_def;
        // f_k := 3/(2*dt) * (4/3 * M * u_{k-1} - 1/3 M * u_{k-2}
        //      = -1/(2*dt) * M * (u_{k-2} - 4*u_{k-1})
        vec_temp.axpy(_prev_sol1, _prev_sol2, -DataType(4));
        the_system_level.local_velo_mass_matrix.apply(
          vec_rhs_time_stepping.local().template at<0>(),
          vec_temp.local().template at<0>());
        vec_rhs_time_stepping.axpy(vec_rhs_time_stepping, _vec_rhs, -DataType(0.5) / delta_t);
      }

      //call the inner solver function
      Solver::Status status = _inner_solver.solve_navier(_vec_sol, vec_rhs_time_stepping, orient_space, second_moment, fourth_moment, verbose);

      //project into our output solution
      vec_sol_orient.format();
      Assembly::FEInterpolator<SpaceOrientationType, SpaceVeloType>::interpolate(vec_sol_orient, _vec_sol.local().template at<0>(), orient_space, _inner_solver._domain_ptr->front()->space_velo);


      ++time_step;

      watch_total_run.stop();

      return status;
    }

    //automatically writes vtk into filename with intern time_step number
    void write_vtk(String filename)
    {
      DomainLevelType& the_domain_level = *(_inner_solver._domain_ptr->front());
        watch_vtk_write.start();
        String now_vtk_name = filename + "." + stringify(time_step).pad_front(5, '0');
        String line = " > Writing VTK file: '" + now_vtk_name + "'";
        if(verbose)
          _comm.print(line);
        {
          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

          // project velocity and pressure
          exporter.add_vertex_vector("velocity", _vec_sol.local().template at<0>());
//           exporter.add_vertex_vector("velo_der", vec_def.local().template at<0>());

          // project pressure
          Cubature::DynamicFactory cub("gauss-legendre:2");
          LAFEM::DenseVector<Mem::Main, DataType, Index> vtx_p, der_p;
          Assembly::DiscreteCellProjector::project(vtx_p, _vec_sol.local().template at<1>(), the_domain_level.space_pres, cub);
//           Assembly::DiscreteCellProjector::project(der_p, vec_def.local().template at<1>(), the_domain_level.space_pres, cub);

          // write pressure
          exporter.add_cell_scalar("pressure", vtx_p.elements());
//           exporter.add_cell_scalar("pres_der", der_p.elements());

          // finally, write the VTK file
          exporter.write(now_vtk_name, _comm);
        }
        watch_vtk_write.stop();
    }


    void compile_statistics()
    {
      _inner_solver.compile_statistics();
      statistics.times[Times::total_run] = watch_total_run.elapsed() + init_time;
      {
        MemoryUsage mi;
        statistics.bytes[Bytes::peak_p] = mi.get_peak_physical();
        statistics.bytes[Bytes::peak_v] = mi.get_peak_virtual();
      }

      statistics.sync(_comm);

      _comm.print(String("\n") + String(80, '=') + "\n");
      //     comm.print(summary.format());
      _comm.print(statistics.format());

      // print multigrid timings
      if(_comm.rank() == 0)
      {
        _comm.print("Multigrid Timings:");
        _comm.print("              Defect /   Smoother /   Transfer /     Coarse");
        _comm.print("Overall : " +
        stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_defect(), 3, 10) + " / " +
        stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_smooth(), 3, 10) + " / " +
        stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_transfer(), 3, 10) + " / " +
        stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_coarse(), 3, 10));
        for(int i(0); i < int(_inner_solver._multigrid_hierarchy->size_physical()); ++i)
        {
          _comm.print("Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_defect(i), 3, 10) + " / " +
          stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_smooth(i), 3, 10) + " / " +
          stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_transfer(i), 3, 10) + " / " +
          stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_coarse(i), 3, 10));
        }
      }

    }

    inline void set_ns(DataType n_s)
    {
      _inner_solver.n_s = n_s;
    }

    inline void set_np(DataType n_p)
    {
      _inner_solver.n_p = n_p;
    }





  }; //class CCNDUnsteadyInterface

  template<typename DomainLevel_>
  class CCNDSteadyInterface
  {
  public:
    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainLevelType::SpaceVeloType SpaceVeloType;
    typedef typename DomainLevelType::SpacePresType SpacePresType;
//     typedef Space::Lagrange1::Element<typename DomainLevelType::TrafoType> SpaceOrientationType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim*(dim+1)/2> DenseVectorBlocked2ndMoment;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim*(dim+1)*(dim+2)*(dim+3)/24> DenseVectorBlocked4thMoment;

    // define our system level
    typedef ModNavierStokesBlockedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    // get our global solver system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalSystemTransfer GlobalSystemTransfer;
  public:
    TimeStamp stamp_init_time;
    double init_time;
    Cubature::DynamicFactory cubature;
    const String InflowFacetName;
    const String OutflowFacetName;
    const Dist::Comm& _comm;
    BenchmarkStats statistics;
    CCNDFiberSolver<DomainLevel_> _inner_solver;
    GlobalSystemVector _vec_sol;
    GlobalSystemVector _vec_rhs;
    //statistics for Solver data


    StopWatch watch_total_run;




  public:
    template<typename DirichletFunctionType, typename RHSFunctionType, typename NeumannFunctionType, typename NoSlipFunctionType>
    CCNDSteadyInterface(const Dist::Comm& comm, SimpleArgParser& args, std::unique_ptr<Control::Domain::PartiDomainControl<DomainLevel_>>&& domain_ptr,
                        DirichletFunctionType& diri_func, RHSFunctionType& rhs_func, NeumannFunctionType& neumann_func, NoSlipFunctionType& no_slip_func, String cubature_name = "auto-degree:6",
                        String inflow_name = "bnd:io:big", String outflow_name = "bnd:io:small")
    : cubature(cubature_name),
    InflowFacetName(inflow_name),
    OutflowFacetName(outflow_name),
    _comm(comm),
    _inner_solver(comm, args, std::move(domain_ptr), diri_func, no_slip_func, InflowFacetName, OutflowFacetName, cubature_name, statistics)
    {
      init_time = stamp_init_time.elapsed_now();
      watch_total_run.start();
      // fetch our finest levels
      DomainLevelType& the_domain_level = *(_inner_solver._domain_ptr->front());
      SystemLevelType& the_system_level = *(_inner_solver.system_levels.front());
      SpaceVeloType& velo_space = the_domain_level.space_velo;
//       SpacePresType& pres_space = the_domain_level.space_pres;

      //init _vec_sol and _vec_rhs
      _vec_sol = the_system_level.matrix_sys.create_vector_r();
      _vec_rhs = the_system_level.matrix_sys.create_vector_r();
      auto vec_temp = the_system_level.matrix_sys.create_vector_r();

      _vec_sol.format();
      _vec_rhs.format();
      vec_temp.format();

      //create a linear functional assembler
      Assembly::Common::ForceFunctional<RHSFunctionType> functional(rhs_func);
      //write this into vec_rhs for now...
      Assembly::LinearFunctionalAssembler::assemble_vector(_vec_rhs.local().template at<0>(), functional, velo_space, cubature);

      //create vector for neumann boundary condition on right side
      //first we get the facet
      auto* mesh_part_bnd_r = the_domain_level.get_mesh_node()->find_mesh_part(OutflowFacetName);

      //and now use the trace assembler to assemble the additional part for our rhs
      Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> neumann_boundary_asm(velo_space.get_trafo());
      if(mesh_part_bnd_r != nullptr)
        neumann_boundary_asm.add_mesh_part(*mesh_part_bnd_r);
      neumann_boundary_asm.compile();

      //create vector of the right size
      //we will just use vec_def for this

      //create the associated lin func assembler
      Assembly::Common::ForceFunctional<NeumannFunctionType> neumann_functional(neumann_func);
      //and now we should be able to compile our vector
      neumann_boundary_asm.assemble_functional_vector(_vec_sol.local().template at<0>(), neumann_functional, velo_space, cubature);

      //and now add this vector onto cor
      vec_temp.axpy(_vec_sol, _vec_rhs, DataType(1.));
      //multiply with rho, save this into rhs
      _vec_rhs.scale(vec_temp, _inner_solver.rho);
      //set help vectors to zero
      _vec_sol.format();


      //sync the rhs vector!
      _vec_rhs.sync_0();



      //filter
      the_system_level.filter_sys.filter_sol(_vec_sol);
      the_system_level.filter_sys.filter_rhs(_vec_rhs);
      watch_total_run.stop();



    }

    template<typename SpaceOrientationType>
    Solver::Status solve_basic_navier(LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim>& vec_sol_orient, SpaceOrientationType& orient_space)
    {
      XASSERTM(vec_sol_orient.size() == orient_space.get_num_dofs(), "ERROR: Space Orientation type does not fit orient vector");
      watch_total_run.start();

      _comm.print("Solving Basic Navier-Stokes system");
      //now simply call the inner solver function
      Solver::Status status = _inner_solver.solve_basic_navier(_vec_sol, _vec_rhs);

      vec_sol_orient.format();
      //and now project the velocity solution into vec_sol_orient
      Assembly::FEInterpolator<SpaceOrientationType, SpaceVeloType>::interpolate(vec_sol_orient, _vec_sol.local().template at<0>(), orient_space, _inner_solver._domain_ptr->front()->space_velo);
      watch_total_run.stop();

      return status;
    }

    template<typename SpaceOrientationType>
    Solver::Status solve_navier(LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim>& vec_sol_orient, const SpaceOrientationType& orient_space, const DenseVectorBlocked2ndMoment& second_moment, const DenseVectorBlocked4thMoment& fourth_moment)
    {
      XASSERTM(vec_sol_orient.size() == orient_space.get_num_dofs(), "ERROR: Space Orientation type does not fit orient vector");
      watch_total_run.start();

      _comm.print("Solving extended Navier-Stokes");
      //call the inner solver function
      Solver::Status status = _inner_solver.solve_navier(_vec_sol, _vec_rhs, orient_space, second_moment, fourth_moment);

      //project into our output solution
      vec_sol_orient.format();
      Assembly::FEInterpolator<SpaceOrientationType, SpaceVeloType>::interpolate(vec_sol_orient, _vec_sol.local().template at<0>(), orient_space, _inner_solver._domain_ptr->front()->space_velo);
      watch_total_run.stop();

      return status;
    }

    void compile_statistics()
    {
      _inner_solver.compile_statistics();
      statistics.times[Times::total_run] = watch_total_run.elapsed() + init_time;
      {
        MemoryUsage mi;
        statistics.bytes[Bytes::peak_p] = mi.get_peak_physical();
        statistics.bytes[Bytes::peak_v] = mi.get_peak_virtual();
      }

      statistics.sync(_comm);

      _comm.print(String("\n") + String(80, '=') + "\n");
      //     comm.print(summary.format());
      _comm.print(statistics.format());

      // print multigrid timings
      if(_comm.rank() == 0)
      {
        _comm.print("Multigrid Timings:");
        _comm.print("              Defect /   Smoother /   Transfer /     Coarse");
        _comm.print("Overall : " +
        stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_defect(), 3, 10) + " / " +
        stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_smooth(), 3, 10) + " / " +
        stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_transfer(), 3, 10) + " / " +
        stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_coarse(), 3, 10));
        for(int i(0); i < int(_inner_solver._multigrid_hierarchy->size_physical()); ++i)
        {
          _comm.print("Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_defect(i), 3, 10) + " / " +
          stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_smooth(i), 3, 10) + " / " +
          stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_transfer(i), 3, 10) + " / " +
          stringify_fp_fix(_inner_solver._multigrid_hierarchy->get_time_coarse(i), 3, 10));
        }
      }

    }

  };

  //this class provides an interface for a single step of a extended ccnd solver for a given PartiDomainControl which has to be created beforehand (see helper function)
  //the system solved itself depends on the outer ControlClass
  template<typename DomainLevel_>
  class CCNDFiberSolver
  {
  public:
    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainLevelType::SpaceVeloType SpaceVeloType;
    typedef typename DomainLevelType::SpacePresType SpacePresType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> LocalVeloVectorType;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim*(dim+1)/2> DenseVectorBlocked2ndMoment;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim*(dim+1)*(dim+2)*(dim+3)/24> DenseVectorBlocked4thMoment;

    // define our system level
    typedef ModNavierStokesBlockedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    // get our global solver system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalSystemTransfer GlobalSystemTransfer;
    typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;

    bool navier = true; //(args.check("stokes") < 0);
    bool newton = true; //(args.check("picard") < 0);
    bool adapt_tol = true; //(args.check("mg-tol-rel") < 0);
    bool testmode = false; //(args.check("test-mode") >= 0);
    bool plot_mg_iter = false; //(args.check("plot-mg-iter") >= 0);

    // diffusion parameter
    DataType mu = DataType(1e-1);
    //convetion speed
    DataType rho = DataType(1.);
    // maximum velocity, default: 2D: 0.3, 3D: 0.45
    DataType v_max = DataType(0.15) * DataType(dim);
    // max. nonlinear solver iterations
    Index max_nl_iter = Index(20);
    // min. multigrid iterations
    Index min_mg_iter = Index(1);
    // max. multigrid iterations
    Index max_mg_iter = Index(25);
    // number of smoothing steps
    Index smooth_steps = Index(8);
    // damping parameter for smoother
    DataType smooth_damp = DataType(0.7);
    // rel. tolerance for linear solver
    DataType mg_tol_rel = DataType(1e-2);
    // absolute tolerance for nonlinear solver
    DataType nl_tol_abs = DataType(1e-8);
    // shear-rate number
    DataType n_s = DataType(0.);
    // particle number
    DataType n_p = DataType(0.);
    // reaction parameter
    DataType alpha = DataType(0.);

    bool umf_cgs_hold = false;
  public:
    const Dist::Comm& _comm;
    std::unique_ptr<Control::Domain::PartiDomainControl<DomainLevel_>> _domain_ptr;
    GlobalSystemVector vec_def;
    GlobalSystemVector vec_cor;

    std::shared_ptr<Solver::IterativeSolver<GlobalSystemVector>> _solver; //can/should be unqiue??!
    std::shared_ptr<Solver::MultiGridHierarchy<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>> _multigrid_hierarchy;

    // array of Vanka pointers - this is only required to collect the memory usage
    // statistics of the Vankas, as we need the memory usage after factorization
    std::deque<
      std::shared_ptr<
        Solver::AmaVanka<
          typename SystemLevelType::LocalSystemMatrix,
          typename SystemLevelType::LocalSystemFilter>>> ama_vankas;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    // cubature for assembly ( degree can be reduced ... but careful, as orientation tensor A does not have to be polynomial)
    Cubature::DynamicFactory cubature;

    //we will directly create a deque for the all orient_mats, which will be constructed due to truncation
    std::deque<DenseVectorBlocked2ndMoment> orient_2_levels;
    std::deque<DenseVectorBlocked4thMoment> orient_4_levels;

    //Reference to the Statistic of the Interface
    BenchmarkStats& statistics;
    //we need a few clocks
    StopWatch watch_tensor_assembly;
    TimeStamp stamp_ass;
    TimeStamp basic_asm;
    // create a few watches
    StopWatch watch_stokes_solve;
    StopWatch watch_nonlin_loop;
    StopWatch watch_nonlin_def_asm;
    StopWatch watch_nonlin_mat_asm;
    StopWatch watch_nonlin_solver_init;
    StopWatch watch_nonlin_solver_apply;




  public:
    template<typename DirichletFunctionType, typename NoSlipFunctionType>
    CCNDFiberSolver(const Dist::Comm& comm, SimpleArgParser& args, std::unique_ptr<Control::Domain::PartiDomainControl<DomainLevel_>>&& domain_ptr, DirichletFunctionType& diri_func, NoSlipFunctionType& no_slip_func, const String& InflowFacetName, const String& OutflowFacetName, const String& cubature_string, BenchmarkStats& interface_statistic)
    : _comm(comm),
     _domain_ptr(std::move(domain_ptr)),
     cubature(cubature_string),
     statistics(interface_statistic)
    {
      navier = (args.check("stokes") < 0);
      newton = (args.check("picard") < 0);
      adapt_tol = (args.check("mg-tol-rel") < 0);
      testmode = (args.check("test-mode") >= 0);
      plot_mg_iter = (args.check("plot-mg-iter") >= 0);

      // diffusion parameter
      mu = parse(args, "mu", DataType(1e-1));
      //convetion speed
      rho = parse(args, "rho", DataType(1.));
      // maximum velocity, default: 2D: 0.3, 3D: 0.45
      v_max = parse(args, "v-max", DataType(dim) * DataType(0.15));
      // max. nonlinear solver iterations
      max_nl_iter = parse(args, "max-nl-iter", Index(20));
      // min. multigrid iterations
      min_mg_iter = parse(args, "min-mg-iter", Index(1));
      // max. multigrid iterations
      max_mg_iter = parse(args, "max-mg-iter", Index(25));
      // number of smoothing steps
      smooth_steps = parse(args, "smooth-steps", Index(8));
      // damping parameter for smoother
      smooth_damp = parse(args, "smooth-damp", DataType(0.7));
      // rel. tolerance for linear solver
      mg_tol_rel = parse(args, "mg-tol-rel", DataType(1E-2));
      // absolute tolerance for nonlinear solver
      nl_tol_abs = parse(args, "nl-tol-abs", DataType(1E-8));
      //shear-rate Parameter
      n_s = parse(args, "n-s", DataType(0.));
      //particle number
      n_p = parse(args, "n-p", DataType(1.));
      alpha = parse(args, "alpha", DataType(0.));



      //get reference to domain
      Control::Domain::PartiDomainControl<DomainLevel_>& domain = *_domain_ptr;

      #ifdef FEAT_HAVE_UMFPACK
      const bool umf_cgs = (domain.back_layer().comm().size() == 1) && (args.check("no-umfpack") < 0);
      #else
      const bool umf_cgs = false;
      #endif

      umf_cgs_hold = umf_cgs;

      const Index num_levels = Index(domain.size_physical());

      // create system levels
      for (Index i(0); i < num_levels; ++i)
      {
        system_levels.push_back(std::make_shared<SystemLevelType>());
      }

      /* --------------------------------------------------------------------------------------------------------------------------------------------------------
      * Symbolic Assembly
      * ----------------------------------------------------------------------------------------------------------------------------------------------------- */

      stamp_ass.stamp();
      // assemble gates, muxers and transfers
      for (Index i(0); i < num_levels; ++i)
      {
        system_levels.at(i)->assemble_gates(domain.at(i));

        if((i+1) < domain.size_virtual())
        {
          system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
          system_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature);
        }
      }


      // assemble velocity truncation operators -- we need those for the assembly of the
      // non-linear burgers operators on the coarser levels
      for (Index i(0); i < num_levels; ++i)
      {
        if(i+1 < num_levels)
          system_levels.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature, system_levels.at(i+1).get());
        else if(i+1 < domain.size_virtual())
          system_levels.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature);
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

      basic_asm.stamp();
      // assemble basic matrices
      for(Index i(0); i < num_levels; ++i)
      {
        system_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
        system_levels.at(i)->assemble_pres_struct(domain.at(i)->space_pres);
        system_levels.at(i)->assemble_velocity_laplace_matrix(domain.at(i)->space_velo, cubature, mu);
        system_levels.at(i)->assemble_velocity_mass_matrix(domain.at(i)->space_velo, cubature);
        system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
        system_levels.at(i)->compile_system_matrix();
      }
      statistics.times[Times::basic_matrix_asm] = basic_asm.elapsed_now();

      // the names of the mesh parts on which to assemble
      std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

      for(Index i(0); i < num_levels; ++i)
      {
        // get our local velocity filter
        auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

        // create unit-filter assembler
        Assembly::UnitFilterAssembler<MeshType> unit_asm_boundary, unit_asm_noflow;

        // loop over all boundary parts
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
            if(name == InflowFacetName)
            {
              // inflow
              unit_asm_boundary.add_mesh_part(*mesh_part);
            }
            else if(name != OutflowFacetName)
              unit_asm_noflow.add_mesh_part(*mesh_part);
          }
        }

        // assemble the filter
        unit_asm_boundary.assemble(fil_loc_v, domain.at(i)->space_velo, diri_func);
        unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo, no_slip_func);

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

      // fetch our finest levels
//       DomainLevelType& the_domain_level = *domain.front();
      SystemLevelType& the_system_level = *system_levels.front();
//       SpaceVeloType& velo_space = the_domain_level.space_velo;
//       SpacePresType& pres_space = the_domain_level.space_pres;

      // get our global solve matrix and filter
      GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
      GlobalSystemFilter& filter = the_system_level.filter_sys;

      // create new vectors
      vec_def = the_system_level.matrix_sys.create_vector_r();
      vec_cor = the_system_level.matrix_sys.create_vector_r();

      // format the vectors
      vec_def.format();
      vec_cor.format();

      {
        // count non-zeros in a and b
        statistics.counts[Counts::nnze_a] = the_system_level.matrix_sys.local().block_a().used_elements();
        statistics.counts[Counts::nnze_b] = the_system_level.matrix_sys.local().block_b().used_elements();
        statistics.counts[Counts::nnze_total] = the_system_level.matrix_sys.local().template used_elements<LAFEM::Perspective::pod>();
      }


      /* ***************************************************************************************** */
      /* ***************************************************************************************** */
      /* ***************************************************************************************** */

      // create a multigrid solver
      _multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>>(domain.size_virtual());

      //     // array of Vanka pointers - this is only required to collect the memory usage
      //     // statistics of the Vankas, as we need the memory usage after factorization
      //     std::deque<
      //       std::shared_ptr<
      //         Solver::AmaVanka<
      //           typename SystemLevelType::LocalSystemMatrix,
      //           typename SystemLevelType::LocalSystemFilter>>> ama_vankas;

      // push levels into multigrid
      for(std::size_t i(0); i < system_levels.size(); ++i)
      {
        SystemLevelType& lvl = *system_levels.at(i);

        if((i+1) < domain.size_virtual())
        {
          auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
          ama_vankas.push_back(vanka);
          auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
//           auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);
          //instead use a bicg solver
          auto smoother = Solver::new_bicgstab(lvl.matrix_sys, lvl.filter_sys, schwarz);
//           smoother->set_min_iter(smooth_steps);
          smoother->set_max_iter(smooth_steps);
//           smoother->skip_defect_calc(true); // skip defect calculation
          _multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
        }
        #ifdef FEAT_HAVE_UMFPACK
        else if(umf_cgs)
        {
          // create UMFPACK coarse grid solver
          auto umfpack = Solver::new_generic_umfpack(lvl.local_matrix_sys);
          auto cgsolver = Solver::new_schwarz_precond(umfpack, lvl.filter_sys);
          _multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
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

          _multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
        }
      }

      // create our multigrid solver
      auto multigrid = Solver::new_multigrid(_multigrid_hierarchy, Solver::MultiGridCycle::V);

      // create our solver
      _solver = Solver::new_richardson(matrix, filter, 1.0, multigrid);

      _solver->set_plot_name("Multigrid");

      _solver->set_min_iter(min_mg_iter);
      _solver->set_max_iter(max_mg_iter);
      _solver->set_tol_rel(mg_tol_rel);
      _solver->set_min_stag_iter(Index(3));

      //init the solver
      watch_nonlin_solver_init.start();
      _multigrid_hierarchy->init_symbolic();
      _solver->init_symbolic();
      watch_nonlin_solver_init.stop();

    }

    ~CCNDFiberSolver()
    {
      _solver->done_symbolic();
      _multigrid_hierarchy->done_symbolic();
    }

    void print_parameters() const
    {
      static constexpr std::size_t pl = 30u;
      static constexpr char pc = '.';
      _comm.print("\nProblem Parameters:");
      _comm.print(String("Mu").pad_back(pl, pc) + ": " + stringify(mu));
      _comm.print(String("Rho").pad_back(pl, pc) + ": " + stringify(rho));
      _comm.print(String("N_s").pad_back(pl, pc) + ": " + stringify(n_s));
      _comm.print(String("N_p").pad_back(pl, pc) + ": " + stringify(n_p));
      _comm.print(String("V-Max").pad_back(pl, pc) + ": " + stringify(v_max));
      _comm.print(String("System").pad_back(pl, pc) + ": " + (navier ? "Navier-Stokes" : "Stokes"));
      _comm.print(String("Tensor").pad_back(pl, pc) + ": " + ("Deformation"));
      _comm.print(String("Nonlinear Solver").pad_back(pl, pc) + ": " + (newton ? "Newton" : "Picard"));
      _comm.print(String("AmaVanka Smoother Steps").pad_back(pl, pc) + ": " + stringify(smooth_steps));
      _comm.print(String("AmaVanka Smoother Damping").pad_back(pl, pc) + ": " + stringify(smooth_damp));
      _comm.print(String("Nonlinear Absolute Tol").pad_back(pl, pc) + ": " + stringify_fp_sci(nl_tol_abs));
      _comm.print(String("Multigrid Relative Tol").pad_back(pl, pc) + ": " + (adapt_tol ? String("adaptive") : stringify_fp_sci(mg_tol_rel)));
      _comm.print(String("Min Multigrid Iterations").pad_back(pl, pc) + ": " + stringify(min_mg_iter));
      _comm.print(String("Max Multigrid Iterations").pad_back(pl, pc) + ": " + stringify(max_mg_iter));
      _comm.print(String("Max Nonlinear Iterations").pad_back(pl, pc) + ": " + stringify(max_nl_iter));
      _comm.print(String("Alpha").pad_back(pl, pc) + ": " + stringify(alpha));
      if(umf_cgs_hold)
        _comm.print(String("Coarse Solver").pad_back(pl, pc) + ": UMFPACK");
      else
        _comm.print(String("Coarse Solver").pad_back(pl, pc) + ": BiCGStab-AmaVanka");
    }

    template<typename Mirror_>
    void assemble_orient_coarse_muxer(Global::Muxer<DenseVectorBlocked2ndMoment, Mirror_>& muxer_2, Global::Muxer<DenseVectorBlocked4thMoment, Mirror_>& muxer_4, const Control::Domain::VirtualLevel<DomainLevel_>& virt_lvl_coarse)
    {
      // assemble muxer parent
      if(virt_lvl_coarse.is_parent())
      {
        XASSERT(virt_lvl_coarse.is_child());

        const auto& layer_c = virt_lvl_coarse.layer_c();
        const DomainLevel_& level_p = virt_lvl_coarse.level_p();

        // loop over all children
        for(Index i(0); i < layer_c.child_count(); ++i)
        {
          const auto* child = level_p.find_patch_part(int(i));
          XASSERT(child != nullptr);
          Mirror_ child_mirror_v;
          Assembly::MirrorAssembler::assemble_mirror(child_mirror_v, level_p.space_velo, *child);
          muxer_2.push_child(child_mirror_v.clone(LAFEM::CloneMode::Deep));
          muxer_4.push_child(std::move(child_mirror_v));
        }
      }

      // assemble muxer child
      if(virt_lvl_coarse.is_child())
      {
        const auto& layer_c = virt_lvl_coarse.layer_c();
        const DomainLevel_& level_c = virt_lvl_coarse.level_c();

        Mirror_ parent_mirror_v;

        // manually set up an identity gather/scatter matrix
        {
          Index n = level_c.space_velo.get_num_dofs();
          parent_mirror_v = Mirror_(n, n);
          auto* idx = parent_mirror_v.indices();
          for(Index i(0); i < n; ++i)
            idx[i] = i;
        }

        // set parent and sibling comms
        muxer_2.set_parent(
          layer_c.sibling_comm_ptr(),
          layer_c.get_parent_rank(),
          parent_mirror_v.clone(LAFEM::CloneMode::Deep)
        );

        muxer_4.set_parent(
          layer_c.sibling_comm_ptr(),
          layer_c.get_parent_rank(),
          std::move(parent_mirror_v)
        );

        // compile muxer
        DenseVectorBlocked2ndMoment tmpl_v2(level_c.space_velo.get_num_dofs());
        DenseVectorBlocked4thMoment tmpl_v4(level_c.space_velo.get_num_dofs());
        muxer_2.compile(tmpl_v2);
        muxer_4.compile(tmpl_v4);
      }
    }

    void assemble_orient_burgers_mat(const Assembly::FullFiberOrientationTensorBurgersAssembler<DataType, IndexType, dim>& burgers_mat, const LocalVeloVectorType& vec_conv)
    {
      XASSERTM(!orient_2_levels.empty() && !orient_4_levels.empty(), "Orientation Matrices are not compiled!");

      SystemLevelType& the_system_level = *system_levels.front();
      //get reference to domain
      Control::Domain::PartiDomainControl<DomainLevel_>& domain = *_domain_ptr;

      typename SystemLevelType::GlobalVeloVector vec_conv_clone(
        &the_system_level.gate_velo, vec_conv.clone(LAFEM::CloneMode::Deep));

      // loop over all system levels
      for(std::size_t i(0); i < system_levels.size(); ++i)
      {
        // assemble our system matrix
        auto& loc_mat_a = system_levels.at(i)->matrix_sys.local().block_a();
        loc_mat_a.format();
        //the first vec_first_velo is just a placeholder variable, since we have to use it anyway, just use it for both...
        burgers_mat.assemble_matrix(loc_mat_a, vec_conv_clone.local(), orient_2_levels.at(i), orient_4_levels.at(i), domain.at(i)->space_velo, cubature);
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
          system_levels.at(i)->transfer_velo.trunc(vec_conv_clone, vec_crs);

          // the coarse vector is our next convection vector
          vec_conv_clone = std::move(vec_crs);
        }
        else
        {
          // this process is a child, so send truncation to parent
          system_levels.at(i)->transfer_velo.trunc_send(vec_conv_clone);
        }
      }
    }

    Solver::Status solve_basic_navier(GlobalSystemVector& vec_sol, const GlobalSystemVector& vec_rhs, bool verbose = false)
    {

      //this function should only be called directly after initialization, since we expect the matrix to be assembled with Stokes data...

      //get reference to domain
      Control::Domain::PartiDomainControl<DomainLevel_>& domain = *_domain_ptr;
      // fetch our finest levels
      DomainLevelType& the_domain_level = *domain.front();
      SystemLevelType& the_system_level = *system_levels.front();
      // --------------------------------------------------------------------------------------------
      // SOLVE STANDARD STOKES
      // --------------------------------------------------------------------------------------------

      // get our global solve matrix and filter
      GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
      GlobalSystemFilter& filter = the_system_level.filter_sys;

      _solver->set_plot_mode(Solver::PlotMode::iter);
      watch_nonlin_solver_init.start();
      _multigrid_hierarchy->init_numeric();
      _solver->init_numeric();
      watch_nonlin_solver_init.stop();

      // accumulate vanka sizes
      statistics.counts[Counts::vanka_data] = ama_vankas.size() > 0 ? Index(ama_vankas.front()->data_size()) : Index(0);
      statistics.bytes[Bytes::vanka] = 0ull;
      for(auto& v : ama_vankas)
        statistics.bytes[Bytes::vanka] += v->bytes();

      watch_stokes_solve.start();
      Solver::Status stokes_status = Solver::solve(*_solver, vec_sol, vec_rhs, matrix, filter);
      watch_stokes_solve.stop();

      statistics.counts[Counts::linsol_iter] = _solver->get_num_iter();

      _solver->done_numeric();
      _multigrid_hierarchy->done_numeric();

      FEAT::Statistics::compress_solver_expressions();

      if(!Solver::status_success(stokes_status))
      {
        _comm.print("\nERROR: BASIC STOKES LINEAR SOLVER BREAKDOWN\n");
        if(testmode)
          _comm.print("Test-Mode: FAILED");
        return stokes_status;
      }

      Solver::Status status = stokes_status;
      // --------------------------------------------------------------------------------------------
      // SOLVE NAVIER-STOKES (IF DESIRED)
      // --------------------------------------------------------------------------------------------
      if(navier)
      {
        if(verbose)
          _comm.print("\nSolving Navier-Stokes system...");

        if(!plot_mg_iter)
          _solver->set_plot_mode(Solver::PlotMode::none);

        // setup burgers assembler for matrix
        Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_mat;
        burgers_mat.nu = mu;
        burgers_mat.beta = rho;
        burgers_mat.frechet_beta = DataType(newton ? 1 : 0);
        burgers_mat.theta = DataType(0);
        burgers_mat.deformation = true;

        // setup burgers assembler for defect vector
        Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_def;
        burgers_def.nu = mu;
        burgers_def.beta = rho;
        burgers_def.theta = DataType(0);
        burgers_def.deformation = true;

        // vector of all non-linear defect norms
        std::vector<DataType> nl_defs;

        watch_nonlin_loop.start();


        // nonlinear loop
        for(Index nl_step(0); nl_step <= max_nl_iter; ++nl_step)
        {
          statistics.counts[Counts::nonlin_iter] = nl_step;

          // assemble nonlinear defect vector
          watch_nonlin_def_asm.start();
          vec_def.format();
          //i think full copy is not needed... for now we will take no chances
          vec_def.local().template at<0>().scale(vec_rhs.local().template at<0>(), DataType(1.));
          //since vec_rhs is a type 0 vector here, vec_def is also...
//           vec_def.from_1_to_0();
          // assemble burgers operator defect
          burgers_def.assemble_vector(vec_def.local().template at<0>(), vec_sol.local().template at<0>(),
                                      vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature, -1.0);
          // compute remainder of defect vector
          the_system_level.matrix_sys.local().block_b().apply(
            vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
          the_system_level.matrix_sys.local().block_d().apply(
            vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);
          // sync and filter
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
            if(verbose)
            {
              _comm.print(line);
              _comm.print("\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
            }
            if(testmode)
              _comm.print("Test-Mode: FAILED");
            return status;
          }
          else if(def_nl < nl_tol_abs)
          {
            if(verbose)
            {
              _comm.print(line);
              _comm.print("\nNonlinear solver converged!");
            }
            break;
          }
          else if(nl_step >= max_nl_iter)
          {
            if(verbose)
            {
              _comm.print(line);
              _comm.print("\nMaximum iterations reached!");
            }
            break;
          }
          else if((nl_step >= 3) && (DataType(0.95)*def_prev < def_nl))
          {
            if(verbose)
            {
              _comm.print(line);
              _comm.print("\nNonlinear solver stagnated!");
            }
            break;
          }

          // assemble burgers matrices on all levels
          watch_nonlin_mat_asm.start();
          {
            // get a clone of the global velocity vector
            typename SystemLevelType::GlobalVeloVector vec_conv(
              &the_system_level.gate_velo, vec_sol.local().template at<0>().clone());

            // loop over all system levels
            for(std::size_t i(0); i < system_levels.size(); ++i)
            {
              // assemble our system matrix
              auto& loc_mat_a = system_levels.at(i)->matrix_sys.local().block_a();
              loc_mat_a.format();
              burgers_mat.assemble_matrix(loc_mat_a, vec_conv.local(), domain.at(i)->space_velo, cubature);
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
          _multigrid_hierarchy->init_numeric();
          _solver->init_numeric();
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
              _solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
              // Also make sure that we gain at least 2 digits.
              _solver->set_tol_rel(1E-2);
            }
            else
            {
              // In the case if Picard iteration, we only expect linear convergence,
              // which (in analogy to Newton) leads us to the following estimate:
              //
              //      def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})
              //
              DataType abs_tol = def_nl * def_improve * DataType(0.1);
              _solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
              _solver->set_tol_rel(1E-2);
          }
        }
        else
        {
          // disable absolute tolerance
          _solver->set_tol_abs(1E+10);
          _solver->set_tol_rel(mg_tol_rel);
        }

        // solve linear system
        watch_nonlin_solver_apply.start();
        status = _solver->apply(vec_cor, vec_def);
        watch_nonlin_solver_apply.stop();

        statistics.counts[Counts::linsol_iter] += _solver->get_num_iter();


        line += String(" | ") + stringify(_solver->get_num_iter()).pad_front(3) + ": "
        + stringify_fp_sci(_solver->get_def_final(), 4) + " / "
        + stringify_fp_sci(_solver->get_def_final() / _solver->get_def_initial(), 4);
        if(adapt_tol && (nl_step > Index(0)))
          line += String(" [") + stringify_fp_sci(_solver->get_tol_abs(), 4) + "]";
        _comm.print(line);

        // release linear solver
        _solver->done_numeric();
        _multigrid_hierarchy->done_numeric();

        if(!Solver::status_success(status))
        {
          _comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
          if(testmode)
            _comm.print("Test-Mode: FAILED");
          return status;
        }

        // update solution
        vec_sol.axpy(vec_cor, vec_sol, DataType(1));

        FEAT::Statistics::compress_solver_expressions();
        // next non-linear iteration
      }

      watch_nonlin_loop.stop();

      // end of Navier-Stokes solve
    }

    return status;

  }

    template<typename SpaceOrientationType>
    void truncate_orientation_moment_vectors(const SpaceOrientationType& orient_space, const DenseVectorBlocked2ndMoment& second_moment, const DenseVectorBlocked4thMoment& fourth_moment)
    {
      //check the tensor sizes
      XASSERTM(second_moment.size() == orient_space.get_num_dofs(), "Vector size does not match dof number");
      XASSERTM(fourth_moment.size() == orient_space.get_num_dofs(), "Vector size does not match dof number");
      XASSERTM(orient_2_levels.empty() && orient_4_levels.empty(), "Orient deques are not empty");

      //get reference to domain
      Control::Domain::PartiDomainControl<DomainLevel_>& domain = *_domain_ptr;


      //push_back a vector of the size for VeloSpace
      orient_2_levels.push_back(DenseVectorBlocked2ndMoment(domain.front()->space_velo.get_num_dofs()));
      orient_4_levels.push_back(DenseVectorBlocked4thMoment(domain.front()->space_velo.get_num_dofs()));
      //and now project into those
      Assembly::FEInterpolator<SpaceVeloType, SpaceOrientationType>::interpolate(orient_2_levels.front(), second_moment, domain.front()->space_velo, orient_space);
      Assembly::FEInterpolator<SpaceVeloType, SpaceOrientationType>::interpolate(orient_4_levels.front(), fourth_moment, domain.front()->space_velo, orient_space);

      //create gates for the right Vectortypes
      //we need to assemble all gates beforhand...
      typedef typename SystemLevelType::BaseClass::VeloMirror VeloMirror;
      std::deque<Global::Gate<DenseVectorBlocked2ndMoment, VeloMirror>> orient_2_gate_deque;
      std::deque<Global::Gate<DenseVectorBlocked4thMoment, VeloMirror>> orient_4_gate_deque;

      // loop over all system levels
      for(std::size_t i(0); i < system_levels.size(); ++i)
      {
        //we create temporary clones of gates, mirrors and muxers to truncate our orientation vectors
        //first of all create gates for our OrienationVectors
        Global::Gate<DenseVectorBlocked2ndMoment, VeloMirror> orient_2_gate(_comm);
        Global::Gate<DenseVectorBlocked4thMoment, VeloMirror> orient_4_gate(_comm);

        //now push copies of ranks and mirrors
        orient_2_gate._ranks = system_levels.at(i)->gate_velo._ranks;
        orient_4_gate._ranks = system_levels.at(i)->gate_velo._ranks;

        //push mirrors
        for(auto& velo_mirrors_i : system_levels.at(i)->gate_velo._mirrors)
        {
          orient_2_gate._mirrors.push_back(velo_mirrors_i.clone());
          orient_4_gate._mirrors.push_back(velo_mirrors_i.clone());
        }

        //now compile _freqs, for this, create a dummy vectors of right vector size
        {
          DenseVectorBlocked2ndMoment orient_2_temp(domain.at(i)->space_velo.get_num_dofs());
          DenseVectorBlocked4thMoment orient_4_temp(domain.at(i)->space_velo.get_num_dofs());

          orient_2_gate.compile(std::move(orient_2_temp));
          orient_4_gate.compile(std::move(orient_4_temp));
        }

        //and now push them into our deque
        orient_2_gate_deque.push_back(std::move(orient_2_gate));
        orient_4_gate_deque.push_back(std::move(orient_4_gate));
      }


      //get a clone of the interpolated vectors and wrap them into a global vector
      Global::Vector<DenseVectorBlocked2ndMoment, VeloMirror> global_second_moment(&orient_2_gate_deque.front(), orient_2_levels.front().clone(LAFEM::CloneMode::Deep));
      Global::Vector<DenseVectorBlocked4thMoment, VeloMirror> global_fourth_moment(&orient_4_gate_deque.front(), orient_4_levels.front().clone(LAFEM::CloneMode::Deep));

      // loop over all system levels
      for(std::size_t i(0); i < system_levels.size(); ++i)
      {
        // restrict our orientation vectors
        //if there are no more "local" prozesses, no need to do anything
        if((i+1) >= domain.size_virtual())
          break;

        //now in any case (i+1) < size_virtual, so assemble muxers and transfers
        // define muxers
        typename Global::Muxer<DenseVectorBlocked2ndMoment, VeloMirror> muxer_2;
        typename Global::Muxer<DenseVectorBlocked4thMoment, VeloMirror> muxer_4;


        //assemble muxers
        assemble_orient_coarse_muxer(muxer_2, muxer_4, domain.at(i+1));


        //some transfer typedefs...
        typedef LAFEM::SparseMatrixBWrappedCSR<MemType, DataType, IndexType, dim*(dim+1)/2> LocalOrient2TransferMatrix;
        typedef LAFEM::SparseMatrixBWrappedCSR<MemType, DataType, IndexType, dim*(dim+1)*(dim+2)*(dim+3)/24> LocalOrient4TransferMatrix;
        typedef LAFEM::Transfer<LocalOrient2TransferMatrix> LocalOrient2Transfer;
        typedef LAFEM::Transfer<LocalOrient4TransferMatrix> LocalOrient4Transfer;

        //Init our global transfers
        Global::Transfer<LocalOrient2Transfer, VeloMirror> orient_2_transfer(&muxer_2);
        Global::Transfer<LocalOrient4Transfer, VeloMirror> orient_4_transfer(&muxer_4);


        //Now the idea is to get the local "unwrapped" transfer matrices and overwrite them with the ones of our system
        const auto& local_velo_trunc_mat_unwrapped = system_levels.at(i)->transfer_velo.local().get_mat_trunc().unwrap();  // i+1 seems to be the right size, but why?

        auto& local_orient_2_truncmat_unwrapped = orient_2_transfer.local().get_mat_trunc().unwrap();
        auto& local_orient_4_truncmat_unwrapped = orient_4_transfer.local().get_mat_trunc().unwrap();

        //and get shallow?! copies
        local_orient_2_truncmat_unwrapped.clone(local_velo_trunc_mat_unwrapped, LAFEM::CloneMode::Deep);
        local_orient_4_truncmat_unwrapped.clone(local_velo_trunc_mat_unwrapped, LAFEM::CloneMode::Deep);

        //important: overwrite the tmp vectors in our transfer operator...
        orient_2_transfer._vec_tmp = orient_2_transfer.get_mat_trunc().create_vector_l();
        orient_4_transfer._vec_tmp = orient_4_transfer.get_mat_trunc().create_vector_l();


        // does this process have another system level?
        if((i+1) < system_levels.size())
        {
          // create a coarse mesh orientation vectors
          Global::Vector<DenseVectorBlocked2ndMoment, VeloMirror> vec_crs_2(&orient_2_gate_deque.at(i+1), DenseVectorBlocked2ndMoment(domain.at(i+1)->space_velo.get_num_dofs()));
          Global::Vector<DenseVectorBlocked4thMoment, VeloMirror> vec_crs_4(&orient_4_gate_deque.at(i+1), DenseVectorBlocked4thMoment(domain.at(i+1)->space_velo.get_num_dofs()));


          // truncate fine mesh orientation vectors
          orient_2_transfer.trunc(global_second_moment, vec_crs_2);
          orient_4_transfer.trunc(global_fourth_moment, vec_crs_4);

          //now push local copies into our orientation vector deques
          orient_2_levels.push_back(vec_crs_2.local().clone(LAFEM::CloneMode::Deep));
          orient_4_levels.push_back(vec_crs_4.local().clone(LAFEM::CloneMode::Deep));

          // the coarse vector is our next convection vector
          global_second_moment = std::move(vec_crs_2);
          global_fourth_moment = std::move(vec_crs_4);
        }
        else
        {
          // this process is a child, so send truncation to parent
          orient_2_transfer.trunc_send(global_second_moment);
          orient_4_transfer.trunc_send(global_fourth_moment);
        }
      }


    }


    template<typename SpaceOrientationType>
    Solver::Status solve_navier(GlobalSystemVector& vec_sol, const GlobalSystemVector& vec_rhs, const SpaceOrientationType& orient_space, const DenseVectorBlocked2ndMoment& second_moment, const DenseVectorBlocked4thMoment& fourth_moment, bool verbose = false)
    {
//       print_parameters();
      //get reference to domain
      Control::Domain::PartiDomainControl<DomainLevel_>& domain = *_domain_ptr;
      // fetch our finest levels
      DomainLevelType& the_domain_level = *domain.front();
      SystemLevelType& the_system_level = *system_levels.front();

      // get our global solve matrix and filter
//       GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
      GlobalSystemFilter& filter = the_system_level.filter_sys;

      //       const Index num_levels = Index(domain.size_physical());

      //project our tensors from orient_space to the velocity_space and truncate through the levels into orient_2_levels and orient_4_levels
      watch_tensor_assembly.start();
      truncate_orientation_moment_vectors(orient_space, second_moment, fourth_moment);
      watch_tensor_assembly.stop();

      //now we can solve
      Solver::Status status{Solver::Status::undefined};


//       //solve full Stokes system
//       {
//         //create BurgersAssembler
//         Assembly::FullFiberOrientationTensorBurgersAssembler<DataType, IndexType, dim> burgers_mat_stokes;
//         burgers_mat_stokes.nu = mu;
//         burgers_mat_stokes.beta = DataType(0.);
//         burgers_mat_stokes.N_s = n_s;
//         burgers_mat_stokes.N_p = n_p;
//         burgers_mat_stokes.theta = alpha;
//
//         //format our solution vector
//         watch_nonlin_mat_asm.start();
//         vec_sol.format();
//         //and filter
//         the_system_level.filter_sys.filter_sol(vec_sol);
//         //now assemble our matrix
//         assemble_orient_burgers_mat(burgers_mat_stokes, vec_sol.local().template at<0>().clone());
//         watch_nonlin_mat_asm.stop();
//
//         // accumulate vanka sizes
//         statistics.counts[Counts::vanka_data] = ama_vankas.size() > 0 ? Index(ama_vankas.front()->data_size()) : Index(0);
//         statistics.bytes[Bytes::vanka] = 0ull;
//         for(auto& v : ama_vankas)
//           statistics.bytes[Bytes::vanka] += v->bytes();
//
//
//         //init our solver
//         watch_nonlin_solver_init.start();
//         _multigrid_hierarchy->init_numeric();
//         _solver->init_numeric();
//         watch_nonlin_solver_init.stop();
//         _solver->set_plot_mode(Solver::PlotMode::iter);
//
//         std::cout << "Starting to solve" << std::endl;
//
//         watch_nonlin_solver_apply.start();
//         status = Solver::solve(*_solver, vec_sol, vec_rhs, matrix, filter);
//         watch_nonlin_solver_apply.stop();
//
//         // release linear solver
//         _solver->done_numeric();
//         _multigrid_hierarchy->done_numeric();
//
//         FEAT::Statistics::compress_solver_expressions();
//
//         if(!Solver::status_success(status))
//         {
//           _comm.print("\nERROR: LINEAR SOLVER BREAKDOWN FOR FULL STOKES\n");
//           if(testmode)
//             _comm.print("Test-Mode: FAILED");
//           return status;
//         }
//       }

      // --------------------------------------------------------------------------------------------
      // SOLVE NAVIER-STOKES (IF DESIRED)
      // --------------------------------------------------------------------------------------------
      if(verbose)
        _comm.print("\nSolving Navier-Stokes system...");

      if(!plot_mg_iter)
        _solver->set_plot_mode(Solver::PlotMode::none);

      // setup burgers assembler for matrix
      Assembly::FullFiberOrientationTensorBurgersAssembler<DataType, IndexType, dim> burgers_mat;
      burgers_mat.nu = mu;
      burgers_mat.beta = rho;
      burgers_mat.frechet_beta = DataType(newton ? 1 : 0);
      burgers_mat.N_s = n_s;
      burgers_mat.N_p = n_p;
      burgers_mat.theta = alpha;

      // setup burgers assembler for defect vector
      Assembly::FullFiberOrientationTensorBurgersAssembler<DataType, IndexType, dim> burgers_def;
      burgers_def.nu = mu;
      burgers_def.beta = rho;
      burgers_def.N_s = n_s;
      burgers_def.N_p = n_p;
      burgers_def.theta = alpha;

      // vector of all non-linear defect norms
      std::vector<DataType> nl_defs;

      watch_nonlin_loop.start();

      // nonlinear loop
      for(Index nl_step(0); nl_step <= max_nl_iter; ++nl_step)
      {
        statistics.counts[Counts::nonlin_iter] = nl_step;

        // assemble nonlinear defect vector
        watch_nonlin_def_asm.start();
        vec_def.format();
        //i think full copy is not needed... for now we will take no chances
        vec_def.local().template at<0>().scale(vec_rhs.local().template at<0>(), DataType(1.));
        //since vec_rhs is a type 0 vector here, vec_def is also...
//         vec_def.from_1_to_0();
        // assemble burgers operator defect
        burgers_def.assemble_vector(vec_def.local().template at<0>(), vec_sol.local().template at<0>(), orient_2_levels.front(), orient_4_levels.front(), vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature, -1.0);
        // compute remainder of defect vector
        the_system_level.matrix_sys.local().block_b().apply(
          vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
        the_system_level.matrix_sys.local().block_d().apply(
          vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);
        // sync and filter
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
          if(verbose)
          {
            _comm.print(line);
            _comm.print("\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
            if(testmode)
              _comm.print("Test-Mode: FAILED");
          }
          return status;
        }
        else if(def_nl < nl_tol_abs)
        {
          if(verbose)
          {
            _comm.print(line);
            _comm.print("\nNonlinear solver converged!");
          }
          break;
        }
        else if(nl_step >= max_nl_iter)
        {
          if(verbose)
          {
            _comm.print(line);
            _comm.print("\nMaximum iterations reached!");
          }
          break;
        }
        else if((nl_step >= 3) && (DataType(0.95)*def_prev < def_nl))
        {
          if(verbose)
          {
            _comm.print(line);
            _comm.print("\nNonlinear solver stagnated!");
          }
          break;
        }

        // assemble burgers matrices on all levels
        watch_nonlin_mat_asm.start();
        assemble_orient_burgers_mat(burgers_mat, vec_sol.local().template at <0>());
        watch_nonlin_mat_asm.stop();

        // initialize linear solver
        watch_nonlin_solver_init.start();
        _multigrid_hierarchy->init_numeric();
        _solver->init_numeric();
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
            _solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
            // Also make sure that we gain at least 2 digits.
            _solver->set_tol_rel(1E-2);
          }
          else
          {
            // In the case if Picard iteration, we only expect linear convergence,
            // which (in analogy to Newton) leads us to the following estimate:
            //
            //      def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})
            //
            DataType abs_tol = def_nl * def_improve * DataType(0.1);
            _solver->set_tol_abs(Math::max(abs_tol, nl_tol_abs * DataType(0.01)));
            _solver->set_tol_rel(1E-2);
          }
        }
        else
        {
          // disable absolute tolerance
          _solver->set_tol_abs(1E+10);
          _solver->set_tol_rel(mg_tol_rel);
        }

        // solve linear system
        watch_nonlin_solver_apply.start();
        status = _solver->apply(vec_cor, vec_def);
        watch_nonlin_solver_apply.stop();

        statistics.counts[Counts::linsol_iter] += _solver->get_num_iter();


        line += String(" | ") + stringify(_solver->get_num_iter()).pad_front(3) + ": "
        + stringify_fp_sci(_solver->get_def_final(), 4) + " / "
        + stringify_fp_sci(_solver->get_def_final() / _solver->get_def_initial(), 4);
        if(adapt_tol && (nl_step > Index(0)))
          line += String(" [") + stringify_fp_sci(_solver->get_tol_abs(), 4) + "]";
        if(verbose)
          _comm.print(line);

        // release linear solver
        _solver->done_numeric();
        _multigrid_hierarchy->done_numeric();

        if(!Solver::status_success(status))
        {
          _comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
          if(testmode)
            _comm.print("Test-Mode: FAILED");
          return status;
        }

        // update solution
        vec_sol.axpy(vec_cor, vec_sol, DataType(1));

        FEAT::Statistics::compress_solver_expressions();
        // next non-linear iteration
      }

      watch_nonlin_loop.stop();

      // end of Navier-Stokes solve




    //we want to clean up our orient deques... depending on the problem, this is not always necessary, but the allociation cost should be negliable
    orient_2_levels.clear();
    orient_4_levels.clear();

    return status;

  }

  void compile_statistics()
  {
    // save timings
    statistics.times[Times::stokes_solve] = watch_stokes_solve.elapsed();
    statistics.times[Times::nonlin_total] = watch_nonlin_loop.elapsed();
    statistics.times[Times::nonlin_asm_def] = watch_nonlin_def_asm.elapsed();
    statistics.times[Times::nonlin_asm_mat] = watch_nonlin_mat_asm.elapsed();
    statistics.times[Times::linsol_init] = watch_nonlin_solver_init.elapsed();
    statistics.times[Times::linsol_apply] = watch_nonlin_solver_apply.elapsed();

    // get multigrid timings
    statistics.times[Times::mg_defect] = _multigrid_hierarchy->get_time_defect();
    statistics.times[Times::mg_smooth] = _multigrid_hierarchy->get_time_smooth();
    statistics.times[Times::mg_coarse] = _multigrid_hierarchy->get_time_coarse();
    statistics.times[Times::mg_transfer] = _multigrid_hierarchy->get_time_transfer();

    // accumulate vanka timings
    for(auto& v : ama_vankas)
    {
      statistics.times[Times::vanka_init_sym] += v->time_init_symbolic();
      statistics.times[Times::vanka_init_num] += v->time_init_numeric();
      statistics.times[Times::vanka_apply] += v->time_apply();
    }
  }

  //only update inflow boundary... this does not work as expected...
  template<typename DirichletFunctionType>
  void update_filters(DirichletFunctionType& diri_func, const String& InflowFacetName)
  {
    //get reference to domain
    Control::Domain::PartiDomainControl<DomainLevel_>& domain = *_domain_ptr;

    for(Index i(0); i < Index(domain.size_physical()); ++i)
    {
      // get our local velocity filter
      auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // first, clone the noflow filter
      fil_loc_v.clone(system_levels.at(i)->local_velo_filter_noflow);

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow;

      // try to fetch the inflow mesh part node
      auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(InflowFacetName);
      XASSERT(mesh_part_node != nullptr);

      auto* mesh_part = mesh_part_node->get_mesh();
      if (mesh_part != nullptr)
      {
        unit_asm_inflow.add_mesh_part(*mesh_part);
        unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, diri_func);
      }

      // finally, compile the system filter
      system_levels.at(i)->compile_system_filter();
    }

  }

  //update all non Neumann boundary
  template<typename DirichletFunctionType, typename NoSlipFunctionType>
  void update_filters(DirichletFunctionType& diri_func, NoSlipFunctionType& noslip_func, const String& InflowFacetName, const String& OutflowFacetName)
  {
    //get reference to domain
    Control::Domain::PartiDomainControl<DomainLevel_>& domain = *_domain_ptr;
    // the names of the mesh parts on which to assemble
    std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);
    const Index num_levels = Index(domain.size_physical());

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
          // bench3 starts with homogeneous inflow
          if((name == InflowFacetName))
          {
            // inflow
            unit_asm_inflow.add_mesh_part(*mesh_part);
          }
          else if(name != OutflowFacetName)
          {
            // outflow
            unit_asm_noflow.add_mesh_part(*mesh_part);
          }
        }
      }

      // assemble the no-flow filter first
      unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo, noslip_func);

      unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, diri_func);

      // finally, compile the system filter
      system_levels.at(i)->compile_system_filter();
    }

  }


  };


//   Control::Domain::PartiDomainControl<Control::Domain::StokesDomainLevel<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>,
//                                                      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>,
//                                                      Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>>,
//                                                      Space::Discontinuous::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>,
//                                                      Space::Discontinuous::Variant::StdPolyP<1>>>


// //   Control::Domain::PartiDomainControl<Control::Domain::StokesDomainLevel<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>, Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>, Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>>, Space::Discontinuous::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>, Space::Discontinuous::Variant::StdPolyP<1>>>>
//  //will auto just work fine?
//  template<int dim_>
// //  std::unique_ptr<Control::Domain::PartiDomainControl<Control::Domain::StokesDomainLevel<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>,
// //                                                      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>,
// //                                                      Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>>,
// //                                                      Space::Discontinuous::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>,
// //                                                      Space::Discontinuous::Variant::StdPolyP<1>>>>
//   Control::Domain::PartiDomainControl<Control::Domain::StokesDomainLevel<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>, Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>, Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>>, Space::Discontinuous::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>, Shape::Hypercube<dim_>::dimension, DataType>>, Space::Discontinuous::Variant::StdPolyP<1>>>>
//   create_hypercube_domain_control(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
//   {
//     // define our mesh type
//     typedef Shape::Hypercube<dim_> ShapeType;
//     typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, DataType> MeshType;
//     typedef Trafo::Standard::Mapping<MeshType> TrafoType;
//     typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
//     typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;
//
//     // create our domain control
//     typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
//     Control::Domain::PartiDomainControl<DomainLevelType>* domain_ptr = new Control::Domain::PartiDomainControl<DomainLevelType>(comm, true);
//
//     domain_ptr->parse_args(args);
//     domain_ptr->set_desired_levels(args.query("level")->second);
//     domain_ptr->create(mesh_reader);
//
//     return std::unique_ptr<Control::Domain::PartiDomainControl<DomainLevelType>>(domain_ptr);
//
//   }


}





#endif