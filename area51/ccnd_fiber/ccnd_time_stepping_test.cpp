#include <area51/ccnd_fiber/ccnd_fiber_solver.hpp>
#include <area51/ccnd_fiber/ccnd_fiber_common.hpp>
#include <area51/ccnd_fiber/ccnd_test_functions.hpp>


namespace CCND_FIBER
{


  template<int dim_, template<typename> class SpaceOrientTemp>
  void run_dim(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
//     auto domain_ptr = create_hypercube_domain_control<dim_>(args, comm, mesh_reader);

    // define our mesh type
    typedef Shape::Hypercube<dim_> ShapeType;
//     typedef typename TrafoType::ShapeType ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, DataType> MeshType;
//     typedef typename TrafoType::MeshType MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    typedef SpaceOrientTemp<typename DomainLevelType::TrafoType> SpaceOrientationType;
    //create unique ptr handling our domain
    std::unique_ptr<Control::Domain::PartiDomainControl<DomainLevelType>> domain_ptr(new Control::Domain::PartiDomainControl<DomainLevelType>(comm, true));

    domain_ptr->parse_args(args);
    domain_ptr->set_desired_levels(args.query("level")->second);
    domain_ptr->create(mesh_reader);

    //we need the constants for our function_wrapper:
    DataType v_max = parse(args, "v-max", DataType(dim_) * DataType(0.5));
    DataType mu = parse(args, "mu", DataType(0.1));
    DataType rho = parse(args, "rho", DataType(1.));
    DataType n_s = parse(args, "n-s", DataType(3.9));
    DataType n_p = parse(args, "n-p", DataType(7.5));
    DataType epsilon = parse(args, "epsilon", DataType(1.));
    DataType period = DataType(Math::pi<DataType>()) * parse(args, "period", DataType(2.));
    DataType alpha = parse(args, "expo_alpha", DataType(-0.7));
    DataType delta_t = parse(args, "delta-t", DataType(0.01));
    DataType t_max = parse(args, "t-max", DataType(8));
    IndexType testing_steps = parse(args, "testing-steps", IndexType(10));
//     bool conv = (args.check("stokes") < 0);
    const bool conv = true;
    const bool verbose = (args.check("verbose") >= 0);


//     {
//       static constexpr std::size_t pl = 30u;
//       static constexpr char pc = '.';
//       comm.print("\nProblem Parameters:");
//       comm.print(String("Mu").pad_back(pl, pc) + ": " + stringify(mu));
//       comm.print(String("Rho").pad_back(pl, pc) + ": " + stringify(rho));
//       comm.print(String("N_s").pad_back(pl, pc) + ": " + stringify(n_s));
//       comm.print(String("N_p").pad_back(pl, pc) + ": " + stringify(n_p));
//       comm.print(String("Conv") + ": " + stringify(conv));
//       comm.print(String("expo_alpha").pad_back(pl, pc) + ": " + stringify(alpha));
//
//     }

    //our DirichletFunction inflow function
//     ContractInflowFunction<dim_> inflow_func;
//     ContractInflowFunction other_func;
    //we create fitting boundary functions by first initializing the wrapper class
    FullTensorTimeExpo<dim_> function_wrapper(v_max, mu, rho, n_s, n_p, epsilon, period, alpha, conv);

    //now we can define our functions
    typename FullTensorTimeExpo<dim_>::Velo inflow_func(function_wrapper);
    typename FullTensorTimeExpo<dim_>::Robin outflow_func(function_wrapper);
    typename FullTensorTimeExpo<dim_>::RHS rhs_func(function_wrapper);

    //set time
    DataType cur_time = DataType(0.);

    inflow_func.set_time(cur_time);
    outflow_func.set_time(cur_time);
    rhs_func.set_time(cur_time);


    //zero boundary/rhs function
    Analytic::Common::ConstantVectorFunction<dim_> const_func(DataType(0));

    String cubature_string("auto-degree:6");
    String inflow_facet("bnd:l");
    String outflow_facet("bnd:r");

    //create our UnsteadySolver interface
    CCNDUnsteadyInterface<DomainLevelType> unsteady_solver(comm, args, std::move(domain_ptr), inflow_func, rhs_func, outflow_func, inflow_func, cubature_string, inflow_facet, outflow_facet);
//     CCNDUnsteadyInterface<DomainLevelType> steady_solver(comm, args, std::move(domain_ptr), inflow_func, const_func, const_func, const_func, cubature_string, inflow_facet, outflow_facet);


    //create solution space and vector
    SpaceOrientationType orient_space(unsteady_solver._inner_solver._domain_ptr->front()->trafo);
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim_> orient_sol(orient_space.get_num_dofs());
    orient_sol.format();

    //get our starting solution, solves the basic, steady Navier-Stokes system
//     Solver::Status status;// = unsteady_solver.solve_basic_navier(orient_sol, orient_space);
    //get our starting solution by projecting in our space
    //for all purposes, lets format and filter...
    unsteady_solver._vec_sol.format();
    Assembly::Interpolator::project(unsteady_solver._vec_sol.local().template at<0>(), inflow_func, unsteady_solver._inner_solver._domain_ptr->front()->space_velo);

    //interpolate into our OrientationSpace
    Assembly::FEInterpolator<SpaceOrientationType, SpaceVeloType>::interpolate(orient_sol, unsteady_solver._vec_sol.local().template at<0>(), orient_space, unsteady_solver._inner_solver._domain_ptr->front()->space_velo);
    //copy vec_sol_1 anyway...
    unsteady_solver._prev_sol1.copy(unsteady_solver._vec_sol);
    unsteady_solver._prev_sol2.copy(unsteady_solver._vec_sol);

    DataType h0_project_error{DataType(0)};
    DataType h1_project_error{DataType(0)};
    //now calculate the error
    {
        // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        // Post-Processing: Computing L2/H1-Errors

//         comm.print("\nComputing errors against reference solution...");
        Cubature::DynamicFactory cubature_postproc("auto-degree:6");

        Assembly::VectorErrorInfo<DataType, dim_> velo_errors = Assembly::VectorErrorComputer<1>::compute(
          orient_sol, inflow_func, orient_space, cubature_postproc);

        // As usual, we want to compute the L2- and H1-errors against our reference solution.
        // For this, we first compute the errors on our patch by analysing our local solution vector:
        //       Assembly::ScalarErrorInfo<DataType> pres_errors = Assembly::ScalarErrorComputer<0>::compute(
        //         vec_sol.local().template at<1>(), pres_sol, pres_space, cubature_postproc);

        // And then we need to synchronize the errors over our communicator to sum up the errors of
        // each patch to obtain the errors over the whole domain:
        velo_errors.synchronize(comm);
        h0_project_error = velo_errors.norm_h0;
        h1_project_error = velo_errors.norm_h1;


        //     if(velo_errors.norm_h0 >= DataType(4e-2)/DataType(the_domain_level.get_level_index()) || velo_errors.norm_h1 >= DataType(1e-1)/DataType(the_domain_level.get_level_index()))
        //       comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");

        // And let's print the errors to the console; we need to use the "format_string" function here,
        // as the "print" function accepts only String objects as input:
//         comm.print("\nErrors for velocity:");
//         comm.print(velo_errors.format_string());
      }
//     std::cout << "?" << std::endl;
//     if(!Solver::status_success(status))
//     {
//       comm.print("\nERROR: OUTER SOLVER BREAKDOWN\n");
//       return;
//     }

    //create dummy Orientation Tensors
    typename CCNDSteadyInterface<DomainLevelType>::DenseVectorBlocked2ndMoment second_moment_vector(orient_space.get_num_dofs(), DataType(0));
    typename CCNDSteadyInterface<DomainLevelType>::DenseVectorBlocked4thMoment fourth_moment_vector(orient_space.get_num_dofs(), DataType(0));

    //create our tensor functions and project onto the orientation space
    //function do not map to the vectors...
    typename FullTensorTimeExpo<dim_>::Orient second_moment_func(function_wrapper);
    typename FullTensorTimeExpo<dim_>::Fourth_Moment fourth_moment_func(function_wrapper);

    IndexType time_step(0u);

    //lets start our time_step
    while((cur_time += delta_t) < t_max && time_step <= testing_steps)
    {
      //next time_step
      ++time_step;

      //set time
      inflow_func.set_time(cur_time);
      outflow_func.set_time(cur_time);
      rhs_func.set_time(cur_time);
      second_moment_func.set_time(cur_time);
      fourth_moment_func.set_time(cur_time);

      Assembly::Interpolator::project(second_moment_vector, second_moment_func, orient_space);
      Assembly::Interpolator::project(fourth_moment_vector, fourth_moment_func, orient_space);

      //reinitaliaze filter and right handside
      unsteady_solver.update_filters(inflow_func, inflow_func);
      unsteady_solver.update_rhs(rhs_func, outflow_func);


      //call our solver
      Solver::Status status = unsteady_solver.solve_time_step(orient_sol, orient_space, second_moment_vector, fourth_moment_vector);

      if(status != Solver::Status::success)
      {
        comm.print("Test FAILED : Solver did not reach convergence criterium!");
        return;
      }

      if(verbose)
      {
        // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        // Post-Processing: Computing L2/H1-Errors

//         comm.print("\nComputing errors against reference solution...");
        Cubature::DynamicFactory cubature_postproc("auto-degree:6");

        Assembly::VectorErrorInfo<DataType, dim_> velo_errors = Assembly::VectorErrorComputer<1>::compute(
          orient_sol, inflow_func, orient_space, cubature_postproc);

        // As usual, we want to compute the L2- and H1-errors against our reference solution.
        // For this, we first compute the errors on our patch by analysing our local solution vector:
        //       Assembly::ScalarErrorInfo<DataType> pres_errors = Assembly::ScalarErrorComputer<0>::compute(
        //         vec_sol.local().template at<1>(), pres_sol, pres_space, cubature_postproc);

        // And then we need to synchronize the errors over our communicator to sum up the errors of
        // each patch to obtain the errors over the whole domain:
        velo_errors.synchronize(comm);



        //     if(velo_errors.norm_h0 >= DataType(4e-2)/DataType(the_domain_level.get_level_index()) || velo_errors.norm_h1 >= DataType(1e-1)/DataType(the_domain_level.get_level_index()))
        //       comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");

        // And let's print the errors to the console; we need to use the "format_string" function here,
        // as the "print" function accepts only String objects as input:
        comm.print("\nErrors for velocity:");
        comm.print(velo_errors.format_string());
      }


    }

    //now calculate the error
      {
        // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        // Post-Processing: Computing L2/H1-Errors

//         comm.print("\nComputing errors against reference solution...");
        Cubature::DynamicFactory cubature_postproc("auto-degree:6");

        Assembly::VectorErrorInfo<DataType, dim_> velo_errors = Assembly::VectorErrorComputer<1>::compute(
          orient_sol, inflow_func, orient_space, cubature_postproc);

        // As usual, we want to compute the L2- and H1-errors against our reference solution.
        // For this, we first compute the errors on our patch by analysing our local solution vector:
        //       Assembly::ScalarErrorInfo<DataType> pres_errors = Assembly::ScalarErrorComputer<0>::compute(
        //         vec_sol.local().template at<1>(), pres_sol, pres_space, cubature_postproc);

        // And then we need to synchronize the errors over our communicator to sum up the errors of
        // each patch to obtain the errors over the whole domain:
        velo_errors.synchronize(comm);

        if(velo_errors.norm_h0 >= h0_project_error * DataType(1.8))
        {
          comm.print("FAILED Test for orient space " + SpaceOrientationType::name() + " failed with h0_project_error = " + stringify(h0_project_error) + " and h0_velocity_error = " + stringify(velo_errors.norm_h0));
        }
        if(velo_errors.norm_h1 >= h1_project_error * DataType(1.8))
        {
          comm.print("FAILED Test for orient space " + SpaceOrientationType::name() + " failed with h1_project_error = " + stringify(h1_project_error) + " and h1_velocity_error = " + stringify(velo_errors.norm_h1));
        }
        if(velo_errors.norm_h0 < h0_project_error * DataType(1.8) && velo_errors.norm_h1 < h1_project_error * DataType(1.8))
        {
          comm.print("SUCCESS Test for orient space " + SpaceOrientationType::name());
        }



        //     if(velo_errors.norm_h0 >= DataType(4e-2)/DataType(the_domain_level.get_level_index()) || velo_errors.norm_h1 >= DataType(1e-1)/DataType(the_domain_level.get_level_index()))
        //       comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");

        // And let's print the errors to the console; we need to use the "format_string" function here,
        // as the "print" function accepts only String objects as input:
//         comm.print("\nErrors for velocity:");
//         comm.print(velo_errors.format_string());
      }
    //print solver statistics
//     unsteady_solver.compile_statistics();




  }


  void main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
//     comm.print("Number of Processes: " + stringify(comm.size()));
//     comm.print("Floating Point Type: " + String(fp_typename) + " precision");

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("level");
    args.support("vtk");
    args.support("mesh");
    args.support("mu");
    args.support("rho");
    args.support("n-s");
    args.support("n-p");
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
    args.support("save-sol");
    args.support("epsilon");
    args.support("euler");
    args.support("delta-t");
    args.support("expo_alpha");
    args.support("t-max");
    args.support("t-expo");
    args.support("testing-steps");
    args.support("orient-space-2");
    args.support("verbose");
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
//     if(args.check("mat-files") < 2)
//     {
//       comm.print(std::cerr, "ERROR: Mandatory option '--mat-files <2ndOrderMat-file> <4thOrderMat-file>' is missing!");
//       FEAT::Runtime::abort();
//     }

    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run 2D or 3D ?
    //if(args.check("isoparam") >= 0)
    const bool orient_lagrange_2 = args.check("orient-space-2") >= 0;
    {
      if(mesh_type == "conformal:hypercube:2:2")
      {
        if(orient_lagrange_2)
        {
          run_dim<2, Space::Lagrange2::Element>(args, comm, mesh_reader);
        }
        else
        {
          run_dim<2, Space::Lagrange1::Element>(args, comm, mesh_reader);
        }
      }
      else if(mesh_type == "conformal:hypercube:3:3")
      {
        if(orient_lagrange_2)
        {
          run_dim<3, Space::Lagrange2::Element>(args, comm, mesh_reader);
        }
        else
        {
          run_dim<3, Space::Lagrange1::Element>(args, comm, mesh_reader);
        }
      }
      else
      {
        comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
        FEAT::Runtime::abort();
      }
    }
  }
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialize(argc, argv);
  try
  {
    CCND_FIBER::main(argc, argv);
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