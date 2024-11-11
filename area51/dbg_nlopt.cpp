// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/assembly/analytic_projector.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/none_filter.hpp>

#include <kernel/solver/fixed_step_linesearch.hpp>
#include <kernel/solver/newton_raphson_linesearch.hpp>
#include <kernel/solver/secant_linesearch.hpp>
#include <kernel/solver/mqc_linesearch.hpp>

#include <kernel/solver/alglib_wrapper.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/solver/nlsd.hpp>
#include <kernel/solver/hessian_precond.hpp>

#include <kernel/solver/test_aux/analytic_function_operator.hpp>
#include <kernel/solver/test_aux/function_traits.hpp>

#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/runtime.hpp>

using namespace FEAT;
using namespace FEAT::Solver;

/**
 * \brief Tool for debugging of nonlinear optimizers
 *
 * The tool can be used for testing nonlinear optimizers on a few well-known functions. It has been used for
 * nonlinear operators generated from AnalyticFunctions using the AnalyticFunctionOperator wrapper class. Have a look
 * at the Analytic::Common namespace candidates. There must be an implementation of a helper traits class in
 * kernel/solver/test_aux/function_traits.hpp specifying the real minima and a starting point.
 *
 * Some part of the solver configuration can be done on the command line (solvertype, direction update,
 * preconditioner), while other options must be specified at compile time (the linesearch, memory architecture and
 * floating point precision). The solver by default logs all iterates and produces a vtk output of the domain of
 * interest and a polyline of all iterates so one can study the convergence behavior of the solver.
 *
 * This can easily be extended to use nonlinear operators from other sources, but the whole plotting procedures
 * require that we optimize a scalar function in two or three variables.
 *
 */

template<typename Solver_, typename Operator_>
int run(Solver_& solver, Operator_& op)
{
  typedef Operator_ OperatorType;
  typedef typename OperatorType::DataType DataType;
  typedef typename OperatorType::IndexType IndexType;

  typedef typename OperatorType::FunctionType FunctionType;
  typedef typename OperatorType::PointType PointType;
  static constexpr int dim = PointType::n;

  typedef LAFEM::DenseVector<DataType, IndexType> VertexVectorType;
  typedef LAFEM::DenseVectorBlocked<DataType, IndexType, dim> VertexFieldVectorType;

  typedef OptimizationTestTraits<DataType, FunctionType> TestTraitsType;

  // Get an initial guess from the Traits class for the given function
  PointType starting_point(DataType(0));
  TestTraitsType::get_starting_point(starting_point);

  // From the traits class, get the set of minimal points
  std::deque<PointType> min_points;
  TestTraitsType::get_minimal_points(min_points);

  // Get the default domain for post processing
  Tiny::Matrix<DataType, dim, 2> domain_bounds(DataType(0));
  for(int d(0); d < dim; ++d)
  {
    TestTraitsType::get_domain_bounds(domain_bounds[d], d);
  }

  std::cout << "dbg-nlopt: Minimising " << TestTraitsType::name();
  std::cout << ", starting point is " << starting_point << "\n";
  std::cout << "Known minima:";
  const auto& jt = min_points.end();
  auto it = min_points.begin();
  for(; it != jt; ++it)
  {
    std::cout << " " << *it;
  }
  std::cout << "\n";

  FunctionType my_function;

  solver->init();
  solver->set_max_iter(10000);
  solver->set_tol_fval(DataType(0));
  solver->set_tol_step(Math::eps<DataType>());
  //solver->set_tol_abs(Math::sqrt(Math::eps<DataType>()));
  //solver->set_tol_rel(Math::sqrt(Math::eps<DataType>()));
  solver->set_tol_abs(Math::eps<DataType>());
  solver->set_tol_rel(Math::eps<DataType>());
  solver->set_plot_mode(Solver::PlotMode::all);

  //PropertyMap config;
  //solver->write_config(&config, "dbg-nlopt-solver");
  //config.dump(std::cout);

  // This will hold the solution
  auto sol = op.create_vector_r();
  sol(0,starting_point);
  // We need a dummy rhs
  auto rhs = op.create_vector_r();

  // Solve the optimization problem
  solver->correct(sol, rhs);

  // Check the distance between solution and minimal points
  DataType min_dist(Math::Limits<DataType>::max());

  it = min_points.begin();
  for(; it != jt; ++it)
  {
    DataType dist((sol(0) - *it).norm_euclid());
    if(dist  < min_dist)
      min_dist = dist;
  }


  String filename("");

  std::cout << "Found solution " << sol << ", distance to nearest minimum: " << stringify_fp_sci(min_dist) << "\n";

  // Write the iterates to file if we can
  if(solver->iterates != nullptr)
  {
    filename = "dbg_nlopt_"+TestTraitsType::name()+"_iterates";
    std::cout << "Writing iterates to file " << filename << ".vtu" <<"\n";

    // This will hold the Tiny::Vectors representing the iterates, which are contained in LAFEM::DenseVectorBlocked
    // containers of length 1 because of the solver interface
    std::deque<typename VertexFieldVectorType::ValueType> points;
    const auto& jt2(solver->iterates->end());
    auto it2(solver->iterates->begin());
    for(; it2 != jt2; ++it2)
    {
      // Get the iterate from the DenseVectorBlocked
      auto tmp2 = (*it2)(0);
      points.push_back(tmp2);

      // Check if we need to adjust the default domain
      for(int d(0); d < dim; ++d)
      {
        domain_bounds[d](0) = Math::min(domain_bounds[d](0), tmp2[d]);
        domain_bounds[d](1) = Math::max(domain_bounds[d](1), tmp2[d]);
      }
    }

    // This Factory will create the polyline mesh
    Geometry::PolylineFactory<dim, DataType> pl_factory(points);
    typedef Geometry::ConformalMesh<Shape::Hypercube<1>, dim, DataType> PolylineMesh;
    PolylineMesh polyline(pl_factory);

    // Write this to file
    Geometry::ExportVTK<PolylineMesh> polyline_writer(polyline);
    polyline_writer.write(filename);

  }

  // Enlarge the plot domain by 10 percent in each direction
  PointType scaling(DataType(0));
  for(int d(0); d < dim; ++d)
    scaling(d) = DataType(1.1)*(domain_bounds[d](1) - domain_bounds[d](0));

  filename = "dbg_nlopt_"+TestTraitsType::name()+"_domain";
  std::cout << "Writing domain to file " << filename << ".vtu" << "\n";

  // Shape for the mesh for plotting
  typedef Shape::Hypercube<dim> ShapeType;
  // The mesh tupe
  typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, DataType> MeshType;
  // We need a transformation so we can project the analytic function
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  // Create the mesh as refined unitcube. Refinement level 6 is usually enough
  Geometry::RefineFactory<MeshType, Geometry::UnitCubeFactory> mesh_factory(6);
  MeshType* mesh(new MeshType(mesh_factory));
  // Trafo for the interpolation of the analytic function
  TrafoType trafo(*mesh);

  auto& vtx = mesh->get_vertex_set();
  // Transform the mesh
  for(IndexType i(0); i < mesh->get_num_entities(0); ++i)
  {
    for(int d(0); d < dim; ++d)
      vtx[i](d) = scaling(d)*vtx[i](d) + (domain_bounds(d,0) - scaling(d)*DataType(0.05));
  }

  // This will hold the function values
  VertexVectorType vtx_vec(mesh->get_num_entities(0));
  // This will hold the gradient
  VertexFieldVectorType vtx_grad(mesh->get_num_entities(0));

  // Project the analytic function
  Assembly::AnalyticVertexProjector::project(vtx_vec, my_function, trafo);
  // Get the gradient and project it
  Analytic::Gradient<FunctionType> my_grad(my_function);
  Assembly::AnalyticVertexProjector::project(vtx_grad, my_grad, trafo);

  // Write the mesh
  Geometry::ExportVTK<MeshType> writer(*mesh);
  writer.add_vertex_scalar("f", vtx_vec.elements());
  writer.add_vertex_vector("grad", vtx_grad);
  writer.write(filename);

  // Clean up
  delete mesh;

  // Finish the solver
  solver->done();

  return 0;

}

static void display_help()
{
  std::cout << "dbg-nlopt: Nonlinear optimizer debugging tool." << "\n";
  std::cout << "Required arguments:" << "\n";
  std::cout << " --solver [String]: Available solvers are NLCG, NLSD" << "\n";
#ifdef FEAT_HAVE_ALGLIB
  std::cout << "                    and ALGLIBMinLBFGS, ALGLIBMinCG" << "\n";
#else
  std::cout << "                    Compiling with the alglib token in the buildid "
  << "\n" <<
  "                    will enable ALGLIBMinLBFGS and ALGLIBMinCG" << "\n";
#endif // FEAT_HAVE_ALGLIB
  std::cout << "Optional arguments:" << "\n";
  std::cout << " --help: Displays this text" << "\n";
  std::cout << " --linesearch [String]: Available linesearches for NLCG and NLSD are NewtonRaphsonLinesearch,"
  << "\n";
  std::cout << "                    SecantLinesearch and MQCLinesearch (default)"
  << "\n";
  std::cout <<
    " --precon [String]: Available preconditioners for NLCG and NLSD are Hessian, ApproximateHessian, none(default)"
    << "\n" <<
    " --direction_update [String]: Available search direction updates for NLCG are DaiYuan, DYHSybrid,"
    << "\n" <<
    "                              FletcherReeves, HestenesStiefel, HagerZhang and PolakRibiere (default)."
    << "\n";
#ifdef FEAT_HAVE_ALGLIB
  std::cout <<
    "                              Available search direction updates for ALGLIBMinCG are DaiYuan,"
    << "\n" <<
    "                              and DYHSHybrid."
    << "\n";
#endif // FEAT_HAVE_ALGLIB
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  typedef double DataType;
  typedef unsigned int IndexType;

  // The analytic function we want to minimize. Look at the Analytic::Common namespace for other candidates.
  // There must be an implementation of a helper traits class in kernel/solver/test_aux/function_traits.hpp
  // specifying the real minima and a starting point.
  typedef Analytic::Common::RosenbrockFunction AnalyticFunctionType;
  typedef AnalyticFunctionOperator<DataType, IndexType, AnalyticFunctionType> OperatorType;
  typedef typename OperatorType::PointType PointType;
  static constexpr int dim = PointType::n;

  typedef LAFEM::NoneFilterBlocked<DataType, IndexType, dim> FilterType;

  // The analytic function
  AnalyticFunctionType my_function;
  // Create the (nonlinear) operator
  OperatorType my_op(my_function);
  // The filter
  FilterType my_filter;

  // Ugly way to get a preconditioner, or not
  std::shared_ptr<NLOptPrecond<typename OperatorType::VectorTypeL, FilterType> > my_precond(nullptr);
  // We might need a linesearch
  std::shared_ptr<Linesearch<OperatorType, FilterType> > my_linesearch(nullptr);

  // create an argument parser
  SimpleArgParser args(argc, argv);

  // add all supported options
  args.support("direction_update");
  args.support("help");
  args.support("linesearch");
  args.support("precon");
  args.support("solver");
  args.support("steplength");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() || args.check("help") > -1)
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << "\n";
    display_help();

    Runtime::abort();
    return 1;
  }

  String solver_name("");
  auto* solver_pair(args.query("solver"));
  if(solver_pair == nullptr)
  {
    display_help();
    return 1;
  }
  else
  {
    solver_name = solver_pair->second.front();
  }

  if(solver_name == "NLCG")
  {
    // The default linesearch is MQCLinesearch
    String linesearch_name("MQCLinesearch");
    auto* linesearch_pair(args.query("linesearch"));
    if(linesearch_pair != nullptr)
      linesearch_name = linesearch_pair->second.front();

    if(linesearch_name== "NewtonRaphsonLinesearch")
    {
      my_linesearch = new_newton_raphson_linesearch(my_op, my_filter);
    }
    else if(linesearch_name== "SecantLinesearch")
    {
      my_linesearch = new_secant_linesearch(my_op, my_filter);
    }
    else if(linesearch_name== "MQCLinesearch")
    {
      my_linesearch = new_mqc_linesearch(my_op, my_filter);
    }
    else if(linesearch_name== "FixedStepLinesearch")
    {
      DataType steplength(1);
      auto* steplength_pair(args.query("steplength"));

      if(steplength_pair != nullptr)
      {
        steplength = DataType(std::stod(steplength_pair->second.front()));
      }

      my_linesearch = new_fixed_step_linesearch(my_op, my_filter, steplength);
    }
    my_linesearch->set_max_iter(20);
    //my_linesearch->set_plot_mode(Solver::PlotMode::all);

    // The default is no preconditioner
    String precon_name("none");
    // Check if any preconditioner was specified on the command line
    auto* precon_pair(args.query("precon"));
    if(precon_pair != nullptr)
      precon_name = precon_pair->second.front();

    if(precon_name== "Hessian")
      my_precond = new_hessian_precond(my_op, my_filter);
    else if(precon_name == "ApproximateHessian")
      my_precond = new_approximate_hessian_precond(my_op, my_filter);
    else if(precon_name != "none")
      XABORTM("Got invalid precon_name: "+precon_name);

    NLCGDirectionUpdate my_direction_update(NLCGDirectionUpdate::DYHSHybrid);
    auto update_pair(args.query("direction_update"));
    if(update_pair != nullptr)
      my_direction_update << update_pair->second.front();

    auto my_solver = new_nlcg(my_op, my_filter, my_linesearch, my_direction_update, true, my_precond);

    int ret = run(my_solver, my_op);

    // Clean up
    my_solver.reset();
    return ret;
  }
  else if(solver_name == "NLSD")
  {
    // The default linesearch is MQCLinesearch
    String linesearch_name("MQCLinesearch");
    auto* linesearch_pair(args.query("linesearch"));
    if(linesearch_pair != nullptr)
      linesearch_name = linesearch_pair->second.front();

    if(linesearch_name== "NewtonRaphsonLinesearch")
      my_linesearch = new_newton_raphson_linesearch(my_op, my_filter);
    else if(linesearch_name== "SecantLinesearch")
      my_linesearch = new_secant_linesearch(my_op, my_filter);
    else if(linesearch_name== "MQCLinesearch")
      my_linesearch = new_mqc_linesearch(my_op, my_filter);
    else if(linesearch_name== "FixedStepLinesearch")
    {
      DataType steplength(1);
      auto* steplength_pair(args.query("steplength"));

      if(steplength_pair != nullptr)
        steplength = DataType(std::stod(steplength_pair->second.front()));

      my_linesearch = new_fixed_step_linesearch(my_op, my_filter, steplength);
    }

    // The default is no preconditioner
    String precon_name("none");
    // Check if any preconditioner was specified on the command line
    auto* precon_pair(args.query("precon"));
    if(precon_pair != nullptr)
      precon_name = precon_pair->second.front();

    if(precon_name== "Hessian")
      my_precond = new_hessian_precond(my_op, my_filter);
    else if(precon_name == "ApproximateHessian")
      my_precond = new_approximate_hessian_precond(my_op, my_filter);
    else if(precon_name != "none")
      XABORTM("Got invalid precon_name: "+precon_name);

    auto* update_pair(args.query("direction_update"));
    if(update_pair != nullptr)
      std::cout << "Parameter direction_update specified for NLSD solver, ignoring it." << "\n";

    auto my_solver = new_nlsd(my_op, my_filter, my_linesearch, true, my_precond);

    int ret = run(my_solver, my_op);

    // Clean up
    my_solver.reset();
    return ret;
  }
#ifdef FEAT_HAVE_ALGLIB
  else if (solver_name == "ALGLIBMinLBFGS")
  {
    auto my_solver = new_alglib_minlbfgs(my_op, my_filter);

    int ret = run(my_solver, my_op);

    my_solver.reset();
    return ret;
  }
  else if (solver_name == "ALGLIBMinCG")
  {
    NLCGDirectionUpdate my_direction_update(NLCGDirectionUpdate::DYHSHybrid);
    auto* update_pair(args.query("direction_update"));
    if(update_pair != nullptr)
    {
      String update_name(update_pair->second.front());
      if(update_name == "DaiYuan")
        my_direction_update = NLCGDirectionUpdate::DaiYuan;
      else if(update_name == "DYHSHybrid")
        my_direction_update = NLCGDirectionUpdate::DYHSHybrid;
      else
        XABORTM("Got invalid ALGLIBMinCG direction update: "+update_name);
    }

    auto my_solver = new_alglib_mincg(my_op, my_filter, my_direction_update, true);

    int ret = run(my_solver, my_op);

    my_solver.reset();
    return ret;
  }
#endif // FEAT_HAVE_ALGLIB
  else
  {
    XABORTM("dbg-nlopt got invalid solver name: "+solver_name);
  }
}
