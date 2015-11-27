#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
//#include <kernel/analytic/auto_derive.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/solver/hessian_precond.hpp>
#include <kernel/solver/test_aux/analytic_function_operator.hpp>
#include <kernel/solver/test_aux/function_traits.hpp>
#include <kernel/util/simple_arg_parser.hpp>
// DEBUG
#include <kernel/assembly/analytic_projector.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/trafo/standard/mapping.hpp>

using namespace FEAST;
using namespace FEAST::Solver;

template
<
  typename Mem_, typename DT_, typename IT_,
  typename Function_,
  template<typename, typename> class Linesearch_
>
int run(String precon_type)
{
  typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
  typedef typename OperatorType::PointType PointType;

  static constexpr int dim = PointType::n;

  typedef LAFEM::NoneFilterBlocked<Mem_, DT_, IT_, dim> FilterType;
  typedef Linesearch_<OperatorType, FilterType> LinesearchType;

  typedef LAFEM::DenseVector<Mem_, DT_, IT_> VertexVectorType;
  typedef LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, dim> VertexFieldVectorType;

  typedef OptimisationTestTraits<DT_, Function_> TestTraitsType;

  // Get an initial guess from the Traits class for the given function
  PointType starting_point(DT_(0));
  TestTraitsType::get_starting_point(starting_point);

  // From the traits class, get the set of minimal points
  std::deque<PointType> min_points;
  TestTraitsType::get_minimal_points(min_points);

  // Get the default domain for post processing
  Tiny::Matrix<DT_, dim, 2> domain_bounds(DT_(0));
  for(int d(0); d < dim; ++d)
  {
    TestTraitsType::get_domain_bounds(domain_bounds[d], d);
  }


  std::cout << "dbg-nlopt: Minimising " << TestTraitsType::name();
  std::cout << ", starting point is " << starting_point << std::endl;
  std::cout << "Known minima:";
  const auto& jt = min_points.end();
  auto it = min_points.begin();
  for(; it != jt; ++it)
  {
    std::cout << " " << *it;
  }
  std::cout << std::endl;

  // The analytic function
  Function_ my_function;
  // Create the (nonlinear) operator
  OperatorType my_op(my_function);
  // The filter
  FilterType my_filter;

  // Create the linesearch
  LinesearchType my_linesearch(my_op, my_filter);
  my_linesearch.set_plot(false);

  // Ugly way to get a preconditioner, or not
  std::shared_ptr<SolverBase<typename OperatorType::VectorTypeL> > my_precond(nullptr);

  if(precon_type == "Hessian")
    my_precond = new_hessian_precond(my_op, my_filter);
  else if(precon_type == "ApproximateHessian")
    my_precond = new_approximate_hessian_precond(my_op, my_filter);
  else if(precon_type != "none")
    throw InternalError("Got invalid precon_type: "+precon_type);

  auto solver = new_nlcg(my_op, my_filter, my_linesearch, true, my_precond);
  solver->init();
  solver->set_tol_rel(Math::eps<DT_>());
  solver->set_plot(true);
  solver->set_max_iter(250);
  std::cout << "Using solver " << solver->get_formated_solver_tree() << std::endl;

  // This will hold the solution
  auto sol = my_op.create_vector_r();
  sol(0,starting_point);
  // We need a dummy rhs
  auto rhs = my_op.create_vector_r();

  my_op.prepare(sol);

  // Solve the optimisation problem
  Status st = solver->correct(sol, rhs);

  // Check the distance betwen solution and minimal points
  DT_ min_dist(Math::Limits<DT_>::max());

  it = min_points.begin();
  for(; it != jt; ++it)
  {
    DT_ dist((sol(0) - *it).norm_euclid());
    if(dist  < min_dist)
      min_dist = dist;
  }
  std::cout << "Found solution " << sol << ", distance to nearest minimum: " << stringify_fp_sci(min_dist) << std::endl;

  // Print solver summary
  std::cout << solver->get_plot_name() << ": " << st << ", " << solver->get_num_iter();
  std::cout << " its, defect initial/final: " << stringify_fp_sci(solver->get_def_initial());
  std::cout << " / " << stringify_fp_sci(solver->get_def_final()) << std::endl;
  std::cout << "Needed evaluations: " << my_op.num_func_evals << " (func) / " << my_op.num_grad_evals;
  std::cout <<  " (grad) / " << my_op.num_hess_evals << " (hess)" << std::endl;

  // Finish the solver
  solver->done();


  solver->done();

  String filename("");

  // Write the iterates to file if we can
  if(solver->iterates != nullptr)
  {
    filename = "dbg_nlopt_"+TestTraitsType::name()+"_iterates";
    std::cout << "Writing iterates to file " << filename << ".vtu" <<std::endl;

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

      // Check if we needd to adjust the default domain
      for(int d(0); d < dim; ++d)
      {
        domain_bounds[d](0) = Math::min(domain_bounds[d](0), tmp2[d]);
        domain_bounds[d](1) = Math::max(domain_bounds[d](1), tmp2[d]);
      }
    }

    // This Factory will create the polyline mesh
    Geometry::PolylineFactory<dim, dim, DT_> pl_factory(points);
    typedef Geometry::ConformalMesh<Shape::Hypercube<1>, dim, dim, DT_> PolylineMesh;
    PolylineMesh polyline(pl_factory);

    // Write this to file
    Geometry::ExportVTK<PolylineMesh> polyline_writer(polyline);
    polyline_writer.write(filename);

  }

  // Enlarge the plot domain by 10 percent in each direction
  PointType scaling(DT_(0));
  for(int d(0); d < dim; ++d)
    scaling(d) = DT_(1.1)*(domain_bounds[d](1) - domain_bounds[d](0));

  filename = "dbg_nlopt_"+TestTraitsType::name()+"_domain";
  std::cout << "Writing domain to file " << filename << ".vtu" << std::endl;

  // Shape for the mesh for plotting
  typedef Shape::Hypercube<dim> ShapeType;
  // The mesh tupe
  typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DT_> MeshType;
  // We need a transformation so we can project the analytic function
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  // Create the mesh as refined unitcube. Refinement level 6 is usually enough
  Geometry::RefineFactory<MeshType, Geometry::UnitCubeFactory> mesh_factory(6);
  MeshType* mesh(new MeshType(mesh_factory));
  // Trafo for the interpolation of the analytic function
  TrafoType trafo(*mesh);

  auto& vtx = mesh->get_vertex_set();
  // Transform the mesh
  for(Index i(0); i < mesh->get_num_entities(0); ++i)
  {
    for(int d(0); d < dim; ++d)
      vtx[i](d) = scaling(d)*vtx[i](d) + (domain_bounds(d,0) - scaling(d)*DT_(0.05));
  }

  // This will hold the function values
  VertexVectorType vtx_vec(mesh->get_num_entities(0));
  // This will hold the gradient
  VertexFieldVectorType vtx_grad(mesh->get_num_entities(0));

  // Project the analytic function
  Assembly::AnalyticVertexProjector::project(vtx_vec, my_function, trafo);
  // Get the gradient and project it
  Analytic::Gradient<Function_> my_grad(my_function);
  Assembly::AnalyticVertexProjector::project(vtx_grad, my_grad, trafo);

  // Write the mesh
  Geometry::ExportVTK<MeshType> writer(*mesh);
  writer.add_scalar_vertex("f", vtx_vec.elements());
  writer.add_field_vertex_blocked_vector("grad", vtx_grad);
  writer.write(filename);

  // Clean up
  delete mesh;

  return 0;

}

int main(int argc, char* argv[])
{
  // The analytic function we want to minimise. Look at the Analytic::Common namespace for other candidates
  typedef Analytic::Common::HimmelblauFunction AnalyticFunctionType;

  // create an argument parser
  SimpleArgParser args(argc, argv);

  // add all supported options
  args.support("precon");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
    return 1;
  }

  String precon_name("none");
  auto* precon_pair(args.query("precon"));
  if(precon_pair != nullptr)
  {
    precon_name = precon_pair->second.front();
  }

  // Canditates for the linesearch type (last template parameter) are SecantLinesearch and NRLinesearch
  return run<Mem::Main, double, Index, AnalyticFunctionType, Solver::SecantLinesearch> (precon_name);

}
