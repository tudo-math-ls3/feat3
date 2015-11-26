#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/analytic/auto_derive.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/solver/hessian_precond.hpp>
#include <kernel/solver/test_aux/analytic_function_operator.hpp>
#include <kernel/solver/test_aux/function_traits.hpp>
#include <kernel/util/math.hpp>
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
int run(DT_ tol, String precon_type)
{
  typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
  typedef typename OperatorType::PointType PointType;
  typedef OptimisationTestTraits<DT_, Function_> TestTraitsType;

  typedef LAFEM::NoneFilterBlocked<Mem_, DT_, IT_, 2> FilterType;
  typedef Linesearch_<OperatorType, FilterType> LinesearchType;

  typedef LAFEM::DenseVector<Mem_, DT_, IT_> VertexVectorType;
  typedef LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, 2> VertexFieldVectorType;

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

  if(precon_type == "ApproximateHessian")
    my_precond = new_approximate_hessian_precond(my_op, my_filter);
  else if(precon_type != "none")
    throw InternalError("Got invalid precon_type: "+precon_type);

  auto solver = new_nlcg(my_op, my_filter, my_linesearch, false, my_precond);
  solver->init();
  solver->set_tol_rel(Math::eps<DT_>());
  solver->set_plot(false);
  solver->set_max_iter(250);
  //std::cout << solver->get_formated_solver_tree() << std::endl;

  // This will hold the solution
  auto sol = my_op.create_vector_r();
  // We need a dummy rhs
  auto rhs = my_op.create_vector_r();

  // Get an initial guess from the Traits class for the given function
  PointType starting_point(DT_(0));
  TestTraitsType::get_starting_point(starting_point);
  sol(0,starting_point);

  my_op.prepare(sol);

  // Solve the optimisation problem
  solver->correct(sol, rhs);

  solver->done();

  // From the traits class, get the set of minimal points
  std::deque<PointType> min_points;
  TestTraitsType::get_minimal_points(min_points);

  // Check the distance betwen solution and minimal points
  DT_ min_dist(Math::Limits<DT_>::max());

  const auto& jt = min_points.end();
  auto it = min_points.begin();
  for(; it != jt; ++it)
  {
    DT_ dist((sol(0) - *it).norm_euclid());
    if(dist  < min_dist)
      min_dist = dist;
  }
  // Check if we found a valid minimum
  std::cout << "Distance from solution to nearest minimum: " << stringify_fp_sci(min_dist) << ", tolerance was "
  << stringify_fp_sci(tol) << std::endl;

  solver->done();

  String filename("");
  if(solver->iterates != nullptr)
  {
    filename = "dbg_nlop_"+TestTraitsType::name()+"_iterates";
    std::deque<typename VertexFieldVectorType::ValueType> points;
    const auto& jt2(solver->iterates->end());
    auto it2(solver->iterates->begin());
    for(; it2 != jt2; ++it2)
    {
      auto tmp2 = (*it2)(0);
      points.push_back(tmp2);
    }
    Geometry::PolylineFactory<2,2, DT_> pl_factory(points);
    typedef Geometry::ConformalMesh<Shape::Hypercube<1>,2,2,DT_> PolylineMesh;
    PolylineMesh polyline(pl_factory);

    Geometry::ExportVTK<PolylineMesh> polyline_writer(polyline);
    polyline_writer.write(filename);

  }

  // Mesh
  filename = "dbg_nlop_"+TestTraitsType::name()+"_domain";
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DT_> MeshType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  Geometry::RefineFactory<MeshType, Geometry::UnitCubeFactory> mesh_factory(6);
  MeshType* mesh(new MeshType(mesh_factory));
  TrafoType trafo(*mesh);

  Tiny::Matrix<DT_, MeshType::world_dim, 2> domain_bounds(DT_(0));
  PointType scaling(DT_(0));

  for(int d(0); d < MeshType::world_dim; ++d)
  {
    TestTraitsType::get_domain_bounds(domain_bounds[d], d);
    scaling(d) = domain_bounds[d](1) - domain_bounds[d](0);
  }

  auto& vtx = mesh->get_vertex_set();
  for(Index i(0); i < mesh->get_num_entities(0); ++i)
  {
    for(int d(0); d < MeshType::world_dim; ++d)
      vtx[i](d) = scaling(d)*vtx[i](d)+domain_bounds(d,0);
  }

  VertexVectorType vtx_vec(mesh->get_num_entities(0));
  VertexFieldVectorType vtx_grad(mesh->get_num_entities(0));

  Assembly::AnalyticVertexProjector::project(vtx_vec, my_function, trafo);
  Analytic::Gradient<Function_> my_grad(my_function);
  Assembly::AnalyticVertexProjector::project(vtx_grad, my_grad, trafo);

  Geometry::ExportVTK<MeshType> writer(*mesh);
  writer.add_scalar_vertex("f", vtx_vec.elements());
  writer.add_field_vertex_blocked_vector("grad", vtx_grad);
  writer.write(filename);

  delete mesh;

  return 0;

}

//int main(int argc, char* argv[])
int main()
{
  double tol_d = Math::pow(Math::eps<double>(), double(0.6));

  return run<Mem::Main, double, Index, Analytic::Common::HimmelblauFunction, Solver::SecantLinesearch> (tol_d,"none");

}
