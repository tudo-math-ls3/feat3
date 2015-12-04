#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/solver/hessian_precond.hpp>
#include <kernel/solver/test_aux/analytic_function_operator.hpp>
#include <kernel/solver/test_aux/function_traits.hpp>

using namespace FEAST;
using namespace FEAST::Solver;
using namespace FEAST::TestSystem;

/**
 * \brief Test class template for the nonlinear CG optimiser
 *
 * This has the analytic function and the linesearch as template parameters because not all linesearches are suitable
 * for all functions. Same with preconditioners.
 *
 */
template
<
  typename Mem_, typename DT_, typename IT_,
  typename Function_,
  template<typename, typename> class Linesearch_
>
class NLCGTest:
  public FullTaggedTest<Mem_, DT_, IT_>
{
  public:
    typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
    typedef typename OperatorType::PointType PointType;
    typedef OptimisationTestTraits<DT_, Function_> TestTraitsType;

    typedef LAFEM::NoneFilterBlocked<Mem_, DT_, IT_, 2> FilterType;
    typedef Linesearch_<OperatorType, FilterType> LinesearchType;

  private:
    DT_ _tol;
    String _precon_type;
    NLCGDirectionUpdate _update;

  public:
    NLCGTest(DT_ exponent_, String precon_type_, NLCGDirectionUpdate update_) :
      FullTaggedTest<Mem_, DT_, IT_>("NLCGTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _precon_type(precon_type_),
      _update(update_)
    {
    }

    void run() const
    {
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

      if(_precon_type == "ApproximateHessian")
        my_precond = new_approximate_hessian_precond(my_op, my_filter);
      else if(_precon_type == "Hessian")
        my_precond = new_hessian_precond(my_op, my_filter);
      else if(_precon_type != "none")
        throw InternalError("Got invalid precon_type: "+_precon_type);

      //auto my_precond = nullptr;
      auto solver = new_nlcg(my_op, my_filter, my_linesearch, false, my_precond);
      solver->init();
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot(false);
      solver->set_max_iter(250);
      solver->set_direction_update(_update);
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
      TEST_CHECK_MSG(min_dist < _tol,"min_dist = "+stringify_fp_sci(min_dist)+" > "+stringify_fp_sci(_tol)+" = tol");
    }
};

// The Himmelblau function is not too difficult, so low tolerance
NLCGTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction, Solver::NRLinesearch>
opt_hb_f(float(0.8),"ApproximateHessian", NLCGDirectionUpdate::PolakRibiere);

// The same with Secant linesearch, without preconditioner and with Fletcher-Reeves update
NLCGTest<Mem::Main, double, Index, Analytic::Common::HimmelblauFunction, Solver::SecantLinesearch> opt_hb_d(double(0.8),"none", NLCGDirectionUpdate::FletcherReeves);

// The Rosenbrock function's steep valley is bad for secant linesearch, so use Newton Raphson
NLCGTest<Mem::Main, float, unsigned int, Analytic::Common::RosenbrockFunction, Solver::NRLinesearch>
opt_rb_d(float(0.8),"Hessian", NLCGDirectionUpdate::PolakRibiere);

// The Hessian of the Bazaraa/Shetty function is singular at the optimal point, so Newton Raphson linesearch does not
// work very well, so just use the secant linesearch
NLCGTest<Mem::Main, double, unsigned int, Analytic::Common::BazaraaShettyFunction, Solver::SecantLinesearch>
opt_bs_d(double(0.3),"none", NLCGDirectionUpdate::PolakRibiere);

// Rosenbrock with secant linesearch, preconditioning and Polak-Ribière in quad precision
#ifdef FEAST_HAVE_QUADMATH
NLCGTest<Mem::Main, __float128, Index, Analytic::Common::RosenbrockFunction, Solver::SecantLinesearch>
opt_rb_q(__float128(0.9),"Hessian", NLCGDirectionUpdate::PolakRibiere);
#endif

// Running this in CUDA is really nonsensical because all operator evaluations use Tiny::Vectors which reside in
// Mem::Main anyway, so apart from the occasional axpy nothing is done on the GPU. It should work nonetheless.
#ifdef FEAST_BACKENDS_CUDA
NLCGTest<Mem::CUDA, double, unsigned int, Analytic::Common::BazaraaShettyFunction, Solver::SecantLinesearch>
opt_bs_f_cuda(double(0.25),"none", NLCGDirectionUpdate::FletcherReeves);
#endif
