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

template<typename Mem_, typename DT_, typename IT_, typename Function_, template<typename, typename> class Linesearch_>
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

  public:
    NLCGTest(DT_ tol_, String precon_type_) :
      FullTaggedTest<Mem_, DT_, IT_>("NLCGTest"),
      _tol(tol_),
      _precon_type(precon_type_)
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
      else if(_precon_type != "none")
        throw InternalError("Got invalid precon_type: "+_precon_type);

      //auto my_precond = nullptr;
      auto solver = new_nlcg(my_op, my_filter, my_linesearch, false, my_precond);
      solver->init();
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot(false);
      solver->set_max_iter(250);
      //std::cout << solver->get_formated_solver_tree() << std::endl;

      auto sol = my_op.create_vector_r();
      auto rhs = my_op.create_vector_r();

      PointType starting_point(DT_(0));
      TestTraitsType::get_starting_point(starting_point);
      sol(0,starting_point);

      my_op.prepare(sol);

      solver->correct(sol, rhs);
      //std::cout << "sol = " << sol << std::endl;
      //std::cout << "f(sol) = " << my_op.compute_fval() << std::endl;

      solver->done();

      std::deque<PointType> min_points;
      TestTraitsType::get_minimal_points(min_points);

      DT_ min_dist(Math::Limits<DT_>::max());

      const auto& jt = min_points.end();
      auto it = min_points.begin();
      for(; it != jt; ++it)
      {
        DT_ dist((sol(0) - *it).norm_euclid());
        //std::cout << "dist = " << stringify_fp_sci(dist) << std::endl;
        if(dist  < min_dist)
          min_dist = dist;
      }
      //std::cout << "min_dist " << stringify_fp_sci(min_dist) << std::endl;
      TEST_CHECK_MSG(min_dist < _tol,"min_dist = "+stringify_fp_sci(min_dist)+" > "+stringify_fp_sci(_tol)+" = tol");

    }
};

float tol_f = Math::pow(Math::eps<float>(), float(0.8));
NLCGTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction, Solver::NRLinesearch>
opt_hb_f(tol_f,"ApproximateHessian");

double tol_d = Math::pow(Math::eps<double>(), double(0.6));
NLCGTest<Mem::Main, double, Index, Analytic::Common::HimmelblauFunction, Solver::SecantLinesearch> opt_hb_d(tol_d,"none");

double tol_d2 = Math::eps<double>();
NLCGTest<Mem::Main, double, unsigned int, Analytic::Common::RosenbrockFunction, Solver::NRLinesearch>
opt_rb_d(tol_d2,"none");

double tol_d3 = Math::pow(Math::eps<double>(), double(0.25));
NLCGTest<Mem::Main, double, unsigned int, Analytic::Common::BazaraaShettyFunction, Solver::SecantLinesearch>
opt_bs_d(tol_d3,"ApproximateHessian");

//OptimiserTest<Mem::Main, double, Index, Analytic::Common::RosenbrockFunction> opt_rb(tol_d);
//OptimiserTest<Mem::Main, float, unsigned int, Analytic::Common::BazaraaShettyFunction> opt_bs(tol_f);
