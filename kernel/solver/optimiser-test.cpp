#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/solver/mqc_linesearch.hpp>
#include <kernel/solver/newton_raphson_linesearch.hpp>
#include <kernel/solver/secant_linesearch.hpp>
#include <kernel/solver/alglib_wrapper.hpp>
#include <kernel/solver/nlcg.hpp>
#include <kernel/solver/nlsd.hpp>
#include <kernel/solver/hessian_precond.hpp>
#include <kernel/solver/nlopt_precond.hpp>
#include <kernel/solver/test_aux/analytic_function_operator.hpp>
#include <kernel/solver/test_aux/function_traits.hpp>

using namespace FEAT;
using namespace FEAT::Solver;
using namespace FEAT::TestSystem;

/**
 * \brief Test class template for the nonlinear CG optimiser
 *
 * This has the analytic function and the linesearch as template parameters because not all linesearches are suitable
 * for all functions. Same with preconditioners.
 *
 */
template<typename Mem_, typename DT_, typename IT_, typename Function_>
class NLCGTest:
  public FullTaggedTest<Mem_, DT_, IT_>
{
  public:
    typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
    typedef typename OperatorType::PointType PointType;
    typedef OptimisationTestTraits<DT_, Function_> TestTraitsType;

    typedef LAFEM::NoneFilterBlocked<Mem_, DT_, IT_, 2> FilterType;

  private:
    const DT_ _tol;
    const Index _max_iter;
    const Index _max_func_evals;
    const String _linesearch_type;
    const String _precon_type;
    NLCGDirectionUpdate _update;

  public:
    explicit NLCGTest(const DT_ exponent_, const Index max_iter_, const Index max_func_evals_, const String& linesearch_type_, const String& precon_type_,
    NLCGDirectionUpdate update_) :
      FullTaggedTest<Mem_, DT_, IT_>("NLCGTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _max_iter(max_iter_),
      _max_func_evals(max_func_evals_),
      _linesearch_type(linesearch_type_),
      _precon_type(precon_type_),
      _update(update_)
    {
    }

    virtual ~NLCGTest()
    {
    }

    void run() const override
    {
      int failed_checks(0);
      // The analytic function
      Function_ my_function;
      // Create the (nonlinear) operator
      OperatorType my_functional(my_function);
      // The filter
      FilterType my_filter;

      // Create the linesearch
      std::shared_ptr<Solver::Linesearch<OperatorType, FilterType>> my_linesearch;
      if(_linesearch_type == "NewtonRaphsonLinesearch")
      {
        my_linesearch = new_newton_raphson_linesearch(my_functional, my_filter, false);
      }
      else if(_linesearch_type == "SecantLinesearch")
      {
        my_linesearch = new_secant_linesearch(my_functional, my_filter, DT_(1e-2), false);
        // We need to make the line search more exact in this case
        my_linesearch->set_tol_decrease(DT_(1e-4));
        my_linesearch->set_tol_curvature(DT_(1e-1));
      }
      else if(_linesearch_type == "MQCLinesearch")
      {
        my_linesearch = new_mqc_linesearch(my_functional, my_filter, false);
      }
      else
      {
        throw InternalError(__func__, __FILE__, __LINE__, "Got invalid linesearch_type: "+_linesearch_type);
      }

      my_linesearch->set_max_iter(20);
      my_linesearch->set_plot_mode(Solver::PlotMode::none);

      // Ugly way to get a preconditioner, or not
      std::shared_ptr<NLOptPrecond<typename OperatorType::VectorTypeL, FilterType>> my_precond(nullptr);

      if(_precon_type == "ApproximateHessian")
      {
        my_precond = new_approximate_hessian_precond(my_functional, my_filter);
      }
      else if(_precon_type == "Hessian")
      {
        my_precond = new_hessian_precond(my_functional, my_filter);
      }
      else if(_precon_type != "none")
      {
        throw InternalError(__func__, __FILE__, __LINE__, "Got invalid precon_type: "+_precon_type);
      }

      auto solver = new_nlcg(my_functional, my_filter, my_linesearch, _update, false, my_precond);

      solver->init();
      solver->set_tol_abs(Math::eps<DT_>());
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot_mode(Solver::PlotMode::summary);
      solver->set_max_iter(100);

      // This will hold the solution
      auto sol = my_functional.create_vector_r();
      // We need a dummy rhs
      auto rhs = my_functional.create_vector_r();

      // Get an initial guess from the Traits class for the given function
      PointType starting_point(DT_(0));
      TestTraitsType::get_starting_point(starting_point);
      sol(0,starting_point);

      // Solve the optimisation problem
      solver->correct(sol, rhs);

      // From the traits class, get the set of minimal points
      std::deque<PointType> min_points;
      TestTraitsType::get_minimal_points(min_points);

      // Check the distance between solution and minimal points
      DT_ min_dist(Math::Limits<DT_>::max());

      const auto& jt = min_points.end();
      auto it = min_points.begin();
      for(; it != jt; ++it)
      {
        DT_ dist((sol(0) - *it).norm_euclid());
        if(dist  < min_dist)
        {
          min_dist = dist;
        }
      }

      // Check if we stayed in the iteration number bound
      if(solver->get_num_iter() > _max_iter)
      {
        ++failed_checks;
        std::cout << "num_iter = "+stringify(solver->get_num_iter())+" > " +stringify(_max_iter)+" = max_iter" << std::endl;
      }

      // Check if we stayed in the functional evaluation number bound
      if(my_functional.get_num_func_evals() > _max_func_evals)
      {
        ++failed_checks;
        std::cout << "num_func_evals = "
          +stringify(my_functional.get_num_func_evals())+" > "+stringify(_max_func_evals)+" = max_func_evals" << std::endl;
      }

      // Check if we found a valid minimum
      if(min_dist > _tol)
      {
        ++failed_checks;
        std::cout << "min_dist = "+stringify_fp_sci(min_dist)+" > " +stringify_fp_sci(_tol)+" = tol" << std::endl;
      }

      solver->done();

      TEST_CHECK_MSG(failed_checks == 0,TestTraitsType::name()+" "+solver->name()+" "+my_linesearch->name()+": "
          +stringify(failed_checks)+ " failed checks");

    }
};

// The first three are for comparison with ALGLIBMinCG
NLCGTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
nlcg_sw_hb_f(float(0.6),Index(14), Index(39),"MQCLinesearch","none", NLCGDirectionUpdate::DaiYuan);

NLCGTest<Mem::Main, double, Index, Analytic::Common::RosenbrockFunction>
nlcg_sw_rb_d(double(0.8),Index(40), Index(118),"MQCLinesearch","none", NLCGDirectionUpdate::DYHSHybrid);

NLCGTest<Mem::Main, double, Index, Analytic::Common::BazaraaShettyFunction>
nlcg_sw_bs_d(double(0.33),Index(64), Index(240),"MQCLinesearch","none", NLCGDirectionUpdate::DaiYuan);

// This is the weird Hager-Zhang update
NLCGTest<Mem::Main, double, Index, Analytic::Common::HimmelblauFunction>
nlcg_s_hb_d(double(0.4),Index(27) ,Index(127),
"NewtonRaphsonLinesearch","ApproximateHessian", NLCGDirectionUpdate::HagerZhang);

NLCGTest<Mem::Main, double, unsigned int, Analytic::Common::BazaraaShettyFunction>
nlcg_s_bs_d(double(0.19), Index(31), Index(97), "SecantLinesearch", "none", NLCGDirectionUpdate::HestenesStiefel);

NLCGTest<Mem::Main, double, unsigned int, Analytic::Common::RosenbrockFunction>
nlcg_nr_rb_d(double(0.5), Index(47), Index(231),"NewtonRaphsonLinesearch","ApproximateHessian", NLCGDirectionUpdate::FletcherReeves);

NLCGTest<Mem::Main, double, Index, Analytic::Common::RosenbrockFunction>
nlcg_sw_hessian_rb_d(double(0.95), Index(25), Index(42), "MQCLinesearch","Hessian", NLCGDirectionUpdate::DYHSHybrid);

NLCGTest<Mem::Main, double, Index, Analytic::Common::GoldsteinPriceFunction>
nlcg_sw_hessian_gp_d(double(0.5), Index(16), Index(99), "MQCLinesearch","Hessian", NLCGDirectionUpdate::PolakRibiere);

#ifdef FEAT_HAVE_QUADMATH
NLCGTest<Mem::Main, __float128, Index, Analytic::Common::RosenbrockFunction>
nlcg_nr_rb_q(__float128(0.55), Index(33), Index(137), "NewtonRaphsonLinesearch", "Hessian",
NLCGDirectionUpdate::PolakRibiere);

NLCGTest<Mem::Main, __float128, Index, Analytic::Common::HimmelblauFunction>
nlcg_sw_bs_q(__float128(1), Index(23), Index(74), "MQCLinesearch", "none", NLCGDirectionUpdate::HestenesStiefel);
#endif

// Running this in CUDA is really nonsensical because all operator evaluations use Tiny::Vectors which reside in
// Mem::Main anyway, so apart from the occasional axpy nothing is done on the GPU. It should work nonetheless.
#ifdef FEAT_HAVE_CUDA
NLCGTest<Mem::CUDA, float, unsigned int, Analytic::Common::HimmelblauFunction>
nlcg_sw_hb_f_cuda(float(0.9), Index(11), Index(23), "MQCLinesearch", "Hessian", NLCGDirectionUpdate::FletcherReeves);

NLCGTest<Mem::CUDA, double, unsigned int, Analytic::Common::BazaraaShettyFunction>
nlcg_s_bs_d_cuda(double(0.18), Index(77), Index(261), "SecantLinesearch", "none", NLCGDirectionUpdate::PolakRibiere);
#endif

/**
 * \brief Test class template for the nonlinear Steepest Descent optimiser
 *
 * This has the analytic function and the linesearch as template parameters because not all linesearches are suitable
 * for all functions. Same with preconditioners.
 *
 */
template < typename Mem_, typename DT_, typename IT_, typename Function_ >
class NLSDTest:
  public FullTaggedTest<Mem_, DT_, IT_>
{
  public:
    typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
    typedef typename OperatorType::PointType PointType;
    typedef OptimisationTestTraits<DT_, Function_> TestTraitsType;

    typedef LAFEM::NoneFilterBlocked<Mem_, DT_, IT_, 2> FilterType;

  private:
    DT_ _tol;
    const Index _max_iter;
    const Index _max_func_evals;
    const String _linesearch_type;
    const String _precon_type;

  public:
    explicit NLSDTest(const DT_ exponent_, const Index max_iter_, const Index max_func_evals_,
    const String& linesearch_type_, const String& precon_type_) :
      FullTaggedTest<Mem_, DT_, IT_>("NLSDTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _max_iter(max_iter_),
      _max_func_evals(max_func_evals_),
      _linesearch_type(linesearch_type_),
      _precon_type(precon_type_)
    {
    }

    virtual ~NLSDTest()
    {
    }

    void run() const override
    {
      int failed_checks(0);
      // The analytic function
      Function_ my_function;
      // Create the (nonlinear) operator
      OperatorType my_functional(my_function);
      // The filter
      FilterType my_filter;

      // Create the linesearch
      std::shared_ptr<Solver::Linesearch<OperatorType, FilterType>> my_linesearch(nullptr);
      if(_linesearch_type == "NewtonRaphsonLinesearch")
      {
        my_linesearch = new_newton_raphson_linesearch(my_functional, my_filter, false);
        // We need to make the line search more exact in this case
        my_linesearch->set_tol_decrease(DT_(1e-4));
        my_linesearch->set_tol_curvature(DT_(1e-1));
      }
      else if(_linesearch_type == "SecantLinesearch")
      {
        my_linesearch = new_secant_linesearch(my_functional, my_filter, DT_(1e-2), false);
        // We need to make the line search more exact in this case
        my_linesearch->set_tol_decrease(DT_(1e-4));
        my_linesearch->set_tol_curvature(DT_(1e-1));
      }
      else if(_linesearch_type == "MQCLinesearch")
      {
        my_linesearch = new_mqc_linesearch(my_functional, my_filter, false);
      }
      else
      {
        throw InternalError(__func__, __FILE__, __LINE__, "Got invalid linesearch_type: "+_linesearch_type);
      }

      my_linesearch->set_max_iter(20);
      my_linesearch->set_plot_mode(Solver::PlotMode::none);

      // Ugly way to get a preconditioner, or not
      std::shared_ptr<NLOptPrecond<typename OperatorType::VectorTypeL, FilterType>> my_precond(nullptr);

      if(_precon_type == "ApproximateHessian")
      {
        my_precond = new_approximate_hessian_precond(my_functional, my_filter);
      }
      else if(_precon_type == "Hessian")
      {
        my_precond = new_hessian_precond(my_functional, my_filter);
      }
      else if(_precon_type != "none")
      {
        throw InternalError(__func__, __FILE__, __LINE__, "Got invalid precon_type: "+_precon_type);
      }

      auto solver = new_nlsd(my_functional, my_filter, my_linesearch, false, my_precond);

      solver->init();
      solver->set_tol_abs(Math::eps<DT_>());
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot_mode(Solver::PlotMode::summary);
      solver->set_max_iter(100);

      // This will hold the solution
      auto sol = my_functional.create_vector_r();
      // We need a dummy rhs
      auto rhs = my_functional.create_vector_r();

      // Get an initial guess from the Traits class for the given function
      PointType starting_point(DT_(0));
      TestTraitsType::get_starting_point(starting_point);
      sol(0,starting_point);

      // Solve the optimisation problem
      solver->correct(sol, rhs);

      solver->done();

      // From the traits class, get the set of minimal points
      std::deque<PointType> min_points;
      TestTraitsType::get_minimal_points(min_points);

      // Check the distance between solution and minimal points
      DT_ min_dist(Math::Limits<DT_>::max());

      const auto& jt = min_points.end();
      auto it = min_points.begin();
      for(; it != jt; ++it)
      {
        DT_ dist((sol(0) - *it).norm_euclid());
        if(dist  < min_dist)
        {
          min_dist = dist;
        }
      }

      // Check if we stayed in the iteration number bound
      if(solver->get_num_iter() > _max_iter)
      {
        ++failed_checks;
        std::cout << "num_iter = "+stringify(solver->get_num_iter())+" > " +stringify(_max_iter)+" = max_iter" << std::endl;
      }

      // Check if we stayed in the functional evaluation number bound
      if(my_functional.get_num_func_evals() > _max_func_evals)
      {
        ++failed_checks;
        std::cout << "num_func_evals = "
          +stringify(my_functional.get_num_func_evals())+" > "+stringify(_max_func_evals)+" = max_func_evals" << std::endl;
      }

      // Check if we found a valid minimum
      if(min_dist > _tol)
      {
        ++failed_checks;
        std::cout << "min_dist = "+stringify_fp_sci(min_dist)+" > " +stringify_fp_sci(_tol)+" = tol" << std::endl;
      }

      solver->done();

      TEST_CHECK_MSG(failed_checks == 0,TestTraitsType::name()+" "+solver->name()+" "+my_linesearch->name()+": "
          +stringify(failed_checks)+ " failed checks");

    }
};

NLSDTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
nlsd_hb_f(float(0.5), Index(19), Index(56), "SecantLinesearch", "ApproximateHessian");

NLSDTest<Mem::Main, double, unsigned int, Analytic::Common::RosenbrockFunction>
nlsd_rb_d(double(0.75), Index(20), Index(87), "MQCLinesearch", "Hessian");

NLSDTest<Mem::Main, double, Index, Analytic::Common::HimmelblauFunction>
nlsd_rb_d_sw(double(0.6), Index(17), Index(77), "NewtonRaphsonLinesearch", "none");

#ifdef FEAT_HAVE_QUADMATH
NLSDTest<Mem::Main, __float128, Index, Analytic::Common::RosenbrockFunction>
nlsd_rb_q(__float128(0.55), Index(96), Index(158), "SecantLinesearch", "Hessian");
#endif

// Running this in CUDA is really nonsensical because all operator evaluations use Tiny::Vectors which reside in
// Mem::Main anyway, so apart from the occasional axpy nothing is done on the GPU. It should work nonetheless.
#ifdef FEAT_HAVE_CUDA
NLSDTest<Mem::CUDA, float, unsigned int, Analytic::Common::HimmelblauFunction>
nlsd_hb_f_cuda(float(0.75), Index(8), Index(35), "MQCLinesearch", "Hessian");
#endif

#ifdef FEAT_HAVE_ALGLIB
/**
 * \brief Test class template for ALGLIB's lBFGS optimiser
 *
 */
template < typename Mem_, typename DT_, typename IT_, typename Function_ >
class ALGLIBMinLBFGSTest:
  public FullTaggedTest<Mem_, DT_, IT_>
{
  public:
    typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
    typedef typename OperatorType::PointType PointType;
    typedef OptimisationTestTraits<DT_, Function_> TestTraitsType;

    typedef LAFEM::NoneFilterBlocked<Mem_, DT_, IT_, 2> FilterType;

  private:
    const DT_ _tol;
    const Index _max_iter;
    const Index _max_func_evals;

  public:
    explicit ALGLIBMinLBFGSTest(const DT_ exponent_, const Index max_iter_, const Index max_func_evals_) :
      FullTaggedTest<Mem_, DT_, IT_>("ALGLIBMinLBFGSTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _max_iter(max_iter_),
      _max_func_evals(max_func_evals_)
    {
    }

    virtual ~ALGLIBMinLBFGSTest()
    {
    }

    void run() const override
    {
      int failed_checks(0);
      // The analytic function
      Function_ my_function;
      // Create the (nonlinear) operator
      OperatorType my_functional(my_function);
      // The filter
      FilterType my_filter;

      //auto my_precond = nullptr;
      auto solver = new_alglib_minlbfgs(my_functional, my_filter);
      solver->init();
      solver->set_tol_abs(Math::eps<DT_>());
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot_mode(Solver::PlotMode::summary);
      solver->set_max_iter(100);

      // This will hold the solution
      auto sol = my_functional.create_vector_r();
      // We need a dummy rhs
      auto rhs = my_functional.create_vector_r();

      // Get an initial guess from the Traits class for the given function
      PointType starting_point(DT_(0));
      TestTraitsType::get_starting_point(starting_point);
      sol(0,starting_point);

      // Solve the optimisation problem
      solver->correct(sol, rhs);

      solver->done();

      // From the traits class, get the set of minimal points
      std::deque<PointType> min_points;
      TestTraitsType::get_minimal_points(min_points);

      // Check the distance between solution and minimal points
      DT_ min_dist(Math::Limits<DT_>::max());

      const auto& jt = min_points.end();
      auto it = min_points.begin();
      for(; it != jt; ++it)
      {
        DT_ dist((sol(0) - *it).norm_euclid());
        if(dist  < min_dist)
        {
          min_dist = dist;
        }
      }

      // Check if we stayed in the iteration number bound
      if(solver->get_num_iter() > _max_iter)
      {
        ++failed_checks;
        std::cout << "num_iter = "+stringify(solver->get_num_iter())+" > " +stringify(_max_iter)+" = max_iter" << std::endl;
      }

      // Check if we stayed in the functional evaluation number bound
      if(my_functional.get_num_func_evals() > _max_func_evals)
      {
        ++failed_checks;
        std::cout << "num_func_evals = "
          +stringify(my_functional.get_num_func_evals())+" > "+stringify(_max_func_evals)+" = max_func_evals" << std::endl;
      }

      // Check if we found a valid minimum
      if(min_dist > _tol)
      {
        ++failed_checks;
        std::cout << "min_dist = "+stringify_fp_sci(min_dist)+" > " +stringify_fp_sci(_tol)+" = tol" << std::endl;
      }

      solver->done();

      TEST_CHECK_MSG(failed_checks == 0,TestTraitsType::name()+" "+solver->name()+": "
          +stringify(failed_checks)+ " failed checks");

    }
};

ALGLIBMinLBFGSTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
alg_lbfgs_hb_f(float(0.9), Index(12), Index(17));

ALGLIBMinLBFGSTest<Mem::Main, double, unsigned int, Analytic::Common::RosenbrockFunction>
alg_lbfgs_rb_d(double(1), Index(36), Index(72));

ALGLIBMinLBFGSTest<Mem::Main, double, unsigned int, Analytic::Common::BazaraaShettyFunction>
alg_lbfgs_bs_d(double(0.48), Index(34), Index(87));

ALGLIBMinLBFGSTest<Mem::Main, double, unsigned int, Analytic::Common::GoldsteinPriceFunction>
alg_lbfgs_gp_d(double(0.65), Index(15), Index(59));

/**
 * \brief Test class template for ALGLIB's mincg optimiser
 *
 */
template < typename Mem_, typename DT_, typename IT_, typename Function_ >
class ALGLIBMinCGTest:
  public FullTaggedTest<Mem_, DT_, IT_>
{
  public:
    typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
    typedef typename OperatorType::PointType PointType;
    typedef OptimisationTestTraits<DT_, Function_> TestTraitsType;

    typedef LAFEM::NoneFilterBlocked<Mem_, DT_, IT_, 2> FilterType;

  private:
    DT_ _tol;
    const Index _max_iter;
    const Index _max_func_evals;
    NLCGDirectionUpdate _direction_update;

  public:
    explicit ALGLIBMinCGTest(const DT_ exponent_, const Index max_iter_, const Index max_func_evals_,
    NLCGDirectionUpdate direction_update_) :
      FullTaggedTest<Mem_, DT_, IT_>("ALGLIBMinCGTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _max_iter(max_iter_),
      _max_func_evals(max_func_evals_),
      _direction_update(direction_update_)
    {
    }

    virtual ~ALGLIBMinCGTest()
    {
    }

    void run() const override
    {
      int failed_checks(0);
      // The analytic function
      Function_ my_function;
      // Create the (nonlinear) operator
      OperatorType my_functional(my_function);
      // The filter
      FilterType my_filter;

      auto solver = new_alglib_mincg(my_functional, my_filter, _direction_update, false);

      solver->init();
      solver->set_tol_abs(Math::eps<DT_>());
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot_mode(Solver::PlotMode::summary);
      solver->set_max_iter(100);

      // This will hold the solution
      auto sol = my_functional.create_vector_r();
      // We need a dummy rhs
      auto rhs = my_functional.create_vector_r();

      // Get an initial guess from the Traits class for the given function
      PointType starting_point(DT_(0));
      TestTraitsType::get_starting_point(starting_point);
      sol(0,starting_point);

      // Solve the optimisation problem
      solver->correct(sol, rhs);

      solver->done();

      // From the traits class, get the set of minimal points
      std::deque<PointType> min_points;
      TestTraitsType::get_minimal_points(min_points);

      // Check the distance between solution and minimal points
      DT_ min_dist(Math::Limits<DT_>::max());

      const auto& jt = min_points.end();
      auto it = min_points.begin();
      for(; it != jt; ++it)
      {
        DT_ dist((sol(0) - *it).norm_euclid());
        if(dist  < min_dist)
        {
          min_dist = dist;
        }
      }

      // Check if we stayed in the iteration number bound
      if(solver->get_num_iter() > _max_iter)
      {
        ++failed_checks;
        std::cout << "num_iter = "+stringify(solver->get_num_iter())+" > " +stringify(_max_iter)+" = max_iter" << std::endl;
      }

      // Check if we stayed in the functional evaluation number bound
      if(my_functional.get_num_func_evals() > _max_func_evals)
      {
        ++failed_checks;
        std::cout << "num_func_evals = "
          +stringify(my_functional.get_num_func_evals())+" > "+stringify(_max_func_evals)+" = max_func_evals" << std::endl;
      }

      // Check if we found a valid minimum
      if(min_dist > _tol)
      {
        ++failed_checks;
        std::cout << "min_dist = "+stringify_fp_sci(min_dist)+" > " +stringify_fp_sci(_tol)+" = tol" << std::endl;
      }

      solver->done();

      TEST_CHECK_MSG(failed_checks == 0,TestTraitsType::name()+" "+solver->name()+": "
          +stringify(failed_checks)+ " failed checks");

    }
};

ALGLIBMinCGTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
alg_mincg_hb_f(float(0.6), Index(14), Index(60), NLCGDirectionUpdate::DaiYuan);

ALGLIBMinCGTest<Mem::Main, double, unsigned int, Analytic::Common::RosenbrockFunction>
alg_mincg_rb_d(double(0.8), Index(41), Index(118), NLCGDirectionUpdate::DYHSHybrid);

ALGLIBMinCGTest<Mem::Main, double, unsigned int, Analytic::Common::BazaraaShettyFunction>
alg_mincg_bs_d(double(0.33), Index(58), Index(194), NLCGDirectionUpdate::DaiYuan);
#endif // FEAT_HAVE_ALGLIB
