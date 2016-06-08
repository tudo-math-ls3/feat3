#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/lafem/none_filter.hpp>
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
    const String _linesearch_type;
    const String _precon_type;
    NLCGDirectionUpdate _update;

  public:
    NLCGTest(DT_ exponent_, Index max_iter_, const String& linesearch_type_, const String& precon_type_,
    NLCGDirectionUpdate update_) :
      FullTaggedTest<Mem_, DT_, IT_>("NLCGTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _max_iter(max_iter_),
      _linesearch_type(linesearch_type_),
      _precon_type(precon_type_),
      _update(update_)
    {
    }

    void run() const override
    {
      // The analytic function
      Function_ my_function;
      // Create the (nonlinear) operator
      OperatorType my_op(my_function);
      // The filter
      FilterType my_filter;

      // Create the linesearch
      std::shared_ptr<Solver::Linesearch<OperatorType, FilterType>> my_linesearch;
      if(_linesearch_type == "NewtonRaphsonLinesearch")
        my_linesearch = new_newton_raphson_linesearch(my_op, my_filter, false);
      else if(_linesearch_type == "SecantLinesearch")
        my_linesearch = new_secant_linesearch(my_op, my_filter, DT_(1e-2), false);
      else if(_linesearch_type == "StrongWolfeLinesearch")
        my_linesearch = new_strong_wolfe_linesearch(my_op, my_filter, false);
      else
        throw InternalError(__func__, __FILE__, __LINE__, "Got invalid linesearch_type: "+_linesearch_type);

      my_linesearch->set_plot(false);

      // Ugly way to get a preconditioner, or not
      std::shared_ptr<NLOptPrecond<typename OperatorType::VectorTypeL, FilterType>> my_precond(nullptr);

      if(_precon_type == "ApproximateHessian")
        my_precond = new_approximate_hessian_precond(my_op, my_filter);
      else if(_precon_type == "Hessian")
        my_precond = new_hessian_precond(my_op, my_filter);
      else if(_precon_type != "none")
        throw InternalError(__func__, __FILE__, __LINE__, "Got invalid precon_type: "+_precon_type);

      std::shared_ptr<Solver::IterativeSolver<typename OperatorType::VectorTypeR>> solver;
      solver = new_nlcg(my_op, my_filter, my_linesearch, _update, false, my_precond);

      solver->init();
      solver->set_tol_abs(Math::sqrt(Math::eps<DT_>()));
      solver->set_tol_rel(Math::sqrt(Math::eps<DT_>()));
      solver->set_plot(false);
      solver->set_max_iter(100);

      // This will hold the solution
      auto sol = my_op.create_vector_r();
      // We need a dummy rhs
      auto rhs = my_op.create_vector_r();

      // Get an initial guess from the Traits class for the given function
      PointType starting_point(DT_(0));
      TestTraitsType::get_starting_point(starting_point);
      sol(0,starting_point);

      // Solve the optimisation problem
      solver->correct(sol, rhs);

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

      // Check if we stayed in the iteration number bound
      TEST_CHECK_MSG(solver->get_num_iter() <= _max_iter, solver->name()+"_"+_precon_type+": num_iter = "+stringify(solver->get_num_iter())+" > "+stringify(_max_iter)+" = max_iter");
      // Check if we found a valid minimum
      TEST_CHECK_MSG(min_dist < _tol, solver->name()+"_"+_precon_type+": min_dist = "+stringify_fp_sci(min_dist)+" > "+stringify_fp_sci(_tol)+" = tol");

      solver->done();
    }
};

NLCGTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
nlcg_sw_hb_f(float(0.6),Index(12),"StrongWolfeLinesearch","none", NLCGDirectionUpdate::DaiYuan);

NLCGTest<Mem::Main, double, Index, Analytic::Common::RosenbrockFunction>
nlcg_sw_rb_d(double(0.6),Index(40),"StrongWolfeLinesearch","none", NLCGDirectionUpdate::DYHSHybrid);

NLCGTest<Mem::Main, double, Index, Analytic::Common::BazaraaShettyFunction>
nlcg_sw_bs_d(double(0.15),Index(25),"StrongWolfeLinesearch","none", NLCGDirectionUpdate::DaiYuan);

NLCGTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
nlcg_s_hb_d(float(0.5),Index(10),"SecantLinesearch","none", NLCGDirectionUpdate::FletcherReeves);

NLCGTest<Mem::Main, double, unsigned int, Analytic::Common::BazaraaShettyFunction>
nlcg_s_bs_d(double(0.175), Index(13), "SecantLinesearch", "none", NLCGDirectionUpdate::HagerZhang);

NLCGTest<Mem::Main, float, unsigned int, Analytic::Common::RosenbrockFunction>
nlcg_nr_rb_d(float(0.6), Index(32),"NewtonRaphsonLinesearch","Hessian", NLCGDirectionUpdate::HestenesStiefel);

NLCGTest<Mem::Main, double, Index, Analytic::Common::RosenbrockFunction>
nlcg_sw_hessian_rb_d(double(0.7), Index(25),"StrongWolfeLinesearch","Hessian", NLCGDirectionUpdate::DYHSHybrid);

NLCGTest<Mem::Main, double, Index, Analytic::Common::GoldsteinPriceFunction>
nlcg_sw_hessian_gp_d(double(0.6), Index(25),"StrongWolfeLinesearch","Hessian", NLCGDirectionUpdate::PolakRibiere);

#ifdef FEAT_HAVE_QUADMATH
NLCGTest<Mem::Main, __float128, Index, Analytic::Common::RosenbrockFunction>
nlcg_nr_rb_q(__float128(0.6), Index(33), "NewtonRaphsonLinesearch", "ApproximateHessian",
NLCGDirectionUpdate::PolakRibiere);

NLCGTest<Mem::Main, __float128, Index, Analytic::Common::BazaraaShettyFunction>
nlcg_sw_bs_q(__float128(0.175), Index(30), "StrongWolfeLinesearch", "none", NLCGDirectionUpdate::HestenesStiefel);
#endif

// Running this in CUDA is really nonsensical because all operator evaluations use Tiny::Vectors which reside in
// Mem::Main anyway, so apart from the occasional axpy nothing is done on the GPU. It should work nonetheless.
#ifdef FEAT_HAVE_CUDA
NLCGTest<Mem::CUDA, float, unsigned int, Analytic::Common::HimmelblauFunction>
nlcg_sw_hb_f_cuda(float(0.7), Index(10), "StrongWolfeLinesearch", "Hessian", NLCGDirectionUpdate::FletcherReeves);

NLCGTest<Mem::CUDA, double, unsigned int, Analytic::Common::BazaraaShettyFunction>
nlcg_s_bs_d_cuda(double(0.175), Index(13), "SecantLinesearch", "none", NLCGDirectionUpdate::HagerZhang);
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
    const String _linesearch_type;
    const String _precon_type;

  public:
    NLSDTest(DT_ exponent_, Index max_iter_, const String& linesearch_type_, const String& precon_type_) :
      FullTaggedTest<Mem_, DT_, IT_>("NLSDTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _max_iter(max_iter_),
      _linesearch_type(linesearch_type_),
      _precon_type(precon_type_)
    {
    }

    void run() const override
    {
      // The analytic function
      Function_ my_function;
      // Create the (nonlinear) operator
      OperatorType my_op(my_function);
      // The filter
      FilterType my_filter;

      // Create the linesearch
      std::shared_ptr<Solver::Linesearch<OperatorType, FilterType>> my_linesearch(nullptr);
      if(_linesearch_type == "NewtonRaphsonLinesearch")
        my_linesearch = new_newton_raphson_linesearch(my_op, my_filter, false);
      else if(_linesearch_type == "SecantLinesearch")
        my_linesearch = new_secant_linesearch(my_op, my_filter, DT_(1e-2), false);
      else if(_linesearch_type == "StrongWolfeLinesearch")
        my_linesearch = new_strong_wolfe_linesearch(my_op, my_filter, false);
      else
        throw InternalError(__func__, __FILE__, __LINE__, "Got invalid linesearch_type: "+_linesearch_type);

      my_linesearch->set_plot(false);

      // Ugly way to get a preconditioner, or not
      std::shared_ptr<NLOptPrecond<typename OperatorType::VectorTypeL, FilterType>> my_precond(nullptr);

      if(_precon_type == "ApproximateHessian")
        my_precond = new_approximate_hessian_precond(my_op, my_filter);
      else if(_precon_type == "Hessian")
        my_precond = new_hessian_precond(my_op, my_filter);
      else if(_precon_type != "none")
        throw InternalError("Got invalid precon_type: "+_precon_type);

      std::shared_ptr<Solver::IterativeSolver<typename OperatorType::VectorTypeR>> solver;
      solver = new_nlsd(my_op, my_filter, my_linesearch, false, my_precond);

      solver->init();
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot(false);
      solver->set_max_iter(100);

      // This will hold the solution
      auto sol = my_op.create_vector_r();
      // We need a dummy rhs
      auto rhs = my_op.create_vector_r();

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

      // Check if we stayed in the iteration number bound
      TEST_CHECK_MSG(solver->get_num_iter() <= _max_iter, solver->name()+"_"+_precon_type+": num_iter = "+stringify(solver->get_num_iter())+" > "+stringify(_max_iter)+" = max_iter");
      // Check if we found a valid minimum
      TEST_CHECK_MSG(min_dist < _tol, solver->name()+"_"+_precon_type+": min_dist = "+stringify_fp_sci(min_dist)+" > "+stringify_fp_sci(_tol)+" = tol");

    }
};

NLSDTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
nlsd_hb_f(float(0.5), Index(13), "SecantLinesearch", "ApproximateHessian");

NLSDTest<Mem::Main, double, unsigned int, Analytic::Common::RosenbrockFunction>
nlsd_rb_d(double(0.75), Index(20), "StrongWolfeLinesearch", "Hessian");

NLSDTest<Mem::Main, double, Index, Analytic::Common::HimmelblauFunction>
nlsd_rb_d_sw(double(0.6), Index(10), "NewtonRaphsonLinesearch", "none");

#ifdef FEAT_HAVE_QUADMATH
NLSDTest<Mem::Main, __float128, Index, Analytic::Common::RosenbrockFunction>
nlsd_rb_q(__float128(1), Index(19), "SecantLinesearch", "Hessian");
#endif

// Running this in CUDA is really nonsensical because all operator evaluations use Tiny::Vectors which reside in
// Mem::Main anyway, so apart from the occasional axpy nothing is done on the GPU. It should work nonetheless.
#ifdef FEAT_HAVE_CUDA
NLSDTest<Mem::CUDA, float, unsigned int, Analytic::Common::HimmelblauFunction>
nlsd_hb_f_cuda(float(0.9), Index(8), "StrongWolfeLinesearch", "Hessian");
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

  public:
    ALGLIBMinLBFGSTest(DT_ exponent_, Index max_iter_) :
      FullTaggedTest<Mem_, DT_, IT_>("ALGLIBMinLBFGSTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _max_iter(max_iter_)
    {
    }

    virtual ~ALGLIBMinLBFGSTest()
    {
    }

    void run() const override
    {
      // The analytic function
      Function_ my_function;
      // Create the (nonlinear) operator
      OperatorType my_op(my_function);
      // The filter
      FilterType my_filter;

      //auto my_precond = nullptr;
      auto solver = new_alglib_minlbfgs(my_op, my_filter);
      solver->init();
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot(false);
      solver->set_max_iter(100);

      // This will hold the solution
      auto sol = my_op.create_vector_r();
      // We need a dummy rhs
      auto rhs = my_op.create_vector_r();

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

      // Check if we stayed in the iteration number bound
      TEST_CHECK_MSG(solver->get_num_iter() <= _max_iter, solver->name()+": num_iter = "+stringify(solver->get_num_iter())+" > "+stringify(_max_iter)+" = max_iter");
      // Check if we found a valid minimum
      TEST_CHECK_MSG(min_dist < _tol, solver->name()+": min_dist = "+stringify_fp_sci(min_dist)+" > "+stringify_fp_sci(_tol)+" = tol");

    }
};

ALGLIBMinLBFGSTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
alg_lbfgs_hb_f(float(0.9), Index(12));

ALGLIBMinLBFGSTest<Mem::Main, double, unsigned int, Analytic::Common::RosenbrockFunction>
alg_lbfgs_rb_d(double(0.7), Index(36));

ALGLIBMinLBFGSTest<Mem::Main, double, unsigned int, Analytic::Common::BazaraaShettyFunction>
alg_lbfgs_bs_d(double(0.35), Index(29));

ALGLIBMinLBFGSTest<Mem::Main, double, unsigned int, Analytic::Common::GoldsteinPriceFunction>
alg_lbfgs_gp_d(double(0.6), Index(15));

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
    Index _max_iter;
    NLCGDirectionUpdate _direction_update;

  public:
    ALGLIBMinCGTest(DT_ exponent_, Index max_iter_, NLCGDirectionUpdate direction_update_) :
      FullTaggedTest<Mem_, DT_, IT_>("ALGLIBMinCGTest"),
      _tol(Math::pow(Math::eps<DT_>(), exponent_)),
      _max_iter(max_iter_),
      _direction_update(direction_update_)
    {
    }

    virtual ~ALGLIBMinCGTest()
    {
    }

    void run() const override
    {
      // The analytic function
      Function_ my_function;
      // Create the (nonlinear) operator
      OperatorType my_op(my_function);
      // The filter
      FilterType my_filter;

      std::shared_ptr<Solver::IterativeSolver<typename OperatorType::VectorTypeR>> solver;
      solver = new_alglib_mincg(my_op, my_filter, _direction_update, false);

      solver->init();
      solver->set_tol_rel(Math::eps<DT_>());
      solver->set_plot(false);
      solver->set_max_iter(100);

      // This will hold the solution
      auto sol = my_op.create_vector_r();
      // We need a dummy rhs
      auto rhs = my_op.create_vector_r();

      // Get an initial guess from the Traits class for the given function
      PointType starting_point(DT_(0));
      TestTraitsType::get_starting_point(starting_point);
      sol(0,starting_point);

      // Solve the optimisation problem
      std::cout << solver->correct(sol, rhs) << std::endl;

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
      TEST_CHECK_MSG(solver->get_num_iter() <= _max_iter, solver->name()+": num_iter = "+stringify(solver->get_num_iter())+" > "+stringify(_max_iter)+" = max_iter");
      // Check if we found a valid minimum
      TEST_CHECK_MSG(min_dist < _tol,solver->name()+": min_dist = "+stringify_fp_sci(min_dist)+" > "+stringify_fp_sci(_tol)+" = tol");
    }
};

ALGLIBMinCGTest<Mem::Main, float, Index, Analytic::Common::HimmelblauFunction>
alg_mincg_hb_f(float(0.6), Index(12), NLCGDirectionUpdate::DaiYuan);

ALGLIBMinCGTest<Mem::Main, double, unsigned int, Analytic::Common::RosenbrockFunction>
alg_mincg_rb_d(double(0.6), Index(40), NLCGDirectionUpdate::DYHSHybrid);

ALGLIBMinCGTest<Mem::Main, double, unsigned int, Analytic::Common::BazaraaShettyFunction>
alg_mincg_bs_d(double(0.15), Index(25), NLCGDirectionUpdate::DaiYuan);
#endif // FEAT_HAVE_ALGLIB
