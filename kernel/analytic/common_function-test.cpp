// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/common.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Analytic;

template<typename DT_>
class CommonFunctionTest :
  public FullTaggedTest<Mem::Main, DT_, Index>
{
public:
  CommonFunctionTest() :
    FullTaggedTest<Mem::Main, DT_, Index>("CommonFunctionTest")
  {
  }

  void test_par_profile_scalar() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create scalar parabolic profile function
    Analytic::Common::ParProfileScalar<DT_> pprof;
    TEST_CHECK(pprof.parse("(1 2, 3 4, 5)"));

    // evaluate endpoints and midpoint
    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(pprof, DT_(1), DT_(2)), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(pprof, DT_(3), DT_(4)), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(pprof, DT_(2), DT_(3)), DT_(5), tol);
  }

  void test_par_profile_vector() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create vector parabolic profile function
    Analytic::Common::ParProfileVector<DT_> pprof;
    TEST_CHECK(pprof.parse("(1 2, 3 4, 5)"));

    // evaluate endpoints and midpoint
    auto v_0 = Analytic::eval_value_x(pprof, DT_(1), DT_(2));
    auto v_1 = Analytic::eval_value_x(pprof, DT_(3), DT_(4));
    auto v_c = Analytic::eval_value_x(pprof, DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(v_0[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v_0[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v_1[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v_1[1], DT_(0), tol);
    const DT_ x_c = DT_(5) * Math::sqrt(DT_(0.5));
    TEST_CHECK_EQUAL_WITHIN_EPS(v_c[0], +x_c, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(v_c[1], -x_c, tol);
  }

  void test_distance_function() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 2d DistanceFunction with origin in (0,2).
    typename Analytic::Common::DistanceFunction<2, DT_>::PointType orig2d;
    orig2d[0] = DT_(0);
    orig2d[1] = DT_(2);
    Analytic::Common::DistanceFunction<2, DT_> dfunc2d(orig2d);

    // evaluate in (4, -1)
    auto val2d = Analytic::eval_value_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(val2d, DT_(5), tol);

    auto grad2d = Analytic::eval_gradient_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad2d[0], DT_( 4)/DT_(5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad2d[1], DT_(-3)/DT_(5), tol);

    auto hess2d = Analytic::eval_hessian_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[0][0], DT_( 9)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[0][1], DT_(12)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[1][0], DT_(12)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[1][1], DT_(16)/DT_(125), tol);

    // create 3d DistanceFunction with origin in (2,2,1).
    typename Analytic::Common::DistanceFunction<3, DT_>::PointType orig3d;
    orig3d[0] = DT_(2);
    orig3d[1] = DT_(2);
    orig3d[2] = DT_(1);
    Analytic::Common::DistanceFunction<3, DT_> dfunc3d(orig3d);

    // evaluate in (4, 0, 2)
    auto val3d = Analytic::eval_value_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val3d, DT_(3), tol);

    auto grad3d = Analytic::eval_gradient_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[0], DT_( 2)/DT_(3), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[1], DT_(-2)/DT_(3), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[2], DT_( 1)/DT_(3), tol);

    auto hess3d = Analytic::eval_hessian_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][0], DT_( 5)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][1], DT_( 4)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][2], DT_(-2)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][0], DT_( 4)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][1], DT_( 5)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][2], DT_( 2)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][0], DT_(-2)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][1], DT_( 2)/DT_(27), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][2], DT_( 8)/DT_(27), tol);
  }

  void test_distance_function_sd() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 2d DistanceFunction with origin in (0,2).
    typename Analytic::Common::DistanceFunctionSD<2, DT_>::PointType orig2d;
    orig2d[0] = DT_(0);
    orig2d[1] = DT_(2);
    Analytic::Common::DistanceFunctionSD<2, DT_> dfunc2d(orig2d, DT_(0.25), 10);

    // evaluate in (4, -1)
    auto val2d = Analytic::eval_value_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(val2d, DT_(50.25), tol);

    auto grad2d = Analytic::eval_gradient_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad2d[0], DT_( 8), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad2d[1], DT_(-6), tol);

    auto hess2d = Analytic::eval_hessian_x(dfunc2d, DT_(4), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[0][0], DT_( 90)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[0][1], DT_(120)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[1][0], DT_(120)/DT_(125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess2d[1][1], DT_(160)/DT_(125), tol);

    // create 3d DistanceFunction with origin in (2,2,1).
    typename Analytic::Common::DistanceFunctionSD<3, DT_>::PointType orig3d;
    orig3d[0] = DT_(2);
    orig3d[1] = DT_(2);
    orig3d[2] = DT_(1);
    Analytic::Common::DistanceFunctionSD<3, DT_> dfunc3d(orig3d, DT_(4), DT_(3));

    // evaluate in (4, 0, 2)
    auto val3d = Analytic::eval_value_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val3d, DT_(13), tol);

    auto grad3d = Analytic::eval_gradient_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[0], DT_( 2), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[1], DT_(-2), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[2], DT_( 1), tol);

    auto hess3d = Analytic::eval_hessian_x(dfunc3d, DT_(4), DT_(0), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][0], DT_( 5)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][1], DT_( 4)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[0][2], DT_(-2)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][0], DT_( 4)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][1], DT_( 5)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[1][2], DT_( 2)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][0], DT_(-2)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][1], DT_( 2)/DT_(9), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess3d[2][2], DT_( 8)/DT_(9), tol);
  }


  void test_plain_distance_function_sd() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 3d DistanceFunction with origin in (2,2,1).
    typename Analytic::Common::PlaneDistanceFunctionSD<1, 3, DT_>::PointType orig3d;
    orig3d[0] = DT_(2);
    orig3d[1] = DT_(2);
    orig3d[2] = DT_(1);
    Analytic::Common::PlaneDistanceFunctionSD<1, 3, DT_> dfunc3d(orig3d, DT_(3));

    auto val3d = Analytic::eval_value_x(dfunc3d, DT_(4), DT_(7), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val3d, DT_(15), tol);

    auto grad3d = Analytic::eval_gradient_x(dfunc3d, DT_(4), DT_(-3), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[1], DT_(-3), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad3d[2], DT_(0), tol);
  }


  void test_min_function() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 3d DistanceFunction with origin in (2,2,1).
    typename Analytic::Common::DistanceFunction<3, DT_>::PointType orig3d;
    orig3d[0] = DT_(2);
    orig3d[1] = DT_(2);
    orig3d[2] = DT_(1);
    Analytic::Common::DistanceFunction<3, DT_> func1(orig3d);
    Analytic::Common::ConstantFunction<3, DT_> func2(DT_(5));
    Analytic::Common::MinOfTwoFunctions<decltype(func1), decltype(func2)> min_func(func1, func2);

    // evaluate in (4, 0, 2) -> func1
    auto val  = Analytic::eval_value_x   (min_func, DT_(4), DT_(0), DT_(2));
    auto grad = Analytic::eval_gradient_x(min_func, DT_(4), DT_(0), DT_(2));
    auto hess = Analytic::eval_hessian_x (min_func, DT_(4), DT_(0), DT_(2));
    auto ref_val  = Analytic::eval_value_x   (func1, DT_(4), DT_(0), DT_(2));
    auto ref_grad = Analytic::eval_gradient_x(func1, DT_(4), DT_(0), DT_(2));
    auto ref_hess = Analytic::eval_hessian_x (func1, DT_(4), DT_(0), DT_(2));
    auto diff_val  = val - ref_val;
    auto diff_grad = (grad - ref_grad).norm_euclid();
    auto diff_hess = (hess - ref_hess).norm_frobenius();

    TEST_CHECK_EQUAL_WITHIN_EPS(diff_val, DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(diff_grad, DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(diff_hess, DT_(0), tol);

    // evaluate in (8, 8, 10) -> func2
    val  = Analytic::eval_value_x   (min_func, DT_(8), DT_(8), DT_(10));
    grad = Analytic::eval_gradient_x(min_func, DT_(8), DT_(8), DT_(10));
    hess = Analytic::eval_hessian_x (min_func, DT_(8), DT_(8), DT_(10));
    ref_val  = Analytic::eval_value_x   (func2, DT_(8), DT_(8), DT_(10));
    ref_grad = Analytic::eval_gradient_x(func2, DT_(8), DT_(8), DT_(10));
    ref_hess = Analytic::eval_hessian_x (func2, DT_(8), DT_(8), DT_(10));
    diff_val  = val - ref_val;
    diff_grad = (grad - ref_grad).norm_euclid();
    diff_hess = (hess - ref_hess).norm_frobenius();

    TEST_CHECK_EQUAL_WITHIN_EPS(diff_val, DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(diff_grad, DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(diff_hess, DT_(0), tol);
  }


  virtual void run() const override
  {
    test_par_profile_scalar();
    test_par_profile_vector();
    test_distance_function();
    test_distance_function_sd();
    test_plain_distance_function_sd();
    test_min_function();
  }
};

CommonFunctionTest<double> common_function_test_double;
