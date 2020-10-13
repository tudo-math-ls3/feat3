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

  void test_distance_function_2d() const
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
  }

  void test_distance_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

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

  void test_distance_function_sd_2d() const
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
  }

  void test_distance_function_sd_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

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

  void test_plane_distance_function_sd() const
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

  void test_sine_bubble_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // some useful constants
    const DT_ pi = Math::pi<DT_>();
    const DT_ s4 = Math::sin(DT_(0.25 )*pi); // = sin(pi/4)
    const DT_ s8 = Math::sin(DT_(0.125)*pi); // = sin(pi/8)
    const DT_ c4 = Math::cos(DT_(0.25 )*pi); // = cos(pi/4)
    const DT_ c8 = Math::cos(DT_(0.125)*pi); // = cos(pi/8)

    // create sine-bubble-function object
    Analytic::Common::SineBubbleFunction<2> func;

    // evaluate function value in point (1/4, 1/8)
    DT_ val = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, s4*s8, tol);

    // evaluate gradient in point (1/4, 1/8)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], pi*c4*s8, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], pi*s4*c8, tol);

    // evaluate hessian in point (1/4, 1/8)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], -pi*pi*s4*s8, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], +pi*pi*c4*c8, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], +pi*pi*c4*c8, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], -pi*pi*s4*s8, tol);
  }

  void test_sine_bubble_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // some useful constants
    const DT_ pi  = Math::pi<DT_>();
    const DT_ s4  = Math::sin(DT_(0.25  ) * pi); // = sin(pi/4)
    const DT_ s8  = Math::sin(DT_(0.125 ) * pi); // = sin(pi/8)
    const DT_ s16 = Math::sin(DT_(0.0625) * pi); // = sin(pi/16)
    const DT_ c4  = Math::cos(DT_(0.25  ) * pi); // = cos(pi/4)
    const DT_ c8  = Math::cos(DT_(0.125 ) * pi); // = cos(pi/8)
    const DT_ c16 = Math::cos(DT_(0.0625) * pi); // = cos(pi/16)

    // create sine-bubble-function object
    Analytic::Common::SineBubbleFunction<3> func;

    // evaluate function value in point (1/4, 1/8, 1/16)
    DT_ val = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125), DT_(0.0625));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, s4 * s8 * s16, tol);

    // evaluate gradient in point (1/4, 1/8, 1/16)
    Tiny::Vector<DT_, 3> grad = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125), DT_(0.0625));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], pi * c4 * s8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], pi * s4 * c8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2], pi * s4 * s8 * c16, tol);

    // evaluate hessian in point (1/4, 1/8, 1/16)
    Tiny::Matrix<DT_, 3, 3> hess = Analytic::eval_hessian_x(func, DT_(0.25), DT_(0.125), DT_(0.0625));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], -pi * pi * s4 * s8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], +pi * pi * c4 * c8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][2], +pi * pi * c4 * s8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], +pi * pi * c4 * c8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], -pi * pi * s4 * s8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][2], +pi * pi * s4 * c8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][0], +pi * pi * c4 * s8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][1], +pi * pi * s4 * c8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][2], -pi * pi * s4 * s8 * s16, tol);
  }

  void test_cosine_wave_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // some useful constants
    const DT_ pi = Math::pi<DT_>();
    const DT_ s4 = Math::sin(DT_(0.25 ) * pi); // = sin(pi/4)
    const DT_ s8 = Math::sin(DT_(0.125) * pi); // = sin(pi/8)
    const DT_ c4 = Math::cos(DT_(0.25 ) * pi); // = cos(pi/4)
    const DT_ c8 = Math::cos(DT_(0.125) * pi); // = cos(pi/8)

    // create cosine-wave-function object
    Analytic::Common::CosineWaveFunction<2> func;

    // evaluate function value in point (1/4, 1/8)
    DT_ val = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, c4 * c8, tol);

    // evaluate gradient in point (1/4, 1/8)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], -pi * s4 * c8, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], -pi * c4 * s8, tol);

    // evaluate hessian in point (1/4, 1/8)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], -pi * pi * c4 * c8, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], +pi * pi * s4 * s8, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], +pi * pi * s4 * s8, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], -pi * pi * c4 * c8, tol);
  }

  void test_cosine_wave_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // some useful constants
    const DT_ pi  = Math::pi<DT_>();
    const DT_ s4  = Math::sin(DT_(0.25  ) * pi); // = sin(pi/4)
    const DT_ s8  = Math::sin(DT_(0.125 ) * pi); // = sin(pi/8)
    const DT_ s16 = Math::sin(DT_(0.0625) * pi); // = sin(pi/16)
    const DT_ c4  = Math::cos(DT_(0.25  ) * pi); // = cos(pi/4)
    const DT_ c8  = Math::cos(DT_(0.125 ) * pi); // = cos(pi/8)
    const DT_ c16 = Math::cos(DT_(0.0625) * pi); // = cos(pi/16)

    // create cosine-wave-function object
    Analytic::Common::CosineWaveFunction<3> func;

    // evaluate function value in point (1/4, 1/8, 1/16)
    DT_ val = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125), DT_(0.0625));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, c4 * c8 * c16, tol);

    // evaluate gradient in point (1/4, 1/8, 1/16)
    Tiny::Vector<DT_, 3> grad = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125), DT_(0.0625));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], -pi * s4 * c8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], -pi * c4 * s8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2], -pi * c4 * c8 * s16, tol);

    // evaluate hessian in point (1/4, 1/8, 1/16)
    Tiny::Matrix<DT_, 3, 3> hess = Analytic::eval_hessian_x(func, DT_(0.25), DT_(0.125), DT_(0.0625));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], -pi * pi * c4 * c8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], +pi * pi * s4 * s8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][2], +pi * pi * s4 * c8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], +pi * pi * s4 * s8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], -pi * pi * c4 * c8 * c16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][2], +pi * pi * c4 * s8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][0], +pi * pi * s4 * c8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][1], +pi * pi * c4 * s8 * s16, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][2], -pi * pi * c4 * c8 * c16, tol);
  }

  void test_q2_bubble_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create q2-bubble-function object
    Analytic::Common::Q2BubbleFunction<2> func;

    //name the constants
    const DT_ x = DT_(0.25 );
    const DT_ y = DT_(0.125);

    // evaluate function value in point (1/4, 1/8)
    DT_ val = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(16) * x * (DT_(1) - x) * y * (DT_(1) - y) , tol);

    // evaluate gradient in point (1/4, 1/8)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(16) * (DT_(2) * x - DT_(1)) * (y - DT_(1)) * y, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], DT_(16) * (x - DT_(1)) * x * (DT_(2) * y - DT_(1)), tol);

    // evaluate hessian in point (1/4, 1/8)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(32) * (y - DT_(1)) * y, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_(16) * (DT_(2) * x - DT_(1)) * (DT_(2) * y - DT_(1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_(16) * (DT_(2) * x - DT_(1)) * (DT_(2) * y - DT_(1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_(32) * (x - DT_(1)) * x, tol);
  }

  void test_q2_bubble_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create q2-bubble-function object
    Analytic::Common::Q2BubbleFunction<3> func;

    // name the constants
    const DT_ x = DT_(0.5  );
    const DT_ y = DT_(0.25 );
    const DT_ z = DT_(0.125);

    // evaluate function value in point (1/2, 1/4, 1/8)
    DT_ val = Analytic::eval_value_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(64) * x * (DT_(1) - x) * y * (DT_(1) - y) * z * (DT_(1) - z), tol);

    // evaluate gradient in point (1/2, 1/4, 1/8)
    Tiny::Vector<DT_, 3> grad = Analytic::eval_gradient_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], -DT_(64) * (DT_(2) * x - DT_(1)) * (y - DT_(1)) * y * (z - DT_(1)) * z, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], -DT_(64) * (x - DT_(1)) * x * (DT_(2) * y - DT_(1)) * (z - DT_(1)) * z, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2], -DT_(64) * (x - DT_(1)) * x * (y - DT_(1)) * y * (DT_(2) * z - DT_(1)), tol);

    // evaluate hessian in point (1/2, 1/4, 1/8)
    Tiny::Matrix<DT_, 3, 3> hess = Analytic::eval_hessian_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], -DT_(128) * (y - DT_(1)) * y * (z - DT_(1)) * z, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_(64) * (DT_(2) * x - DT_(1)) * (DT_(1) - DT_(2) * y) * (z - DT_(1)) * z, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][2], DT_(64) * (DT_(2) * x - DT_(1)) * (y - DT_(1)) * y * (DT_(1) - DT_(2) * z), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_(64) * (DT_(1) - DT_(2) * x) * (DT_(2) * y - DT_(1)) * (z - DT_(1)) * z, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], -DT_(128) * (x - DT_(1)) * x * (z - DT_(1)) * z, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][2], DT_(64) * (x - DT_(1)) * x * (DT_(2) * y - DT_(1)) * (DT_(1) - DT_(2) * z), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][0], DT_(64) * (DT_(1) - DT_(2) * x) * (y - DT_(1)) * y * (DT_(2) * z - DT_(1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][1], DT_(64) * (x - DT_(1)) * x * (DT_(1) - DT_(2) * y) * (DT_(2) * z - DT_(1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][2], -DT_(128) * (x - DT_(1)) * x * (y - DT_(1)) * y, tol);
  }

  void test_exp_bubble_function_1d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create exp-bubble-function object
    Analytic::Common::ExpBubbleFunction<1> func;

    // some useful constants
    const DT_ x = DT_(0.25);
    const DT_ e_1     = Math::exp( DT_(1)); // = exp( 1)
    const DT_ e_neg_1 = Math::exp(-DT_(1)); // = exp(-1)
    const DT_ e_x = Math::exp(-DT_(4) * (x - DT_(1)) * x); // = exp(-4(x-1)x)

    // function value, gradient and hessian in 1/4
    const DT_ u_x = (Math::exp(-Math::pow(DT_(2) * x - DT_(1), DT_(2))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(x)
    const DT_ u_x_grad = -(Math::exp(-DT_(4) * (x - DT_(1)) * x) * (DT_(8) * x - DT_(4))) / (e_1 - DT_(1)); // = u'(x)
    const DT_ u_x_hess = (DT_(8) * e_x * (DT_(8) * DT_(0.0625) - DT_(8) * x + DT_(1))) / (e_1 - DT_(1)); // = u''(x)

    // evaluate function value in point (1/4)
    DT_ val = Analytic::eval_value_x(func, DT_(0.25));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, u_x, tol);

    // evaluate gradient in point (1/4)
    Tiny::Vector<DT_, 1> grad = Analytic::eval_gradient_x(func, DT_(0.25));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], u_x_grad, tol);

    // evaluate hessian in point (1/4)
    Tiny::Matrix<DT_, 1, 1> hess = Analytic::eval_hessian_x(func, DT_(0.25));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], u_x_hess, tol);
  }

  void test_exp_bubble_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create exp-bubble-function object
    Analytic::Common::ExpBubbleFunction<2> func;

    // some useful constants
    const DT_ x = DT_(0.25 );
    const DT_ y = DT_(0.125);
    const DT_ e_1     = Math::exp( DT_(1)); // = exp( 1)
    const DT_ e_neg_1 = Math::exp(-DT_(1)); // = exp(-1)
    const DT_ e_x = Math::exp(-DT_(4) * (x - DT_(1)) * x); // = exp(-4(x-1)x)
    const DT_ e_y = Math::exp(-DT_(4) * (y - DT_(1)) * y); // = exp(-4(y-1)y)

    // function value, gradient and hessian in 1/4 and 1/8
    const DT_ u_x = (Math::exp(-Math::pow(DT_(2) * x - DT_(1), DT_(2))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(x)
    const DT_ u_y = (Math::exp(-Math::pow(DT_(2) * y - DT_(1), DT_(2))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(y)
    const DT_ u_x_grad = -(Math::exp(-DT_(4) * (x - DT_(1)) * x) * (DT_(8) * x - DT_(4))) / (e_1 - DT_(1)); // = u'(x)
    const DT_ u_y_grad = -(Math::exp(-DT_(4) * (y - DT_(1)) * y) * (DT_(8) * y - DT_(4))) / (e_1 - DT_(1)); // = u'(y)
    const DT_ u_x_hess = (DT_(8) * e_x * (DT_(8) * DT_(0.0625  ) - DT_(8) * x + DT_(1))) / (e_1 - DT_(1)); // = u''(x)
    const DT_ u_y_hess = (DT_(8) * e_y * (DT_(8) * DT_(0.015625) - DT_(8) * y + DT_(1))) / (e_1 - DT_(1)); // = u''(y)

    // evaluate function value in point (1/4, 1/8)
    DT_ val = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, u_x * u_y , tol);

    // evaluate gradient in point (1/4, 1/8)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], u_x_grad * u_y      , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], u_x      * u_y_grad , tol);

    // evaluate hessian in point (1/4, 1/8)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], u_x_hess * u_y     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], u_x_grad * u_y_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], u_x_grad * u_y_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], u_x      * u_y_hess, tol);
  }

  void test_exp_bubble_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create exp-bubble-function object
    Analytic::Common::ExpBubbleFunction<3> func;

    // some useful constants
    const DT_ x = DT_(0.5  );
    const DT_ y = DT_(0.25 );
    const DT_ z = DT_(0.125);
    const DT_ e_1     = Math::exp( DT_(1)); // = exp( 1)
    const DT_ e_neg_1 = Math::exp(-DT_(1)); // = exp(-1)
    const DT_ e_x = Math::exp(-DT_(4) * (x - DT_(1)) * x); // = exp(-4(x-1)x) for the calculation of the hessian
    const DT_ e_y = Math::exp(-DT_(4) * (y - DT_(1)) * y); // = exp(-4(y-1)y)
    const DT_ e_z = Math::exp(-DT_(4) * (z - DT_(1)) * z); // = exp(-4(z-1)z)

    // function value, gradient and hessian in 1/2, 1/4 and 1/8
    const DT_ u_x = (Math::exp(-Math::pow(DT_(2) * x - DT_(1), DT_(2))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(x)
    const DT_ u_y = (Math::exp(-Math::pow(DT_(2) * y - DT_(1), DT_(2))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(y)
    const DT_ u_z = (Math::exp(-Math::pow(DT_(2) * z - DT_(1), DT_(2))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(z)
    const DT_ u_x_grad = -(Math::exp(-DT_(4) * (x - DT_(1)) * x) * (DT_(8) * x - DT_(4))) / (e_1 - DT_(1)); // = u'(x)
    const DT_ u_y_grad = -(Math::exp(-DT_(4) * (y - DT_(1)) * y) * (DT_(8) * y - DT_(4))) / (e_1 - DT_(1)); // = u'(y)
    const DT_ u_z_grad = -(Math::exp(-DT_(4) * (z - DT_(1)) * z) * (DT_(8) * z - DT_(4))) / (e_1 - DT_(1)); // = u'(z)
    const DT_ u_x_hess = (DT_(8) * e_x * (DT_(8) * DT_(0.25    ) - DT_(8) * x + DT_(1))) / (e_1 - DT_(1)); // = u''(x)
    const DT_ u_y_hess = (DT_(8) * e_y * (DT_(8) * DT_(0.0625  ) - DT_(8) * y + DT_(1))) / (e_1 - DT_(1)); // = u''(y)
    const DT_ u_z_hess = (DT_(8) * e_z * (DT_(8) * DT_(0.015625) - DT_(8) * z + DT_(1))) / (e_1 - DT_(1)); // = u''(z)

    // evaluate function value in point (1/2, 1/4, 1/8)
    DT_ val = Analytic::eval_value_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, u_x * u_y * u_z, tol);

    // evaluate gradient in point (1/2, 1/4, 1/8)
    Tiny::Vector<DT_, 3> grad = Analytic::eval_gradient_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], u_x_grad * u_y      * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], u_x      * u_y_grad * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2], u_x      * u_y      * u_z_grad, tol);

    // evaluate hessian in point (1/2, 1/4, 1/8)
    Tiny::Matrix<DT_, 3, 3> hess = Analytic::eval_hessian_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], u_x_hess * u_y      * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], u_x_grad * u_y_grad * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][2], u_x_grad * u_y      * u_z_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], u_x_grad * u_y_grad * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], u_x      * u_y_hess * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][2], u_x      * u_y_grad * u_z_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][0], u_x_grad * u_y      * u_z_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][1], u_x      * u_y_grad * u_z_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][2], u_x      * u_y      * u_z_hess, tol);
  }

  void test_exp_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create exp-function object
    Analytic::Common::ExpFunction<2> func;

    // some useful constants
    const DT_ x = DT_(0.25 );
    const DT_ y = DT_(0.125);
    const DT_ e_10 = Math::exp(DT_(10)); // = exp(10)
    const DT_ e_x = Math::exp(DT_(10) * Math::sqr(x)); // = exp(10*x^2)
    const DT_ e_y = Math::exp(DT_(10) * Math::sqr(y)); // = exp(10*y^2)

    // function value, gradient and hessian in 1/4 and 1/8
    const DT_ u_x = (e_10 - e_x) / (e_10 - DT_(1)); // = u(x)
    const DT_ u_y = (e_10 - e_y) / (e_10 - DT_(1)); // = u(y)
    const DT_ u_x_grad = -(DT_(20) * e_x * x) / (e_10 - DT_(1)); // = u'(x)
    const DT_ u_y_grad = -(DT_(20) * e_y * y) / (e_10 - DT_(1)); // = u'(y)
    const DT_ u_x_hess = -(DT_(20) * e_x * (DT_(20) * Math::sqr(x) + 1)) / (e_10 - DT_(1)); // = u''(x)
    const DT_ u_y_hess = -(DT_(20) * e_y * (DT_(20) * Math::sqr(y) + 1)) / (e_10 - DT_(1)); // = u''(y)

    // evaluate function value in point (1/4, 1/8)
    DT_ val = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, u_x * u_y, tol);

    // evaluate gradient in point (1/4, 1/8)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], u_x_grad * u_y     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], u_x      * u_y_grad, tol);

    // evaluate hessian in point (1/4, 1/8)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], u_x_hess * u_y     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], u_x_grad * u_y_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], u_x_grad * u_y_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], u_x      * u_y_hess, tol);
  }

  void test_exp_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create exp-function object
    Analytic::Common::ExpFunction<3> func;

    // some useful constants
    const DT_ x = DT_(0.5  );
    const DT_ y = DT_(0.25 );
    const DT_ z = DT_(0.125);
    const DT_ e_10 = Math::exp(DT_(10)); // = exp(10)
    const DT_ e_x = Math::exp(DT_(10) * Math::sqr(x)); // = exp(10*x^2)
    const DT_ e_y = Math::exp(DT_(10) * Math::sqr(y)); // = exp(10*y^2)
    const DT_ e_z = Math::exp(DT_(10) * Math::sqr(z)); // = exp(10*z^2)

    // function value, gradient and hessian in 1/2, 1/4 and 1/8
    const DT_ u_x = (e_10 - e_x) / (e_10 - DT_(1)); // = u(x)
    const DT_ u_y = (e_10 - e_y) / (e_10 - DT_(1)); // = u(y)
    const DT_ u_z = (e_10 - e_z) / (e_10 - DT_(1)); // = u(z)
    const DT_ u_x_grad = -(DT_(20) * e_x * x) / (e_10 - DT_(1)); // = u'(x)
    const DT_ u_y_grad = -(DT_(20) * e_y * y) / (e_10 - DT_(1)); // = u'(y)
    const DT_ u_z_grad = -(DT_(20) * e_z * z) / (e_10 - DT_(1)); // = u'(z)
    const DT_ u_x_hess = -(DT_(20) * e_x * (DT_(20) * Math::sqr(x) + 1)) / (e_10 - DT_(1)); // = u''(x)
    const DT_ u_y_hess = -(DT_(20) * e_y * (DT_(20) * Math::sqr(y) + 1)) / (e_10 - DT_(1)); // = u''(y)
    const DT_ u_z_hess = -(DT_(20) * e_z * (DT_(20) * Math::sqr(z) + 1)) / (e_10 - DT_(1)); // = u''(z)

    // evaluate function value in point (1/2, 1/4, 1/8)
    DT_ val = Analytic::eval_value_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, u_x * u_y * u_z, tol);

    // evaluate gradient in point (1/2, 1/4, 1/8)
    Tiny::Vector<DT_, 3> grad = Analytic::eval_gradient_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], u_x_grad * u_y      * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], u_x      * u_y_grad * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2], u_x      * u_y      * u_z_grad, tol);

    // evaluate hessian in point (1/2, 1/4, 1/8)
    Tiny::Matrix<DT_, 3, 3> hess = Analytic::eval_hessian_x(func, DT_(0.5), DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], u_x_hess * u_y      * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], u_x_grad * u_y_grad * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][2], u_x_grad * u_y      * u_z_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], u_x_grad * u_y_grad * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], u_x      * u_y_hess * u_z     , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][2], u_x      * u_y_grad * u_z_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][0], u_x_grad * u_y      * u_z_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][1], u_x      * u_y_grad * u_z_grad, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][2], u_x      * u_y      * u_z_hess, tol);
  }

  void test_heaviside_function_1d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create heaviside-function object
    Analytic::Common::HeavisideFunction<1> func;

    // evaluate function value in points (2), (0), (-2)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(1), tol);

    DT_ val_2 = Analytic::eval_value_x(func, DT_(0));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(1), tol);

    DT_ val_3 = Analytic::eval_value_x(func, DT_(-2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_3, DT_(0), tol);
  }

  void test_heaviside_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create heaviside-function object
    Analytic::Common::HeavisideFunction<2> func;

    // evaluate function value in points (2,2), (-2,2), (-2,-2)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(2), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(1), tol);

    DT_ val_2 = Analytic::eval_value_x(func, DT_(-2), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(0), tol);

    DT_ val_3 = Analytic::eval_value_x(func, DT_(-2), DT_(-2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_3, DT_(0), tol);
  }

  void test_heaviside_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create heaviside-function object
    Analytic::Common::HeavisideFunction<3> func;

    // evaluate function value in points (2,2,2), (-2,2,2), (-2,-2,-2)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(2), DT_(2), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(1), tol);

    DT_ val_2 = Analytic::eval_value_x(func, DT_(-2), DT_(2), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(0), tol);

    DT_ val_3 = Analytic::eval_value_x(func, DT_(-2), DT_(-2), DT_(-2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_3, DT_(0), tol);
  }

  void test_heaviside_reg_function_1d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create heaviside-reg-function object
    Analytic::Common::HeavisideRegFunction<1> func;

    // some useful constants
    const DT_ ch2 = Math::cosh(DT_(2));
    const DT_ ch0 = Math::cosh(DT_(0));

    // evaluate function value in points (2), (0), (-2)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(2)*(ch2-DT_(1)), tol);

    DT_ val_2 = Analytic::eval_value_x(func, DT_(0));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(2) * (ch0 - DT_(1)), tol);

    DT_ val_3 = Analytic::eval_value_x(func, DT_(-2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_3, DT_(0), tol);
  }

  void test_heaviside_reg_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create heaviside-reg-function object
    Analytic::Common::HeavisideRegFunction<2> func;

    // some useful constants
    const DT_ ch2 = Math::cosh(DT_(2));

    // evaluate function value in points (2,2), (2,-2), (-2,-2)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(2), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(2) * (ch2 - DT_(1)) * DT_(2) * (ch2 - DT_(1)), tol);

    DT_ val_2 = Analytic::eval_value_x(func, DT_(2), DT_(-2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(0), tol);

    DT_ val_3 = Analytic::eval_value_x(func, DT_(-2), DT_(-2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_3, DT_(0), tol);
  }

  void test_heaviside_reg_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create heaviside-reg-function object
    Analytic::Common::HeavisideRegFunction<3> func;

    // some useful constants
    const DT_ ch2 = Math::cosh(DT_(2));

    // evaluate function value in points (2,2,2), (2,-2,2), (-2,-2,-2)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(2), DT_(2), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(8) * (ch2 - DT_(1)) * (ch2 - DT_(1)) * (ch2 - DT_(1)), tol);

    DT_ val_2 = Analytic::eval_value_x(func, DT_(2), DT_(-2), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(0), tol);

    DT_ val_3 = Analytic::eval_value_x(func, DT_(-2), DT_(-2), DT_(-2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_3, DT_(0), tol);
  }

  void test_goldstein_price_function() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create goldstein-price-function object
    Analytic::Common::GoldsteinPriceFunction func;

    // evaluate function value in point (1,2)
    DT_ val = Analytic::eval_value_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(137150), tol);

    // evaluate gradient in point (1,2)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(-15840), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], DT_(530160), tol);

    // evaluate hessian in point (1,2)
    Tiny::Matrix<DT_, 2,2> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(- 31680), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_( 127320), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_( 127320), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_(1904820), tol);

    // global minimum in (0,-1) with func(0,-1)=3
    DT_ min = Analytic::eval_value_x(func, DT_(0), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(min, DT_(3), tol);

    // evaluate gradient in global minimum
    Tiny::Vector<DT_, 2> grad_min = Analytic::eval_gradient_x(func, DT_(0), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_min[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_min[1], DT_(0), tol);

    // evaluate hessian in global minimum
    Tiny::Matrix<DT_, 2, 2> hess_min = Analytic::eval_hessian_x(func, DT_(0), DT_(-1));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_min[0][0], DT_( 504), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_min[0][1], DT_(-216), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_min[1][0], DT_(-216), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_min[1][1], DT_( 864), tol);
  }

  void test_bazaraa_shetty_function() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create bazaraa-shetty-function object
    Analytic::Common::BazaraaShettyFunction func;

    // evaluate function value in point (1,2)
    DT_ val = Analytic::eval_value_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(10), tol);

    // evaluate gradient in point (1,2)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(-10), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], DT_( 12), tol);

    // evaluate hessian in point (1,2)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_( 14), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_(-4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_(-4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_( 8), tol);

    // global minimum in (2,1) with func(2,1)=0
    DT_ min = Analytic::eval_value_x(func, DT_(2), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(min, DT_(0), tol);

    // evaluate gradient in global minimum
    Tiny::Vector<DT_, 2> grad_min = Analytic::eval_gradient_x(func, DT_(2), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_min[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_min[1], DT_(0), tol);

    // evaluate hessian in global minimum, because it should be singular
    Tiny::Matrix<DT_, 2, 2> hess_min = Analytic::eval_hessian_x(func, DT_(2), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_min[0][0], DT_( 2), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_min[0][1], DT_(-4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_min[1][0], DT_(-4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_min[1][1], DT_( 8), tol);
  }

  void test_himmelblau_function() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create himmelblau-function object
    Analytic::Common::HimmelblauFunction func;

    // evaluate function value in point (1,2)
    DT_ val = Analytic::eval_value_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(68), tol);

    // evaluate gradient in point (1,2)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(-36), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], DT_(-32), tol);

    // evaluate hessian in point (1,2)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(-22), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_( 12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_( 12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_( 26), tol);

    // test local minima
    DT_ min_1 = Analytic::eval_value_x(func, DT_(-3.77931025337774689189076584129), DT_(-3.28318599128616941226600051437));
    TEST_CHECK_EQUAL_WITHIN_EPS(min_1, DT_(0), tol);

    DT_ min_2 = Analytic::eval_value_x(func, DT_(-2.80511808695274485305357239809), DT_(3.13131251825057296580430072341));
    TEST_CHECK_EQUAL_WITHIN_EPS(min_2, DT_(0), tol);

    DT_ min_3 = Analytic::eval_value_x(func, DT_(3), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(min_3, DT_(0), tol);

    DT_ min_4 = Analytic::eval_value_x(func, DT_(3.58442834033049174494433823938), DT_(-1.84812652696440355353830020904));
    TEST_CHECK_EQUAL_WITHIN_EPS(min_4, DT_(0), tol);

    Tiny::Vector<DT_, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(-3.77931025337774689189076584129), DT_(-3.28318599128616941226600051437));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(0), tol);

    Tiny::Vector<DT_, 2> grad_2 = Analytic::eval_gradient_x(func, DT_(-2.80511808695274485305357239809), DT_(3.13131251825057296580430072341));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1], DT_(0), tol);

    Tiny::Vector<DT_, 2> grad_3 = Analytic::eval_gradient_x(func, DT_(3), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_3[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_3[1], DT_(0), tol);

    Tiny::Vector<DT_, 2> grad_4 = Analytic::eval_gradient_x(func, DT_(3.58442834033049174494433823938), DT_(-1.84812652696440355353830020904));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_4[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_4[1], DT_(0), tol);
  }

  void test_rosenbrock_function() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create rosenbrock-function object
    Analytic::Common::RosenbrockFunction func;

    // evaluate function value in point (1,2)
    DT_ val = Analytic::eval_value_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(100), tol);

    // evaluate gradient in point (1,2)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(-400), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], DT_( 200), tol);

    // evaluate hessian in point (1,2)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_( 402), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_(-400), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_(-400), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_( 200), tol);

    // test local minimum
    DT_ min = Analytic::eval_value_x(func, DT_(1), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(min, DT_(0), tol);

    Tiny::Vector<DT_, 2> grad_min = Analytic::eval_gradient_x(func, DT_(1), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_min[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_min[1], DT_(0), tol);
  }

  void test_constant_function_1d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create constant-function object
    Analytic::Common::ConstantFunction<1,DT_> func(3);

    // evaluate function value in point (1)
    DT_ val = Analytic::eval_value_x(func, DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(3), tol);

    // evaluate gradient in point (1)
    Tiny::Vector<DT_, 1> grad = Analytic::eval_gradient_x(func, DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(0), tol);

    // evaluate hessian in point (1)
    Tiny::Matrix<DT_, 1,1> hess = Analytic::eval_hessian_x(func, DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(0), tol);
  }

  void test_constant_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create constant-function object
    Analytic::Common::ConstantFunction<2,DT_> func(4);

    // evaluate function value in point (1,1)
    DT_ val = Analytic::eval_value_x(func, DT_(1), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(4), tol);

    // evaluate gradient in point (1,1)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], DT_(0), tol);

    // evaluate hessian in point (1,1)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_(0), tol);
  }

  void test_xy_plane_rotation() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    //create origin vector
    typename Analytic::Common::XYPlaneRotation<DT_, 3>::PointType origin;
    origin[0] = DT_(1);
    origin[1] = DT_(1);
    origin[2] = DT_(1);

    // create xy-plane-rotation object
    Analytic::Common::XYPlaneRotation<DT_, 3> func(0.5, origin);

    // evaluate function value in point (1,2,3)
    Tiny::Vector<DT_, 3> val = Analytic::eval_value_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val[0], DT_(0.5)*-(DT_(2)-DT_(1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[1], DT_(0.5) * (DT_(1) - 1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[2], DT_(0), tol);

    // evaluate gradient in point (1,2,3)
    Tiny::Matrix<DT_, 3,3> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][0], DT_( 0  ), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][1], DT_(-0.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][2], DT_( 0  ), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][0], DT_( 0.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][1], DT_( 0  ), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][2], DT_( 0  ), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][0], DT_( 0  ), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][1], DT_( 0  ), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][2], DT_( 0  ), tol);

    // evaluate hessian in point (1,2,3)
    Tiny::Tensor3<DT_, 3,3,3,3,3,3> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess(0, 0, 0), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess(0, 0, 1), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess(0, 1, 0), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess(0, 1, 1), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess(1, 0, 0), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess(1, 0, 1), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess(1, 1, 0), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess(1, 1, 1), DT_(0), tol);
  }

  void test_yz_plane_parabolic_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    //create vectors
    typename Analytic::Common::YZPlaneParabolic<DT_, 2>::PointType zeros_y;
    zeros_y[0] = DT_(0);
    zeros_y[1] = DT_(4);

    // create yz-plane-parabolic object
    Analytic::Common::YZPlaneParabolic<DT_, 2> func(0.5, zeros_y);

    // evaluate function value in point (2,3)
    Tiny::Vector<DT_, 2> val = Analytic::eval_value_x(func, DT_(2), DT_(3));
    //TEST_CHECK_EQUAL_WITHIN_EPS(val[0], DT_(0.375), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[1], DT_(0), tol);

    // check the zeros of the function
    Tiny::Vector<DT_, 2> zer_2 = Analytic::eval_value_x(func, DT_(2), DT_(0));
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_2[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_2[1], DT_(0), tol);

    Tiny::Vector<DT_, 2> zer_3 = Analytic::eval_value_x(func, DT_(3), DT_(4));
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_3[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_3[1], DT_(0), tol);

    // check the maximum of the function
    Tiny::Vector<DT_, 2> max = Analytic::eval_value_x(func, DT_(4), DT_(2));
    //TEST_CHECK_EQUAL_WITHIN_EPS(max[0], DT_(0.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(max[1], DT_(0), tol);

    // evaluate gradient in maximum
    Tiny::Matrix<DT_, 2, 2> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][1], DT_(0), tol);
  }

  void test_yz_plane_parabolic_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    //create vectors
    typename Analytic::Common::YZPlaneParabolic<DT_, 2>::PointType zeros_y;
    zeros_y[0] = DT_(0);
    zeros_y[1] = DT_(4);
    typename Analytic::Common::YZPlaneParabolic<DT_, 2>::PointType zeros_z;
    zeros_z[0] = DT_(1);
    zeros_z[1] = DT_(5);

    // create yz-plane-parabolic object
    Analytic::Common::YZPlaneParabolic<DT_, 3> func(DT_(0.5), zeros_y, zeros_z);

    // check the zeros of the function
    Tiny::Vector<DT_, 3> zer_1 = Analytic::eval_value_x(func, DT_(2), DT_(0), DT_(1));
    //TEST_CHECK_EQUAL_WITHIN_EPS(zer_1[0], DT_(0), tol);
    //TEST_CHECK_EQUAL_WITHIN_EPS(zer_1[1], DT_(0), tol);
  }

  void test_sin_yt0_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constant
    const DT_ pi = Math::pi<DT_>();

    // create sin-yt0 object
    Analytic::Common::SinYT0<DT_, 2> func(pi);

    // evaluate function value in point (1,2)
    Tiny::Vector<DT_, 2> val_1 = Analytic::eval_value_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(0), tol);

    // evaluate gradient in point (1,2)
    Tiny::Matrix<DT_, 2,2> grad_1 = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], pi    , tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(0), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate function value in point (1,2)
    Tiny::Vector<DT_, 2> val_2 = Analytic::eval_value_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[1], DT_(0), tol);

    // evaluate gradient in point (1,2)
    Tiny::Matrix<DT_, 2, 2> grad_2 = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][0],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][1], -DT_(0.5) * pi, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][0],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][1],  DT_(0), tol);
  }

  void test_sin_yt0_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constant
    const DT_ pi = Math::pi<DT_>();

    // create sin-yt0 object
    Analytic::Common::SinYT0<DT_, 3> func(DT_(0.25)*pi);

    // evaluate function value in point (1,2,3)
    Tiny::Vector<DT_, 3> val_1 = Analytic::eval_value_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[2], DT_(0), tol);

    // evaluate gradient in point (1,2,3)
    Tiny::Matrix<DT_, 3, 3> grad_1 = Analytic::eval_gradient_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][2], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][2], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][2], DT_(0), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate function value in point (1,2,3)
    Tiny::Vector<DT_, 3> val_2 = Analytic::eval_value_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[2], DT_(0), tol);

    // evaluate gradient in point (1,2,3)
    Tiny::Matrix<DT_, 3, 3> grad_2 = Analytic::eval_gradient_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][0],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][1], -DT_(0.5) * pi, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][2],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][0],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][1],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][2],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[2][0],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[2][1],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[2][2],  DT_(0), tol);
  }

  void test_sin_yt0_stokes_rhs_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constant
    const DT_ pi = Math::pi<DT_>();

    // create sin-yt0-stokes-rhs object
    Analytic::Common::SinYT0StokesRhs<DT_, 2> func(DT_(2300),pi);

    // evaluate function value in point (1,2)
    Tiny::Vector<DT_, 2> val_1 = Analytic::eval_value_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(2), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(0), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate function value in point (0,1)
    Tiny::Vector<DT_, 2> val_2 = Analytic::eval_value_x(func, DT_(0), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[0], DT_(1) / DT_(2300) * DT_(0.5) * pi * DT_(0.5) * pi, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[1], DT_(0), tol);
  }

  void test_sin_yt0_stokes_rhs_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constant
    const DT_ pi = Math::pi<DT_>();

    // create sin-yt0-stokes-rhs object with t=0
    Analytic::Common::SinYT0StokesRhs<DT_, 3> func(DT_(2300));

    // evaluate function value in point (1,2,3)
    Tiny::Vector<DT_, 3> val_1 = Analytic::eval_value_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(2), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[2], DT_(0), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate function value in point (1,2,3)
    Tiny::Vector<DT_, 3> val_2 = Analytic::eval_value_x(func, DT_(1), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[0], DT_(-2), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[1], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[2], DT_( 0), tol);
  }

  void test_guermond_stokes_sol_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constant
    const DT_ pi = Math::pi<DT_>();

    // create guermond-stokes-sol object
    Analytic::Common::GuermondStokesSol<DT_, 2> func(pi);

    // evaluate function value in point (pi, 2*pi)
    Tiny::Vector<DT_, 2> val = Analytic::eval_value_x(func, pi, DT_(2)*pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[0], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[1], DT_(-1), tol);

    // evaluate gradient in point (pi, 2*pi)
    Tiny::Matrix<DT_, 2, 2> grad_1 = Analytic::eval_gradient_x(func, pi, DT_(2)*pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(0), tol);

    // change time t
    func.set_time(DT_(0.5)*pi);

    // evaluate function value in point (0.5*pi, pi)
    Tiny::Vector<DT_, 2> val_2 = Analytic::eval_value_x(func, DT_(0.5) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[1], DT_(0), tol);

    // evaluate gradient in point (0.5*pi, pi)
    Tiny::Matrix<DT_, 2,2> grad_2 = Analytic::eval_gradient_x(func, DT_(0.5) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][0], DT_( 1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][1], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][0], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][1], DT_(-1), tol);
  }

  void test_guermond_stokes_sol_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constant
    const DT_ pi = Math::pi<DT_>();

    // create guermond-stokes-sol object
    Analytic::Common::GuermondStokesSol<DT_, 3> func(pi);

    // evaluate function value in point (pi, 2*pi, pi)
    Tiny::Vector<DT_, 3> val = Analytic::eval_value_x(func, pi, DT_(2) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[0], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[1], DT_(-1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[2], DT_( 0), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate gradient in point (0.5 * pi, pi, pi)
    Tiny::Matrix<DT_, 3, 3> grad = Analytic::eval_gradient_x(func, DT_(0.5) * pi, pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][0], DT_( 1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][1], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][2], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][0], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][1], DT_(-1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][2], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][0], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][1], DT_( 0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][2], DT_( 0), tol);
  }

  void test_guermond_stokes_sol_pressure_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constants
    const DT_ pi = Math::pi<DT_>();
    const DT_ s1 = Math::sin(DT_(0.5) * pi + DT_(1));

    // create guermond-stokes-sol-pressure object
    Analytic::Common::GuermondStokesSolPressure<DT_, 2> func(pi);

    // evaluate function value in point (0.5 * pi, pi)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(0.5) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(1), tol);

    // evaluate gradient in point (0.5 * pi, pi)
    Tiny::Vector<DT_, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(0.5) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(0), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate function value in point (0.5 * pi, pi))
    DT_ val_2 = Analytic::eval_value_x(func, DT_(0.5) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, - (DT_(2) - DT_(2) * s1), tol);

    // evaluate gradient in point (0.5 * pi, pi))
    Tiny::Vector<DT_, 2> grad_2 = Analytic::eval_gradient_x(func, DT_(0.5) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0], DT_( 1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1], DT_(-1), tol);
  }

  void test_guermond_stokes_sol_pressure_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constants
    const DT_ pi = Math::pi<DT_>();
    const DT_ s1 = Math::sin(DT_(0.5) * pi + DT_(1));

    // create guermond-stokes-sol-pressure object
    Analytic::Common::GuermondStokesSolPressure<DT_, 3> func(pi);

    // evaluate function value in point (0.5 * pi, pi, pi)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(0.5) * pi, pi,pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(1), tol);

    // evaluate gradient in point (0.5 * pi, pi, pi)
    Tiny::Vector<DT_, 3> grad_1 = Analytic::eval_gradient_x(func, DT_(0.5) * pi, pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2], DT_(0), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate function value in point (0.5 * pi, pi, pi)
    DT_ val_2 = Analytic::eval_value_x(func, DT_(0.5) * pi, pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, -(DT_(2) - DT_(2) * s1), tol);

    // evaluate gradient in point (0.5 * pi, pi, pi)
    Tiny::Vector<DT_, 3> grad_2 = Analytic::eval_gradient_x(func, DT_(0.5) * pi, pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0], DT_( 1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1], DT_(-1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[2], DT_( 0), tol);
  }

  void test_guermond_stokes_sol_rhs_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constant
    const DT_ pi = Math::pi<DT_>();

    // create guermond-stokes-sol-rhs object
    Analytic::Common::GuermondStokesSolRhs<DT_, 2> func(DT_(2300),pi);

    // evaluate function value in point (pi, 2*pi)
    Tiny::Vector<DT_, 2>  val_1 = Analytic::eval_value_x(func, pi, DT_(2) * pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0],  DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], -DT_(2) / DT_(2300) - DT_(1), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate function value in point (pi, 2*pi)
    Tiny::Vector<DT_, 2> val_2 = Analytic::eval_value_x(func, pi, DT_(2) * pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[0], -DT_(2) / DT_(2300), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[1],  DT_(0), tol);
  }

  void test_guermond_stokes_sol_rhs_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // useful constant
    const DT_ pi = Math::pi<DT_>();

    // create guermond-stokes-sol-rhs object
    Analytic::Common::GuermondStokesSolRhs<DT_, 3> func(DT_(2300), pi);

    // evaluate function value in point (pi, 2*pi, pi)
    Tiny::Vector<DT_, 3>  val_1 = Analytic::eval_value_x(func, pi, DT_(2) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0],  DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], -DT_(2) / DT_(2300) - DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[2],  DT_(0), tol);

    // change time t
    func.set_time(DT_(0.5) * pi);

    // evaluate function value in point (pi, 2*pi, pi)
    Tiny::Vector<DT_, 3> val_2 = Analytic::eval_value_x(func, pi, DT_(2) * pi, pi);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[0], -DT_(2) / DT_(2300), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[1],  DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[2],  DT_(0), tol);
  }

  virtual void run() const override
  {
    test_par_profile_scalar();
    test_par_profile_vector();
    test_distance_function_2d();
    test_distance_function_3d();
    test_distance_function_sd_2d();
    test_distance_function_sd_3d();
    test_plane_distance_function_sd();
    test_min_function();
    test_sine_bubble_function_2d();
    test_sine_bubble_function_3d();
    test_cosine_wave_function_2d();
    test_cosine_wave_function_3d();
    test_q2_bubble_function_2d();
    test_q2_bubble_function_3d();
    test_exp_bubble_function_1d();
    test_exp_bubble_function_2d();
    test_exp_bubble_function_3d();
    test_exp_function_2d();
    test_exp_function_3d();
    test_heaviside_function_1d();
    test_heaviside_function_2d();
    test_heaviside_function_3d();
    test_heaviside_reg_function_1d();
    test_heaviside_reg_function_2d();
    test_heaviside_reg_function_3d();
    test_goldstein_price_function();
    test_bazaraa_shetty_function();
    test_himmelblau_function();
    test_rosenbrock_function();
    test_constant_function_1d();
    test_constant_function_2d();
    test_xy_plane_rotation();
    test_yz_plane_parabolic_2d();
    //test_yz_plane_parabolic_3d();
    test_sin_yt0_2d();
    test_sin_yt0_3d();
    test_sin_yt0_stokes_rhs_2d();
    test_sin_yt0_stokes_rhs_3d();
    test_guermond_stokes_sol_2d();
    test_guermond_stokes_sol_3d();
    test_guermond_stokes_sol_pressure_2d();
    test_guermond_stokes_sol_pressure_3d();
    test_guermond_stokes_sol_rhs_2d();
    test_guermond_stokes_sol_rhs_3d();
  }
};

CommonFunctionTest<double> common_function_test_double;
