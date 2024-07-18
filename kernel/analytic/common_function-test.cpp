// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/distance_function.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Analytic;

template<typename DT_>
class CommonFunctionTest :
  public UnitTest
{
public:
  CommonFunctionTest() :
    UnitTest("CommonFunctionTest", Type::Traits<DT_>::name(), "none", PreferredBackend::generic)
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

  void test_min_function() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create 3d DistanceFunction with origin in (2,2,1).
    typename Analytic::Distance::DistanceFunction<3, DT_>::PointType orig3d;
    orig3d[0] = DT_(2);
    orig3d[1] = DT_(2);
    orig3d[2] = DT_(1);
    Analytic::Distance::DistanceFunction<3, DT_> func1(orig3d);
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
    const DT_ u_x = (Math::exp(-(DT_(2) * x - DT_(1))*(DT_(2) * x - DT_(1))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(x)
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
    const DT_ u_x = (Math::exp(-(DT_(2) * x - DT_(1))*(DT_(2) * x - DT_(1))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(x)
    const DT_ u_y = (Math::exp(-(DT_(2) * y - DT_(1))*(DT_(2) * y - DT_(1))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(y)
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
    const DT_ u_x = (Math::exp(-(DT_(2) * x - DT_(1))*(DT_(2) * x - DT_(1))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(x)
    const DT_ u_y = (Math::exp(-(DT_(2) * y - DT_(1))*(DT_(2) * y - DT_(1))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(y)
    const DT_ u_z = (Math::exp(-(DT_(2) * z - DT_(1))*(DT_(2) * z - DT_(1))) - e_neg_1) / (DT_(1) - e_neg_1); // = u(z)
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
    const DT_ u_x_hess = -(DT_(20) * e_x * (DT_(20) * Math::sqr(x) + DT_(1))) / (e_10 - DT_(1)); // = u''(x)
    const DT_ u_y_hess = -(DT_(20) * e_y * (DT_(20) * Math::sqr(y) + DT_(1))) / (e_10 - DT_(1)); // = u''(y)

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
    const DT_ u_x_hess = -(DT_(20) * e_x * (DT_(20) * Math::sqr(x) + DT_(1))) / (e_10 - DT_(1)); // = u''(x)
    const DT_ u_y_hess = -(DT_(20) * e_y * (DT_(20) * Math::sqr(y) + DT_(1))) / (e_10 - DT_(1)); // = u''(y)
    const DT_ u_z_hess = -(DT_(20) * e_z * (DT_(20) * Math::sqr(z) + DT_(1))) / (e_10 - DT_(1)); // = u''(z)

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
    DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));

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
    DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.6));

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
    DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));

    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

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

  void test_constant_vector_function_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create a value vector
    Tiny::Vector<DT_, 2> values;
    for(int i(0); i < 2; ++i)
      values[i] = DT_(i) + DT_(0.731);

    // create constant function object
    Analytic::Common::ConstantVectorFunction<2,DT_> func(4);
    //create constant
    Analytic::Common::ConstantVectorFunction<2,DT_> func2(values);

    // evaluate function value in point (1,1)
    Tiny::Vector<DT_, 2> val = Analytic::eval_value_x(func, DT_(1), DT_(1));
    Tiny::Vector<DT_, 2> val2 = Analytic::eval_value_x(func2, DT_(1), DT_(1));
    for(int i(0); i < 2; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(val[i], DT_(4), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(val2[i], DT_(i) + DT_(0.731), tol);
    }

    // evaluate gradient in point (1,1)
    Tiny::Matrix<DT_, 2, 2> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(1));
    Tiny::Matrix<DT_, 2, 2> grad2 = Analytic::eval_gradient_x(func2, DT_(1), DT_(1));

    for(int i(0); i < 2; ++i)
    {
      for(int j(0); j < 2; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(grad(i,j), DT_(0), tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(grad2(i,j), DT_(0), tol);
      }
    }

    // evaluate hessian in point (1,1)
    Tiny::Tensor3<DT_, 2, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(1));
    Tiny::Tensor3<DT_, 2, 2, 2> hess2 = Analytic::eval_hessian_x(func2, DT_(1), DT_(1));

    for(int i(0); i < 2; ++i)
    {
      for(int j(0); j < 2; ++j)
      {
        for(int k(0); k < 2; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(hess(i,j,k), DT_(0), tol);
          TEST_CHECK_EQUAL_WITHIN_EPS(hess2(i,j,k), DT_(0), tol);
        }
      }
    }
  }

  void test_constant_vector_function_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create a value vector
    Tiny::Vector<DT_, 3> values;
    for(int i(0); i < 3; ++i)
      values[i] = DT_(i) + DT_(0.731);

    // create constant function object
    Analytic::Common::ConstantVectorFunction<3,DT_> func(4);
    //create constant
    Analytic::Common::ConstantVectorFunction<3,DT_> func2(values);

    // evaluate function value in point (1,1)
    Tiny::Vector<DT_, 3> val = Analytic::eval_value_x(func, DT_(1), DT_(1), DT_(1));
    Tiny::Vector<DT_, 3> val2 = Analytic::eval_value_x(func2, DT_(1), DT_(1), DT_(1));
    for(int i(0); i < 3; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(val[i], DT_(4), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(val2[i], DT_(i) + DT_(0.731), tol);
    }

    // evaluate gradient in point (1,1)
    Tiny::Matrix<DT_, 3, 3> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(1), DT_(1));
    Tiny::Matrix<DT_, 3, 3> grad2 = Analytic::eval_gradient_x(func2, DT_(1), DT_(1), DT_(1));

    for(int i(0); i < 3; ++i)
    {
      for(int j(0); j < 3; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(grad(i,j), DT_(0), tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(grad2(i,j), DT_(0), tol);
      }
    }

    // evaluate hessian in point (1,1)
    Tiny::Tensor3<DT_, 3, 3, 3> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(1), DT_(1));
    Tiny::Tensor3<DT_, 3, 3, 3> hess2 = Analytic::eval_hessian_x(func2, DT_(1), DT_(1), DT_(1));

    for(int i(0); i < 3; ++i)
    {
      for(int j(0); j < 3; ++j)
      {
        for(int k(0); k < 3; ++k)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(hess(i,j,k), DT_(0), tol);
          TEST_CHECK_EQUAL_WITHIN_EPS(hess2(i,j,k), DT_(0), tol);
        }
      }
    }
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
    TEST_CHECK_EQUAL_WITHIN_EPS(val[1], DT_(0.5) * (DT_(1) - DT_(1)), tol);
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
    typename Analytic::Common::YZPlaneParabolic<DT_, 2>::RangeType range_y;
    range_y[0] = DT_(0);
    range_y[1] = DT_(4);

    // create yz-plane-parabolic object
    Analytic::Common::YZPlaneParabolic<DT_, 2> func(0.5, range_y);

    // evaluate function value in point (2,3)
    Tiny::Vector<DT_, 2> val = Analytic::eval_value_x(func, DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val[0], DT_(0.375), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[1], DT_(0), tol);

    // evaluate gradient in point (2,3)
    Tiny::Matrix<DT_, 2, 2> grad = Analytic::eval_gradient_x(func, DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][0], DT_(-0.25), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][1], DT_(0), tol);

    // check the zeros of the function
    Tiny::Vector<DT_, 2> zer_2 = Analytic::eval_value_x(func, DT_(2), DT_(0));
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_2[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_2[1], DT_(0), tol);

    Tiny::Vector<DT_, 2> zer_3 = Analytic::eval_value_x(func, DT_(3), DT_(4));
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_3[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_3[1], DT_(0), tol);

    // check the maximum of the function
    Tiny::Vector<DT_, 2> max = Analytic::eval_value_x(func, DT_(4), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(max[0], DT_(0.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(max[1], DT_(0), tol);

    // evaluate gradient in maximum
    Tiny::Matrix<DT_, 2, 2> grad_max = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[0][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[1][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[1][1], DT_(0), tol);
  }

  void test_yz_plane_parabolic_3d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    //create vectors
    typename Analytic::Common::YZPlaneParabolic<DT_, 2>::RangeType range_y;
    range_y[0] = DT_(0); // a_1
    range_y[1] = DT_(4); // b_1
    typename Analytic::Common::YZPlaneParabolic<DT_, 2>::RangeType range_z;
    range_z[0] = DT_(1); // a_2
    range_z[1] = DT_(5); //b_2

    // create yz-plane-parabolic object
    Analytic::Common::YZPlaneParabolic<DT_, 3> func(DT_(0.5), range_y, range_z);

    // evaluate function value in point (2,3,4)
    Tiny::Vector<DT_, 3> val = Analytic::eval_value_x(func, DT_(2), DT_(3), DT_(4));
    TEST_CHECK_EQUAL_WITHIN_EPS(val[0], DT_(0.28125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val[2], DT_(0), tol);

    // evaluate gradient in point (2,3,4)
    Tiny::Matrix<DT_, 3, 3> grad = Analytic::eval_gradient_x(func, DT_(2), DT_(3), DT_(4));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0][2], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][0], DT_(-0.1875), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1][2], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][0], DT_(-0.1875), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2][2], DT_(0), tol);

    // check the zeros of the function in point (2,0,1)
    Tiny::Vector<DT_, 3> zer_1 = Analytic::eval_value_x(func, DT_(2), DT_(0), DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_1[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_1[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_1[2], DT_(0), tol);

    // check the zeros of the function in point (2,3,5)
    Tiny::Vector<DT_, 3> zer_2 = Analytic::eval_value_x(func, DT_(2), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_2[0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_2[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(zer_2[2], DT_(0), tol);

    // check the maximum of the function in point (2,2,3)
    Tiny::Vector<DT_, 3> max = Analytic::eval_value_x(func, DT_(2), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(max[0], DT_(0.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(max[1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(max[2], DT_(0), tol);

    // evaluate gradient in maximum
    Tiny::Matrix<DT_, 3, 3> grad_max = Analytic::eval_gradient_x(func, DT_(2), DT_(2), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[0][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[0][2], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[1][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[1][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[1][2], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[2][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[2][1], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_max[2][2], DT_(0), tol);
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

  void test_polynomial_function_1d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create polynomial-function-1d object
    Analytic::Common::PolynomialFunction1D<DT_> p;
    p.set_coeff(0, 5.0); // 5*x^0
    p.set_coeff(1, -7.0); // -7*x^1
    p.set_coeff(2, 3.0); // 3*x^2

    // evaluate function value in point (13)
    DT_ val_1 = Analytic::eval_value_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(421), tol);

    // evaluate gradient in point (13)
    Tiny::Vector<DT_, 1> grad = Analytic::eval_gradient_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(71), tol);

    // evaluate hessian in point (13)
    Tiny::Matrix<DT_, 1,1> hess = Analytic::eval_hessian_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(6), tol);

    // change coefficient
    p.set_coeff(2, 0.0); // 0*x^2

    // evaluate function value in point (13)
    DT_ val_2 = Analytic::eval_value_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(-86), tol);

    // evaluate gradient in point (13)
    Tiny::Vector<DT_, 1> grad_2 = Analytic::eval_gradient_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0], DT_(-7), tol);

    // evaluate hessian in point (13)
    Tiny::Matrix<DT_, 1, 1> hess_2 = Analytic::eval_hessian_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_2[0][0], DT_(0), tol);

    // set coefficient
    p.set_coeff(3, 3.0); // 3*x^2

    // evaluate function value in point (13)
    DT_ val_3 = Analytic::eval_value_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_3, DT_(6505), tol);

    // evaluate gradient in point (13)
    Tiny::Vector<DT_, 1> grad_3 = Analytic::eval_gradient_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_3[0], DT_(1514), tol);

    // evaluate hessian in point (13)
    Tiny::Matrix<DT_, 1, 1> hess_3 = Analytic::eval_hessian_x(p, DT_(13));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_3[0][0], DT_(234), tol);
  }

  void test_standing_vortex_function_2d() const
  {
    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create standing-vortex function object
    Analytic::Common::StandingVortexFunction2D func;

    // evaluate in point (0.6,0.575)
    Tiny::Vector<DT_, 2> val_1 = Analytic::eval_value_x(func, DT_(0.6), DT_(0.575));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(0.375), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(-0.5), tol);

    // evaluate in p0int (0.7, 0.8)
    Tiny::Vector<DT_, 2> val_2 = Analytic::eval_value_x(func, DT_(0.7), DT_(0.8));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[0], DT_(0.1641005886756876), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2[1], DT_(-0.1094003924504584), tol);

    // evaluate gradient in point (0.6,0.575)
    Tiny::Matrix<DT_, 2, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(0.6), DT_(0.575));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(-5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(0), tol);

    // evaluate gradient in point (0.7, 0.8)
    Tiny::Matrix<DT_, 2, 2> grad_2 = Analytic::eval_gradient_x(func, DT_(0.7), DT_(0.8));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][0], DT_(-2.560154751808750), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[0][1], DT_(-3.293230165460834), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][0], DT_(1.159767872286875), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_2[1][1], DT_(2.560154751808750), tol);
  }

  void test_taylor_green_vortex_velo_2d() const
  {
    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    Analytic::Common::TaylorGreenVortexVelo2D<DT_> func(DT_(0.01), DT_(0.5));

    // evaluate function
    Tiny::Vector<DT_, 2> val_1 = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_( 0.59188481860155275663), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(-0.24516671922750232133), tol);

    // evaluate gradient
    Tiny::Matrix<DT_, 2, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_( 1.8594609978899655386), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(-0.77021396402983280152), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_( 0.77021396402983280152), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(-1.8594609978899655386), tol);

    // evaluate hessian
    Tiny::Tensor3<DT_, 2, 2, 2> hess_1 = Analytic::eval_hessian_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(-5.8416690106078617623), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(-2.4196985310883959903), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(-2.4196985310883959903), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(-5.8416690106078617623), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_( 2.4196985310883959903), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_( 5.8416690106078617623), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_( 5.8416690106078617623), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_( 2.4196985310883959903), tol);
  }

  void test_taylor_green_vortex_pres_2d() const
  {
    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    Analytic::Common::TaylorGreenVortexPres2D<DT_> func(DT_(0.01), DT_(0.5));

    // evaluate function
    DT_ val_1 = Analytic::eval_value_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(0.14511045913710802783), tol);

    // evaluate gradient
    Tiny::Vector<DT_, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(0.25), DT_(0.125));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(-1.2894175660971681167), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(-0.91175590476836093505), tol);
  }

  void test_poiseuille_pipe_flow_3d() const
  {
    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));
    const DT_ tol2 = Math::pow(Math::eps<DT_>(), DT_(0.6));

    Tiny::Vector<DT_, 3> ori{DT_(0.1), DT_(0.2), DT_(0.3)};
    Tiny::Vector<DT_, 3> axis{DT_(2),DT_(3),DT_(4)};
    axis.normalize();

    Analytic::Common::PoiseuillePipeFlow<DT_, 3> func(ori, axis, DT_(2.3), DT_(1.7));

    // evaluate function
    Tiny::Vector<DT_, 3> val_1 = Analytic::eval_value_x(func, DT_(0.7), DT_(1.7), DT_(2.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(0.59580593160049614744), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(0.89370889740074422114), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[2], DT_(1.1916118632009922948), tol);

    // evaluate gradient
    Tiny::Matrix<DT_, 3, 3> grad_1 = Analytic::eval_gradient_x(func, DT_(0.7), DT_(1.7), DT_(2.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(0.10865011117118946116), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(0.019754565667488992928), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][2], DT_(-0.069140979836211475279), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(0.16297516675678419174), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(0.029631848501233489392), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][2], DT_(-0.10371146975431721292), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][0], DT_(0.21730022234237892233), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][1], DT_(0.039509131334977985857), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][2], DT_(-0.13828195967242295055), tol);

    // evaluate hessian
    Tiny::Tensor3<DT_, 3, 3, 3> hess_1 = Analytic::eval_hessian_x(func, DT_(0.7), DT_(1.7), DT_(2.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(-0.20577672567667092222), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(0.049386414162401021332), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][2], DT_(0.065848552216534695113), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(0.049386414162401021332), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(-0.16462138054133673778), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][2], DT_(0.098772828324802042661), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][0], DT_(0.065848552216534695113), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][1], DT_(0.098772828324802042661), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][2], DT_(-0.10700389735186887955), tol2);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_(-0.30866508851500638331), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_(0.074079621243601531993), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][2], DT_(0.098772828324802042661), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_(0.074079621243601531993), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_(-0.24693207081200510665), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][2], DT_(0.14815924248720306399), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][0], DT_(0.098772828324802042661), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][1], DT_(0.14815924248720306399), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][2], DT_(-0.16050584602780331932), tol2);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][0], DT_(-0.41155345135334184442), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][1], DT_(0.098772828324802042661), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][2], DT_(0.13169710443306939021), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][0], DT_(0.098772828324802042661), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][1], DT_(-0.32924276108267347554), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][2], DT_(0.19754565664960408532), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][0], DT_(0.13169710443306939021), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][1], DT_(0.19754565664960408532), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][2], DT_(-0.21400779470373775910), tol2);
  }


  void test_sine_ring_vortex_velo_2d() const
  {
    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));

    Analytic::Common::SineRingVortexVelo2D<DT_> func;

    // evaluate function
    Tiny::Vector<DT_, 2> val_1 = Analytic::eval_value_x(func, DT_(0.9), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_( 0.47348043213632444483), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(-0.60876055560384571479), tol);

    // evaluate gradient
    Tiny::Matrix<DT_, 2, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(0.9), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_( 1.6105289973850738860), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_( 1.9290342819704415944), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(-2.7470807568327013460), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(-1.6105289973850738860), tol);

    // evaluate hessian
    Tiny::Tensor3<DT_, 2, 2, 2> hess_1 = Analytic::eval_hessian_x(func, DT_(0.9), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(-13.202182918630301466), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(-9.3594239646432791943), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(-9.3594239646432791943), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(-3.7005986449779418496), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_(12.372723759995890782), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_(13.202182918630301466), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_(13.202182918630301466), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_(9.3594239646432791943), tol);
  }

  void test_sine_ring_vortex_pres_2d() const
  {
    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));

    Analytic::Common::SineRingVortexPres2D<DT_> func;

    // evaluate function
    DT_ val_1 = Analytic::eval_value_x(func, DT_(0.9), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(-0.22406392059917757436), tol);

    // evaluate gradient
    Tiny::Vector<DT_, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(0.9), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(-0.50248135191051632366), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(-0.39081882926373491839), tol);
  }

  void test_sine_ring_vortex_rhs_2d() const
  {
    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));

    Analytic::Common::SineRingVortexRHS2D<DT_> func(DT_(1.0), DT_(1.2), DT_(1.5), DT_(1.7));

    // evaluate function
    Tiny::Vector<DT_, 2> val_1 = Analytic::eval_value_x(func, DT_(0.9), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_( 16.264664694819381659), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(-23.693995515706431689), tol);
  }

  void test_ball_cap_function_2d() const
  {
    // skip this test for quad precision
    if(sizeof(DT_) > 8u)
      return;

    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));

    Analytic::Common::BallCapFunction2D func;

    // evaluate function
    DT_ val_1 = Analytic::eval_value_x(func, DT_(0.6), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(0.02138656368050375964), tol);

    // evaluate gradient
    Tiny::Vector<DT_, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(0.6), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(-0.16903085094570331550), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(-0.19720265943665386808), tol);

    // evaluate hessian
    Tiny::Matrix<DT_, 2, 2> hess_1 = Analytic::eval_hessian_x(func, DT_(0.6), DT_(0.7));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(-0.31391443747059187164), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1], DT_(-0.037562411321267403442), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0], DT_(-0.037562411321267403442), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1], DT_(-0.32554089811765082985), tol);
  }

  void test_corner_singularity_2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));
    {
      Common::CornerSingularity2DRadial<DT_> corner_nat(Math::pi<DT_>()*DT_(0.5));
      Common::CornerSingluratity2DSimple<DT_> corner_sing(corner_nat);

      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(corner_sing, DT_(1.0), DT_(0.)), DT_(0.), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(corner_sing, DT_(0.0), DT_(1.)), DT_(0.), tol);
      // std::cout << "Val is " << Analytic::eval_hessian_x(corner_sing, DT_(4.0), DT_(4.0)).trace() << "\n";

      Tiny::Vector<DT_, 2> corner{0., 0.}, r_vec{1.0, 0.}, l_vec{0., 1.0};
      Common::CornerSingularity2D<DT_> n_corner_sing(corner, r_vec, l_vec);
      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(n_corner_sing, DT_(1.0), DT_(0.)), DT_(0.), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(n_corner_sing, DT_(0.0), DT_(1.)), DT_(0.), tol);
      {
        DT_ x = DT_(0.4), y = DT_(3.1);
        DT_ val = Analytic::eval_value_x(corner_sing, x, y) - Analytic::eval_value_x(n_corner_sing, x, y);
        TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(0.), tol);
      }
      {
        DT_ x = DT_(0.4), y = DT_(3.1);
        auto val_t = Analytic::eval_gradient_x(corner_sing, x, y) - Analytic::eval_gradient_x(n_corner_sing, x, y);
        TEST_CHECK_EQUAL_WITHIN_EPS(val_t.norm_euclid(), DT_(0.), tol);
      }
      {
        DT_ x = DT_(0.4), y = DT_(3.1);
        auto val_t = Analytic::eval_hessian_x(corner_sing, x, y) - Analytic::eval_hessian_x(n_corner_sing, x, y);
        TEST_CHECK_EQUAL_WITHIN_EPS(val_t.norm_frobenius(), DT_(0.), tol);
      }

      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(corner_sing, DT_(0.1), DT_(0.6)), Analytic::eval_value_x(n_corner_sing, DT_(0.1), DT_(0.6)), tol);

      // std::cout << "Val is " << Analytic::eval_hessian_x(corner_sing, DT_(4.0), DT_(4.0)).trace() << "\n";
    }
    //shifted and rotated coordinate system
    {
      DT_ alpha = DT_(FEAT_F128C(0.4));
      Tiny::Vector<DT_, 2> corner{DT_(FEAT_F128C(0.0)), DT_(FEAT_F128C(1.0))}, r_vec{DT_(FEAT_F128C(0.0)), DT_(FEAT_F128C(1.0))}, l_vec{DT_(FEAT_F128C(1.4)), DT_(FEAT_F128C(-1.1))};
      DT_ theta = Math::calc_opening_angle(r_vec[0], r_vec[1], l_vec[0], l_vec[1]);
      DT_ offset = Math::calc_opening_angle(DT_(FEAT_F128C(1.0)), DT_(FEAT_F128C(0.0)), r_vec[0], r_vec[1]);
      Common::CornerSingularity2DRadial<DT_> corner_nat(theta, alpha);
      Common::CornerSingularity2DRadial<DT_> corner_alt(r_vec, l_vec, alpha);
      Common::CornerSingluratity2DSimple<DT_> corner_sing(corner_nat, corner, offset);
      Common::CornerSingluratity2DSimple<DT_> corner_sing_alt(corner_alt, corner, r_vec);

      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(corner_sing, DT_(FEAT_F128C(0.0)), DT_(FEAT_F128C(1.3))), DT_(FEAT_F128C(0.)), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(corner_sing, DT_(FEAT_F128C(1.4)), DT_(FEAT_F128C(-0.1))), DT_(FEAT_F128C(0.)), tol);
      // std::cout << "Val is " << Analytic::eval_hessian_x(corner_sing, DT_(4.0), DT_(4.0)).trace() << "\n";

      Common::CornerSingularity2D<DT_> n_corner_sing(corner, r_vec, l_vec, alpha);
      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(n_corner_sing, DT_(FEAT_F128C(0.0)), DT_(FEAT_F128C(1.3))), DT_(FEAT_F128C(0.)), tol);
      TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(n_corner_sing, DT_(FEAT_F128C(1.4)), DT_(FEAT_F128C(-0.1))), DT_(FEAT_F128C(0.)), tol);
      {
        DT_ x = DT_(-0.4), y = DT_(1.1);
        DT_ val = Analytic::eval_value_x(corner_sing, x, y) - Analytic::eval_value_x(n_corner_sing, x, y);
        TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(0.), tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_value_x(corner_sing, x, y), Analytic::eval_value_x(corner_sing_alt, x, y), tol);
      }
      {
        DT_ x = DT_(1.2), y = DT_(-2.1);
        auto val_t = Analytic::eval_gradient_x(corner_sing, x, y)- Analytic::eval_gradient_x(n_corner_sing, x, y);
        TEST_CHECK_EQUAL_WITHIN_EPS(val_t.norm_euclid(), DT_(0.), tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_gradient_x(corner_sing, x, y)[1], Analytic::eval_gradient_x(corner_sing_alt, x, y)[1], tol);
      }
      {
        DT_ x = DT_(0.4), y = DT_(-3.1);
        auto val_t = Analytic::eval_hessian_x(corner_sing, x, y)- Analytic::eval_hessian_x(n_corner_sing, x, y);
        TEST_CHECK_EQUAL_WITHIN_EPS(val_t.norm_frobenius(), DT_(0.), tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_hessian_x(corner_sing, x, y)[1][0], Analytic::eval_hessian_x(corner_sing_alt, x, y)[0][1], tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(Analytic::eval_hessian_x(n_corner_sing, x, y).trace(), DT_(0), tol);
      }

    }
  }

  void test_frankes_function() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    // create goldstein-price-function object
    Analytic::Common::FrankesFunction<DT_> func;

    // evaluate function value in point (1,2)
    DT_ val = Analytic::eval_value_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(FEAT_F128C(0.014574258847493253357381166783366)), tol);

    // evaluate gradient in point (1,2)
    Tiny::Vector<DT_, 2> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(FEAT_F128C(-0.053538093725485420496502610127097)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], DT_(FEAT_F128C(-0.013116832962743928021647611368389)), tol);

    // evaluate hessian in point (1,2)
    Tiny::Matrix<DT_, 2,2> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(FEAT_F128C(0.14848626402639731929540320127819)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_(FEAT_F128C(0.048184284352936878446893400485549)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_(FEAT_F128C(0.048184284352936878446893400485549)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_(FEAT_F128C(0.011805149666469535219787961767313)), tol);

  }

  void test_frankes_3d_variant_function() const
  {
    #ifdef FEAT_HAVE_QUADMATH
      const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));
    #else
      const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));
    #endif

    // create goldstein-price-function object
    Analytic::Common::Frankes3DVariantFunction<DT_> func;

    // evaluate function value in point (1,2, 0.5)
    DT_ val = Analytic::eval_value_x(func, DT_(1), DT_(2), DT_(0.5));
    TEST_CHECK_EQUAL_WITHIN_EPS(val, DT_(FEAT_F128C(0.010168115367927312326141809033119)), tol);

    // evaluate gradient in point (1,2, 0.5)
    Tiny::Vector<DT_, 3> grad = Analytic::eval_gradient_x(func, DT_(1), DT_(2), DT_(0.5));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[0], DT_(FEAT_F128C(-0.037352260535243188136847716266949)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[1], DT_(FEAT_F128C(-0.00915130383113458109353081041527)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad[2], DT_(FEAT_F128C(-0.024403476883025549582740341679484)), tol);

    // evaluate hessian in point (1,2, 0.5)
    Tiny::Matrix<DT_, 3,3> hess = Analytic::eval_hessian_x(func, DT_(1), DT_(2), DT_(0.5));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(FEAT_F128C(0.10359535115794998505709471473385)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_(FEAT_F128C(0.033617034481718869323191585210066)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][2], DT_(FEAT_F128C(0.089645425284583651528434519040676)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_(FEAT_F128C(0.033617034481718869323191585210066)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_(FEAT_F128C(0.0082361734480211229843905984690556)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][2], DT_(FEAT_F128C(0.021963129194722994624473944996648)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][0], DT_(FEAT_F128C(0.089645425284583651528434519040676)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][1], DT_(FEAT_F128C(0.021963129194722994624473944996648)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[2][2], DT_(FEAT_F128C(-0.022776578424157179610557652234185)), tol);

  }

  virtual void run() const override
  {
    test_par_profile_scalar();
    test_par_profile_vector();
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
    #ifndef FEAT_HAVE_HALFMATH
    test_goldstein_price_function();
    #endif
    test_bazaraa_shetty_function();
    test_himmelblau_function();
    test_rosenbrock_function();
    test_constant_function_1d();
    test_constant_function_2d();
    test_constant_vector_function_2d();
    test_constant_vector_function_3d();
    test_xy_plane_rotation();
    test_yz_plane_parabolic_2d();
    test_yz_plane_parabolic_3d();
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
    test_polynomial_function_1d();
    test_standing_vortex_function_2d();
    test_taylor_green_vortex_velo_2d();
    test_taylor_green_vortex_pres_2d();
    test_poiseuille_pipe_flow_3d();
    test_sine_ring_vortex_velo_2d();
    test_sine_ring_vortex_pres_2d();
    test_sine_ring_vortex_rhs_2d();
    test_ball_cap_function_2d();
    test_corner_singularity_2d();
    test_frankes_function();
    test_frankes_3d_variant_function();
  }
};

CommonFunctionTest <double> common_function_test_double;
CommonFunctionTest <float> common_function_test_float;
#ifdef FEAT_HAVE_QUADMATH
CommonFunctionTest <__float128> common_function_test_float128;
#endif
#ifdef FEAT_HAVE_HALFMATH
CommonFunctionTest <Half> common_function_test_half;
#endif
#ifdef FEAT_HAVE_CUDA
CommonFunctionTest <float> cuda_common_function_test_float;
CommonFunctionTest <double> cuda_common_function_test_double;
#endif
