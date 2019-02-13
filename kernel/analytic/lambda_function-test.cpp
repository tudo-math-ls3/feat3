// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2022 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/lambda_function.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Analytic;

template<typename DT_, typename IT_>
class LambdaFunctionTest :
  public UnitTest
{
public:
  LambdaFunctionTest(PreferredBackend backend) :
    UnitTest("LambdaFunctionTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }


  void test_scalar_1d_a() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2
    auto func = create_lambda_function_scalar_1d(
      [](DT_ x) -> DT_ {return DT_(2)*x*x;},
      [](DT_ x) -> DT_ {return DT_(4)*x;},
      [](DT_  ) -> DT_ {return DT_(4);}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(18), tol);

    Tiny::Vector<DT_,1> grad_1 = Analytic::eval_gradient_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(12), tol);

    Tiny::Matrix<DT_,1,1> hess_1 = Analytic::eval_hessian_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(4), tol);
  }

  void test_scalar_1d_b() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2
    auto func = create_lambda_function_scalar_1d(
      [](DT_ x) -> DT_ {return DT_(2)*x*x;},
      [](DT_ x) -> DT_ {return DT_(4)*x;}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(18), tol);

    Tiny::Vector<DT_,1> grad_1 = Analytic::eval_gradient_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(12), tol);

    Tiny::Matrix<DT_,1,1> hess_1 = Analytic::eval_hessian_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(4), tol);
  }

  void test_scalar_1d_c() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2
    auto func = create_lambda_function_scalar_1d(
      [](DT_ x) -> DT_ {return DT_(2)*x*x;}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(18), tol);

    Tiny::Vector<DT_,1> grad_1 = Analytic::eval_gradient_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(12), tol);

    Tiny::Matrix<DT_,1,1> hess_1 = Analytic::eval_hessian_x(func, DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(4), tol);
  }


  void test_scalar_2d_1() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.6));

    // local variable to capture by lambda
    DT_ a(2.0);

    // create lambda function: a*x - y^2
    auto func = create_lambda_function_scalar_2d([&](DT_ x, DT_ y) -> DT_ {return a*x - y*y;});

    // 2*3 - 2*2 = 6 - 4 = 2
    DT_ val_1 = Analytic::eval_value_x(func, DT_(3), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(2), tol);

    // change captured variable
    a = DT_(5.0);

    // 5*3 - 2*2 = 15 - 4 = 11
    DT_ val_2 = Analytic::eval_value_x(func, DT_(3), DT_(2));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_2, DT_(11), tol);
  }

  void test_scalar_2d_a() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2 + x*y^2
    auto func = create_lambda_function_scalar_2d(
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*x*x + y*y*x;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(4)*x + y*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*y*x;},
      [](DT_  , DT_  ) -> DT_ {return DT_(4);},
      [](DT_ x, DT_  ) -> DT_ {return DT_(2)*x;},
      [](DT_  , DT_ y) -> DT_ {return DT_(2)*y;}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(68), tol);

    Tiny::Vector<DT_,2> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(25), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(24), tol);

    Tiny::Matrix<DT_,2,2> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1], DT_(6), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0], DT_(6), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1], DT_(8), tol);
  }

  void test_scalar_2d_b() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2 + x*y^2
    auto func = create_lambda_function_scalar_2d(
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*x*x + y*y*x;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(4)*x + y*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*y*x;}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(68), tol);

    Tiny::Vector<DT_,2> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(25), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(24), tol);

    Tiny::Matrix<DT_,2,2> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1], DT_(6), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0], DT_(6), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1], DT_(8), tol);
  }

  void test_scalar_2d_c() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2 + x*y^2
    auto func = create_lambda_function_scalar_2d(
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*x*x + y*y*x;}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(68), tol);

    Tiny::Vector<DT_,2> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(25), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(24), tol);

    Tiny::Matrix<DT_,2,2> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1], DT_(6), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0], DT_(6), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1], DT_(8), tol);
  }

  void test_vector_2d_a() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*y + 3*x*y^2
    auto func = create_lambda_function_vector_2d(
      // value
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*x*x*y + DT_(3)*x*y*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(3)*y*x*x + DT_(5)*x*y*y;},
      // dx
      [](DT_ x, DT_ y) -> DT_ {return DT_(4)*x*y + DT_(3)*y*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(6)*y*x + DT_(5)*y*y;},
      // dy
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*x*x + DT_(6)*x*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(3)*x*x + DT_(10)*x*y;},
      // dxx
      [](DT_  , DT_ y) -> DT_ {return DT_(4)*y;},
      [](DT_  , DT_ y) -> DT_ {return DT_(6)*y;},
      // dyy
      [](DT_ x, DT_  ) -> DT_ {return DT_(6)*x;},
      [](DT_ x, DT_  ) -> DT_ {return DT_(10)*x;},
      // dxy
      [](DT_ x, DT_ y) -> DT_ {return DT_(4)*x + DT_(6)*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(6)*x + DT_(10)*y;}
      );

    Tiny::Vector<DT_,2> val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(204), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(324), tol);

    Tiny::Matrix<DT_,2,2> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(75), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(104), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(117), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(168), tol);

    Tiny::Tensor3<DT_,2,2,2> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(34), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(34), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(24), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_(54), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_(54), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_(40), tol);
  }

  void test_vector_2d_b() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*y + 3*x*y^2
    auto func = create_lambda_function_vector_2d(
      // value
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*x*x*y + DT_(3)*x*y*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(3)*y*x*x + DT_(5)*x*y*y;},
      // dx
      [](DT_ x, DT_ y) -> DT_ {return DT_(4)*x*y + DT_(3)*y*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(6)*y*x + DT_(5)*y*y;},
      // dy
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*x*x + DT_(6)*x*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(3)*x*x + DT_(10)*x*y;}
    );

    Tiny::Vector<DT_,2> val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(204), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(324), tol);

    Tiny::Matrix<DT_,2,2> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(75), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(104), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(117), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(168), tol);

    Tiny::Tensor3<DT_,2,2,2> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(34), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(34), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(24), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_(54), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_(54), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_(40), tol);
  }

  void test_vector_2d_c() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*y + 3*x*y^2
    auto func = create_lambda_function_vector_2d(
      // value
      [](DT_ x, DT_ y) -> DT_ {return DT_(2)*x*x*y + DT_(3)*x*y*y;},
      [](DT_ x, DT_ y) -> DT_ {return DT_(3)*y*x*x + DT_(5)*x*y*y;}
    );

    Tiny::Vector<DT_,2> val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(204), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(324), tol);

    Tiny::Matrix<DT_,2,2> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(75), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(104), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(117), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(168), tol);

    Tiny::Tensor3<DT_,2,2,2> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(34), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(34), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(24), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_(54), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_(54), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_(40), tol);
  }


  void test_scalar_3d_a() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*z - 3*y*z^2 + x*y^2
    auto func = create_lambda_function_scalar_3d(
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x*z - DT_(3)*y*z*z + y*y*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*x*z + y*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*y - DT_(3)*z*z;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x - DT_(6)*y*z;},
      [](DT_  , DT_  , DT_ z) -> DT_ {return DT_(4)*z;},
      [](DT_ x, DT_  , DT_  ) -> DT_ {return DT_(2)*x;},
      [](DT_  , DT_ y, DT_  ) -> DT_ {return -DT_(6)*y;},
      [](DT_  , DT_ y, DT_  ) -> DT_ {return DT_(2)*y;},
      [](DT_  , DT_  , DT_ z) -> DT_ {return -DT_(6)*z;},
      [](DT_ x, DT_  , DT_  ) -> DT_ {return DT_(4)*x;}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(26), tol);

    Tiny::Vector<DT_,3> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(84), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], -DT_(59), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2], -DT_(28), tol);

    Tiny::Matrix<DT_,3,3> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1], DT_(8), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2], -DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1], -DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2], -DT_(30), tol);
  }

  void test_scalar_3d_b() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*z - 3*y*z^2 + x*y^2
    auto func = create_lambda_function_scalar_3d(
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x*z - DT_(3)*y*z*z + y*y*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*x*z + y*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*y - DT_(3)*z*z;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x - DT_(6)*y*z;}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(26), tol);

    Tiny::Vector<DT_,3> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(84), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], -DT_(59), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2], -DT_(28), tol);

    Tiny::Matrix<DT_,3,3> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1], DT_(8), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2], -DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1], -DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2], -DT_(30), tol);
  }

  void test_scalar_3d_c() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*z - 3*y*z^2 + x*y^2
    auto func = create_lambda_function_scalar_3d(
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x*z - DT_(3)*y*z*z + y*y*x;}
    );

    DT_ val_1 = Analytic::eval_value_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(26), tol);

    Tiny::Vector<DT_,3> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_(84), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], -DT_(59), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2], -DT_(28), tol);

    Tiny::Matrix<DT_,3,3> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(2), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1], DT_(8), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2], -DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1], -DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2], -DT_(30), tol);
  }

  void test_vector_3d_a() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*z - 3*y*z^2 + x*y^2
    auto func = create_lambda_function_vector_3d(
      // value
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x*z + DT_(3)*y*z*z + DT_(4)*y*y*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(3)*y*y*x + DT_(4)*z*x*x + DT_(2)*z*z*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*z*z*y + DT_(2)*x*y*y + DT_(3)*x*x*z;},
      // dx
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*x*z + DT_(4)*y*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(3)*y*y + DT_(8)*z*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*y*y + DT_(6)*x*z;},
      // dy
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(3)*z*z + DT_(8)*y*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(6)*y*x + DT_(2)*z*z;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*z*z + DT_(4)*x*y;},
      // dz
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x + DT_(6)*y*z;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*x*x + DT_(4)*z*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(8)*z*y + DT_(3)*x*x;},
      // dxx
      [](DT_  , DT_  , DT_ z) -> DT_ {return DT_(4)*z;},
      [](DT_  , DT_  , DT_ z) -> DT_ {return DT_(8)*z;},
      [](DT_  , DT_  , DT_ z) -> DT_ {return DT_(6)*z;},
      // dyy
      [](DT_ x, DT_  , DT_  ) -> DT_ {return DT_(8)*x;},
      [](DT_ x, DT_  , DT_  ) -> DT_ {return DT_(6)*x;},
      [](DT_ x, DT_  , DT_  ) -> DT_ {return DT_(4)*x;},
      // dzz
      [](DT_  , DT_ y, DT_  ) -> DT_ {return DT_(6)*y;},
      [](DT_  , DT_ y, DT_  ) -> DT_ {return DT_(4)*y;},
      [](DT_  , DT_ y, DT_  ) -> DT_ {return DT_(8)*y;},
      // dxy
      [](DT_  , DT_ y, DT_  ) -> DT_ {return DT_(8)*y;},
      [](DT_  , DT_ y, DT_  ) -> DT_ {return DT_(6)*y;},
      [](DT_  , DT_ y, DT_  ) -> DT_ {return DT_(4)*y;},
      // dyz
      [](DT_  , DT_  , DT_ z) -> DT_ {return DT_(6)*z;},
      [](DT_  , DT_  , DT_ z) -> DT_ {return DT_(4)*z;},
      [](DT_  , DT_  , DT_ z) -> DT_ {return DT_(8)*z;},
      // dzx
      [](DT_ x, DT_  , DT_  ) -> DT_ {return DT_(4)*x;},
      [](DT_ x, DT_  , DT_  ) -> DT_ {return DT_(8)*x;},
      [](DT_ x, DT_  , DT_  ) -> DT_ {return DT_(6)*x;}
    );

    Tiny::Vector<DT_,3> val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(529), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(578), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[2], DT_(612), tol);

    Tiny::Matrix<DT_,3,3> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(116), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(171), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][2], DT_(122), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(187), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(122), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][2], DT_(124), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][0], DT_(138), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][1], DT_(148), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][2], DT_(168), tol);

    Tiny::Tensor3<DT_,3,3,3> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][2], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][2], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][0], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][1], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][2], DT_(18), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][2], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][2], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][0], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][1], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][2], DT_(12), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][0], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][1], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][2], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][0], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][1], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][2], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][0], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][1], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][2], DT_(24), tol);
  }

  void test_vector_3d_b() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*z - 3*y*z^2 + x*y^2
    auto func = create_lambda_function_vector_3d(
      // value
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x*z + DT_(3)*y*z*z + DT_(4)*y*y*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(3)*y*y*x + DT_(4)*z*x*x + DT_(2)*z*z*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*z*z*y + DT_(2)*x*y*y + DT_(3)*x*x*z;},
      // dx
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*x*z + DT_(4)*y*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(3)*y*y + DT_(8)*z*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*y*y + DT_(6)*x*z;},
      // dy
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(3)*z*z + DT_(8)*y*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(6)*y*x + DT_(2)*z*z;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*z*z + DT_(4)*x*y;},
      // dz
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x + DT_(6)*y*z;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*x*x + DT_(4)*z*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(8)*z*y + DT_(3)*x*x;}
    );

    Tiny::Vector<DT_,3> val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(529), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(578), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[2], DT_(612), tol);

    Tiny::Matrix<DT_,3,3> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(116), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(171), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][2], DT_(122), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(187), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(122), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][2], DT_(124), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][0], DT_(138), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][1], DT_(148), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][2], DT_(168), tol);

    Tiny::Tensor3<DT_,3,3,3> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][2], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][2], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][0], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][1], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][2], DT_(18), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][2], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][2], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][0], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][1], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][2], DT_(12), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][0], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][1], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][2], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][0], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][1], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][2], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][0], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][1], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][2], DT_(24), tol);
  }

  void test_vector_3d_c() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.5));

    // create lambda function: 2*x^2*z - 3*y*z^2 + x*y^2
    auto func = create_lambda_function_vector_3d(
      // value
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(2)*x*x*z + DT_(3)*y*z*z + DT_(4)*y*y*x;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(3)*y*y*x + DT_(4)*z*x*x + DT_(2)*z*z*y;},
      [](DT_ x, DT_ y, DT_ z) -> DT_ {return DT_(4)*z*z*y + DT_(2)*x*y*y + DT_(3)*x*x*z;}
    );

    Tiny::Vector<DT_,3> val_1 = Analytic::eval_value_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(529), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(578), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[2], DT_(612), tol);

    Tiny::Matrix<DT_,3,3> grad_1 = Analytic::eval_gradient_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_(116), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(171), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][2], DT_(122), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(187), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_(122), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][2], DT_(124), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][0], DT_(138), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][1], DT_(148), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[2][2], DT_(168), tol);

    Tiny::Tensor3<DT_,3,3,3> hess_1 = Analytic::eval_hessian_x(func, DT_(4), DT_(3), DT_(5));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][0], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][1], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][0][2], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][0], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][1], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][1][2], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][0], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][1], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[0][2][2], DT_(18), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][0], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][1], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][0][2], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][0], DT_(18), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][1], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][1][2], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][0], DT_(32), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][1], DT_(20), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[1][2][2], DT_(12), tol);

    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][0], DT_(30), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][1], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][0][2], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][0], DT_(12), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][1], DT_(16), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][1][2], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][0], DT_(24), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][1], DT_(40), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess_1[2][2][2], DT_(24), tol);
  }

  virtual void run() const override
  {
    test_scalar_1d_a();
    test_scalar_1d_b();
    test_scalar_1d_c();

    test_scalar_2d_a();
    test_scalar_2d_b();
    test_scalar_2d_c();
    test_scalar_2d_1();

    test_scalar_3d_a();
    test_scalar_3d_b();
    test_scalar_3d_c();

    test_vector_2d_a();
    test_vector_2d_b();
    test_vector_2d_c();

    test_vector_3d_a();
    test_vector_3d_b();
    test_vector_3d_c();
  }
};

LambdaFunctionTest<double, Index> lambda_function_test_double(PreferredBackend::generic);
