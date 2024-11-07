// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/auto_derive.hpp>
#include <kernel/analytic/function.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Analytic;

// 2*x^2 - 3*y^2
class ScalarTestFunction :
  public Analytic::Function
{
public:
  typedef Analytic::Image::Scalar ImageType;
  static constexpr int domain_dim = 2;
  static constexpr bool can_value = true;
  static constexpr bool can_grad = false;
  static constexpr bool can_hess = false;

  template<typename EvalTraits_>
  class Evaluator :
    public Analytic::Function::Evaluator<EvalTraits_>
  {
  public:
    typedef typename EvalTraits_::DataType DataType;
    typedef typename EvalTraits_::PointType PointType;
    typedef typename EvalTraits_::ValueType ValueType;

  public:
    explicit Evaluator(const ScalarTestFunction&)
    {
    }

    ValueType value(const PointType& point)
    {
      return point[0]*point[1] + DataType(2)*Math::sqr(point[0]) - DataType(3)*Math::sqr(point[1]);
    }
  };
}; // class ScalarTestFunction


// [2*x^2 - 3*y^2 + x*y, 2*y^3 - 3*x^3 - x*y]
class VectorTestFunction :
  public Analytic::Function
{
public:
  typedef Analytic::Image::Vector<2> ImageType;
  static constexpr int domain_dim = 2;
  static constexpr bool can_value = true;
  static constexpr bool can_grad = false;
  static constexpr bool can_hess = false;

  template<typename EvalTraits_>
  class Evaluator :
    public Analytic::Function::Evaluator<EvalTraits_>
  {
  public:
    typedef typename EvalTraits_::DataType DataType;
    typedef typename EvalTraits_::PointType PointType;
    typedef typename EvalTraits_::ValueType ValueType;

  public:
    explicit Evaluator(const VectorTestFunction&)
    {
    }

    ValueType value(const PointType& point)
    {
      return ValueType({
        DataType(2)*Math::sqr(point[0]) - DataType(3)*Math::sqr(point[1]) + point[0]*point[1],
        DataType(2)*Math::cub(point[1]) - DataType(3)*Math::cub(point[0]) - point[0]*point[1],
      });
    }
  };
}; // class VectorTestFunction

template<typename DT_, typename IT_>
class AutoDeriveTest :
  public UnitTest
{
public:
  AutoDeriveTest(PreferredBackend backend) :
    UnitTest("AutoDeriveTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  void test_auto_derive_scalar2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.6));

    Analytic::AutoDerive<ScalarTestFunction, DT_> func;

    // evaluate value in point (1.25, 0.75)
    DT_ val_1 = Analytic::eval_value_x(func, DT_(1.25), DT_(0.75));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1, DT_(2.375), tol);

    // evaluate gradient in point (1.25, 0.75)
    Tiny::Vector<DT_, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(1.25), DT_(0.75));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0], DT_( 5.75), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1], DT_(-3.25), tol);

    // evaluate hessian in point (1.25, 0.75)
    Tiny::Matrix<DT_, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(1.25), DT_(0.75));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1], DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0], DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1], DT_(-6), tol);
  }

  void test_auto_derive_vector2d() const
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.6));

    Analytic::AutoDerive<VectorTestFunction, DT_> func;

    // evaluate value in point (1.25, 0.75)
    Tiny::Vector<DT_, 2> val_1 = Analytic::eval_value_x(func, DT_(1.25), DT_(0.75));
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[0], DT_(2.375), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(val_1[1], DT_(-5.953125), tol);

    // evaluate gradient in point (1.25, 0.75)
    Tiny::Matrix<DT_, 2, 2> grad_1 = Analytic::eval_gradient_x(func, DT_(1.25), DT_(0.75));
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][0], DT_( 5.75), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[0][1], DT_(-3.25), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][0], DT_(-14.8125), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(grad_1[1][1], DT_( 2.125), tol);

    // evaluate hessian in point (1.25, 0.75)
    Tiny::Tensor3<DT_, 2, 2, 2> hess = Analytic::eval_hessian_x(func, DT_(1.25), DT_(0.75));
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0][0], DT_(4), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][0][1], DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1][0], DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[0][1][1], DT_(-6), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0][0], DT_(-22.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][0][1], DT_(-1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1][0], DT_(-1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(hess[1][1][1], DT_(9), tol);
  }

  virtual void run() const override
  {
    test_auto_derive_scalar2d();
    test_auto_derive_vector2d();
  }
}; // class AutoDeriveTest<...>

AutoDeriveTest<double, Index> auto_derive_test_double(PreferredBackend::generic);
