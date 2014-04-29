#include <test_system/test_system.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Math;

class BasicMathTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  BasicMathTest() :
    TaggedTest<Archs::None, Archs::None>("BasicMathTest")
  {
  }

  virtual void run() const
  {
    test_factorial();
    test_binomial();
  }

  void test_factorial() const
  {
    // 0! = 1 (by definition)
    TEST_CHECK_EQUAL(factorial(0), 1);
    // 1! = 1
    TEST_CHECK_EQUAL(factorial(1,0), 1);
    // 5! = 120
    TEST_CHECK_EQUAL(factorial(5), 120);
    // 5*6*7 = 210
    TEST_CHECK_EQUAL(factorial(7,5), 210);
  }

  void test_binomial() const
  {
    TEST_CHECK_EQUAL(binomial(0,0), 1);
    TEST_CHECK_EQUAL(binomial(3,0), 1);
    TEST_CHECK_EQUAL(binomial(3,3), 1);
    TEST_CHECK_EQUAL(binomial(5,2), 10);
    TEST_CHECK_EQUAL(binomial(49,6), 13983816);
  }
} basic_math_test;

/**
 * \brief Test class for the Math functions.
 *
 * \test Tests the class templates in the MetaMath namespace.
 *
 * \author Peter Zajac
 */
template<typename DT_>
class MathTest
  : public TaggedTest<Archs::None, DT_>
{
public:
  const DT_ tol;
  MathTest() :
    TaggedTest<Archs::None, DT_>("MathTest"),
    tol(Math::pow(Math::eps<DT_>(), DT_(0.9)))
  {
  }

  virtual void run() const
  {
    test_sqrt();
    test_sin();
    test_cos();
    test_sinh();
    test_cosh();
    test_exp();
    test_log();
    test_log10();
    test_pow();
    test_atan();
    test_atan2();
    test_pi();
    test_asin();
    test_acos();
  }

  void test_sqrt() const
  {
    // we need a weaker tolerance for sqrt
    const DT_ tol2(Math::pow(Math::eps<DT_>(), DT_(0.4)));

    // exact root tests
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(0)), DT_(0), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(1)), DT_(1), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(64)), DT_(8), tol2);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(2)), Math::sqrt(DT_(2)), tol2);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<DT_>(DT_(0.5)), Math::sqrt(DT_(0.5)), tol2);
  }

  void test_sin() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sin<DT_>(DT_(0)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sin<DT_>(DT_(0.7)), Math::sin(DT_(0.7)), tol);
  }

  void test_cos() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cos<DT_>(DT_(0)), DT_(1), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cos<DT_>(DT_(0.7)), Math::cos(DT_(0.7)), tol);
  }

  void test_sinh() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sinh<DT_>(DT_(0)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sinh<DT_>(DT_(0.7)), Math::sinh(DT_(0.7)), tol);
  }

  void test_cosh() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cosh<DT_>(DT_(0)), DT_(1), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cosh<DT_>(DT_(0.7)), Math::cosh(DT_(0.7)), tol);
  }

  void test_exp() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::exp<DT_>(DT_(0)), DT_(1), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::exp<DT_>( DT_(2)), Math::exp( DT_(2)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::exp<DT_>(-DT_(2)), Math::exp(-DT_(2)), tol);
  }

  void test_log() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log<DT_>(DT_(1)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log<DT_>(DT_(0.5)), Math::log(DT_(0.5)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log<DT_>(DT_(2)), Math::log(DT_(2)), tol);
  }

  void test_log10() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<DT_>(DT_(1)), DT_(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<DT_>(DT_(10)), DT_(1), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<DT_>(DT_(0.5)), Math::log10(DT_(0.5)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<DT_>(DT_(2)), Math::log10(DT_(2)), tol);
  }

  void test_pow() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(1), DT_(1)), DT_(1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(2), DT_(3)), DT_(8), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(4), DT_(0.5)), DT_(2), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(0.7), DT_(3.1)), Math::pow(DT_(0.7), DT_(3.1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<DT_>(DT_(3.1), DT_(0.7)), Math::pow(DT_(3.1), DT_(0.7)), tol);
  }

  void test_atan() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<DT_>(DT_(0)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<DT_>( DT_(1.5)), Math::atan( DT_(1.5)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<DT_>(-DT_(1.5)), Math::atan(-DT_(1.5)), tol);
  }

  void test_atan2() const
  {
    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan2<DT_>(DT_(0.7), DT_(3.1)), Math::atan2(DT_(0.7), DT_(3.1)), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan2<DT_>(DT_(3.1), DT_(0.7)), Math::atan2(DT_(3.1), DT_(0.7)), tol);
  }

  void test_pi() const
  {
    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pi<DT_>(), DT_(2)*Math::acos(DT_(0)), tol);
  }

  void test_asin() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::asin<DT_>(DT_(0)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::asin<DT_>(DT_(0.7)), Math::asin(DT_(0.7)), tol);
  }

  void test_acos() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::acos<DT_>(DT_(1)), DT_(0), tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::acos<DT_>(DT_(0.7)), Math::acos(DT_(0.7)), tol);
  }
};

MathTest<double> math_test_double;
#ifdef FEAST_HAVE_QUADMATH
MathTest<__float128> math_test_float128;
#endif // FEAST_HAVE_QUADMATH
