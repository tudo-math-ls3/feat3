#include <test_system/test_system.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::Math;

/**
 * \brief Test class for the Math functions.
 *
 * \test Tests the class templates in the MetaMath namespace.
 *
 * \author Peter Zajac
 */
class MathTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  const double tol;
  MathTest() :
    TaggedTest<Archs::None, Archs::None>("MathTest"),
    tol(Math::pow(Math::eps<double>(), 0.9))
  {
  }

  virtual void run() const
  {
    test_factorial();
    test_binomial();
    test_sqrt();
    test_sin();
    test_cos();
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

  void test_sqrt() const
  {
    // we need a weaker tolerance for sqrt
    const double tol(Math::pow(Math::eps<double>(), 0.4));

    // exact root tests
    TEST_CHECK_EQUAL(Math::sqrt<double>(0.0), 0.0);
    TEST_CHECK_EQUAL(Math::sqrt<double>(1.0), 1.0);
    TEST_CHECK_EQUAL(Math::sqrt<double>(64.0), 8.0);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<double>(2.0), std::sqrt(2.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sqrt<double>(0.5), std::sqrt(0.5), tol);
  }

  void test_sin() const
  {
    // test exact
    TEST_CHECK_EQUAL(Math::sin<double>(0.0), 0.0);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::sin<double>(0.7), std::sin(0.7), tol);
  }

  void test_cos() const
  {
    // test exact
    TEST_CHECK_EQUAL(Math::cos<double>(0.0), 1.0);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::cos<double>(0.7), std::cos(0.7), tol);
  }

  void test_exp() const
  {
    // test exact
    TEST_CHECK_EQUAL(Math::exp<double>(0.0), 1.0);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::exp<double>( 2.0), std::exp( 2.0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::exp<double>(-2.0), std::exp(-2.0), tol);
  }

  void test_log() const
  {
    // test exact
    TEST_CHECK_EQUAL(Math::log<double>(1.0), 0.0);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log<double>(0.5), std::log(0.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log<double>(2.0), std::log(2.0), tol);
  }

  void test_log10() const
  {
    // test exact
    TEST_CHECK_EQUAL(Math::log10<double>(1.0), 0.0);
    TEST_CHECK_EQUAL(Math::log10<double>(10.0), 1.0);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<double>(0.5), std::log10(0.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::log10<double>(2.0), std::log10(2.0), tol);
  }

  void test_pow() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<double>(1.0, 1.0), 1.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<double>(2.0, 3.0), 8.0, tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<double>(4.0, 0.5), 2.0, tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<double>(0.7, 3.1), std::pow(0.7, 3.1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pow<double>(3.1, 0.7), std::pow(3.1, 0.7), tol);
  }

  void test_atan() const
  {
    // test exact
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<double>(0.0), 0.0, tol);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<double>( 1.5), std::atan( 1.5), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan<double>(-1.5), std::atan(-1.5), tol);
  }

  void test_atan2() const
  {
    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan2<double>(0.7, 3.1), std::atan2(0.7, 3.1), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::atan2<double>(3.1, 0.7), std::atan2(3.1, 0.7), tol);
  }

  void test_pi() const
  {
    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::pi<double>(), 2.0*std::acos(0.0), tol);
  }

  void test_asin() const
  {
    // test exact
    TEST_CHECK_EQUAL(Math::asin<double>(0.0), 0.0);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::asin<double>(0.7), std::asin(0.7), tol);
  }

  void test_acos() const
  {
    // test exact
    TEST_CHECK_EQUAL(Math::acos<double>(1.0), 0.0);

    // test against std
    TEST_CHECK_EQUAL_WITHIN_EPS(Math::acos<double>(0.7), std::acos(0.7), tol);
  }

} math_test;
