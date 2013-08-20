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

  MathTest() :
    TaggedTest<Archs::None, Archs::None>("MathTest")
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
} math_test;
