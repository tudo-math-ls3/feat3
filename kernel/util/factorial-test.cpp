#include <kernel/base_header.hpp>
#include <kernel/util/factorial.hpp>
#include <test_system/test_system.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for the Factorial class template.
 *
 * \test Tests the Factorial class template.
 *
 * \author Peter Zajac
 */
class FactorialTest
  : public TaggedTest<Nil, Nil>
{
public:

  FactorialTest() :
    TaggedTest<Nil, Nil>("factorial_test")
  {
  }

  virtual void run() const
  {
    // Note:
    // We have to store the result of the compile-time constant factorials in an int variable, as the
    // template magic within the TEST_CHECK_* macros gets confused otherwise. Besides, the compiler might
    // warn for compile-time constant if-expressions otherwise.

    // 0! = 1 (by definition)
    int fact_0_0 = Factorial<0>::value;
    TEST_CHECK_EQUAL(fact_0_0, 1);

    // 1! = 1
    int fact_1_0 = Factorial<1,0>::value;
    TEST_CHECK_EQUAL(fact_1_0, 1);

    // 5! = 120
    int fact_5 = Factorial<5>::value;
    TEST_CHECK_EQUAL(fact_5, 120);

    // 5*6*7 = 210
    int fact_7_5 = Factorial<7,5>::value;
    TEST_CHECK_EQUAL(fact_7_5, 210);
  }
} factorial_test;

/**
 * \brief Test class for the Binomial class template.
 *
 * \test Tests the Binomial class template.
 *
 * \author Peter Zajac
 */
class BinomialTest
  : public TaggedTest<Nil, Nil>
{
public:

  BinomialTest() :
    TaggedTest<Nil, Nil>("binomial_test")
  {
  }

  virtual void run() const
  {
    // Note:
    // We have to store the result of the compile-time constant binomials in an int variable, as the
    // template magic within the TEST_CHECK_* macros gets confused otherwise. Besides, the compiler might
    // warn for compile-time constant if-expressions otherwise.

    int bin_0_0 = Binomial<0, 0>::value;
    TEST_CHECK_EQUAL(bin_0_0, 1);

    int bin_3_0 = Binomial<3, 0>::value;
    TEST_CHECK_EQUAL(bin_3_0, 1);

    int bin_3_3 = Binomial<3, 3>::value;
    TEST_CHECK_EQUAL(bin_3_3, 1);

    int bin_5_2 = Binomial<5, 2>::value;
    TEST_CHECK_EQUAL(bin_5_2, 10);

    int bin_49_6 = Binomial<49, 6>::value;
    TEST_CHECK_EQUAL(bin_49_6, 13983816);
  }
} binomial_test;
