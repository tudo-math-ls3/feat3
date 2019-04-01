// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/util/meta_math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::MetaMath;

/**
 * \brief Test class for the MetaMath class templates.
 *
 * \test Tests the class templates in the MetaMath namespace.
 *
 * \author Peter Zajac
 */
class MetaMathTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:

  MetaMathTest() :
    TaggedTest<Archs::None, Archs::None>("MetaMathTest")
  {
  }

  virtual ~MetaMathTest()
  {
  }

  virtual void run() const override
  {
    test_factorial();
    test_binomial();
  }

  void test_factorial() const
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

  void test_binomial() const
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
} meta_math_test;
