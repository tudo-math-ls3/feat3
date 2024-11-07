// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

class A
{
  public:
  bool flag;

  explicit A(bool f) :
    flag(f)
  {
  }

  A(const A & other) = delete;
};

bool operator==(const A & a, const A & b)
{
  return a.flag == b.flag;
}

std::ostream & operator<< (std::ostream & lhs, const A & b)
{
  lhs << "A: "<< b.flag << "\n";
  return lhs;
}



/**
*
* \brief test class for the unittest framework itself
*
* \test test description missing
*
* \author Dirk Ribbrock
*/
class TestUnitTest
  : public UnitTest
{
public:
  /// Constructor
  explicit TestUnitTest()
    : UnitTest("TestUnitTest")
  {
  }

  virtual ~TestUnitTest()
  {
  }

  /// runs the tests
  virtual void run() const override
  {
    TEST_CHECK(true);
    TEST_CHECK_EQUAL(1, 1);
    TEST_CHECK_NOT_EQUAL(1.0, 0.5);
    TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
    TEST_CHECK_EQUAL_WITHIN_EPS(25., 23., 2.2);
    TEST_CHECK_EQUAL_WITHIN_EPS(6, 7, 1);
    TEST_CHECK_THROWS(auto a = String("0").at(10);(void)a, std::exception);


    TEST_CHECK_RELATIVE(1235., 1234., 0.001);
    TEST_CHECK_RELATIVE(-5673., -5678., 0.001);

    A a(true);
    A b(true);
    TEST_CHECK_EQUAL(a, b);
    A c(false);
    TEST_CHECK_NOT_EQUAL(a, c);
  }
} testunittest;
