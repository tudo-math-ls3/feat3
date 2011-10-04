// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <test_system/test_system.hpp>

#include <string>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the assertion class.
*
* \test test description missing
*
* \tparam Tag_
* description missing
*
* \tparam DT_
* description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Tag_,
  typename DT_>
class AssertionTest
  : public TaggedTest<Tag_, DT_>
{

public:

  AssertionTest()
    : TaggedTest<Tag_, DT_>("assertion_test")
  {
  }

  virtual void run() const
  {
    CONTEXT("When breaching the surface");
    CONTEXT("When going deeper");
    CONTEXT("When reaching the ground");
    TEST_CHECK_THROWS(ASSERT(false, "Should throw!"), Assertion);

    bool no_exception_thrown(true);
    try
    {
      ASSERT(true, "Shouldn't throw!");
    }
    catch (...)
    {
      no_exception_thrown = false;
    }
    TEST_CHECK(no_exception_thrown);

  }
};
AssertionTest<Nil, Nil> assertion_test;
