// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

// includes, FEAST
#include <kernel/util/inheritance.hpp>
#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/util/assertion.hpp>

using namespace TestSystem;
using namespace FEAST;

class A{};
class B : public A{};
class C{};

/**
* \brief Test class for the inheritance macro.
*
* \test test description missing
*
* \author Dirk Ribbrock
*/
template <typename Tag_, typename DT_>
class InheritanceTest
  : public TaggedTest<Tag_, DT_>
{
public:
  InheritanceTest()
    : TaggedTest<Tag_, DT_>("inheritance_test")
  {
  }

  virtual void run() const
  {
    TEST_CHECK(SUPERSUBCLASS(A, B));
    TEST_CHECK(SUPERSUBCLASS(B, B));
    TEST_CHECK(! SUPERSUBCLASS(B, A));
    TEST_CHECK(! SUPERSUBCLASS(B, A));
    TEST_CHECK(! SUPERSUBCLASS(A, C));
    TEST_CHECK(! SUPERSUBCLASS(C, A));
    STATIC_ASSERT(SUPERSUBCLASS(A, B), wrong_inheritance_detetect_at_compile_time);
  }
};
InheritanceTest<Nil, Nil> inheritance_test;
