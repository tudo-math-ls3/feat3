#include <test_system/test_system.hpp>
#include <kernel/util/string.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for the String class.
 *
 * \test Tests the String class.
 *
 * \author Peter Zajac
 */
class StringTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  StringTest() :
    TaggedTest<Archs::None, Archs::None>("string_test")
  {
  }

  virtual void run() const
  {
    bool b;
    int i;

    // test bool stringify/parse
    TEST_CHECK(stringify(b = true).parse(b));
    TEST_CHECK(b == true);
    TEST_CHECK(stringify(b = false).parse(b));
    TEST_CHECK(b == false);

    // test int stringify/parse
    TEST_CHECK(stringify(i = 42).parse(i));
    TEST_CHECK(i == 42);
    TEST_CHECK(stringify(i = -7).parse(i));
    TEST_CHECK(i == -7);

    // should not parse
    TEST_CHECK(!stringify(27).parse(b));
    TEST_CHECK(!stringify("nonsense").parse(i));
  }
} string_test;
