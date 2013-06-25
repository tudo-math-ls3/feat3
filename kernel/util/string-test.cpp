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
    std::vector<String> words;
    String s;

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
    TEST_CHECK(!String("nonsense").parse(i));

    // test split-by-charset
    String("  0 5 \n 3\tfoo ").split_by_charset(words);
    TEST_CHECK_EQUAL(words.size(), 4u);
    TEST_CHECK_EQUAL(words[0], "0");
    TEST_CHECK_EQUAL(words[1], "5");
    TEST_CHECK_EQUAL(words[2], "3");
    TEST_CHECK_EQUAL(words[3], "foo");

    // test split-by-string
    String("0, 4,,7 , ").split_by_string(words, ",");
    TEST_CHECK_EQUAL(words.size(), 5u);
    TEST_CHECK_EQUAL(words[0], "0");
    TEST_CHECK_EQUAL(words[1], " 4");
    TEST_CHECK_EQUAL(words[2], "");
    TEST_CHECK_EQUAL(words[3], "7 ");
    TEST_CHECK_EQUAL(words[4], " ");

    // test replace_all
    TEST_CHECK_EQUAL((s = "abcde").replace_all("", "xy"), size_t(0));
    TEST_CHECK_EQUAL((s = "abcde").replace_all("bc", "xy"), size_t(1));
    TEST_CHECK_EQUAL(s, "axyde");
    TEST_CHECK_EQUAL((s = "abcde").replace_all("cb", "xy"), size_t(0));
    TEST_CHECK_EQUAL((s = "abcde").replace_all("bc", "xbc"), size_t(1));
    TEST_CHECK_EQUAL(s, "axbcde");
    TEST_CHECK_EQUAL((s = "aaaa").replace_all("aa", "pa"), size_t(2));
    TEST_CHECK_EQUAL(s, "papa");
  }
} string_test;
