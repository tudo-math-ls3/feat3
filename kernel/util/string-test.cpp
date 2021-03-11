// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/util/string.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

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

  virtual ~StringTest()
  {
  }

  virtual void run() const override
  {
    bool b;
    int i;
    std::deque<String> words;
    String s;

    // test operator+
    TEST_CHECK_EQUAL(String("ab") + String("cd"), "abcd");

    // test compare_no_case
    TEST_CHECK_EQUAL(String("g1").compare_no_case("G1"), 0);

    // test trim
    TEST_CHECK_EQUAL(String("  ab ").trim_front(), "ab ");
    TEST_CHECK_EQUAL(String("  ab ").trim_back(), "  ab");
    TEST_CHECK_EQUAL(String("  ab ").trim(), "ab");

    // test pad
    TEST_CHECK_EQUAL(String("ab").pad_front(3, 'q'), "qab");
    TEST_CHECK_EQUAL(String("ab").pad_back(4, 'c'), "abcc");
    TEST_CHECK_EQUAL(String("abcd").pad_front(2), "abcd");
    TEST_CHECK_EQUAL(String("abcd").pad_back(4), "abcd");

    // test starts_with(String)
    TEST_CHECK_EQUAL(String("abcde").starts_with("abc"), true);
    TEST_CHECK_EQUAL(String("abcde").starts_with("foo"), false);
    TEST_CHECK_EQUAL(String("abcde").starts_with(""), true);
    TEST_CHECK_EQUAL(String("abcde").starts_with("abcdef"), false);
    TEST_CHECK_EQUAL(String("").starts_with(""), true);
    TEST_CHECK_EQUAL(String("").starts_with("foo"), false);

    // test ends_with(String)
    TEST_CHECK_EQUAL(String("abcde").ends_with("cde"), true);
    TEST_CHECK_EQUAL(String("abcde").ends_with("foo"), false);
    TEST_CHECK_EQUAL(String("abcde").ends_with(""), true);
    TEST_CHECK_EQUAL(String("abcde").ends_with("abcdef"), false);
    TEST_CHECK_EQUAL(String("").ends_with(""), true);
    TEST_CHECK_EQUAL(String("").ends_with("foo"), false);

    // test starts_with(char)
    TEST_CHECK_EQUAL(String("abcde").starts_with('a'), true);
    TEST_CHECK_EQUAL(String("abcde").starts_with('x'), false);
    TEST_CHECK_EQUAL(String("").starts_with('x'), false);

    // test ends_with(char)
    TEST_CHECK_EQUAL(String("abcde").ends_with('e'), true);
    TEST_CHECK_EQUAL(String("abcde").ends_with('x'), false);
    TEST_CHECK_EQUAL(String("").ends_with('x'), false);

    // test trunc
    TEST_CHECK_EQUAL(String("abcde").trunc_front(7u), "abcde");
    TEST_CHECK_EQUAL(String("abcde").trunc_front(3u), "cde");
    TEST_CHECK_EQUAL(String("abcde").trunc_back(7u), "abcde");
    TEST_CHECK_EQUAL(String("abcde").trunc_back(3u), "abc");

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
    words = String("  0 5 \n 3\tfoo ").split_by_whitespaces();
    TEST_CHECK_EQUAL(words.size(), 4u);
    TEST_CHECK_EQUAL(words[0], "0");
    TEST_CHECK_EQUAL(words[1], "5");
    TEST_CHECK_EQUAL(words[2], "3");
    TEST_CHECK_EQUAL(words[3], "foo");

    // test split-by-string
    words = String("0, 4,,7 , ").split_by_string(",");
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
