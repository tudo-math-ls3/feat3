#include <test_system/test_system.hpp>
#include <kernel/util/simple_arg_parser.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

// an enumeration for testing purposes
enum class MyPrecon
{
  none,
  ssor,
  spai,
  iluk
};

// This one's required by the TEST_CHECK macros in the case that an error occurs
inline std::ostream& operator<<(std::ostream& os, MyPrecon p)
{
  switch(p)
  {
  case MyPrecon::none: return os << "none";
  case MyPrecon::ssor: return os << "ssor";
  case MyPrecon::spai: return os << "spai";
  case MyPrecon::iluk: return os << "iluk";
  default: return os;
  }
}

/**
 * \brief Test class for the SimpleArgParser class.
 *
 * \test Tests the SimpleArgParser class.
 *
 * \author Peter Zajac
 */
class SimpleArgParserTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  SimpleArgParserTest() :
    TaggedTest<Archs::None, Archs::None>("SimpleArgParserTest")
  {
  }

  virtual void run() const override
  {
    // create a string map for our precon
    std::map<String, MyPrecon> precon_map;
    precon_map.insert(std::make_pair("none", MyPrecon::none));
    precon_map.insert(std::make_pair("ssor", MyPrecon::ssor));
    precon_map.insert(std::make_pair("spai", MyPrecon::spai));
    precon_map.insert(std::make_pair("iluk", MyPrecon::iluk));

    const char* argv[] =
    {
      "--ignored",            // 0; this is a placeholder for the application command
      "skip_me",              // 1; this will be skipped by the parser
      "--lvlmin", "1",        // 2, 3
      "--quiet", "-noisy",    // 4, 5
      "--foobar",             // 6
      "--range", "-2", "+7",  // 7, 8, 9
      "--precon1", "spai",    // 10, 11
      "--precon2", "mamba",   // 12, 13
      "--lvlmin", "17"        // 14, 15; overrides arguments #2, #3
    };

    // compute number of arguments
    const int argc = int(sizeof(argv) / sizeof(const char*));

    // create the argument parser
    SimpleArgParser args(argc, argv);

    // push a set of supported options
    args.support("lvlmin");
    args.support("range");
    args.support("precon1");
    args.support("quiet");
    // we leave 'foobar' and 'precon2' unsupported here

    // check for unsupported options
    std::deque<std::pair<int,String> > unsupp = args.query_unsupported();
    TEST_CHECK_EQUAL(unsupp.size(), std::size_t(2));
    TEST_CHECK_EQUAL(unsupp.front().first, 6);
    TEST_CHECK_EQUAL(unsupp.front().second, "foobar");
    TEST_CHECK_EQUAL(unsupp.back().first, 12);
    TEST_CHECK_EQUAL(unsupp.back().second, "precon2");

    // check number of (skipped) arguments
    TEST_CHECK_EQUAL(args.num_args(), argc);
    TEST_CHECK_EQUAL(args.num_skipped_args(), 2);

    // get skipped arguments
    std::deque<String> skipped(args.skipped_args());
    TEST_CHECK_EQUAL(skipped.size(), std::size_t(2));
    TEST_CHECK_EQUAL(skipped.at(0), "--ignored");
    TEST_CHECK_EQUAL(skipped.at(1), "skip_me");

    // check for existing options
    TEST_CHECK_EQUAL(args.check("quiet"), 1);
    TEST_CHECK_EQUAL(args.check("lvlmin"), 1);
    TEST_CHECK_EQUAL(args.check("foobar"), 0);
    TEST_CHECK_EQUAL(args.check("range"), 2);
    TEST_CHECK_EQUAL(args.check("precon1"), 1);
    TEST_CHECK_EQUAL(args.check("precon2"), 1);

    // check for missing options
    TEST_CHECK_EQUAL(args.check("ignored"), -1);
    TEST_CHECK_EQUAL(args.check("minerva"), -1);
    TEST_CHECK_EQUAL(args.check("athena"), -1);

    // try to query 'quiet'
    const auto* isd_quiet = args.query("quiet");
    TEST_CHECK_NOT_EQUAL(isd_quiet, nullptr);
    TEST_CHECK_EQUAL(isd_quiet->first, 4);
    TEST_CHECK_EQUAL(isd_quiet->second.size(), std::size_t(1));
    TEST_CHECK_EQUAL(isd_quiet->second.front(), "-noisy");

    // try to query a missing option
    const auto* isd_juno = args.query("juno");
    TEST_CHECK_EQUAL(isd_juno, nullptr);

    // try to parse lvlmin
    Index lvlmin(0);
    TEST_CHECK_EQUAL(args.parse("lvlmin", lvlmin), 1);
    TEST_CHECK_EQUAL(lvlmin, Index(17));

    // try to parse range 0 and 1; range 2 is missing
    int range0(0), range1(0), range2(0);
    TEST_CHECK_EQUAL(args.parse("range", range0, range1, range2), 2);
    TEST_CHECK_EQUAL(range0, -2);
    TEST_CHECK_EQUAL(range1,  7);
    TEST_CHECK_EQUAL(range2,  0);

    // try to parse precon1
    MyPrecon precon1(MyPrecon::none);
    TEST_CHECK_EQUAL(args.parse("precon1", string_mapped(precon1, precon_map)), 1);
    TEST_CHECK_EQUAL(precon1, MyPrecon::spai);

    // try to parse argument #5 "-noisy" as an int
    int noisy(0);
    TEST_CHECK_EQUAL(args.parse("quiet", noisy), -5);

    // try to parse argument #13 "mamba" as a precon
    MyPrecon precon2(MyPrecon::none);
    TEST_CHECK_EQUAL(args.parse("precon2", string_mapped(precon2, precon_map)), -13);

  }
} simple_arg_parser_test;
