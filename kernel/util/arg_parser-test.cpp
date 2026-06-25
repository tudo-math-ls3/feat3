// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <deque>
#include <kernel/util/arg_parser.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <stdexcept>
#include <test_system/test_system.hpp>

#include <array>
#include <cstdint>
#include <cstdlib>

#ifdef _WIN32
#include <stdlib.h>
#endif

using namespace FEAT;
using namespace FEAT::TestSystem;

class ArgParserTest : public TestSystem::UnitTest
{
public:
  ArgParserTest() : TestSystem::UnitTest("ArgParserTest")
  {
  }

  ~ArgParserTest() override = default;

  struct TestArgParser : public ArgParser
  {
    Parameter<bool> verbose = parameter(false)
                                .long_flag("--verbose")
                                .short_flag("-v")
                                .env("FEAT_VERBOSE")
                                .property("verbose")
                                .help_string("Enable additional output");

    Parameter<std::int32_t> count =
      parameter(10)
        .long_flag("--count")
        .short_flag("-c")
        .env("FEAT_COUNT")
        .property("count")
        .help_string(
          "Set count for smoother steps.\nIf the used solver is a multigrid solver, this sets the number of smoothing "
          "steps on each layer.\nIf the used solver is a GMRES, this sets the number of iterations.")
        .placeholder("solver-iters");

    Parameter<std::uint16_t> n = parameter<std::uint16_t>().short_flag("-n").env("FEAT_N").property("n");

    Parameter<bool> no_env = parameter(false).property("no-env");

    Parameter<std::int32_t> baz = parameter(0).property("foo/bar/baz");

    Parameter<double> scale = parameter(1.0).name("scale");
  };

  void run() const override
  {
    test_metaprogramming();
    test_flag_definition();
    test_parsing();
    test_overwrite();
    test_priority();
    test_multi_option();
    test_map_parsing();
    test_validation();
    test_custom_parser();
    test_required();
    test_duplicate_parameter();
    test_needs();
    test_derived_params();
    test_vector_parsing();
    test_matrix_parsing();
    test_compose();
    test_unknown_flags();
    test_adding_to_property_map();
  }

  void test_metaprogramming() const
  {
    TEST_CHECK((Intern::can_default_parse_v<Real>));
    TEST_CHECK((Intern::can_default_parse_v<String>));
    TEST_CHECK((Intern::can_default_parse_v<bool>));
    TEST_CHECK((Intern::can_default_parse_v<std::uint8_t>));

    struct Foo
    {
    };

    TEST_CHECK((!Intern::can_default_parse_v<Foo>));
  }

  void test_flag_definition() const
  {
    // Check that both --foo and foo (as well as -f and f) are accepted as flags
    struct Arg : public ArgParser
    {
      Parameter<int> a = parameter(0).long_flag("--foo").short_flag("-f");
      Parameter<int> b = parameter(0).long_flag("bar").short_flag("b");
    };

    // Definition did not abort
    TEST_CHECK(true);
  }

  void test_parsing() const
  {
    std::array<const char*, 4> args = {"/path/to/program", "--verbose", "--count", "20"};

    PropertyMap map;
    map.add_entry("no-env", "true");
    map.add_section("foo")->add_section("bar")->add_entry("baz", "42");

    set_env("FEAT_N", "20", 1);

    TestArgParser parser;
    if(!parser.parse(int(args.size()), args.data(), &map))
    {
      for(const auto& error : parser.errors())
      {
        std::cout << error << "\n";
      }
    }

    TEST_CHECK_EQUAL(*parser.verbose, true);
    TEST_CHECK_EQUAL(*parser.count, 20);
    TEST_CHECK_EQUAL(*parser.n, 20);
    TEST_CHECK_EQUAL(*parser.baz, 42);
    TEST_CHECK_EQUAL(*parser.no_env, true);
  }

  /// Reusing the same parser multiple times should behave just like creating a new parser
  void test_overwrite() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::int32_t> arg = parameter(0).long_flag("--arg").short_flag("-a").env("FEAT_ARG").property("arg");
    };

    std::array<const char*, 3> argv = {"/path/to/program", "--arg", "1"};

    Args args;

    TEST_CHECK(args.parse(int(argv.size()), argv.data()));
    TEST_CHECK_EQUAL(*args.arg, 1);

    // Parsing without any arguments should set defaults again
    TEST_CHECK(args.parse(0, nullptr));
    TEST_CHECK_EQUAL(*args.arg, 0);

    set_env("FEAT_ARG", "2", 1);

    TEST_CHECK(args.parse(0, nullptr));
    TEST_CHECK_EQUAL(*args.arg, 2);
  }

  void test_priority() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::int32_t> arg = parameter(0).long_flag("--arg").short_flag("-a").env("FEAT_ARG").property("arg");
    };

    std::array<const char*, 3> argv = {
      "path/to/program",
      "-a",
      "3",
    };

    set_env("FEAT_ARG", "1", 1);

    PropertyMap map;
    map.add_entry("arg", "2");

    Args args;

    // Environment variable should have priority over default value
    TEST_CHECK(args.parse(0, nullptr));
    TEST_CHECK_EQUAL(*args.arg, 1);

    // Property map should have priority over environment variable
    TEST_CHECK(args.parse(0, nullptr, &map));
    TEST_CHECK_EQUAL(*args.arg, 2);

    // Command line arg should have priority over property map
    TEST_CHECK(args.parse(int(argv.size()), argv.data(), &map));
    TEST_CHECK_EQUAL(*args.arg, 3);
  }

  void test_multi_option() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::deque<double>> arg =
        parameter(std::deque<double>{}).long_flag("--arg").short_flag("-a").env("FEAT_ARG").property("arg");
      Parameter<std::int32_t> arg2 = parameter(0).long_flag("--b").short_flag("-b").env("FEAT_B").property("b");
    };

    std::array<const char*, 8> argv = {
      "path/to/program",
      "-a",
      "1",
      "2",
      "3",
      "4",
      "-b",
      "10",
    };

    Args args;
    TEST_CHECK(args.parse(int(argv.size()), argv.data()));

    // -a flag should have collected the 4 following values
    TEST_CHECK_EQUAL(args.arg->size(), 4);
    TEST_CHECK_EQUAL((*args.arg)[0], 1.0);
    TEST_CHECK_EQUAL((*args.arg)[1], 2.0);
    TEST_CHECK_EQUAL((*args.arg)[2], 3.0);
    TEST_CHECK_EQUAL((*args.arg)[3], 4.0);

    // -b flag should have stopped collection of arguments for -a flag
    TEST_CHECK_EQUAL(args.arg2.value(), 10);
  }

  void test_map_parsing() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::map<String, int>> string_int_map = parameter(std::map<String, int>{}).long_flag("--string-int-map");
      Parameter<std::map<int, String>> int_string_map = parameter(std::map<int, String>{}).long_flag("--int-string-map");
    };

    std::array<const char*, 15> argv = {
      "path/to/program",
      "--string-int-map",
      "Foo",
      "1",
      "Bar",
      "2",
      "Baz",
      "3",
      "--int-string-map",
      "1",
      "Foo",
      "2",
      "Bar",
      "3",
      "Baz"
    };

    Args args;
    TEST_CHECK(args.parse(int(argv.size()), argv.data()));

    TEST_CHECK_EQUAL(args.string_int_map->size(), 3);
    TEST_CHECK_EQUAL((*args.string_int_map).at("Foo"), 1);
    TEST_CHECK_EQUAL((*args.string_int_map).at("Bar"), 2);
    TEST_CHECK_EQUAL((*args.string_int_map).at("Baz"), 3);

    TEST_CHECK_EQUAL(args.int_string_map->size(), 3);
    TEST_CHECK_EQUAL((*args.int_string_map).at(1), String("Foo"));
    TEST_CHECK_EQUAL((*args.int_string_map).at(2), String("Bar"));
    TEST_CHECK_EQUAL((*args.int_string_map).at(3), String("Baz"));
  }

  void test_validation() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::int32_t> arg =
        parameter(std::int32_t(0)).long_flag("--arg").validator([](const auto& i) { return i < 10; });
      Parameter<String> arg2 = parameter(String("foo"))
                                 .long_flag("--arg2")
                                 .validator(
                                   [](const auto& s)
                                   {
                                     if(s != "foo")
                                     {
                                       throw std::invalid_argument("Bad string");
                                     }
                                   });
    };

    Args args;

    // Should be good
    {
      std::array<const char*, 3> argv_good = {
        "path/to/program",
        "--arg",
        "1",
      };
      TEST_CHECK(args.parse(int(argv_good.size()), argv_good.data()));
      TEST_CHECK_EQUAL(args.arg.value(), 1);
    }

    // Should fail due to unknown flag
    {
      std::array<const char*, 5> argv_unknown = {"path/to/program", "--unknown", "foo", "--arg", "1"};
      TEST_CHECK(!args.parse(int(argv_unknown.size()), argv_unknown.data()));
      // But properly given arg should be fine
      TEST_CHECK_EQUAL(args.arg.value(), 1);
    }

    // Should fail due to parse error
    {
      std::array<const char*, 3> argv_parseerror = {
        "path/to/program",
        "--arg",
        "foo",
      };
      TEST_CHECK(!args.parse(int(argv_parseerror.size()), argv_parseerror.data()));
      // Parameter should still be default initialized
      TEST_CHECK_EQUAL(args.arg.value(), 0);
    }

    // Should fail due to validator
    {
      std::array<const char*, 3> argv_badvalidator = {
        "path/to/program",
        "--arg",
        "20",
      };
      TEST_CHECK(!args.parse(int(argv_badvalidator.size()), argv_badvalidator.data()));
      TEST_CHECK_EQUAL(args.arg.value(), 20);
    }

    // Should fail due to validator
    {
      std::array<const char*, 3> argv_exception = {"path/to/program", "--arg2", "baz"};
      TEST_CHECK(!args.parse(int(argv_exception.size()), argv_exception.data()));
    }

    // Setting parameter multiple times is an error
    {
      std::array<const char*, 5> argv = {"path/to/program", "--arg", "1", "--arg", "2"};

      TEST_CHECK(!args.parse(int(argv.size()), argv.data()));
    }
  }

  void test_custom_parser() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::int32_t> arg = parameter(std::int32_t(0))
                                      .long_flag("--arg")
                                      .parser(+[](const String& s)
                                              {
                                                if(s == "foo")
                                                {
                                                  return 10;
                                                }
                                                else
                                                {
                                                  return 1;
                                                }
                                              });
    };

    Args args;

    std::array<const char*, 3> argv_foo = {"path/to/program", "--arg", "foo"};

    std::array<const char*, 3> argv = {"path/to/program", "--arg", "baz"};

    TEST_CHECK(args.parse(int(argv_foo.size()), argv_foo.data()));
    TEST_CHECK_EQUAL(*args.arg, 10);

    TEST_CHECK(args.parse(int(argv.size()), argv.data()));
    TEST_CHECK_EQUAL(*args.arg, 1);
  }

  void test_required() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::int32_t> arg = parameter(std::int32_t(0)).long_flag("--arg").required();
    };

    std::array<const char*, 1> argv_bad = {
      "path/to/program",
    };

    std::array<const char*, 3> argv_good = {"path/to/program", "--arg", "15"};

    Args args;
    TEST_CHECK(!args.parse(int(argv_bad.size()), argv_bad.data()));

    TEST_CHECK(args.parse(int(argv_good.size()), argv_good.data()));
    TEST_CHECK_EQUAL(*args.arg, 15);
  }

  void test_duplicate_parameter() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::int32_t> arg = parameter(std::int32_t(0)).long_flag("--arg").required();
      Parameter<std::int32_t> arg2 = parameter(std::int32_t(0)).long_flag("--arg").required();
    };

    Args args;
    TEST_CHECK(!args.parse(0, nullptr));
  }

  void test_needs() const
  {
    struct Args : public ArgParser
    {
      Parameter<std::int32_t> a = parameter(std::int32_t(0)).short_flag("-a");
      Parameter<std::int32_t> b = parameter(std::int32_t(1)).short_flag("-b").needs(a);
      Parameter<std::int32_t> c = parameter(std::int32_t(2))
                                    .short_flag("-c")
                                    .needs_if(a, [](const auto& self) { return self == 5; })
                                    .needs_if(b, [](const auto& self) { return self == 10; });
    };

    Args args;

    // Setting just a is ok
    {
      std::array<const char*, 3> argv = {"path/to/program", "-a", "10"};
      TEST_CHECK(args.parse(int(argv.size()), argv.data()));
    }

    // Just setting b is forbidden
    {
      std::array<const char*, 3> argv = {"path/to/program", "-b", "10"};

      TEST_CHECK(!args.parse(int(argv.size()), argv.data()));
    }

    // Setting a and b is ok
    {
      std::array<const char*, 5> argv = {"path/to/program", "-a", "10", "-b", "10"};
      TEST_CHECK(args.parse(int(argv.size()), argv.data()));
    }

    // Setting c to 5 without setting a is forbidden
    {
      std::array<const char*, 3> argv = {"path/to/program", "-c", "5"};
      TEST_CHECK(!args.parse(int(argv.size()), argv.data()));
    }

    // Setting c to 5 with setting a is ok
    {
      std::array<const char*, 5> argv = {"path/to/program", "-a", "10", "-c", "5"};
      TEST_CHECK(args.parse(int(argv.size()), argv.data()));
    }

    // Setting c to 10 without setting b is forbidden
    {
      std::array<const char*, 3> argv = {"path/to/program", "-c", "10"};
      TEST_CHECK(!args.parse(int(argv.size()), argv.data()));
    }

    // Setting c to 10 with setting b is ok
    // Note that b then requires a, thus all parameters have to be given
    {
      std::array<const char*, 7> argv = {"path/to/program", "-a", "10", "-b", "10", "-c", "10"};
      TEST_CHECK(args.parse(int(argv.size()), argv.data()));
    }
  }

  void test_derived_params() const
  {
    struct Args : public ArgParser
    {
      Parameter<bool> flag = parameter(false).long_flag("--flag");
      Parameter<bool> not_flag = parameter(true).derives_from(flag, [](const bool& bb) { return !bb; });

      Parameter<int> a = parameter(0).long_flag("--a");
      Parameter<int> b = parameter(0).long_flag("--b");
      Parameter<int> c = parameter(0).long_flag("--c").derives_from(a, b, [](const int& foo, const int& bar) { return foo + bar; });
      Parameter<int> d = parameter(0).derives_from(flag, a, b, [](const bool& sub, const int& foo, const int& bar) { return sub ? foo + bar : foo - bar; });
    };

    // Test that derived values are correct
    {
      Args args;
      std::array<const char*, 6> argv = {
        "path/to/program",
        "--flag",
        "--a",
        "5",
        "--b",
        "10",
      };

      TEST_CHECK(args.parse(int(argv.size()), argv.data()));

      TEST_CHECK_EQUAL(*args.not_flag, false);
      TEST_CHECK_EQUAL(*args.c, 15);
      TEST_CHECK_EQUAL(*args.d, 15);
    }

    // Test that derived values are not overwriting values set by user
    {
      Args args;
      std::array<const char*, 8> argv = {
        "path/to/program",
        "--flag",
        "--a",
        "5",
        "--b",
        "10",
        "--c",
        "50",
      };

      TEST_CHECK(args.parse(int(argv.size()), argv.data()));
      TEST_CHECK_EQUAL(*args.c, 50);
    }
  }

  void test_vector_parsing() const
  {
    struct Args : public ArgParser
    {
      Parameter<Tiny::Vector<double, 3>> a = parameter(Tiny::Vector<double, 3>{0.0, 0.0, 0.0}).short_flag("-a");
    };

    Args args;

    std::array<const char*, 3> argv = {"path/to/program", "-a", "[ 10 20 30 ]"};

    TEST_CHECK(args.parse(int(argv.size()), argv.data()));

    TEST_CHECK_EQUAL((*args.a)(0), 10.0);
    TEST_CHECK_EQUAL((*args.a)(1), 20.0);
    TEST_CHECK_EQUAL((*args.a)(2), 30.0);
  }

  void test_matrix_parsing() const
  {
    struct Args : public ArgParser
    {
      Parameter<Tiny::Matrix<double, 2, 2>> a = parameter(Tiny::Matrix<double, 2, 2>{}).short_flag("-a");
    };

    Args args;

    std::array<const char*, 3> argv = {"path/to/program", "-a", "[ 10 20 ] [ 30 40 ]"};

    TEST_CHECK(args.parse(int(argv.size()), argv.data()));

    TEST_CHECK_EQUAL((*args.a)(0, 0), 10.0);
    TEST_CHECK_EQUAL((*args.a)(0, 1), 20.0);
    TEST_CHECK_EQUAL((*args.a)(1, 0), 30.0);
    TEST_CHECK_EQUAL((*args.a)(1, 1), 40.0);
  }

  void test_compose() const
  {
    struct A : public ArgParser
    {
      using ArgParser::ArgParser;

      Parameter<int> a = parameter(0).long_flag("--option");
      Parameter<int> b = parameter(0).property("prop");
      Parameter<int> c = parameter(0).env("ENV");
    };

    struct B : public ArgParser
    {
      using ArgParser::ArgParser;

      A a = parser<A>("A", "a", "a");
    };

    struct C : public ArgParser
    {
      using ArgParser::ArgParser;

      B b = parser<B>("B", "b", "b");
    };

    {
      C args_cli;
      std::array<const char*, 3> argv = {"path/to/program", "--b-a-option", "5"};
      TEST_CHECK(args_cli.parse(int(argv.size()), argv.data()));
      TEST_CHECK_EQUAL((*args_cli.b.a.a), 5);
    }

    {
      C args_prop;
      PropertyMap pmap;
      auto* b = pmap.add_section("b");
      auto* a = b->add_section("a");
      a->add_entry("prop", "5");
      TEST_CHECK(args_prop.parse(0, nullptr, &pmap));
      TEST_CHECK_EQUAL((*args_prop.b.a.b), 5);
    }

    {
      C args_env;
      set_env("B_A_ENV", "5", 1);
      TEST_CHECK(args_env.parse(0, nullptr));
      TEST_CHECK_EQUAL((*args_env.b.a.c), 5);
    }
  }

  void test_unknown_flags() const
  {
    struct Args : public ArgParser
    {
      Parameter<int> foo = parameter(0).long_flag("--foo");
    };

    // Parse with unknown flag should fail
    {
      Args args;
      std::array<const char*, 5> argv = {"path/to/program", "--foo", "1", "--bar" , "5"};

      TEST_CHECK(!args.parse(int(argv.size()), argv.data()));
    }

    // Parse with unknown flag shouldn't fail, if allow_unknown_flags is set
    {
      Args args;
      std::array<const char*, 5> argv = {"path/to/program", "--foo", "1", "--bar" , "5"};

      args.allow_unknown_cli_flags();
      TEST_CHECK(args.parse(int(argv.size()), argv.data()));
    }
  }

  void test_adding_to_property_map() const
  {
    struct Args : public ArgParser
    {
      Parameter<int> a = parameter(0).property("foo/bar/a").long_flag("--a");
      Parameter<int> b = parameter(1).property("foo/bar/b").long_flag("--b");
      Parameter<int> c = parameter(2).property("foo/c").long_flag("--c");
    };

    {
      Args args;
      std::array<const char*, 7> argv = {"path/to/program", "--a", "5", "--b" , "10", "--c", "15"};

      TEST_CHECK(args.parse(int(argv.size()), argv.data()));

      PropertyMap default_map;
      PropertyMap value_map;

      args.add_values_to_property_map(value_map);
      args.add_defaults_to_property_map(default_map);

      TEST_CHECK_NOT_EQUAL(default_map.query_section("foo"), nullptr);
      TEST_CHECK_NOT_EQUAL(default_map.query_section("foo/bar"), nullptr);
      TEST_CHECK_EQUAL(default_map.query("foo/bar/a").first, "0");
      TEST_CHECK_EQUAL(default_map.query("foo/bar/b").first, "1");
      TEST_CHECK_EQUAL(default_map.query("foo/c").first, "2");

      TEST_CHECK_NOT_EQUAL(value_map.query_section("foo"), nullptr);
      TEST_CHECK_NOT_EQUAL(value_map.query_section("foo/bar"), nullptr);
      TEST_CHECK_EQUAL(value_map.query("foo/bar/a").first, "5");
      TEST_CHECK_EQUAL(value_map.query("foo/bar/b").first, "10");
      TEST_CHECK_EQUAL(value_map.query("foo/c").first, "15");
    }
  }
private:
  static int set_env(const char* name, const char* value, int overwrite)
  {
#ifdef _WIN32
    (void)overwrite; // Silence unused variable warning
    return _putenv_s(name, value);
#else
    return setenv(name, value, overwrite);
#endif
  }
};

static const ArgParserTest arg_parser_test;
