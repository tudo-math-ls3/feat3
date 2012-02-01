#include <kernel/util/param_section.hpp>
#include <test_system/test_system.hpp>
#include <sstream>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for the ParamSection class.
 *
 * \test Tests the ParamSection class.
 *
 * \author Peter Zajac
 */
class ParamSectionTest
  : public TaggedTest<Nil, Nil>
{
public:
  ParamSectionTest() :
    TaggedTest<Nil, Nil>("param_section_test")
  {
  }

  void test_0() const
  {
    using namespace std;
    stringstream ioss;

    // let's write an INI file into the stream
    ioss << "# This is a comment" << endl;
    ioss << "Hello = World!" << endl;
    ioss << "E = m*c^2 # = h*f" << endl;
    ioss << "Key1 = # This key has no value" << endl;
    ioss << "Key2 = I once had a another value..." << endl;
    ioss << "KEY2 = ...but then it was overwritten" << endl;
    ioss << "[YourSection]" << endl;
    ioss << "{" << endl;
    ioss << "  # nest another section into [YourSection]" << endl;
    ioss << "  [Constants]" << endl;
    ioss << "  pi = 3.1415926535& " << endl;
    ioss << "         8979323846&" << endl;
    ioss << "         2643383279..." << endl;
    ioss << "} # End of [YourSection]" << endl;
    ioss << "# Now we're back in the root section" << endl;
    ioss << /* This is an empty line */ endl;
    ioss << "[MySection] # This is a section..." << endl;
    ioss << "# ...without braces" << endl;
    ioss << "Key3 = This entry & # ruthlessly" << endl;
    ioss << "       makes use of &" << endl;
    ioss << "     # comments within" << endl;
    ioss << "       line continuation." << endl;

    // parse the stream
    ParamSection parsec;
    parsec.parse(ioss);

    // test the root section entries
    TEST_CHECK(parsec.get_entry("hello").first == "World!");
    TEST_CHECK(parsec.get_entry("E").first == "m*c^2");
    TEST_CHECK(parsec.get_entry("KEY1").second);      // must be true as the key exists...
    TEST_CHECK(parsec.get_entry("KEY1").first == ""); // ...but it has no value
    TEST_CHECK(parsec.get_entry("key2").first == "...but then it was overwritten");

    // okay, we survived the root section, so let's go for YourSection
    ParamSection* yoursec = parsec.get_section("YourSection");
    TEST_CHECK(yoursec != nullptr);

    // let's go for pi
    ParamSection* consec = yoursec->get_section("CONSTANTS");
    TEST_CHECK(consec != nullptr);
    TEST_CHECK(consec->get_entry("pI").first == "3.141592653589793238462643383279...");

    // we're back in the root section; try MySection now
    ParamSection* mysec = parsec.get_section("MySeCtIoN");
    TEST_CHECK(mysec != nullptr);
    TEST_CHECK(mysec->get_entry("Key3").first == "This entry makes use of line continuation.");

    // okay, test passed
  }

  virtual void run() const
  {
    // run test #0
    test_0();

    // there are more to come...
  }
} param_section_test;
