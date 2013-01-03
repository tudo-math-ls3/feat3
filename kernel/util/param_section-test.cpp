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
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  ParamSectionTest() :
    TaggedTest<Archs::None, Archs::None>("param_section_test")
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
    ioss << "  sqrt4 = 2" << endl;
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
    TEST_CHECK_EQUAL(parsec.get_entry("hello").first, "World!");
    TEST_CHECK_EQUAL(parsec.get_entry("E").first, "m*c^2");
    TEST_CHECK(parsec.get_entry("KEY1").second);          // must be true as the key exists...
    TEST_CHECK_EQUAL(parsec.get_entry("KEY1").first, ""); // ...but it has no value
    TEST_CHECK_EQUAL(parsec.get_entry("key2").first, "...but then it was overwritten");

    // okay, we survived the root section, so let's go for YourSection
    ParamSection* yoursec = parsec.get_section("YourSection");
    TEST_CHECK(yoursec != nullptr);

    // let's go for pi
    ParamSection* consec = yoursec->get_section("CONSTANTS");
    TEST_CHECK(consec != nullptr);
    TEST_CHECK_EQUAL(consec->get_entry("pI").first, "3.141592653589793238462643383279...");

    // we're back in the root section; try MySection now
    ParamSection* mysec = parsec.get_section("MySeCtIoN");
    TEST_CHECK(mysec != nullptr);
    TEST_CHECK_EQUAL(mysec->get_entry("Key3").first, "This entry makes use of line continuation.");

    // test the query functions now
    TEST_CHECK_EQUAL(parsec.query("E", ""), "m*c^2");
    TEST_CHECK_EQUAL(parsec.query("YourSection.sqrt4", ""), "2");
    TEST_CHECK_EQUAL(parsec.query("M.B", "test"), "test"); // does not exist

    // okay, test passed
  } // test_0

  void test_1() const
  {
    using namespace std;
    stringstream ioss;

    // let's write an INI file into the stream
    ioss << "[YourSection]" << endl;
    ioss << "{" << endl;
    ioss << "  # nest another section into [YourSection]" << endl;
    ioss << "  [Constants]" << endl;
    ioss << "  pi = 3.1415926535& " << endl;
    ioss << "         8979323846&" << endl;
    ioss << "         2643383279..." << endl;
    ioss << "# The closing brace is missing here" << endl;

    // does parsing fail?
    ParamSection parsec;
    TEST_CHECK_THROWS(parsec.parse(ioss), SyntaxError);
  } //test_1

  void test_2() const
  {
    using namespace std;
    stringstream ioss;

    // let's write an INI file into the stream
    ioss << "[YourSection]" << endl;
    ioss << "{" << endl;
    ioss << "  # nest another section into [YourSection]" << endl;
    ioss << "  [] # the section name is missing" << endl;
    ioss << "  pi = 3.1415926535& " << endl;
    ioss << "         8979323846&" << endl;
    ioss << "         2643383279..." << endl;
    ioss << "} # End of [YourSection]" << endl;

    // does parsing fail?
    ParamSection parsec;
    TEST_CHECK_THROWS(parsec.parse(ioss), SyntaxError);
  } //test_2

  void test_3() const
  {
    using namespace std;
    stringstream ioss;

    // let's write something senseless
    ioss << "6214125()/1321&%/(%????" << endl;

    // does parsing fail?
    ParamSection parsec;
    TEST_CHECK_THROWS(parsec.parse(ioss), SyntaxError);
  } //test_3

  void test_4() const
  {
    using namespace std;
    stringstream ioss;

    // let's write a key-value pair without any key
    ioss << " = 4" << endl;

    // does parsing fail?
    ParamSection parsec;
    TEST_CHECK_THROWS(parsec.parse(ioss), SyntaxError);
  } //test_4

  void test_5() const
  {
    using namespace std;
    stringstream ioss;

    // let's write an INI file into the stream
    ioss << "# This is a comment" << endl;
    ioss << "Hello = World!" << endl;
    ioss << "[YourSection]" << endl;
    ioss << "{" << endl;
    ioss << "  { # this brace is misplaced" << endl;

    // does parsing fail?
    ParamSection parsec;
    TEST_CHECK_THROWS(parsec.parse(ioss), SyntaxError);
  } //test_5

  void test_6() const
  {
    using namespace std;
    stringstream ioss1;
    stringstream ioss2;

    ParamSection parsec1;
    ParamSection parsec2;

    // let's write an INI file into the first stream
    ioss1 << "# This is a comment" << endl;
    ioss1 << "Key1 = # This key has no value" << endl;
    ioss1 << "Key2 = I once had a another value..." << endl;
    ioss1 << "KEY2 = ...but then it was overwritten" << endl;
    ioss1 << "KEY3 = Something" << endl;
    ioss1 << "[Section1]" << endl;
    ioss1 << "{" << endl;
    ioss1 << "  # nest another section into [Section1]" << endl;
    ioss1 << "  [Subsection1]" << endl;
    ioss1 << "  pi = 3.1415926535& " << endl;
    ioss1 << "         8979323846&" << endl;
    ioss1 << "         2643383279..." << endl;
    ioss1 << "} # End of [Section1]" << endl;
    ioss1 << "# Now we're back in the root section" << endl;
    ioss1 << "[Section2] # This is another section..." << endl;
    ioss1 << "# ...without braces" << endl;
    ioss1 << "Key3 = something else" << endl;

    parsec1.parse(ioss1, true); //parse first stream with replace = YES

    TEST_CHECK(! parsec1.get_entry("key?").second);    // has to be false
    TEST_CHECK(parsec1.get_entry("key2").second);      // must be true as the key exists...
    TEST_CHECK_EQUAL(parsec1.get_entry("key2").first, "...but then it was overwritten");

    // let's write an INI file into the second stream
    ioss2 << "Key1 = 1" << endl;
    ioss2 << "Key2 = I once had a another value..." << endl;
    ioss2 << "KEY2 = ...but then it was overwritten by file 2" << endl;
    ioss2 << "KEY42 = 23" << endl;
    ioss2 << "[Section1] # the same name as above" << endl;
    ioss2 << "{" << endl;
    ioss2 << "  [Subsection2] # another subsection" << endl;
    ioss2 << "  pi = circa 3" << endl;
    ioss2 << "} # End of [Section1]" << endl;
    ioss2 << "# Now we're back in the root section" << endl;
    ioss2 << "[Section3]" << endl;
    ioss2 << "# ...without braces" << endl;
    ioss2 << "Key3 = no idea" << endl;

    parsec2.parse(ioss2, false); //parse second stream with replace = false

    TEST_CHECK(!parsec2.get_entry("key?").second);    // has to be false
    TEST_CHECK(parsec2.get_entry("key2").second);     // must be true as the key exists...
    TEST_CHECK_EQUAL(parsec2.get_entry("key2").first, "I once had a another value...");

    // test of the merge function

    parsec1.merge(parsec2, true); // merge with replace = true

    // check if everything is right

    // overwritten key-value pairs
    TEST_CHECK(parsec1.get_entry("key2").second);
    TEST_CHECK_EQUAL(parsec1.get_entry("key2").first, "I once had a another value...");
    TEST_CHECK(parsec1.get_entry("Key1").second);
    TEST_CHECK_EQUAL(parsec1.get_entry("Key1").first, "1");

    // added key-value
    TEST_CHECK(parsec1.get_entry("KEY42").second);
    TEST_CHECK_EQUAL(parsec1.get_entry("KEY42").first, "23");

    // okay, we survived the root section, so let's go for Section1
    ParamSection* sec1 = parsec1.get_section("Section1");
    TEST_CHECK(sec1 != nullptr);

    // let's go for pi in the first subsection
    ParamSection* subsec1 = sec1->get_section("Subsection1");
    TEST_CHECK(subsec1 != nullptr);
    TEST_CHECK_EQUAL(subsec1->get_entry("pI").first, "3.141592653589793238462643383279...");

    // an now for the other one
    ParamSection* subsec2 = sec1->get_section("Subsection2");
    TEST_CHECK(subsec2 != nullptr);
    TEST_CHECK_EQUAL(subsec2->get_entry("pi").first, "circa 3");

    // let's check if the Section 2 and 3 are there
    ParamSection* sec2 = parsec1.get_section("SECTION2");
    TEST_CHECK(sec2 != nullptr);
    TEST_CHECK_EQUAL(sec2->get_entry("Key3").first, "something else");

    ParamSection* sec3 = parsec1.get_section("SECTION3");
    TEST_CHECK(sec3 != nullptr);
    TEST_CHECK_EQUAL(sec3->get_entry("Key3").first, "no idea");

    //ok, everything is right
  } //test_6

  virtual void run() const
  {
    // run test #0
    test_0();

    // run test #1
    test_1();

    // run test #2
    test_2();

    // run test #3
    test_3();

    // run test #4
    test_4();

    // run test #5
    test_5();

    // run test #6
    test_6();

  }
} param_section_test;