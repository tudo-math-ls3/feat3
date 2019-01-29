#include <kernel/util/property_map.hpp>
#include <test_system/test_system.hpp>
#include <sstream>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the PropertyMap class.
 *
 * \test Tests the PropertyMap class.
 *
 * \author Peter Zajac
 * \author Constantin Christof
 */
class PropertyMapTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  PropertyMapTest() :
    TaggedTest<Archs::None, Archs::None>("PropertyMapTest")
  {
  }

  virtual ~PropertyMapTest()
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
    PropertyMap parsec;
    parsec.read(ioss);

    // test the root section entries
    TEST_CHECK_EQUAL(parsec.get_entry("hello").first, "World!");
    TEST_CHECK_EQUAL(parsec.get_entry("E").first, "m*c^2");
    TEST_CHECK(parsec.get_entry("KEY1").second);          // must be true as the key exists...
    TEST_CHECK_EQUAL(parsec.get_entry("KEY1").first, ""); // ...but it has no value
    TEST_CHECK_EQUAL(parsec.get_entry("key2").first, "...but then it was overwritten");

    // okay, we survived the root section, so let's go for YourSection
    PropertyMap* yoursec = parsec.get_sub_section("YourSection");
    TEST_CHECK(yoursec != nullptr);

    // let's go for pi
    PropertyMap* consec = yoursec->get_sub_section("CONSTANTS");
    TEST_CHECK(consec != nullptr);
    TEST_CHECK_EQUAL(consec->get_entry("pI").first, "3.141592653589793238462643383279...");

    // we're back in the root section; try MySection now
    PropertyMap* mysec = parsec.get_sub_section("MySeCtIoN");
    TEST_CHECK(mysec != nullptr);
    TEST_CHECK_EQUAL(mysec->get_entry("Key3").first, "This entry makes use of line continuation.");

    // test the query_section function
    PropertyMap* consec2 = yoursec->query_section("Constants"); // child
    TEST_CHECK_EQUAL(consec2, consec);
    PropertyMap* yoursec2 = consec->query_section("/~"); // parent
    TEST_CHECK_EQUAL(yoursec2, yoursec);
    PropertyMap* rootsec2 = yoursec->query_section("!/"); // root
    TEST_CHECK_EQUAL(rootsec2, &parsec);
    PropertyMap* consec3 = mysec->query_section("!/MySection/~//YourSection/Constants//");
    TEST_CHECK_EQUAL(consec3, consec);

    // test the query function
    TEST_CHECK_EQUAL(parsec.query("E", ""), "m*c^2");
    TEST_CHECK_EQUAL(parsec.query("YourSection/sqrt4", ""), "2");
    TEST_CHECK_EQUAL(parsec.query("M/B", "test"), "test"); // does not exist
    TEST_CHECK_EQUAL(consec->query("~/sqrt4", ""), "2");
    TEST_CHECK_EQUAL(mysec->query("!/E", ""), "m*c^2");

    // test the parse_entry function
    int ival = 0;
    TEST_CHECK(yoursec->parse_entry("sqrt4", ival, false)); // exists
    TEST_CHECK(yoursec->parse_entry("foobar", ival, true)); // not exists

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
    PropertyMap parsec;
    TEST_CHECK_THROWS(parsec.read(ioss), SyntaxError);
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
    PropertyMap parsec;
    TEST_CHECK_THROWS(parsec.read(ioss), SyntaxError);
  } //test_2

  void test_3() const
  {
    using namespace std;
    stringstream ioss;

    // let's write something senseless
    ioss << "6214125()/1321&%/(%????" << endl;

    // does parsing fail?
    PropertyMap parsec;
    TEST_CHECK_THROWS(parsec.read(ioss), SyntaxError);
  } //test_3

  void test_4() const
  {
    using namespace std;
    stringstream ioss;

    // let's write a key-value pair without any key
    ioss << " = 4" << endl;

    // does parsing fail?
    PropertyMap parsec;
    TEST_CHECK_THROWS(parsec.read(ioss), SyntaxError);
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
    PropertyMap parsec;
    TEST_CHECK_THROWS(parsec.read(ioss), SyntaxError);
  } //test_5

  void test_6() const
  {
    using namespace std;
    stringstream ioss1;
    stringstream ioss2;

    PropertyMap parsec1;
    PropertyMap parsec2;

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

    parsec1.read(ioss1, true); //parse first stream with replace = YES

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

    parsec2.read(ioss2, false); //parse second stream with replace = false

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
    PropertyMap* sec1 = parsec1.get_sub_section("Section1");
    TEST_CHECK(sec1 != nullptr);

    // let's go for pi in the first subsection
    PropertyMap* subsec1 = sec1->get_sub_section("Subsection1");
    TEST_CHECK(subsec1 != nullptr);
    TEST_CHECK_EQUAL(subsec1->get_entry("pI").first, "3.141592653589793238462643383279...");

    // an now for the other one
    PropertyMap* subsec2 = sec1->get_sub_section("Subsection2");
    TEST_CHECK(subsec2 != nullptr);
    TEST_CHECK_EQUAL(subsec2->get_entry("pi").first, "circa 3");

    // let's check if the Section 2 and 3 are there
    PropertyMap* sec2 = parsec1.get_sub_section("SECTION2");
    TEST_CHECK(sec2 != nullptr);
    TEST_CHECK_EQUAL(sec2->get_entry("Key3").first, "something else");

    PropertyMap* sec3 = parsec1.get_sub_section("SECTION3");
    TEST_CHECK(sec3 != nullptr);
    TEST_CHECK_EQUAL(sec3->get_entry("Key3").first, "no idea");

    //ok, everything is right
  } //test_6

  virtual void run() const override
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
} property_map_test;
