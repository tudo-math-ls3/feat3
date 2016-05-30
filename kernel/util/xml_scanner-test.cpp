#include <test_system/test_system.hpp>
#include <kernel/util/xml_scanner.hpp>
#include <sstream>

using namespace FEAT;
using namespace FEAT::TestSystem;

// cannot be closed
class TestParser1 :
  public Xml::DummyParser
{
public:
  virtual void create(int iline, const String& sline, const String&, const std::map<String, String>&, bool closed) override
  {
    if(closed)
      throw Xml::GrammarError(iline, sline, "Invalid closed markup");
  }
};

// cannot have content
class TestParser2 :
  public Xml::DummyParser
{
public:
  virtual bool content(int, const String&) override
  {
    return false; // no content allowed
  }
};

// Test root parser #1: create test parsers
class TestRootParser1 :
  public Xml::DummyParser
{
public:
  virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String& name) override
  {
    if(name == "Test1") return std::make_shared<TestParser1>();
    if(name == "Test2") return std::make_shared<TestParser2>();
    return nullptr;
  }
};

// Test root parser #2: check attributes
class TestRootParser2 :
  public Xml::DummyParser
{
public:
  virtual bool attribs(std::map<String,bool>& attrs) const override
  {
    attrs.emplace("mandat", true);
    attrs.emplace("option", false);
    return true;
  }

  virtual void create(int, const String&, const String&, const std::map<String, String>& attrs, bool) override
  {
    // make sure that we have the mandatory attribute
    // no XML error is thrown here, but an InternalError instead to
    // indicate failure of the code rather than an XML error
    if(attrs.find("mandat") == attrs.end())
      throw InternalError("Missing mandatory attribute!");
  }
};

/**
 * \brief Test class for the Xml::Scanner class.
 *
 * \test Tests the String class.
 *
 * \author Peter Zajac
 */
class XmlScannerTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  XmlScannerTest() :
    TaggedTest<Archs::None, Archs::None>("XmlScannerTest")
  {
  }

  void test_syntax_1() const
  {
    using namespace std;

    // empty root
    {
      stringstream ioss;
      ioss << "<>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // Bogus root
    {
      stringstream ioss;
      ioss << "<Foo bar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // closed root
    {
      stringstream ioss;
      ioss << "<Foobar />" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // terminated root
    {
      stringstream ioss;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // content root
    {
      stringstream ioss;
      ioss << "Foobar" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // invalid root markup
    {
      stringstream ioss;
      ioss << "<Foobar" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // invalid root markup
    {
      stringstream ioss;
      ioss << "Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // invalid markup name character
    {
      stringstream ioss;
      ioss << "<Foob@r>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // invalid markup name (starts with digit)
    {
      stringstream ioss;
      ioss << "<7Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.read_root(), Xml::SyntaxError);
    }
    // valid root
    {
      stringstream ioss;
      ioss << "<F00bar>" << std::endl;
      Xml::Scanner scanner(ioss);
      scanner.read_root();
      TEST_CHECK(scanner.is_cur_markup());
      TEST_CHECK(!scanner.is_cur_closed());
      TEST_CHECK(!scanner.is_cur_termin());
      TEST_CHECK_EQUAL(scanner.get_cur_name(), "F00bar");
    }
  }

  void test_syntax_2() const
  {
    using namespace std;

    // dummy root parser
    auto root_parser = std::make_shared<Xml::DummyParser>();

    // no terminator (end-of-file)
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // bogus terminator
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "</Foobar/>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // bogus terminator
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "</Foobar test=\"fail\">" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // invalid terminator
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "</Deadbeef>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // valid terminator
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "</Foobar>" << std::endl;
      // this should not throw anything
      Xml::Scanner scanner(ioss);
      scanner.scan(root_parser);
    }
    // valid terminator
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "</  Foobar  >" << std::endl;
      // this should not throw anything
      Xml::Scanner scanner(ioss);
      scanner.scan(root_parser);
    }
  }

  void test_syntax_3() const
  {
    using namespace std;

    // dummy root parser
    auto root_parser = std::make_shared<Xml::DummyParser>();

    // invalid markup (bogus data)
    {
      stringstream ioss;
      ioss << "<Foobar bla>" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // invalid markup (no attribute value)
    {
      stringstream ioss;
      ioss << "<Foobar name= >" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // invalid markup (missing quotes)
    {
      stringstream ioss;
      ioss << "<Foobar name=\" >" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // invalid markup (no attribute name)
    {
      stringstream ioss;
      ioss << "<Foobar =\"bla\" >" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // invalid markup (invalid attribute character)
    {
      stringstream ioss;
      ioss << "<Foobar f@t=\"bla\" >" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // invalid markup (invalid attribute character)
    {
      stringstream ioss;
      ioss << "<Foobar 1ft=\"bla\" >" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::SyntaxError);
    }
    // valid markup (two attributes)
    {
      stringstream ioss;
      ioss << "<Foobar dead=\"beef\" cr0w=\"\" >" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      scanner.read_root();
      TEST_CHECK( scanner.is_cur_markup());
      TEST_CHECK(!scanner.is_cur_termin());
      TEST_CHECK(!scanner.is_cur_closed());
      TEST_CHECK_EQUAL(scanner.get_cur_name(), "Foobar");
      const auto& attr = scanner.get_cur_attribs();
      TEST_CHECK_EQUAL(int(attr.size()), 2);
      auto vit = attr.begin();
      TEST_CHECK_EQUAL(vit->first, "cr0w");
      TEST_CHECK      (vit->second.empty());
      ++vit;
      TEST_CHECK_EQUAL(vit->first, "dead");
      TEST_CHECK_EQUAL(vit->second, "beef");
      // this should not throw anything
      scanner.set_root_parser(root_parser);
      scanner.scan();
    }
  }

  void test_grammar_1() const
  {
    using namespace std;

    // test root parser
    auto root_parser = std::make_shared<TestRootParser1>();

    // unsupported markup
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "  <Unknown />" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::GrammarError);
    }
    // invalid closed markup
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "  <Test1 />" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::GrammarError);
    }
    // invalid content in markup
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "  <Test2>" << std::endl;
      ioss << "    content" << std::endl;
      ioss << "  </Test2>" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::GrammarError);
    }
    // valid
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "  <Test1>" << std::endl;
      ioss << "    content" << std::endl;
      ioss << "    <Bla/>" << std::endl;
      ioss << "  </Test1>" << std::endl;
      ioss << "  <Test2/>" << std::endl;
      ioss << "  <Test2>" << std::endl;
      ioss << "    <Bla/>" << std::endl;
      ioss << "  </Test2>" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      // this should not throw anything
      scanner.scan(root_parser);
    }
  }

  void test_grammar_2() const
  {
    using namespace std;

    // test root parser
    auto root_parser = std::make_shared<TestRootParser2>();

    // missing mandatory attribute
    {
      stringstream ioss;
      ioss << "<Foobar>" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::GrammarError);
    }
    // unsupported attribute
    {
      stringstream ioss;
      ioss << "<Foobar mandat=\"\" unknown=\"attribute\">" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      TEST_CHECK_THROWS(scanner.scan(root_parser), Xml::GrammarError);
    }
    // valid attributes
    {
      stringstream ioss;
      ioss << "<Foobar mandat=\"given\">" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      // this should not throw anything
      scanner.scan(root_parser);
    }
    // valid attributes
    {
      stringstream ioss;
      ioss << "<Foobar mandat=\"given\" option=\"also given\">" << std::endl;
      ioss << "</Foobar>" << std::endl;
      Xml::Scanner scanner(ioss);
      // this should not throw anything
      scanner.scan(root_parser);
    }
  }

  virtual void run() const override
  {
    test_syntax_1();
    test_syntax_2();
    test_syntax_3();
    test_grammar_1();
    test_grammar_2();
  }
} xml_scanner_test;
