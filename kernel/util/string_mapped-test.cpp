#include <test_system/test_system.hpp>
#include <kernel/util/string_mapped.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

// an enumeration for testing purposes
enum class MyPrecon
{
  none,
  ssor,
  spai,
  iluk,
  juno
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
  case MyPrecon::juno: return os << "juno";
  default: return os;
  }
}

/**
 * \brief Test class for the StringMapped class template.
 *
 * \test Tests the StringMapped class.
 *
 * \author Peter Zajac
 */
class StringMappedTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  StringMappedTest() :
    TaggedTest<Archs::None, Archs::None>("StringMappedTest")
  {
  }

  virtual ~StringMappedTest()
  {
  }

  virtual void run() const override
  {
    // create a string-map for our preconditioner and
    // insert all supported enumeration values into it
    std::map<String, MyPrecon> my_precon_map;
    my_precon_map.insert(std::make_pair("none", MyPrecon::none));
    my_precon_map.insert(std::make_pair("ssor", MyPrecon::ssor));
    my_precon_map.insert(std::make_pair("spai", MyPrecon::spai));
    my_precon_map.insert(std::make_pair("iluk", MyPrecon::iluk));
    my_precon_map.insert(std::make_pair("juno", MyPrecon::juno));

    // create a few enumeration objects and initialise them to 'none'
    MyPrecon my_precon1(MyPrecon::none);
    MyPrecon my_precon2(MyPrecon::none);
    MyPrecon my_precon3(MyPrecon::none);
    MyPrecon my_precon4(MyPrecon::none);
    MyPrecon my_precon5(MyPrecon::none);
    MyPrecon my_precon6(MyPrecon::none);
    MyPrecon my_precon7(MyPrecon::none);
    MyPrecon my_precon8(MyPrecon::none);

    // Variant 1: lookup 'my_precon1' by using the 'string_maped_lookup' function
    TEST_CHECK(string_mapped_lookup(my_precon1, my_precon_map, "ssor"));
    TEST_CHECK_EQUAL(my_precon1, MyPrecon::ssor);

    // Variant 2: lookup 'my_precon2' by creating a StringMappedValue object
    //            and using its 'lookup' member function
    TEST_CHECK(string_mapped(my_precon2, my_precon_map).lookup("spai"));
    TEST_CHECK_EQUAL(my_precon2, MyPrecon::spai);

    // Variant 3: lookup 'my_precon3' by using the 'parse' member function of a String object
    auto my_mapped_precon3 = string_mapped(my_precon3, my_precon_map);
    TEST_CHECK(String("iluk").parse(my_mapped_precon3));
    TEST_CHECK_EQUAL(my_precon3, MyPrecon::iluk);

    // Variant 4: lookup 'my_precon4' by using the 'operator>>' of a std::istream object
    std::istringstream my_istream1("juno");
    TEST_CHECK(!(my_istream1 >> string_mapped(my_precon4, my_precon_map)).fail());
    TEST_CHECK_EQUAL(my_precon4, MyPrecon::juno);

    // check failures
    TEST_CHECK(!string_mapped_lookup(my_precon5, my_precon_map, "io"));
    TEST_CHECK(!string_mapped(my_precon6, my_precon_map).lookup("europa"));
    auto my_mapped_precon7 = string_mapped(my_precon7, my_precon_map);
    TEST_CHECK(!String("ganimede").parse(my_mapped_precon7));
    std::istringstream my_istream2("callisto");
    TEST_CHECK((my_istream2 >> string_mapped(my_precon8, my_precon_map)).fail());
  }
} string_mapped_test;
