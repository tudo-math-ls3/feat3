#include <test_system/test_system.hpp>

using namespace TestSystem;
using namespace FEAST;

/**
 *
 * \brief Test class for the unittest framework itself
 *
 * \author Dirk Ribbrock
 * \test
 */
class UnitTest :
  public BaseTest
{
  public:
    /// Constructor
    UnitTest(const std::string & id) :
      BaseTest(id)
    {
    }

    /// runs the tests
    virtual void run() const
    {
      TEST_CHECK(true);
      TEST_CHECK_EQUAL(1,1);
      TEST_CHECK_NOT_EQUAL(1, 0.5);
      TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
      TEST_CHECK_EQUAL_WITHIN_EPS(25,23,2.2);
      TEST_CHECK_THROWS(std::string("0").at(10), std::exception);
    }
} unittest("UnitTest-test");

/**
 * \brief Tagged-Test class for the unittest framework itself
 *
 * \author Dirk Ribbrock
 * \test
 */
template <typename Tag_, typename DT_>
class TaggedUnitTest :
  public TaggedTest<Tag_, DT_>
{
  public:
    /// Constructor
    TaggedUnitTest(const std::string & id) :
      TaggedTest<Tag_, DT_>(id)
    {
    }

    /// runs the tests
    virtual void run() const
    {
      TEST_CHECK(true);
      TEST_CHECK_EQUAL(1,1);
      TEST_CHECK_NOT_EQUAL(1, 0.5);
      TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
      TEST_CHECK_EQUAL_WITHIN_EPS(25,23,2.2);
      TEST_CHECK_THROWS(std::string("0").at(10), std::exception);
    }
};
TaggedUnitTest<Nil, float> taggedunittestf ("TaggedUnitTest-test float");
TaggedUnitTest<Nil, double> taggedunittestd ("TaggedUnitTest-test double");
TaggedUnitTest<Nil, unsigned long> taggedunittestul ("TaggedUnitTest-test unsigned long");
TaggedUnitTest<Nil, int> taggedunittesti ("TaggedUnitTest-test int");
