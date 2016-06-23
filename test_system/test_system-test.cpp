// includes, FEAT
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

class A
{
  public:
  bool flag;

  explicit A(bool f) :
    flag(f)
  {
  }

  A(const A & other) = delete;
};

bool operator==(const A & a, const A & b)
{
  return a.flag == b.flag;
}

std::ostream & operator<< (std::ostream & lhs, const A & b)
{
  lhs << "A: "<< b.flag << std::endl;
  return lhs;
}



/**
*
* \brief test class for the unittest framework itself
*
* \test test description missing
*
* \author Dirk Ribbrock
*/
class UnitTest
  : public BaseTest
{
public:
  /// Constructor
  explicit UnitTest(const String & id_in)
    : BaseTest(id_in)
  {
  }

  /// runs the tests
  virtual void run() const override
  {
    TEST_CHECK(true);
    TEST_CHECK_EQUAL(1,1);
    TEST_CHECK_NOT_EQUAL(1, 0.5);
    TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
    TEST_CHECK_EQUAL_WITHIN_EPS(25,23,2.2);
    TEST_CHECK_EQUAL_WITHIN_EPS(6, 7, 1);
    TEST_CHECK_THROWS(String("0").at(10), std::exception);

    A a(true);
    A b(true);
    TEST_CHECK_EQUAL(a, b);
    A c(false);
    TEST_CHECK_NOT_EQUAL(a, c);
  }
} unittest("UnitTest-test");



/**
* \brief tagged-test class for the unittest framework itself
*
* \test
* test description missing
*
* \author Dirk Ribbrock
*/
template <typename Tag_, typename DT_>
class TaggedUnitTest
  : public TaggedTest<Tag_, DT_>
{
public:
  /// Constructor
  explicit TaggedUnitTest(const String & id_in)
    : TaggedTest<Tag_, DT_>(id_in)
  {
  }

  /// runs the tests
  virtual void run() const override
  {
    TEST_CHECK(true);
    TEST_CHECK_EQUAL(1,1);
    TEST_CHECK_NOT_EQUAL(1, 0.5);
    TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
    TEST_CHECK_EQUAL_WITHIN_EPS(25,23,2.2);
    TEST_CHECK_THROWS(String("0").at(10), std::exception);

    A a(true);
    A b(true);
    TEST_CHECK_EQUAL(a, b);
    A c(false);
    TEST_CHECK_NOT_EQUAL(a, c);
  }
};
TaggedUnitTest<Archs::None, float> taggedunittestf ("TaggedUnitTest-test float");
TaggedUnitTest<Archs::None, double> taggedunittestd ("TaggedUnitTest-test double");
TaggedUnitTest<Archs::None, unsigned long> taggedunittestul ("TaggedUnitTest-test unsigned long");
TaggedUnitTest<Archs::None, int> taggedunittesti ("TaggedUnitTest-test int");

/**
* \brief full-tagged-test class for the unittest framework itself
*
* \test
* test description missing
*
* \author Dirk Ribbrock
*/
template <typename Mem_, typename DT_, typename IT_>
class FullTaggedUnitTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  /// Constructor
  explicit FullTaggedUnitTest(const String & id_in)
    : FullTaggedTest<Mem_, DT_, IT_>(id_in)
  {
  }

  /// runs the tests
  virtual void run() const override
  {
    TEST_CHECK(true);
    TEST_CHECK_EQUAL(1,1);
    TEST_CHECK_NOT_EQUAL(1, 0.5);
    TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
    TEST_CHECK_EQUAL_WITHIN_EPS(25,23,2.2);
    TEST_CHECK_THROWS(String("0").at(10), std::exception);

    A a(true);
    A b(true);
    TEST_CHECK_EQUAL(a, b);
    A c(false);
    TEST_CHECK_NOT_EQUAL(a, c);
  }
};
FullTaggedUnitTest<Archs::None, float, Index> fulltaggedunittestf ("FullTaggedUnitTest-test float");
FullTaggedUnitTest<Archs::None, double, Index> fulltaggedunittestd ("FullTaggedUnitTest-test double");
