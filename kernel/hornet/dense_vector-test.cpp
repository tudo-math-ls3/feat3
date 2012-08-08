// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/hornet/dense_vector.hpp>
#include <list>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the dense vector class.
*
* \test test description missing
*
* \tparam Tag_
* description missing
*
* \tparam DT_
* description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Tag_,
  typename DT_>
class DenseVectorTest
  : public TaggedTest<Tag_, DT_>
{

public:

  DenseVectorTest()
    : TaggedTest<Tag_, DT_>("dense_vector_test")
  {
  }

  virtual void run() const
  {
    DenseVector<Tag_, DT_> a(10, 7.);
    DenseVector<Tag_, DT_> b(10, 5.);
    b(7, DT_(42));
    DenseVector<Tag_, DT_> c(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c(7), b(7));
    std::list<DenseVector<Tag_, DT_> > list;
    list.push_back(a);
    list.push_back(b);
    list.push_back(c);
    DenseVector<Tag_, DT_> d = a;
    list.push_back(d);
    DenseVector<Tag_, DT_> e(10, 42.);
    e = a;
    TEST_CHECK_EQUAL(e(5), a(5));
  }
};
DenseVectorTest<Archs::CPU, float> dense_vector_test_float;
DenseVectorTest<Archs::CPU, double> dense_vector_test_double;
