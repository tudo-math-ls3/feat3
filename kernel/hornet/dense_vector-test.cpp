// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/hornet/dense_vector.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the assertion class.
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
    DenseVector<Tag_, DT_> a(10);
    DenseVector<Tag_, DT_> b(10, 5.);
    DenseVector<Tag_, DT_> c(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
  }
};
DenseVectorTest<Nil, double> dense_vector_test_double;
