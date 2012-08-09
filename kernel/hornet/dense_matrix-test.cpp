// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/hornet/dense_matrix.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the dense matrix class.
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
class DenseMatrixTest
  : public TaggedTest<Tag_, DT_>
{

public:

  DenseMatrixTest()
    : TaggedTest<Tag_, DT_>("dense_matrix_test")
  {
  }

  virtual void run() const
  {
    DenseMatrix<Tag_, DT_> a(10, 10);
    DenseMatrix<Tag_, DT_> b(10, 10, 5.);
    b(7, 6, DT_(42));
    DenseMatrix<Tag_, DT_> c(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c.rows(), b.rows());
    TEST_CHECK_EQUAL(c(7,6), b(7,6));
    TEST_CHECK_EQUAL(c, b);

    DenseMatrix<Tag_, DT_> e(11, 12, 5.);
    TEST_CHECK_EQUAL(e.rows(), 11ul);
    TEST_CHECK_EQUAL(e.columns(), 12ul);

    DenseMatrix<Tag_, DT_> f(11, 12, 42.);
    f = e;
    TEST_CHECK_EQUAL(f(7,8), e(7,8));
    TEST_CHECK_EQUAL(f, e);
    std::cout<<f;
  }
};
DenseMatrixTest<Archs::CPU, float> dense_matrix_test_float;
DenseMatrixTest<Archs::CPU, double> dense_matrix_test_double;
