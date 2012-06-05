// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <kernel/base_header.hpp>
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
    DenseMatrix<Tag_, DT_> c(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c.rows(), b.rows());
    DenseMatrix<Tag_, DT_> e(11, 12, 5.);
    TEST_CHECK_EQUAL(e.rows(), 11);
    TEST_CHECK_EQUAL(e.columns(), 12);
  }
};
DenseMatrixTest<Archs::CPU, double> dense_matrix_test_double;
