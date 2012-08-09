#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/hornet/sparse_matrix_coo.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the sparse matrix coo class.
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
class SparseMatrixCOOTest
  : public TaggedTest<Tag_, DT_>
{

public:

  SparseMatrixCOOTest()
    : TaggedTest<Tag_, DT_>("sparse_matrix_coo_test")
  {
  }

  virtual void run() const
  {
    SparseMatrixCOO<Tag_, DT_> a(10, 10);
    a(1,2,7);
    a(5,5,2);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);
    SparseMatrixCOO<Tag_, DT_> b(a);
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(a(1,2), b(1,2));
    TEST_CHECK_EQUAL(a(0,2), b(0,2));
    TEST_CHECK_EQUAL(a, b);

    SparseMatrixCOO<Tag_, DT_> c(10, 10);
    c = b;
    TEST_CHECK_EQUAL(c(0,2), b(0,2));
    TEST_CHECK_EQUAL(c, b);
    TEST_CHECK_EQUAL(c.used_elements(), b.used_elements());
    std::cout<<c;
  }
};
SparseMatrixCOOTest<Archs::CPU, double> sparse_matrix_coo_test_double;
