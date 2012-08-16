#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
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
class SparseMatrixCSRTest
  : public TaggedTest<Tag_, DT_>
{

public:

  SparseMatrixCSRTest()
    : TaggedTest<Tag_, DT_>("sparse_matrix_csr_test")
  {
  }

  virtual void run() const
  {
    SparseMatrixCOO<Tag_, DT_> a(10, 10);
    a(1,2,7);
    a(5,5,2);
    SparseMatrixCSR<Tag_, DT_> b(a);
    TEST_CHECK_EQUAL(b.used_elements(), 2ul);
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(b(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(b(5, 5), a(5, 5));

    SparseMatrixCSR<Tag_, DT_> c;
    c = b;
    TEST_CHECK_EQUAL(c.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(c(0,2), b(0,2));
    TEST_CHECK_EQUAL(c(1,2), b(1,2));
    TEST_CHECK_EQUAL(c, b);

    DenseVector<Tag_, Index> Aj(c.used_elements(), c.Aj());
    DenseVector<Tag_, DT_> Ax(c.used_elements(), c.Ax());
    DenseVector<Tag_, Index> Ar(c.rows() + 1, c.Ar());
    SparseMatrixCSR<Tag_, DT_> d(c.rows(), c.columns(), Aj, Ax, Ar);
    TEST_CHECK_EQUAL(d, c);
  }
};
SparseMatrixCSRTest<Archs::CPU, float> sparse_matrix_csr_test_float;
SparseMatrixCSRTest<Archs::CPU, double> sparse_matrix_csr_test_double;
