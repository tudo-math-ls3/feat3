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
    SparseMatrixCOO<Archs::CPU, DT_> a(10, 10);
    a(1,2,7);
    a.clear();
    a(1,2,7);
    a(5,5,2);
    SparseMatrixCSR<Tag_, DT_> b(a);
    TEST_CHECK_EQUAL(b.used_elements(), 2ul);
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(b(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(b(5, 5), a(5, 5));

    SparseMatrixCSR<Tag_, DT_> z(b);
    TEST_CHECK_EQUAL(z.used_elements(), 2ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(z(5, 5), a(5, 5));

    SparseMatrixCSR<Tag_, DT_> c;
    c = b;
    TEST_CHECK_EQUAL(c.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(c(0,2), b(0,2));
    TEST_CHECK_EQUAL(c(1,2), b(1,2));
    TEST_CHECK_EQUAL(c, b);

    DenseVector<Tag_, Index> col_ind(c.used_elements(), c.col_ind());
    DenseVector<Tag_, DT_> val(c.used_elements(), c.val());
    DenseVector<Tag_, Index> row_ptr(c.rows() + 1, c.row_ptr());
    DenseVector<Tag_, Index> row_ptr_end(c.rows(), c.row_ptr_end());
    SparseMatrixCSR<Tag_, DT_> d(c.rows(), c.columns(), col_ind, val, row_ptr, row_ptr_end);
    TEST_CHECK_EQUAL(d, c);

    SparseMatrixCSR<Archs::CPU, DT_> e(c);
    TEST_CHECK_EQUAL(e, c);
    e = c;
    TEST_CHECK_EQUAL(e, c);
  }
};
SparseMatrixCSRTest<Archs::CPU, float> cpu_sparse_matrix_csr_test_float;
SparseMatrixCSRTest<Archs::CPU, double> cpu_sparse_matrix_csr_test_double;
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCSRTest<Archs::GPU, float> gpu_sparse_matrix_csr_test_float;
SparseMatrixCSRTest<Archs::GPU, double> gpu_sparse_matrix_csr_test_double;
#endif
