#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>

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
class SparseMatrixELLTest
  : public TaggedTest<Tag_, DT_>
{

public:

  SparseMatrixELLTest()
    : TaggedTest<Tag_, DT_>("sparse_matrix_ell_test")
  {
  }

  virtual void run() const
  {
    SparseMatrixCOO<Mem::Main, DT_> a(10, 10);
    a(1,2,7);
    a.clear();
    a(1,2,7);
    a(5,5,2);
    SparseMatrixCSR<Tag_, DT_> a2(a);
    SparseMatrixELL<Tag_, DT_> b(a2);
    TEST_CHECK_EQUAL(b.used_elements(), 2ul);
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(b(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(b(5, 5), a(5, 5));

    SparseMatrixELL<Tag_, DT_> z(b);
    TEST_CHECK_EQUAL(z.used_elements(), 2ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z.stride(), b.stride());
    TEST_CHECK_EQUAL(z.num_cols_per_row(), b.num_cols_per_row());
    TEST_CHECK_EQUAL(z(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(z(5, 5), a(5, 5));
    TEST_CHECK_EQUAL(z(1, 3), a(1, 3));

    SparseMatrixELL<Tag_, DT_> c;
    c = b;
    TEST_CHECK_EQUAL(c.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(c(0,2), b(0,2));
    TEST_CHECK_EQUAL(c(1,2), b(1,2));
    TEST_CHECK_EQUAL(c, b);

    SparseMatrixELL<Mem::Main, DT_> e(c);
    TEST_CHECK_EQUAL(e, c);
    e = c;
    TEST_CHECK_EQUAL(e, c);

    SparseMatrixELL<Tag_, DT_> f("5pt_10x10.ell");
  }
};
SparseMatrixELLTest<Mem::Main, float> cpu_sparse_matrix_ell_test_float;
SparseMatrixELLTest<Mem::Main, double> cpu_sparse_matrix_ell_test_double;
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixELLTest<Mem::CUDA, float> gpu_sparse_matrix_ell_test_float;
SparseMatrixELLTest<Mem::CUDA, double> gpu_sparse_matrix_ell_test_double;
#endif
