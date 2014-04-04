#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_csr_blocked.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the sparse matrix csr blocked class.
*
* \test test description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRBlockedTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  SparseMatrixCSRBlockedTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("SparseMatrixCSRBlockedTest")
  {
  }

  virtual void run() const
  {
    SparseMatrixCSRBlocked<Mem_, DT_, IT_, 2, 3> a;
    SparseMatrixCSRBlocked<Mem_, DT_, IT_, 2, 3> b;
    TEST_CHECK_EQUAL(b, a);

    DenseVector<Mem_, DT_, IT_> dv1(12);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i+1));
    DenseVector<Mem_, IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<Mem_, IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixCSRBlocked<Mem_, DT_, IT_, 2, 3> c(2, 2, dv2, dv1, dv3);
    SparseMatrixCSRBlocked<Mem_, DT_, IT_, 2, 3> d;
    d.convert(c);
    TEST_CHECK_EQUAL(d.rows(), c.rows());
    TEST_CHECK_EQUAL(d.columns(), c.columns());
    TEST_CHECK_EQUAL(d.used_elements(), c.used_elements());
    TEST_CHECK_EQUAL(d, c);
    TEST_CHECK_EQUAL((void*)d.raw_val(), (void*)c.raw_val());
    TEST_CHECK_EQUAL((void*)d.row_ptr(), (void*)c.row_ptr());
    SparseMatrixCSRBlocked<Mem_, DT_, IT_, 2, 3> e;
    e.clone(c);
    TEST_CHECK_EQUAL(e, c);
    TEST_CHECK_NOT_EQUAL((void*)e.raw_val(), (void*)c.raw_val());
    TEST_CHECK_EQUAL((void*)e.row_ptr(), (void*)c.row_ptr());
    e = c.clone(true);
    TEST_CHECK_EQUAL(e, c);
    TEST_CHECK_NOT_EQUAL((void*)e.raw_val(), (void*)c.raw_val());
    TEST_CHECK_NOT_EQUAL((void*)e.row_ptr(), (void*)c.row_ptr());

    SparseMatrixCSRBlocked<Mem_, DT_, IT_, 2, 3> f(c.layout());
    TEST_CHECK_EQUAL(f.rows(), c.rows());
    TEST_CHECK_EQUAL(f.columns(), c.columns());
    TEST_CHECK_EQUAL(f.used_elements(), c.used_elements());
    TEST_CHECK_NOT_EQUAL((void*)f.raw_val(), (void*)c.raw_val());
    TEST_CHECK_EQUAL((void*)f.row_ptr(), (void*)c.row_ptr());
  }
};
SparseMatrixCSRBlockedTest<Mem::Main, NotSet, float, unsigned long> cpu_sparse_matrix_csr_blocked_test_float_ulong;
SparseMatrixCSRBlockedTest<Mem::Main, NotSet, double, unsigned long> cpu_sparse_matrix_csr_blocked_test_double_ulong;
SparseMatrixCSRBlockedTest<Mem::Main, NotSet, float, unsigned int> cpu_sparse_matrix_csr_blocked_test_float_uint;
SparseMatrixCSRBlockedTest<Mem::Main, NotSet, double, unsigned int> cpu_sparse_matrix_csr_blocked_test_double_uint;
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCSRBlockedTest<Mem::CUDA, NotSet, float, unsigned long> cuda_sparse_matrix_csr_blocked_test_float_ulong;
SparseMatrixCSRBlockedTest<Mem::CUDA, NotSet, double, unsigned long> cuda_sparse_matrix_csr_blocked_test_double_ulong;
SparseMatrixCSRBlockedTest<Mem::CUDA, NotSet, float, unsigned int> cuda_sparse_matrix_csr_blocked_test_float_uint;
SparseMatrixCSRBlockedTest<Mem::CUDA, NotSet, double, unsigned int> cuda_sparse_matrix_csr_blocked_test_double_uint;
#endif

/**
* \brief Test class for the sparse matrix csr blocked apply method.
*
* \test test description missing
*
* \author Dirk Ribbrock
*/
template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRBlockedApplyTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  SparseMatrixCSRBlockedApplyTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("SparseMatrixCSRBlockedApplyTest")
  {
  }

  virtual void run() const
  {
    DenseVector<Mem_, DT_, IT_> dv1(12);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i+1));
    DenseVector<Mem_, IT_, IT_> dv2(2);
    dv2(0, IT_(0));
    dv2(1, IT_(1));
    DenseVector<Mem_, IT_, IT_> dv3(3);
    dv3(0, IT_(0));
    dv3(1, IT_(1));
    dv3(2, IT_(2));
    SparseMatrixCSRBlocked<Mem_, DT_, IT_, 2, 3> c(2, 2, dv2, dv1, dv3);

    DenseVector<Mem_, DT_, IT_> x(c.raw_columns());
    DenseVector<Mem_, DT_, IT_> r(c.raw_rows());
    DenseVector<Mem_, DT_, IT_> ref(c.raw_rows());
    for (Index i(0) ; i < x.size() ; ++i)
    {
      x(i, DT_(i));
    }
    for (Index i(0) ; i < r.size() ; ++i)
    {
      r(i, DT_(4711));
      ref(i, DT_(4711));
    }
    DenseVectorBlocked<Mem_, DT_, IT_, 3> xb(x);
    DenseVectorBlocked<Mem_, DT_, IT_, 2> rb(r);

    SparseMatrixCSR<Mem_, DT_, IT_> csr;
    csr.convert(c);
    csr.template apply<Algo_>(ref, x);

    c.template apply<Algo_>(r, x);

    TEST_CHECK_EQUAL(r, ref);

    c.template apply<Algo_>(rb, x);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);

    c.template apply<Algo_>(r, xb);
    TEST_CHECK_EQUAL(r, ref);

    c.template apply<Algo_>(rb, xb);
    r.convert(rb);
    TEST_CHECK_EQUAL(r, ref);
  }
};
SparseMatrixCSRBlockedApplyTest<Mem::Main, Algo::Generic, float, unsigned long> cpu_sparse_matrix_csr_blocked_apply_test_float_ulong;
SparseMatrixCSRBlockedApplyTest<Mem::Main, Algo::Generic, double, unsigned long> cpu_sparse_matrix_csr_blocked_apply_test_double_ulong;
SparseMatrixCSRBlockedApplyTest<Mem::Main, Algo::Generic, float, unsigned int> cpu_sparse_matrix_csr_blocked_apply_test_float_uint;
SparseMatrixCSRBlockedApplyTest<Mem::Main, Algo::Generic, double, unsigned int> cpu_sparse_matrix_csr_blocked_apply_test_double_uint;
