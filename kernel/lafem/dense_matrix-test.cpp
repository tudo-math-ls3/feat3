#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_matrix.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for the dense matrix class.
 *
 * \test test description missing
 *
 * \tparam Mem_
 * description missing
 *
 * \tparam DT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename Mem_,
  typename DT_,
  typename IT_>
class DenseMatrixTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{

public:

  DenseMatrixTest()
    : FullTaggedTest<Mem_, DT_, IT_>("dense_matrix_test")
  {
  }

  virtual void run() const override
  {
    DenseMatrix<Mem_, DT_, IT_> zero1;
    DenseMatrix<Mem::Main, DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseMatrix<Mem_, DT_, IT_> a(10, 10);
    DenseMatrix<Mem_, DT_, IT_> b(10, 10, 5.);
    b(7, 6, DT_(42));
    DenseMatrix<Mem_, DT_, IT_> c;
    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c.rows(), b.rows());
    TEST_CHECK_EQUAL(c(7,6), b(7,6));
    TEST_CHECK_EQUAL(c, b);

    DenseMatrix<Mem_, DT_, IT_> e(11, 12, 5.);
    TEST_CHECK_EQUAL(e.rows(), 11ul);
    TEST_CHECK_EQUAL(e.columns(), 12ul);

    DenseMatrix<Mem_, DT_, IT_> h;
    h.clone(c);
    TEST_CHECK_EQUAL(h, c);
    h(1,2,3);
    TEST_CHECK_NOT_EQUAL(h, c);
    TEST_CHECK_NOT_EQUAL((void *)h.elements(), (void *)c.elements());

    h = c.shared();
    TEST_CHECK_EQUAL(h, c);
    h(1,2,3);
    TEST_CHECK_EQUAL(h, c);
    TEST_CHECK_EQUAL((void *)h.elements(), (void *)c.elements());

    auto kp = b.serialise();
    DenseMatrix<Mem_, DT_, IT_> k(kp);
    TEST_CHECK_EQUAL(k, b);
  }
};
DenseMatrixTest<Mem::Main, float, Index> cpu_dense_matrix_test_float;
DenseMatrixTest<Mem::Main, double, Index> cpu_dense_matrix_test_double;
#ifdef FEAST_BACKENDS_CUDA
DenseMatrixTest<Mem::CUDA, float, Index> cuda_dense_matrix_test_float;
DenseMatrixTest<Mem::CUDA, double, Index> cuda_dense_matrix_test_double;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixDenseApplyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  SparseMatrixDenseApplyTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixDenseApplyTest")
  {
  }

  virtual void run() const override
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseMatrix<Mem::Main, DT_, IT_> a(size, size, DT_(0));
      DenseVector<Mem::Main, DT_, IT_> x_local(size);
      DenseVector<Mem::Main, DT_, IT_> y_local(size);
      DenseVector<Mem::Main, DT_, IT_> ref_local(size);
      DenseVector<Mem_, DT_, IT_> ref(size);
      DenseVector<Mem::Main, DT_, IT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        x_local(i, DT_(i % 100 * DT_(1.234)));
        y_local(i, DT_(2 - DT_(i % 42)));
      }
      DenseVector<Mem_, DT_, IT_> x(size);
      x.copy(x_local);
      DenseVector<Mem_, DT_, IT_> y(size);
      y.copy(y_local);

      Index ue(0);
      for (Index row(0) ; row < a.rows() ; ++row)
      {
        for (Index col(0) ; col < a.columns() ; ++col)
        {
          if(row == col)
          {
            a(row, col, DT_(2));
            ++ue;
          }
          else if((row == col+1) || (row+1 == col))
          {
            a(row, col, DT_(-1));
            ++ue;
          }
        }
      }

      DenseVector<Mem_, DT_, IT_> r(size);

      // apply-test for alpha = 0.0
      a.apply(r, x, y, DT_(0.0));
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(y(i), r(i), 1e-2);

      // apply-test for alpha = -1.0
      a.apply(r, x, y, DT_(-1.0));
      a.apply(ref, x);
      ref.scale(ref, DT_(-1.0));
      ref.axpy(ref, y);

      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(ref(i), r(i), 1e-2);

      // apply-test for alpha = 4711.1
      //r.axpy(s, a, x, y);
      a.apply(r, x, y, s);
      result_local.copy(r);

      //ref.product_matvec(a, x);
      a.apply(ref, x);
      ref.scale(ref, s);
      ref.axpy(ref, y);
      ref_local.copy(ref);

      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), 1e-2);

      a.apply(r, x);
      result_local.copy(r);
      a.apply(ref, x);
      ref_local.copy(ref);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), 1e-2);
    }
  }
};

SparseMatrixDenseApplyTest<Mem::Main, float, unsigned long> sm_dense_apply_test_float_ulong;
SparseMatrixDenseApplyTest<Mem::Main, double, unsigned long> sm_dense_apply_test_double_ulong;
SparseMatrixDenseApplyTest<Mem::Main, float, unsigned int> sm_dense_apply_test_float_uint;
SparseMatrixDenseApplyTest<Mem::Main, double, unsigned int> sm_dense_apply_test_double_uint;
