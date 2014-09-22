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
  typename Algo_,
  typename DT_,
  typename IT_>
class DenseMatrixTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{

public:

  DenseMatrixTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("dense_matrix_test")
  {
  }

  virtual void run() const
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
DenseMatrixTest<Mem::Main, NotSet, float, Index> cpu_dense_matrix_test_float;
DenseMatrixTest<Mem::Main, NotSet, double, Index> cpu_dense_matrix_test_double;
#ifdef FEAST_BACKENDS_CUDA
DenseMatrixTest<Mem::CUDA, NotSet, float, Index> cuda_dense_matrix_test_float;
DenseMatrixTest<Mem::CUDA, NotSet, double, Index> cuda_dense_matrix_test_double;
#endif
