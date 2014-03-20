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
  typename DT_>
class DenseMatrixTest
  : public TaggedTest<Mem_, DT_>
{

public:

  DenseMatrixTest()
    : TaggedTest<Mem_, DT_>("dense_matrix_test")
  {
  }

  virtual void run() const
  {
    DenseMatrix<Mem_, DT_> zero1;
    DenseMatrix<Mem::Main, DT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseMatrix<Mem_, DT_> a(10, 10);
    DenseMatrix<Mem_, DT_> b(10, 10, 5.);
    b(7, 6, DT_(42));
    DenseMatrix<Mem_, DT_> c;
    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c.rows(), b.rows());
    TEST_CHECK_EQUAL(c(7,6), b(7,6));
    TEST_CHECK_EQUAL(c, b);

    DenseMatrix<Mem_, DT_> e(11, 12, 5.);
    TEST_CHECK_EQUAL(e.rows(), 11ul);
    TEST_CHECK_EQUAL(e.columns(), 12ul);

    DenseMatrix<Mem_, DT_> h;
    h.clone(c);
    TEST_CHECK_EQUAL(h, c);
    h(1,2,3);
    TEST_CHECK_NOT_EQUAL(h, c);
    TEST_CHECK_NOT_EQUAL((void *)h.elements(), (void *)c.elements());
  }
};
DenseMatrixTest<Mem::Main, float> cpu_dense_matrix_test_float;
DenseMatrixTest<Mem::Main, double> cpu_dense_matrix_test_double;
#ifdef FEAST_BACKENDS_CUDA
DenseMatrixTest<Mem::CUDA, float> cuda_dense_matrix_test_float;
DenseMatrixTest<Mem::CUDA, double> cuda_dense_matrix_test_double;
#endif
