#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Meta-Vector axpy test class
 *
 * \test The 'axpy' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Peter Zajac
 */
template<typename MemType_, typename DataType_, typename IndexType_>
class MetaVectorAxpyTest
  : public MetaVectorTestBase<MemType_, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<MemType_, DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

  MetaVectorAxpyTest() : BaseClass("MetaVectorAxpyTest") {}

  using BaseClass::fx00;
  using BaseClass::fx01;
  using BaseClass::fx1;
  using BaseClass::fy00;
  using BaseClass::fy01;
  using BaseClass::fy1;

  virtual void run() const override
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.7));

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    MetaVector x(this->gen_vector_x(n00, n01, n1));
    MetaVector y(this->gen_vector_y(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    // test: z <- 0.7*x + y
    // purpose: general test
    z.axpy(x, y, DataType(0.7));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), DataType(0.7)*fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), DataType(0.7)*fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), DataType(0.7)*fx1(i) + fy1(i), tol);

    // test: z <- 0*x + y; z <- 0.7*x + z
    // purpose: case z=y; case alpha = 0
    z.axpy(x, y, DataType(0));
    z.axpy(x, z, DataType(0.7));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), DataType(0.7)*fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), DataType(0.7)*fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), DataType(0.7)*fx1(i) + fy1(i), tol);

    // test: z <- x; z <- 0.7*z + y
    // purpose: case z=x
    z.copy(x);
    z.axpy(z, y, DataType(0.7));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), DataType(0.7)*fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), DataType(0.7)*fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), DataType(0.7)*fx1(i) + fy1(i), tol);

    // test: z <- x + y
    // purpose: alpha = 1
    z.axpy(x, y, DataType(1));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), fx1(i) + fy1(i), tol);

    // test: z <- -x + y
    // purpose: alpha = -1
    z.axpy(x, y, -DataType(1));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), -fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), -fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), -fx1(i) + fy1(i), tol);
  }
};

MetaVectorAxpyTest<Mem::Main, float, Index> meta_vector_axpy_test_generic_float;
MetaVectorAxpyTest<Mem::Main, double, Index> meta_vector_axpy_test_generic_double;
#ifdef FEAT_BACKENDS_CUDA
MetaVectorAxpyTest<Mem::CUDA, float, Index> meta_vector_axpy_test_cuda_float;
MetaVectorAxpyTest<Mem::CUDA, double, Index> meta_vector_axpy_test_cuda_double;
#endif
