#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Meta-Vector scale test class
 *
 * \test The 'scale' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Peter Zajac
 */
template<typename Algo_, typename DataType_>
class MetaVectorScaleTest
  : public MetaVectorTestBase<Algo_, DataType_>
{
public:
  typedef Algo_ AlgoType;
  typedef typename AlgoType::mem_type MemType;
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<Algo_, DataType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

  MetaVectorScaleTest() : BaseClass("MetaVectorScaleTest") {}

  using BaseClass::fx00;
  using BaseClass::fx01;
  using BaseClass::fx1;

  virtual void run() const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.7));

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    MetaVector x(this->gen_vector_x(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    // test: z <- 0.7*x
    // purpose: general test
    z.template scale<AlgoType>(x, DataType(0.7));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), DataType(0.7)*fx00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), DataType(0.7)*fx01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), DataType(0.7)*fx1(i), tol);

    // test: z <- x
    // purpose: alpha = 1
    z.template scale<AlgoType>(x, DataType(1));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), fx00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), fx01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), fx1(i), tol);

    // test: z <- -x
    // purpose: alpha = -1
    z.template scale<AlgoType>(x, -DataType(1));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), -fx00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), -fx01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), -fx1(i), tol);

    // test: z <- 0*x
    // purpose: alpha = 0
    z.template scale<AlgoType>(x, -DataType(0));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), DataType(0), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), DataType(0), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), DataType(0), tol);
  }
};

MetaVectorScaleTest<Algo::Generic, float> meta_vector_scale_test_generic_float;
MetaVectorScaleTest<Algo::Generic, double> meta_vector_scale_test_generic_double;
#ifdef FEAST_BACKENDS_MKL
MetaVectorScaleTest<Algo::MKL, float> meta_vector_scale_test_mkl_float;
MetaVectorScaleTest<Algo::MKL, double> meta_vector_scale_test_mkl_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
MetaVectorScaleTest<Algo::CUDA, float> meta_vector_scale_test_cuda_float;
MetaVectorScaleTest<Algo::CUDA, double> meta_vector_scale_test_cuda_double;
#endif
