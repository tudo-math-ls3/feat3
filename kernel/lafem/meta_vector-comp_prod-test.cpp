#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Meta-Vector component-product test class
 *
 * \test The 'component_product' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Peter Zajac
 */
template<typename Algo_, typename DataType_>
class MetaVectorCompProdTest
  : public MetaVectorTestBase<Algo_, DataType_>
{
public:
  typedef Algo_ AlgoType;
  typedef typename AlgoType::MemType MemType;
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<Algo_, DataType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

  MetaVectorCompProdTest() : BaseClass("MetaVectorCompProdTest") {}

  using BaseClass::fx00;
  using BaseClass::fx01;
  using BaseClass::fx1;
  using BaseClass::fy00;
  using BaseClass::fy01;
  using BaseClass::fy1;

  virtual void run() const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.7));

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    MetaVector x(this->gen_vector_x(n00, n01, n1));
    MetaVector y(this->gen_vector_y(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    // test: z <- x * y
    // purpose: general test
    z.template component_product<AlgoType>(x, y);
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), fx00(i) * fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), fx01(i) * fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), fx1(i) * fy1(i), tol);

    // test: z <- x; z <- z * y
    // purpose: z = x
    z.copy(x);
    z.template component_product<AlgoType>(z, y);
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), fx00(i) * fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), fx01(i) * fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), fx1(i) * fy1(i), tol);

    // test: z <- y; z <- x * z
    // purpose: z = y
    z.copy(y);
    z.template component_product<AlgoType>(x, z);
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), fx00(i) * fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), fx01(i) * fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), fx1(i) * fy1(i), tol);

    // test: z <- 1; z <- x * y + z
    // purpose: general test
    z.format(DataType(1));
    z.template component_product<AlgoType>(x, y, z);
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), fx00(i) * fy00(i) + DataType(1), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), fx01(i) * fy01(i) + DataType(1), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), fx1(i) * fy1(i) + DataType(1), tol);
  }
};

MetaVectorCompProdTest<Algo::Generic, float> meta_vector_comp_prod_test_generic_float;
MetaVectorCompProdTest<Algo::Generic, double> meta_vector_comp_prod_test_generic_double;
#ifdef FEAST_BACKENDS_MKL
MetaVectorCompProdTest<Algo::MKL, float> meta_vector_comp_prod_test_mkl_float;
MetaVectorCompProdTest<Algo::MKL, double> meta_vector_comp_prod_test_mkl_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
MetaVectorCompProdTest<Algo::CUDA, float> meta_vector_comp_prod_test_cuda_float;
MetaVectorCompProdTest<Algo::CUDA, double> meta_vector_comp_prod_test_cuda_double;
#endif
