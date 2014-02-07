#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Meta-Vector 'dot' and 'norm2' test class
 *
 * \test The 'dot' and 'norm2' operations of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Peter Zajac
 */
template<typename Algo_, typename DataType_>
class MetaVectorDotNorm2Test
  : public MetaVectorTestBase<Algo_, DataType_>
{
public:
  typedef Algo_ AlgoType;
  typedef typename AlgoType::mem_type MemType;
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<Algo_, DataType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

  MetaVectorDotNorm2Test() : BaseClass("MetaVectorDotNorm2Test") {}

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

    // compute x*x and x*y
    DataType x_dot_y(DataType(0));
    DataType x_dot_x(DataType(0));

    // compute reference results
    for(Index i(0); i < n00; ++i)
    {
      x_dot_y += fx00(i) * fy00(i);
      x_dot_x += Math::sqr(fx00(i));
    }
    for(Index i(0); i < n01; ++i)
    {
      x_dot_y += fx01(i) * fy01(i);
      x_dot_x += Math::sqr(fx01(i));
    }
    for(Index i(0); i < n1; ++i)
    {
      x_dot_y += fx1(i) * fy1(i);
      x_dot_x += Math::sqr(fx1(i));
    }
    DataType x_norm2(Math::sqrt(x_dot_x));

    // test x*y
    TEST_CHECK_EQUAL_WITHIN_EPS(x.template dot<AlgoType>(y), x_dot_y, tol);

    // test x*x
    TEST_CHECK_EQUAL_WITHIN_EPS(x.template dot<AlgoType>(x), x_dot_x, tol);

    // test norm2(x)
    TEST_CHECK_EQUAL_WITHIN_EPS(x.template norm2<AlgoType>(), x_norm2, tol);
  }
};

MetaVectorDotNorm2Test<Algo::Generic, float> meta_vector_dot_norm2_test_generic_float;
MetaVectorDotNorm2Test<Algo::Generic, double> meta_vector_dot_norm2_test_generic_double;
#ifdef FEAST_BACKENDS_MKL
MetaVectorDotNorm2Test<Algo::MKL, float> meta_vector_dot_norm2_test_mkl_float;
MetaVectorDotNorm2Test<Algo::MKL, double> meta_vector_dot_norm2_test_mkl_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
MetaVectorDotNorm2Test<Algo::CUDA, float> meta_vector_dot_norm2_test_cuda_float;
MetaVectorDotNorm2Test<Algo::CUDA, double> meta_vector_dot_norm2_test_cuda_double;
#endif
