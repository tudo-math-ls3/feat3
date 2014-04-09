#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Meta-Vector component-product test class
 *
 * \test The 'component_invert' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Peter Zajac
 */
template<typename Algo_, typename DataType_, typename IndexType_>
class MetaVectorCompInvertTest
  : public MetaVectorTestBase<Algo_, DataType_, IndexType_>
{
public:
  typedef Algo_ AlgoType;
  typedef typename AlgoType::MemType MemType;
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<Algo_, DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

  MetaVectorCompInvertTest() : BaseClass("MetaVectorCompInvertTest") {}

  using BaseClass::fy00;
  using BaseClass::fy01;
  using BaseClass::fy1;

  virtual void run() const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.7));

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    // Note: All entries of y are positive, so we can safely use it for component inversion
    MetaVector y(this->gen_vector_y(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    // set z <- 1 / y
    z.template component_invert<AlgoType>(y);

    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(0)>()(i), DataType(1) / fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(0)>().template at<Index(1)>()(i), DataType(1) / fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<Index(1)>()(i), DataType(1) / fy1(i), tol);
  }
};

MetaVectorCompInvertTest<Algo::Generic, float, Index> meta_vector_comp_invert_test_generic_float;
MetaVectorCompInvertTest<Algo::Generic, double, Index> meta_vector_comp_invert_test_generic_double;
