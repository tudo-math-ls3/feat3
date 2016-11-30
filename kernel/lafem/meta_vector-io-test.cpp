#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>
#include <kernel/util/binary_stream.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Meta-Vector scale test class
 *
 * \test The 'write/read' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Dirk Ribbrock
 */
template<typename MemType_, typename DataType_, typename IndexType_>
class MetaVectorIOTest
  : public MetaVectorTestBase<MemType_, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<MemType_, DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

   MetaVectorIOTest() :
    BaseClass("MetaVectorIOTest")
  {
  }

  virtual ~MetaVectorIOTest()
  {
  }

  using BaseClass::fx00;
  using BaseClass::fx01;
  using BaseClass::fx1;

  virtual void run() const override
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.7));

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    MetaVector x(this->gen_vector_x(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    BinaryStream bs;
    x.write_out(FileMode::fm_binary, bs);
    bs.seekg(0);
    z.read_from(FileMode::fm_binary, bs);

    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), x.template at<0>().template at<0>()(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), x.template at<0>().template at<1>()(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), x.template at<1>()(i), tol);
  }
};

MetaVectorIOTest<Mem::Main, float, Index> meta_vector_io_test_generic_float;
MetaVectorIOTest<Mem::Main, double, Index> meta_vector_io_test_generic_double;
#ifdef FEAT_HAVE_CUDA
MetaVectorIOTest<Mem::CUDA, float, Index> meta_vector_io_test_cuda_float;
MetaVectorIOTest<Mem::CUDA, double, Index> meta_vector_io_test_cuda_double;
#endif
