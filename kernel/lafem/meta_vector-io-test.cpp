// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
template<
  typename DataType_,
  typename IndexType_>
class MetaVectorIOTest
  : public MetaVectorTestBase<DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

   MetaVectorIOTest(PreferredBackend backend) :
    BaseClass("MetaVectorIOTest", Type::Traits<DataType>::name(), Type::Traits<IndexType_>::name(), backend)
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

MetaVectorIOTest<float, unsigned int> meta_vector_io_test_generic_float_uint(PreferredBackend::generic);
MetaVectorIOTest<double, unsigned int> meta_vector_io_test_generic_double_uint(PreferredBackend::generic);
MetaVectorIOTest<float, unsigned long> meta_vector_io_test_generic_float_ulong(PreferredBackend::generic);
MetaVectorIOTest<double, unsigned long> meta_vector_io_test_generic_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MetaVectorIOTest<float, unsigned long> mkl_meta_vector_io_test_float_ulong(PreferredBackend::mkl);
MetaVectorIOTest<double, unsigned long> mkl_meta_vector_io_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MetaVectorIOTest<__float128, unsigned int> meta_vector_io_test_generic_float128_uint(PreferredBackend::generic);
MetaVectorIOTest<__float128, unsigned long> meta_vector_io_test_generic_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MetaVectorIOTest<Half, unsigned int> meta_vector_io_test_half_uint(PreferredBackend::generic);
MetaVectorIOTest<Half, unsigned long> meta_vector_io_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MetaVectorIOTest<float, unsigned int> meta_vector_io_test_cuda_float_uint(PreferredBackend::cuda);
MetaVectorIOTest<double, unsigned int> meta_vector_io_test_cuda_double_uint(PreferredBackend::cuda);
MetaVectorIOTest<float, unsigned long> meta_vector_io_test_cuda_float_ulong(PreferredBackend::cuda);
MetaVectorIOTest<double, unsigned long> meta_vector_io_test_cuda_double_ulong(PreferredBackend::cuda);
#endif
