// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Meta-Vector component-product test class
 *
 * \test The 'component_invert' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Peter Zajac
 */
template<
  typename DataType_,
  typename IndexType_>
class MetaVectorCompInvertTest
  : public MetaVectorTestBase<DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

   MetaVectorCompInvertTest(PreferredBackend backend) :
    BaseClass("MetaVectorCompInvertTest", Type::Traits<DataType>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~MetaVectorCompInvertTest()
  {
  }

  using BaseClass::fy00;
  using BaseClass::fy01;
  using BaseClass::fy1;

  virtual void run() const override
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.7));

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    // Note: All entries of y are positive, so we can safely use it for component inversion
    MetaVector y(this->gen_vector_y(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    // set z <- 1 / y
    z.component_invert(y);

    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), DataType(1) / fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), DataType(1) / fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), DataType(1) / fy1(i), tol);
  }
};

MetaVectorCompInvertTest <float, std::uint32_t> meta_vector_comp_invert_test_float_uint32(PreferredBackend::generic);
MetaVectorCompInvertTest <double, std::uint32_t> meta_vector_comp_invert_test_double_uint32(PreferredBackend::generic);
MetaVectorCompInvertTest <float, std::uint64_t> meta_vector_comp_invert_test_float_uint64(PreferredBackend::generic);
MetaVectorCompInvertTest <double, std::uint64_t> meta_vector_comp_invert_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MetaVectorCompInvertTest <float, std::uint64_t> mkl_meta_vector_comp_invert_test_float_uint64(PreferredBackend::mkl);
MetaVectorCompInvertTest <double, std::uint64_t> mkl_meta_vector_comp_invert_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MetaVectorCompInvertTest <__float128, std::uint32_t> meta_vector_comp_invert_test_generic_float128_uint32(PreferredBackend::generic);
MetaVectorCompInvertTest <__float128, std::uint64_t> meta_vector_comp_invert_test_generic_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MetaVectorCompInvertTest <Half, std::uint32_t> meta_vector_comp_invert_test_half_uint32(PreferredBackend::generic);
MetaVectorCompInvertTest <Half, std::uint64_t> meta_vector_comp_invert_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MetaVectorCompInvertTest <float, std::uint32_t> meta_vector_comp_invert_test_cuda_float_uint32(PreferredBackend::cuda);
MetaVectorCompInvertTest <double, std::uint32_t> meta_vector_comp_invert_test_cuda_double_uint32(PreferredBackend::cuda);
MetaVectorCompInvertTest <float, std::uint64_t> meta_vector_comp_invert_test_cuda_float_uint64(PreferredBackend::cuda);
MetaVectorCompInvertTest <double, std::uint64_t> meta_vector_comp_invert_test_cuda_double_uint64(PreferredBackend::cuda);
#endif
