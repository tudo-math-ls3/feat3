// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
template<
  typename DataType_,
  typename IndexType_>
class MetaVectorAxpyTest
  : public MetaVectorTestBase<DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

   MetaVectorAxpyTest(PreferredBackend backend) :
    BaseClass("MetaVectorAxpyTest", Type::Traits<DataType>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~MetaVectorAxpyTest()
  {
  }

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
    z.copy(y);
    z.axpy(x, DataType(0.7));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), DataType(0.7)*fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), DataType(0.7)*fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), DataType(0.7)*fx1(i) + fy1(i), tol);

    // test: z <- 0*x + y; z <- 0.7*x + z
    // purpose: case z=y; case alpha = 0
    z.copy(y);
    z.axpy(x, DataType(0));
    z.axpy(x, DataType(0.7));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), DataType(0.7)*fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), DataType(0.7)*fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), DataType(0.7)*fx1(i) + fy1(i), tol);

    // test: z <- x; z <- 0.7*z + y
    // purpose: case z=x
    z.copy(x);
    z.scale(z, DataType(0.7));
    z.axpy(y); /// \todo use axpby here
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), DataType(0.7)*fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), DataType(0.7)*fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), DataType(0.7)*fx1(i) + fy1(i), tol);

    // test: z <- x + y
    // purpose: alpha = 1
    z.copy(y);
    z.axpy(x, DataType(1));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), fx1(i) + fy1(i), tol);

    // test: z <- -x + y
    // purpose: alpha = -1
    z.copy(y);
    z.axpy(x, -DataType(1));
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<0>()(i), -fx00(i) + fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<0>().template at<1>()(i), -fx01(i) + fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(z.template at<1>()(i), -fx1(i) + fy1(i), tol);
  }
};

MetaVectorAxpyTest <float, std::uint32_t> meta_vector_axpy_test_generic_float_uint32(PreferredBackend::generic);
MetaVectorAxpyTest <double, std::uint32_t> meta_vector_axpy_test_generic_double_uint32(PreferredBackend::generic);
MetaVectorAxpyTest <float, std::uint64_t> meta_vector_axpy_test_generic_float_uint64(PreferredBackend::generic);
MetaVectorAxpyTest <double, std::uint64_t> meta_vector_axpy_test_generic_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MetaVectorAxpyTest <float, std::uint64_t> mkl_meta_vector_axpy_test_float_uint64(PreferredBackend::mkl);
MetaVectorAxpyTest <double, std::uint64_t> mkl_meta_vector_axpy_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MetaVectorAxpyTest <__float128, std::uint32_t> meta_vector_axpy_test_generic_float128_uint32(PreferredBackend::generic);
MetaVectorAxpyTest <__float128, std::uint64_t> meta_vector_axpy_test_generic_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MetaVectorAxpyTest <Half, std::uint32_t> meta_vector_axpy_test_half_uint32(PreferredBackend::generic);
MetaVectorAxpyTest <Half, std::uint64_t> meta_vector_axpy_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MetaVectorAxpyTest <float, std::uint32_t> meta_vector_axpy_test_cuda_float_uint32(PreferredBackend::cuda);
MetaVectorAxpyTest <double, std::uint32_t> meta_vector_axpy_test_cuda_double_uint32(PreferredBackend::cuda);
MetaVectorAxpyTest <float, std::uint64_t> meta_vector_axpy_test_cuda_float_uint64(PreferredBackend::cuda);
MetaVectorAxpyTest <double, std::uint64_t> meta_vector_axpy_test_cuda_double_uint64(PreferredBackend::cuda);
#endif
