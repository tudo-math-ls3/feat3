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
 * \brief Meta-Vector scale test class
 *
 * \test The 'scale' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Peter Zajac
 */
template<
  typename DataType_,
  typename IndexType_>
class MetaVectorScaleTest
  : public MetaVectorTestBase<DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

  explicit MetaVectorScaleTest(PreferredBackend backend) :
    BaseClass("MetaVectorScaleTest", Type::Traits<DataType>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  using BaseClass::fx00;
  using BaseClass::fx01;
  using BaseClass::fx1;

  virtual void run() const override
  {
    const DataType tol = TestSystem::tol<DataType>();

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    MetaVector x(this->gen_vector_x(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    // test: z <- 0.7*x
    // purpose: general test
    {
      z.scale(x, DataType(0.7));
      const Memory::TypedView<DataType> vz00(z.template at<0>().template at<0>().elements_view_r());
      const Memory::TypedView<DataType> vz01(z.template at<0>().template at<1>().elements_view_r());
      const Memory::TypedView<DataType> vz1(z.template at<1>().elements_view_r());
      for(Index i(0); i < n00; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz00(i), DataType(0.7)*fx00(i), tol);
      for(Index i(0); i < n01; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz01(i), DataType(0.7)*fx01(i), tol);
      for(Index i(0); i < n1; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz1(i), DataType(0.7)*fx1(i), tol);
    }

    // test: z <- x
    // purpose: alpha = 1
    z.scale(x, DataType(1));
    {
      const Memory::TypedView<DataType> vz00(z.template at<0>().template at<0>().elements_view_r());
      const Memory::TypedView<DataType> vz01(z.template at<0>().template at<1>().elements_view_r());
      const Memory::TypedView<DataType> vz1(z.template at<1>().elements_view_r());
      for(Index i(0); i < n00; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz00(i), fx00(i), tol);
      for(Index i(0); i < n01; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz01(i), fx01(i), tol);
      for(Index i(0); i < n1; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz1(i), fx1(i), tol);
    }

    // test: z <- -x
    // purpose: alpha = -1
    z.scale(x, -DataType(1));
    {
      const Memory::TypedView<DataType> vz00(z.template at<0>().template at<0>().elements_view_r());
      const Memory::TypedView<DataType> vz01(z.template at<0>().template at<1>().elements_view_r());
      const Memory::TypedView<DataType> vz1(z.template at<1>().elements_view_r());
      for(Index i(0); i < n00; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz00(i), -fx00(i), tol);
      for(Index i(0); i < n01; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz01(i), -fx01(i), tol);
      for(Index i(0); i < n1; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz1(i), -fx1(i), tol);
    }

    // test: z <- 0*x
    // purpose: alpha = 0
    z.scale(x, -DataType(0));
    {
      const Memory::TypedView<DataType> vz00(z.template at<0>().template at<0>().elements_view_r());
      const Memory::TypedView<DataType> vz01(z.template at<0>().template at<1>().elements_view_r());
      const Memory::TypedView<DataType> vz1(z.template at<1>().elements_view_r());
      for(Index i(0); i < n00; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz00(i), DataType(0), tol);
      for(Index i(0); i < n01; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz01(i), DataType(0), tol);
      for(Index i(0); i < n1; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(vz1(i), DataType(0), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorScaleTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
