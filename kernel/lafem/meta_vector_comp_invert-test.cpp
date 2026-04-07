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

  explicit MetaVectorCompInvertTest(PreferredBackend backend) :
    BaseClass("MetaVectorCompInvertTest", Type::Traits<DataType>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  using BaseClass::fy00;
  using BaseClass::fy01;
  using BaseClass::fy1;

  virtual void run() const override
  {
    const DataType tol = TestSystem::tol<DataType>();

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    // Note: All entries of y are positive, so we can safely use it for component inversion
    MetaVector y(this->gen_vector_y(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    // set z <- 1 / y
    z.component_invert(y);

    const Memory::TypedView<DataType> vz00(z.template at<0>().template at<0>().elements_view_r());
    const Memory::TypedView<DataType> vz01(z.template at<0>().template at<1>().elements_view_r());
    const Memory::TypedView<DataType> vz1(z.template at<1>().elements_view_r());
    for(Index i(0); i < n00; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(vz00(i), DataType(1) / fy00(i), tol);
    for(Index i(0); i < n01; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(vz01(i), DataType(1) / fy01(i), tol);
    for(Index i(0); i < n1; ++i)
      TEST_CHECK_EQUAL_WITHIN_EPS(vz1(i), DataType(1) / fy1(i), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorCompInvertTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
