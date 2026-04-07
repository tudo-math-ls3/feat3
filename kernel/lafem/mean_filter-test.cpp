// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/mean_filter_blocked.hpp>


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
* \brief Test class for MeanFilter class template
*
* \author Pia Ritter
*/
template<
  typename DT_,
  typename IT_>
class MeanFilterVectorTest
  : public UnitTest
{
public:
  explicit MeanFilterVectorTest(PreferredBackend backend)
    : UnitTest("MeanFilterVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    const DT_ eps = TestSystem::tol<DT_>();
    for(Index size(1); size < Index(1e3); size*=2)
    {
      DT_ tol = eps * DT_(size);

      DenseVector<DT_, IT_> vec_prim(size, DT_(1));
      DenseVector<DT_, IT_> vec_dual(size, DT_(DT_(1)/DT_(size)));

      DT_ sol_mean = DT_(1.563);

      DenseVector<DT_, IT_> vec_test_prim(rng, size, DT_(0), DT_(1));
      DenseVector<DT_, IT_> vec_test_dual(rng, size, DT_(0), DT_(1));

      MeanFilter<DT_, IT_> filter(vec_prim.clone(), vec_dual.clone(), sol_mean);

      filter.filter_def(vec_test_dual);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_dual.dot(vec_prim), DT_(0), tol);

      filter.filter_cor(vec_test_prim);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_prim.dot(vec_dual), DT_(0), tol);

      filter.filter_rhs(vec_test_dual);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_dual.dot(vec_prim), DT_(0), tol);

      filter.filter_sol(vec_test_prim);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_prim.dot(vec_dual), sol_mean, tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MeanFilterVectorTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
