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
* \brief Test class for MeanFilterBlocked class template
*
* \author Pia Ritter
*/

template<
  typename DT_,
  typename IT_,
  int block_size>
class MeanFilterBlockedVectorTest
  : public UnitTest
{
public:
  explicit MeanFilterBlockedVectorTest(PreferredBackend backend)
    : UnitTest("MeanFilterBlockedVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    DT_ eps = TestSystem::tol<DT_>();
    for(Index size(1); size < Index(1e3); size*=2)
    {
      DT_ tol = eps * DT_(size);
      DenseVectorBlocked<DT_, IT_, block_size> vec_prim(size, DT_(1));
      DenseVectorBlocked<DT_, IT_, block_size> vec_dual(size, DT_(DT_(1)/DT_(vec_prim.size())));

      Tiny::Vector<DT_, block_size> sol_mean(DT_(1.563));

      DenseVectorBlocked<DT_, IT_, block_size> vec_test_prim(rng, size, DT_(0), DT_(1));
      DenseVectorBlocked<DT_, IT_, block_size> vec_test_dual(rng, size, DT_(0), DT_(1));

      MeanFilterBlocked<DT_, IT_, block_size> filter(vec_prim.clone(), vec_dual.clone(), sol_mean);

      filter.filter_def(vec_test_dual);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_dual.dot(vec_prim), DT_(0), tol);

      filter.filter_cor(vec_test_prim);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_prim.dot(vec_dual), DT_(0), tol);

      filter.filter_rhs(vec_test_dual);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_dual.dot(vec_prim), DT_(0), tol);

      filter.filter_sol(vec_test_prim);

      auto r_blocked = vec_test_prim.dot_blocked(vec_dual);
      for(int i(0); i< block_size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(r_blocked(i), sol_mean(i), tol);
      }
    }
  }
};

SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(MeanFilterBlockedVectorTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif
