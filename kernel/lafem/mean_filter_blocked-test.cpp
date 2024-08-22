// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
  MeanFilterBlockedVectorTest(PreferredBackend backend)
    : UnitTest("MeanFilterBlockedVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~MeanFilterBlockedVectorTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.7));
    for(Index size(1); size < Index(1e3); size*=2)
    {
      DenseVectorBlocked<DT_, IT_, block_size> vec_prim(size, DT_(1));
      DenseVectorBlocked<DT_, IT_, block_size> vec_dual(size, DT_(DT_(1)/DT_(vec_prim.size())));

      Tiny::Vector<DT_, block_size> sol_mean(DT_(1.563));
      //DT_ sol_mean = DT_(1.563);

      Random rng;
      DenseVectorBlocked<DT_, IT_, block_size> vec_test_prim(rng, size, DT_(0), DT_(1));
      DenseVectorBlocked<DT_, IT_, block_size> vec_test_dual(rng, size, DT_(0), DT_(1));

      MeanFilterBlocked<DT_, IT_, block_size> filter(vec_prim.clone(), vec_dual.clone(), sol_mean);

      filter.filter_def(vec_test_dual);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_dual.dot(vec_prim), DT_(0), eps*DT_(10));

      filter.filter_cor(vec_test_prim);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_prim.dot(vec_dual), DT_(0), eps);

      filter.filter_rhs(vec_test_dual);
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_dual.dot(vec_prim), DT_(0), eps);

      filter.filter_sol(vec_test_prim);
      for(int i(0); i<block_size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(vec_test_prim.dot_blocked(vec_dual)(i), sol_mean(i), eps);
      }
    }
  }
};

MeanFilterBlockedVectorTest <float, std::uint32_t, 2> mean_filter_blocked_vector_test_generic_float_uint32(PreferredBackend::generic);
MeanFilterBlockedVectorTest <double, std::uint32_t, 2> mean_filter_blocked_vector_test_generic_double_uint32(PreferredBackend::generic);
MeanFilterBlockedVectorTest <float, std::uint64_t, 3> mean_filter_blocked_vector_test_generic_float_uint64(PreferredBackend::generic);
MeanFilterBlockedVectorTest <double, std::uint64_t, 3> mean_filter_blocked_vector_test_generic_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MeanFilterBlockedVectorTest <float, std::uint64_t, 2> mkl_mean_filter_blocked_vector_test_float_uint64(PreferredBackend::mkl);
MeanFilterBlockedVectorTest <double, std::uint64_t, 3> mkl_mean_filter_blocked_vector_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MeanFilterBlockedVectorTest <__float128, std::uint64_t, 2> mean_filter_blocked_vector_test_float128_uint64(PreferredBackend::generic);
MeanFilterBlockedVectorTest <__float128, std::uint32_t, 3> mean_filter_blocked_vector_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MeanFilterBlockedVectorTest <Half, std::uint32_t, 2> mean_filter_blockedvector_test_half_uint32(PreferredBackend::generic);
MeanFilterBlockedVectorTest <Half, std::uint64_t, 3> mean_filter_blockedvector_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MeanFilterBlockedVectorTest <float, std::uint32_t, 2> mean_filter_blocked_vector_test_cuda_float_uint32(PreferredBackend::cuda);
MeanFilterBlockedVectorTest <double, std::uint32_t, 2> mean_filter_blocked_vector_test_cuda_double_uint32(PreferredBackend::cuda);
MeanFilterBlockedVectorTest <float, std::uint64_t, 3> mean_filter_blocked_vector_test_cuda_float_uint64(PreferredBackend::cuda);
MeanFilterBlockedVectorTest <double, std::uint64_t, 3> mean_filter_blocked_vector_test_cuda_double_uint64(PreferredBackend::cuda);
#endif