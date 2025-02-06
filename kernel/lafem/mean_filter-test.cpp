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
  MeanFilterVectorTest(PreferredBackend backend)
    : UnitTest("MeanFilterVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~MeanFilterVectorTest()
  {
  }

  virtual void run() const override
  {
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.75));
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

MeanFilterVectorTest <float, std::uint32_t> mean_filter_vector_test_generic_float_uint32(PreferredBackend::generic);
MeanFilterVectorTest <double, std::uint32_t> mean_filter_vector_test_generic_double_uint32(PreferredBackend::generic);
MeanFilterVectorTest <float, std::uint64_t> mean_filter_vector_test_generic_float_uint64(PreferredBackend::generic);
MeanFilterVectorTest <double, std::uint64_t> mean_filter_vector_test_generic_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MeanFilterVectorTest <float, std::uint64_t> mkl_mean_filter_vector_test_float_uint64(PreferredBackend::mkl);
MeanFilterVectorTest <double, std::uint64_t> mkl_mean_filter_vector_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MeanFilterVectorTest <__float128, std::uint64_t> mean_filter_vector_test_float128_uint64(PreferredBackend::generic);
MeanFilterVectorTest <__float128, std::uint32_t> mean_filter_vector_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MeanFilterVectorTest <Half, std::uint32_t> mean_filter_vector_test_half_uint32(PreferredBackend::generic);
MeanFilterVectorTest <Half, std::uint64_t> mean_filter_vector_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MeanFilterVectorTest <float, std::uint32_t> mean_filter_vector_test_cuda_float_uint32(PreferredBackend::cuda);
MeanFilterVectorTest <double, std::uint32_t> mean_filter_vector_test_cuda_double_uint32(PreferredBackend::cuda);
MeanFilterVectorTest <float, std::uint64_t> mean_filter_vector_test_cuda_float_uint64(PreferredBackend::cuda);
MeanFilterVectorTest <double, std::uint64_t> mean_filter_vector_test_cuda_double_uint64(PreferredBackend::cuda);
#endif