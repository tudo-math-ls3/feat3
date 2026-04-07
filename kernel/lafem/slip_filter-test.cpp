// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/slip_filter.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for SlipFilter class template
 *
 * Primitive tests based on add()-ing values to the filter
 *
 * \author Jordi Paul
 *
 */
template<
  typename DT_,
  typename IT_,
  int block_size_>
class SlipFilterVectorTest
  : public UnitTest
{
  typedef Tiny::Vector<DT_, block_size_> ValueType;
  typedef DenseVectorBlocked<DT_, IT_, block_size_> VectorType;
  typedef DenseVectorBlocked<IT_, IT_, block_size_> IVectorType;
  typedef SlipFilter<DT_, IT_, block_size_> FilterType;

public:
  explicit SlipFilterVectorTest(PreferredBackend backend)
    : UnitTest("SlipFilterVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    const Index nn(100), nnz(8);

    IT_ jj[8];

    for(Index i(0); i < nnz; ++i)
      jj[i] = IT_(i*(3 + i));

    FilterType my_filter(nn, nn, nnz);
    {
      Memory::TypedView<IT_> idx = my_filter.get_filter_vector().indices_view_w();
      Memory::TypedView<ValueType> val = my_filter.get_filter_vector().elements_view_w();

      for(Index i(0); i < nnz; ++i)
      {
        idx[i] = jj[i];
        val[i] = Math::sqrt(DT_(jj[i]+1));
        val[i][0] *= DT_(Math::pow(-DT_(0.5), DT_(i)));
        val[i][block_size_-1] = -DT_(i);
      }
    }

    VectorType my_vector(nn, DT_(-2));
    {
      const Memory::TypedView<ValueType> vfil(my_filter.get_filter_vector().elements_view_r());
      Memory::TypedView<ValueType> vvec(my_vector.elements_view_rw());
      // 0: Set element to 0
      vvec[0] = DT_(0);
      // 1: Set element to value of the filter
      vvec[jj[1]] = vfil[1];
      // 2: Set element to negative value of the filter
      vvec[jj[2]] = DT_(-1) * vfil[2];
      // 3: Set element to first unit vector
      vvec[jj[3]] = DT_(0);
      vvec[jj[3]][0] = DT_(1);
      // 4: Set element to -(last unit vector)
      vvec[jj[4]] = DT_(0);
      vvec[jj[4]][block_size_-1] = DT_(-1);
      // 5: Set element to twice the negative value of the filter
      vvec[jj[5]] = DT_(-2) * vfil[5];
    }

    // Filter vector
    my_filter.filter_def(my_vector);

    const Memory::TypedView<ValueType> nu = my_filter.get_filter_vector().elements_view_r();
    const Memory::TypedView<ValueType> v = my_vector.elements_view_r();

    // Check results
    for(IT_ i(0); i < nn; ++i)
    {
      // Check if we have a filtered entry
      int k = 0;
      for(; k < 8; ++k)
      {
        if(i == jj[k])
        {
          break;
        }
      }

      if(k < 8)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(Tiny::dot(nu[Index(k)], v[i]), DT_(0), tol);
      }
      else
      {
        for(int d(0); d < block_size_; ++d)
          TEST_CHECK_EQUAL_WITHIN_EPS(v(i)[d], -DT_(2), tol);
      }
    }
  }
};

SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, Half, std::uint32_t, 2, PreferredBackend::generic);
//SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, float, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, float, std::uint64_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, double, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(SlipFilterVectorTest, double, std::uint64_t, 2, PreferredBackend::cuda);
#endif
