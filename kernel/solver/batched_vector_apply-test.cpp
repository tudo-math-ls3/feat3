// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/solver/krylov_orthogonalizer.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/type_traits.hpp>

#include <cstdio>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the batched dot product kernel.
 *
 * \test Tests whether the standard dot product and the batched dot product arrive at the same result
 *
 * \author Maximilian Esser
 */
template<
  typename DT_,
  typename IT_>
class BatchedDotTest
  : public UnitTest
{
public:
  BatchedDotTest(PreferredBackend backend)
    : UnitTest("BatchedDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~BatchedDotTest()
  {
  }


  virtual void run() const override
  {
    constexpr std::size_t dot_size = 1047;
    constexpr std::size_t num_dots = 13;
    const DT_ tol = Math::pow(Math::Limits<DT_>::epsilon(), DT_(0.5)) * (std::is_same_v<float, DT_> ? DT_(2) : DT_(1));
    {
      LAFEM::DenseVector<DT_, IT_> a(dot_size);
      a.format();
      std::vector<LAFEM::DenseVector<DT_, IT_>> vecs;
      std::vector<Memory::Arbiter const*> arbiter;
      vecs.reserve(num_dots);
      arbiter.reserve(num_dots);
      for(std::size_t k = 0; k < num_dots; ++k)
      {
        vecs.emplace_back(dot_size);
        vecs.back().format();
        arbiter.emplace_back(&vecs.back().elements_arbiter());
      }
      {
        // fill via host write views (raw .elements() access was replaced by the arbiter API)
        auto a_view = a.elements_view_w();
        DT_* a_data = a_view.get_w();
        std::vector<Memory::TypedView<DT_>> vec_views;
        vec_views.reserve(num_dots);
        std::vector<DT_*> vec_data(num_dots);
        for(std::size_t k = 0; k < num_dots; ++k)
        {
          vec_views.emplace_back(vecs[k].elements_view_w());
          vec_data[k] = vec_views[k].get_w();
        }
        for(std::size_t i = 0; i < dot_size; ++i)
        {
          a_data[i] = DT_(i%7)/DT_(5.);
          for(std::size_t k = 0; k < num_dots; ++k)
            vec_data[k][i] = DT_((11*(i+1)) %17)/DT_(k+1);
        }
      }
      std::array<DT_, num_dots> res_a;
      for(std::size_t k = 0; k < num_dots; ++k)
      {
        res_a[k] = a.dot(vecs[k]);
      }
      // the batched dot results are passed as a Memory::Arbiter; zero it, run the dot, read it back
      Memory::Arbiter res_arb(num_dots * sizeof(DT_));
      res_arb.format(Memory::Location::main);
      // stride equals the number of vectors for a single DenseVector leaf (no tuple components)
      FEAT::Solver::Intern::BatchedHelper<LAFEM::DenseVector<DT_, IT_>>::template batched_dot<DT_>(res_arb, a, arbiter.data(), num_dots, num_dots);

      std::array<DT_, num_dots> res_b;
      {
        Memory::TypedView<DT_> rv(res_arb.view(Memory::Location::main, Memory::Access::read));
        const DT_* rp = rv.get_r();
        for(std::size_t k = 0; k < num_dots; ++k)
          res_b[k] = rp[k];
      }

      for(std::size_t k = 0; k < num_dots; ++k)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(res_a[k], res_b[k], tol);
      }
    }
  }
};


BatchedDotTest<double, std::uint32_t> batched_dot_test_double_generic(PreferredBackend::generic);
BatchedDotTest<float, std::uint32_t> batched_dot_test_float_generic(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
BatchedDotTest<double, std::uint32_t> batched_dot_test_double_cuda(PreferredBackend::cuda);
BatchedDotTest<float, std::uint32_t> batched_dot_test_float_cuda(PreferredBackend::cuda);
#endif

/**
 * \brief Test class for the batched dot product kernel.
 *
 * \test Tests whether the standard dot product and the batched dot product arrive at the same result
 *
 * \author Maximilian Esser
 */
template<
  typename DT_,
  typename IT_>
class CumulativeAxpyTest
  : public UnitTest
{
public:
  CumulativeAxpyTest(PreferredBackend backend)
    : UnitTest("CommulativeAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~CumulativeAxpyTest()
  {
  }


  virtual void run() const override
  {
    constexpr std::size_t dot_size = 1047;
    constexpr std::size_t num_dots = 17;
    const DT_ tol = Math::pow(Math::Limits<DT_>::epsilon(), DT_(0.4));
    {
      LAFEM::DenseVector<DT_, IT_> a(dot_size), b(dot_size);
      a.format();
      b.format();
      std::vector<LAFEM::DenseVector<DT_, IT_>> vecs;
      std::vector<Memory::Arbiter const*> arbiter;
      vecs.reserve(num_dots);
      arbiter.reserve(num_dots);
      for(std::size_t k = 0; k < num_dots; ++k)
      {
        vecs.emplace_back(dot_size);
        vecs.back().format();
        arbiter.emplace_back(&vecs.back().elements_arbiter());
        auto v_view = vecs.back().elements_view_w();
        DT_* v_data = v_view.get_w();
        for(std::size_t i = 0; i < dot_size; ++i)
        {
          v_data[i] = DT_((11*(i+1)) %17)/DT_(k+1);
        }
      }
      std::array<DT_, num_dots> res_a;
      for(std::size_t k = 0; k < num_dots; ++k)
      {
        res_a[k] = DT_(k+1)/DT_(1.7);
      }
      for(std::size_t k = 0; k < num_dots; ++k)
      {
        a.axpy(vecs.at(k), res_a[k] * DT_(-1.2));
      }

      // stride equals the number of vectors for a single DenseVector leaf (no tuple components)
      Solver::Intern::BatchedHelper<LAFEM::DenseVector<DT_, IT_>>::cumulative_axpy(b, arbiter.data(), res_a.data(), num_dots, DT_(-1.2), num_dots);

      b.axpy(a, DT_(-1));
      TEST_CHECK_EQUAL_WITHIN_EPS(b.norm2(), DT_(0), tol);
    }
  }
};

CumulativeAxpyTest<double, std::uint32_t> cumulative_axpy_test_double_generic(PreferredBackend::generic);
CumulativeAxpyTest<float, std::uint32_t> cumulative_axpy_test_float_generic(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
CumulativeAxpyTest<double, std::uint32_t> cumulative_axpy_test_double_cuda(PreferredBackend::cuda);
CumulativeAxpyTest<float, std::uint32_t> cumulative_axpy_test_float_cuda(PreferredBackend::cuda);
#endif