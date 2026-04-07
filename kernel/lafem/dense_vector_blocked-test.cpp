// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/binary_stream.hpp>

#include <list>
#include <sstream>
#include <cstdio>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the dense vector blocked class.
 *
 * \test test description missing
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class DenseVectorBlockedBasicTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedBasicTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedBasicTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    DenseVectorBlocked<DT_, IT_, 2> zero;
    TEST_CHECK(zero.empty());
    TEST_CHECK_EQUAL(zero.bytes(), 0);

    DenseVectorBlocked<DT_, IT_, 2> a(10, DT_(7));
    TEST_CHECK(!a.empty());
    TEST_CHECK_EQUAL(a.size(), 10);
    TEST_CHECK_EQUAL(a.size_raw(), 20);
    TEST_CHECK_EQUAL(a.bytes(), 20 * sizeof(DT_));

    DenseVectorBlocked<DT_, IT_, 2> b(10, DT_(5));
    b.elements_view_rw()[7] = DT_(42);

    DenseVectorBlocked<DT_, IT_, 2> b_r(b, 4, 6);
    {
      auto vb = b.elements_view_r();
      auto vb_r = b_r.elements_view_r();
      TEST_CHECK_EQUAL(vb_r(0)[0], vb(0+6)[0]);
      TEST_CHECK_EQUAL(vb_r(0)[1], vb(0+6)[1]);
      TEST_CHECK_EQUAL(vb_r(1)[0], vb(1+6)[0]);
      TEST_CHECK_EQUAL(vb_r(1)[1], vb(1+6)[1]);
    }

    DenseVectorBlocked<DT_, IT_, 2> c(b.clone());
    TEST_CHECK_EQUAL(c.size(), b.size());
    {
      auto cv = c.elements_view_r()(Index(7));
      auto bv = b.elements_view_r()(Index(7));
      for (int i(0) ; i < 2 ; ++i)
      {
        TEST_CHECK_EQUAL(cv[i], bv[i]);
      }
    }
    TEST_CHECK_LESS_THAN(c.max_rel_diff(b), tol);

    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    {
      auto cv = c.elements_view_r()(Index(7));
      auto bv = b.elements_view_r()(Index(7));
      for (int i(0) ; i < 2 ; ++i)
      {
        TEST_CHECK_EQUAL(cv[i], bv[i]);
      }
    }
    TEST_CHECK_LESS_THAN(c.max_rel_diff(b), tol);

    DenseVectorBlocked<float, unsigned int, 2> d;
    d.convert(c);
    DenseVectorBlocked<float, unsigned int, 2> e;
    e.convert(b);
    TEST_CHECK_EQUAL(e.size(), d.size());
    {
      auto dv = d.elements_view_r();
      auto ev = e.elements_view_r();
      for (int i(0) ; i < 2 ; ++i)
      {
        TEST_CHECK_EQUAL(ev(7)[i], dv(7)[i]);
      }
    }
    TEST_CHECK_LESS_THAN(e.max_rel_diff(d), float(tol));

    b.clone(a);
    TEST_CHECK_NOT_EQUAL(b.elements_arbiter(), a.elements_arbiter());
    c.convert(a);
    TEST_CHECK_EQUAL(c.elements_arbiter(), a.elements_arbiter());
    TEST_CHECK_LESS_THAN(b.max_rel_diff(c), tol);
    a.elements_view_rw()[3] = DT_(23);
    TEST_CHECK_LESS_THAN(a.max_rel_diff(c), tol);
    TEST_CHECK_LESS_THAN(tol, a.max_rel_diff(b));

    DenseVector<DT_, IT_> dv(12, DT_(2));
    dv.elements_view_rw()[7] = DT_(3);
    DenseVectorBlocked<DT_, IT_, 3> f(dv);
    auto t3 = f.elements_view_r()(2);
    TEST_CHECK_EQUAL(t3.v[0], DT_(2));
    TEST_CHECK_EQUAL(t3.v[1], DT_(3));
    TEST_CHECK_EQUAL(f.elements_arbiter(), dv.elements_arbiter());
    DenseVector<DT_, IT_> dv2(f);
    TEST_CHECK_LESS_THAN(dv2.max_rel_diff(dv), tol);
    TEST_CHECK_EQUAL(dv2.elements_arbiter(), dv.elements_arbiter());

    DenseVectorBlocked<DT_, IT_, 3> g(f.size(), f.elements_arbiter().attach());
    TEST_CHECK_LESS_THAN(g.max_rel_diff(f), tol);
    TEST_CHECK_EQUAL(g.elements_arbiter(), f.elements_arbiter());

    // random constructor check
    DT_ rnd_range[2];
    IT_ rnd_size = 3*1234;
    rnd_range[0] = DT_(-10);
    rnd_range[1] = DT_(+10);
    DenseVectorBlocked<DT_, IT_, 3> rnd_vec(rng, rnd_size, rnd_range[0], rnd_range[1]);
    TEST_CHECK_EQUAL(rnd_vec.size(), rnd_size);
    DT_ rnd_max = rnd_vec.max_abs_element();
    TEST_CHECK_IN_RANGE(rnd_max, rnd_range[0], rnd_range[1]);
    rnd_vec.scale(rnd_vec, DT_(-1));
    DT_ rnd_min = -rnd_vec.max_abs_element();
    TEST_CHECK_IN_RANGE(rnd_min, rnd_range[0], rnd_range[1]);

  }
};
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedBasicTest, double, std::uint64_t, PreferredBackend::cuda);
#endif


template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMaxRelDiffTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(12); size < Index(500); size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);

      const Index size_raw = a.size_raw();
      const Index off0 = (3*size_raw) / 8;
      const Index off1 = (1*size_raw) / 8;
      const Index off2 = (6*size_raw) / 8;

      // a = i, b = i
      {
        Memory::TypedView<DT_> va = a.elements_view_raw_w();
        Memory::TypedView<DT_> vb = b.elements_view_raw_w();
        for (Index i(0); i < size_raw; ++i)
        {
          va[i] = vb[i] = DT_(int(i) - int(off0)) * DT_(0.123);
        }
      }

      // identical vectors, result should be zero
      TEST_CHECK_LESS_THAN(a.max_rel_diff(b), tol*tol);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(a), tol*tol);

      // two values close to zero
      const DT_ delta_a0(Math::sqrt(Math::eps<DT_>()));
      const DT_ delta_b0(Math::sqr(Math::eps<DT_>()));
      const DT_ ref0 = (delta_a0 + delta_b0) / (DT_(1) + delta_a0 + delta_b0);
      a.elements_view_raw_w()[off0] += delta_a0;
      b.elements_view_raw_w()[off0] -= delta_b0;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref0, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref0, tol);

      const DT_ delta1 = DT_(0.17);
      const DT_ ref1 = delta1 / (DT_(off0 - off1)*DT_(0.246) + delta1 + DT_(1));
      a.elements_view_raw_w()[off1] -= delta1;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref1, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref1, tol);

      const DT_ delta2 = DT_(0.73);
      const DT_ ref2 = delta2 / (DT_(off2 - off0)*DT_(0.246) + delta2 + DT_(1));
      b.elements_view_raw_w()[off2] += delta2;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref2, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref2, tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, double, std::uint64_t, 3, PreferredBackend::generic);
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, float, std::uint64_t, 2, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, double, std::uint64_t, 3, PreferredBackend::mkl);
//#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, Half, std::uint64_t, 3, PreferredBackend::generic);
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, Half, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, Half, std::uint64_t, 3, PreferredBackend::cuda);
//#endif
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxRelDiffTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedSameLayoutTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedSameLayoutTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const Index size = 72;

    DenseVectorBlocked<DT_, IT_, block_size_> a(size, DT_(0));

    // weak copy
    DenseVectorBlocked<DT_, IT_, block_size_> b = a.clone(CloneMode::Weak);
    TEST_CHECK(a.same_layout(b));

    // shallow copy
    DenseVectorBlocked<DT_, IT_, block_size_> c = a.clone(CloneMode::Shallow);
    TEST_CHECK(a.same_layout(c));

    // same size
    DenseVectorBlocked<DT_, IT_, block_size_> d(size, DT_(0));
    TEST_CHECK(a.same_layout(d));

    // different size
    DenseVectorBlocked<DT_, IT_, block_size_> e(size + 2, DT_(0));
    TEST_CHECK(!a.same_layout(e));

    // empty vector
    d.clear();
    TEST_CHECK(!a.same_layout(d));

    // two empty vectors
    e.clear();
    TEST_CHECK(e.same_layout(d));
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, Half, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, Half, std::uint64_t, 3, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, float, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, double, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, float, std::uint64_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedSameLayoutTest, double, std::uint64_t, 3, PreferredBackend::cuda);
#endif


template<
  typename DT_,
  typename IT_>
class DenseVectorBlockedSerializeTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedSerializeTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    const Index io_vector_size = 16*3*5*7;

    DenseVector<DT_, IT_> k(io_vector_size);
    {
      Memory::TypedView<DT_> vk = k.elements_view_w();
      for (Index i(0) ; i < k.size() ; ++i)
        vk[i] = DT_(i) / DT_(12);
    }
    DenseVectorBlocked<DT_, IT_, 3> g(k);

    {
      std::stringstream ioss;
      g.write_out(FileMode::fm_mtx, ioss);
      DenseVectorBlocked<DT_, IT_, 3> l(FileMode::fm_mtx, ioss);
      TEST_CHECK_LESS_THAN(l.max_rel_diff(g), tol);
    }

    {
      std::stringstream ioss;
      g.write_out(FileMode::fm_exp, ioss);
      DenseVectorBlocked<DT_, IT_, 3> m(FileMode::fm_exp, ioss);
      TEST_CHECK_LESS_THAN(m.max_rel_diff(g), tol);
    }

    {
      BinaryStream bs;
      g.write_out(FileMode::fm_dvb, bs);
      bs.seekg(0);
      DenseVectorBlocked<DT_, IT_, 3> n(FileMode::fm_dvb, bs);
      TEST_CHECK_LESS_THAN(n.max_rel_diff(g), tol);
    }

    {
#ifdef FEAT_HAVE_ZLIB
      auto zb = g.serialize(LAFEM::SerialConfig(true,false));
      DenseVectorBlocked<DT_, IT_, 3> z(zb);
      TEST_CHECK_LESS_THAN(z.max_rel_diff(g), tol);
#endif
    }

    {
#ifdef FEAT_HAVE_ZFP
      auto zp = g.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-5)));
      DenseVectorBlocked<DT_, IT_, 3> y(zp);
      TEST_CHECK_LESS_THAN(y.max_rel_diff(g), DT_(2e-5));
#endif
    }
  }
};
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorBlockedSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedAxpyTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedAxpyTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    typedef Tiny::Vector<DT_, block_size_> ValueType;

    DT_ s(DT_(0.4711));
    ValueType bs;
    for(int j(0); j < block_size_; ++j)
      bs[j] = s + DT_(j);

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref_a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref_b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref_c(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref_d(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> vb(b.elements_view_w());
        Memory::TypedView<ValueType> vref_a(ref_a.elements_view_w());
        Memory::TypedView<ValueType> vref_b(ref_b.elements_view_w());

        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
          {
            va[i][j] = DT_(i % 100) * DT_(1.234) + DT_(j);
            vb[i][j] = DT_(2) + DT_(i % Index(42 + j));
          }

          vref_a[i] = s * va(i) + vb(i);
          vref_b[i] = s * vb(i) + vb(i);
        }
      }

      // r != x
      a.scale(a, s);
      a.axpy(b); /// \todo use axpby here
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref_a), tol);

      // r == x
      b.axpy(b, s);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref_b), tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, float, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, double, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, float, std::uint64_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyTest, double, std::uint64_t, 3, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedAxpyBlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedAxpyBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedAxpyBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    typedef Tiny::Vector<DT_, block_size_> ValueType;

    DT_ s(DT_(0.4711));
    ValueType bs;
    for(int j(0); j < block_size_; ++j)
      bs[j] = s + DT_(j);

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref_a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref_b(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> vb(b.elements_view_w());
        Memory::TypedView<ValueType> vref_a(ref_a.elements_view_w());
        Memory::TypedView<ValueType> vref_b(ref_b.elements_view_w());

        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
          {
            va[i][j] = DT_(i % 100) * DT_(1.234) + DT_(j);
            vb[i][j] = DT_(2) + DT_(i % Index(42 + j));
            vref_a[i][j] = va[i][j] * (DT_(1) + bs[j]);
            vref_b[i][j] = vb[i][j] + bs[j] * va[i][j];
          }
        }
      }

      // r != x
      b.axpy_blocked(a, bs);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref_b), tol);

      // r == x
      a.axpy_blocked(a, bs);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref_a), tol);
    }
  }
};

SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedAxpyBlockedTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedDotTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedDotTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> a2(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      const DT_ den(DT_(1) / DT_(size * block_size_));

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> va2(a2.elements_view_w());
        Memory::TypedView<ValueType> vb(b.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
          {
            va[i][j] = (DT_(i) * DT_(block_size_) + DT_(j + 1)) * den;
            va2[i][j] = Math::sqrt(DT_(i))*DT_(j);
            vb[i][j] = DT_(1) / DT_(int(i) * block_size_ + j + 1);
          }
        }
      }

      // a*b = 1
      DT_ ref(DT_(1));
      DT_ c  = a.dot(b);
      TEST_CHECK_RELATIVE(c, ref, tol);
      c = b.dot(a);
      TEST_CHECK_RELATIVE(c, ref, tol);
      c = b.dot(b);
      ref = b.norm2();
      ref *= ref;
      TEST_CHECK_RELATIVE(c, ref, tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, double, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, double, std::uint64_t, 2, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, double, std::uint64_t, 2, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, Half, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, Half, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, float, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, double, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, float, std::uint64_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotTest, double, std::uint64_t, 2, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedDotBlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedDotBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedDotBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> a2(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      const DT_ den(DT_(1) / DT_(size * block_size_));

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> va2(a2.elements_view_w());
        Memory::TypedView<ValueType> vb(b.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
          {
            va[i][j] = (DT_(i) * DT_(block_size_) + DT_(j + 1)) * den;
            va2[i][j] = Math::sqrt(DT_(i))*DT_(j);
            vb[i][j] = DT_(1) / DT_(int(i) * block_size_ + j + 1);
          }
        }
      }

      Tiny::Vector<DT_, block_size_> d(a.dot_blocked(b));
      Tiny::Vector<DT_, block_size_> e(b.dot_blocked(a));
      Tiny::Vector<DT_, block_size_> f(a2.dot_blocked(a2));

      Tiny::Vector<DT_, block_size_> ref_bs;
      Tiny::Vector<DT_, block_size_> ref2_bs;
      for (int j(0) ; j < block_size_ ; ++j)
      {
        ref_bs.v[j] = DT_(size)/DT_(size*block_size_);
        ref2_bs.v[j] = DT_(j*j)*(DT_(size*size)-DT_(size))/DT_(2);
      }

      for (int j(0) ; j < block_size_ ; ++j)
      {
        TEST_CHECK_RELATIVE(d.v[j], ref_bs.v[j], tol);
        TEST_CHECK_RELATIVE(e.v[j], ref_bs.v[j], tol);
        TEST_CHECK_RELATIVE(f.v[j], ref2_bs.v[j], tol);
      }
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, double, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, double, std::uint64_t, 2, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, double, std::uint64_t, 2, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, Half, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, Half, std::uint64_t, 2, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, float, std::uint32_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, double, std::uint32_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, float, std::uint64_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedDotBlockedTest, double, std::uint64_t, 2, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedTripleDotTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedTripleDotTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedTripleDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> c(size);
      Tiny::Vector<DT_, block_size_> ref_bs(0);
      const DT_ den(DT_(1) / DT_(size * block_size_));

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> vb(b.elements_view_w());
        Memory::TypedView<ValueType> vc(c.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
          {
            va[i][j] = DT_(int(i) * block_size_ + j + 1) * den;
            vb[i][j] = DT_(1) / DT_(int(i) * block_size_ + j + 1);
            vc[i][j] =  DT_(int(i) * block_size_ + j + 1);
          }
        }
      }

      DenseVector<DT_, IT_> ref_a;
      ref_a.convert(a);
      DenseVector<DT_, IT_> ref_b;
      ref_b.convert(b);
      DenseVector<DT_, IT_> ref_c;
      ref_c.convert(c);

      DT_ ref(ref_a.triple_dot(ref_b, ref_c));
      DT_ res  = a.triple_dot(b, c);
      TEST_CHECK_RELATIVE(res, ref, tol);
      res = b.triple_dot(a, c);
      TEST_CHECK_RELATIVE(res, ref, tol);
      res = c.triple_dot(b, a);
      TEST_CHECK_RELATIVE(res, ref, tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, double, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, double, std::uint64_t, 2, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, double, std::uint64_t, 2, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, Half, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, Half, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, float, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, double, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, float, std::uint64_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotTest, double, std::uint64_t, 2, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedTripleDotBlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedTripleDotBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedTripleDotBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> c(size);
      Tiny::Vector<DT_, block_size_> ref_bs(0);
      const DT_ den(DT_(1) / DT_(size * block_size_));

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> vb(b.elements_view_w());
        Memory::TypedView<ValueType> vc(c.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
          {
            va[i][j] = DT_(int(i) * block_size_ + j + 1) * den;
            vb[i][j] = DT_(1) / DT_(int(i) * block_size_ + j + 1);
            vc[i][j] =  DT_(int(i) * block_size_ + j + 1);
          }
        }
      }

      for(int j(0); j<block_size_; ++j)
        ref_bs.v[j] = DT_(size-1)/DT_(2)+DT_(j+1)/DT_(block_size_);

      Tiny::Vector<DT_, block_size_> d(a.triple_dot_blocked(b, c));
      Tiny::Vector<DT_, block_size_> e(b.triple_dot_blocked(a, c));
      Tiny::Vector<DT_, block_size_> f(c.triple_dot_blocked(a, b));

      for(int j(0); j < block_size_; ++j)
      {
        TEST_CHECK_RELATIVE(d.v[j], ref_bs.v[j], tol);
        TEST_CHECK_RELATIVE(e.v[j], ref_bs.v[j], tol);
        TEST_CHECK_RELATIVE(f.v[j], ref_bs.v[j], tol);
      }
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, double, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, double, std::uint64_t, 2, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, double, std::uint64_t, 2, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, Half, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, Half, std::uint64_t, 2, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, float, std::uint32_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, double, std::uint32_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, float, std::uint64_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedTripleDotBlockedTest, double, std::uint64_t, 2, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedComponentProductTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedComponentProductTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedComponentProductTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> b(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref2(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> vb(b.elements_view_w());
        Memory::TypedView<ValueType> vref(ref.elements_view_w());
        Memory::TypedView<ValueType> vref2(ref2.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
          {
            va[i][j] = DT_((DT_(j) + DT_(i)) * DT_(1.234) / DT_(size));
            vb[i][j] = DT_(int(size)*2*block_size_ - int(i) + 2*j);
          }
          vref[i] = Tiny::component_product(va(i), vb(i));
          vref2[i] = Tiny::component_product(va(i), va(i));
        }
      }

      DenseVectorBlocked<DT_, IT_, block_size_> a_tmp(a.clone());
      DenseVectorBlocked<DT_, IT_, block_size_> b_tmp(b.clone());
      DenseVectorBlocked<DT_, IT_, block_size_> c(size);

      c.component_product(a, b);
      TEST_CHECK_LESS_THAN(c.max_rel_diff(ref), tol);

      a.component_product(a, b);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), tol);

      a.copy(a_tmp);
      b.component_product(a, b);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref), tol);

      a.component_product(a, a);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref2), tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, double, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, double, std::uint64_t, 2, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, double, std::uint64_t, 2, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, Half, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, Half, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, float, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, double, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, float, std::uint64_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentProductTest, double, std::uint64_t, 2, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedComponentCopyTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedComponentCopyTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedComponentCopyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for(Index size(1); size < Index(1000); size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVector<DT_,IT_> v(size);

      {
        Memory::TypedView<DT_> vv(v.elements_view_w());
        for(Index i(0); i < size; ++i)
          vv[i] = DT_(i + 1);
      }

      a.format();
      a.component_copy(v, 0);
      {
        Memory::TypedView<ValueType> va(a.elements_view_r());
        for(Index i(0); i < size; ++i)
        {
          TEST_CHECK_RELATIVE(va(i)[0], DT_(i + 1), tol);
          for(int j = 1; j < block_size_; ++j)
            TEST_CHECK_EQUAL_WITHIN_EPS(va(i)[j], DT_(0), tol);
        }
      }

      DenseVector<DT_, IT_> w(size);
      w.format();
      a.component_copy_to(w, 0);
      TEST_CHECK_LESS_THAN(w.max_rel_diff(v), tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, double, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, double, std::uint64_t, 2, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, double, std::uint64_t, 2, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, Half, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, Half, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, float, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, double, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, float, std::uint64_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedComponentCopyTest, double, std::uint64_t, 2, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedScaleTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedScaleTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DT_ s(DT_(4.321));

      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref(size);
      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> vref(ref.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j] = DT_(DT_(i) * DT_(block_size_) + DT_(j) * DT_(1.234));
          vref[i] = va(i) * s;
        }
      }

      DenseVectorBlocked<DT_, IT_, block_size_> b(size);

      b.scale(a, s);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref), tol);

      a.scale(a, s);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, double, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, double, std::uint64_t, 2, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, double, std::uint64_t, 2, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, Half, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, Half, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, float, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, double, std::uint32_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, float, std::uint64_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleTest, double, std::uint64_t, 2, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedScaleBlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedScaleBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedScaleBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      Tiny::Vector<DT_, block_size_> bs;
      for(int j(0); j < block_size_; ++j)
        bs[j] = DT_(4.321) + DT_(j)*DT_(0.1);

      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      DenseVectorBlocked<DT_, IT_, block_size_> ref(size);
      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        Memory::TypedView<ValueType> vref(ref.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
          {
            va[i][j] = DT_(DT_(i) * DT_(block_size_) + DT_(j) * DT_(1.234));
            vref[i][j] = va[i][j] * bs[j];
          }
        }
      }

      DenseVectorBlocked<DT_, IT_, block_size_> b(size);

      b.scale_blocked(a, bs);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref), tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, float, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, double, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, float, std::uint64_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, double, std::uint64_t, 2, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, float, std::uint64_t, 3, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, double, std::uint64_t, 2, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, __float128, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, __float128, std::uint64_t, 2, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, Half, std::uint32_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, Half, std::uint64_t, 2, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, float, std::uint32_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, double, std::uint32_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, float, std::uint64_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedScaleBlockedTest, double, std::uint64_t, 2, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedNorm2Test
  : public UnitTest
{
public:
  explicit DenseVectorBlockedNorm2Test(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedNorm2Test", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      Tiny::Vector<DT_, block_size_> ref_bs(0);
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j] = Math::pow(DT_(0.5), DT_(0.5) * DT_(int(i) * block_size_ + j));
        }
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      const DT_ ref(Math::sqrt(DT_(2) - Math::pow(DT_(0.5), DT_(int(size)*block_size_-1))));

      DT_ c = a.norm2();
      TEST_CHECK_RELATIVE(c, ref, tol);

      c = a.norm2sqr();
      TEST_CHECK_RELATIVE(c, ref*ref, tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, float, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, double, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, float, std::uint64_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2Test, double, std::uint64_t, 3, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedNorm2BlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedNorm2BlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedNorm2BlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      Tiny::Vector<DT_, block_size_> ref_bs(0);
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j] = Math::pow(DT_(0.5), DT_(0.5) * DT_(int(i) * block_size_ + j));
        }
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      for (int j(0) ; j < block_size_ ; ++j)
      {
        ref_bs.v[j] = Math::pow(DT_(0.5), DT_(j))*(DT_(1)-Math::pow(DT_(0.5), DT_(block_size_*size)))/(DT_(1)-Math::pow(DT_(0.5), DT_(block_size_)));
        ref_bs.v[j] = Math::sqrt(ref_bs.v[j]);
      }

      Tiny::Vector<DT_, block_size_> d(a.norm2_blocked());
      for(Index j(0); j < block_size_; ++j)
        TEST_CHECK_RELATIVE(d.v[j], ref_bs.v[j], tol);

      Tiny::Vector<DT_, block_size_> e(a.norm2sqr_blocked());
      for(Index j(0); j < block_size_; ++j)
        TEST_CHECK_RELATIVE(e.v[j], ref_bs.v[j]*ref_bs.v[j], tol);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, double, std::uint64_t, 3, PreferredBackend::generic);
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, float, std::uint64_t, 2, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, double, std::uint64_t, 3, PreferredBackend::mkl);
//#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedNorm2BlockedTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMaxAbsElementTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMaxAbsElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMaxAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        }
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);
      DT_ max = a.max_abs_element();

      TEST_CHECK_EQUAL(max, DT_((size*block_size_) -1));
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, float, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, double, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, float, std::uint64_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementTest, double, std::uint64_t, 3, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMaxAbsElementBlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMaxAbsElementBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMaxAbsElementBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        }
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      Tiny::Vector<DT_, block_size_> b(a.max_abs_element_blocked());
      Tiny::Vector<DT_, block_size_> ref_bs;
      for (int j(0) ; j < block_size_ ; ++j)
        ref_bs.v[j] = DT_((int(size)-1) * block_size_  + j);

      for(int j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL(b.v[j], ref_bs.v[j]);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxAbsElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMinAbsElementTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMinAbsElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMinAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        }
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ min = a.min_abs_element();

      TEST_CHECK_EQUAL(min, DT_(0));
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, float, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, double, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, float, std::uint64_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementTest, double, std::uint64_t, 3, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMinAbsElementBlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMinAbsElementBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMinAbsElementBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);

      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        }
      }

      //Adjacency::Permutation prm_rnd(a.size(), rng);
      //a.permute(prm_rnd);

      Tiny::Vector<DT_, block_size_> b(a.min_abs_element_blocked());
      Tiny::Vector<DT_, block_size_> ref_bs;
      for (int j(0) ; j < block_size_ ; ++j)
        ref_bs.v[j] = DT_(j);

      for(int j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL(b.v[j], ref_bs.v[j]);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinAbsElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMaxElementTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMaxElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMaxElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)));
        }
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ max = a.max_element();

      TEST_CHECK_EQUAL(max, DT_((size*block_size_) -1));
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, double, std::uint64_t, 3, PreferredBackend::generic);
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, float, std::uint64_t, 2, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, double, std::uint64_t, 3, PreferredBackend::mkl);
//#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMaxElementBlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMaxElementBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMaxElementBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)));
        }
      }

      //Adjacency::Permutation prm_rnd(a.size(), rng);
      //a.permute(prm_rnd);

      Tiny::Vector<DT_, block_size_> b(a.max_element_blocked());
      Tiny::Vector<DT_, block_size_> ref_bs;
      for (int j(0) ; j < block_size_ ; ++j)
        ref_bs.v[j] = DT_((int(size)-1)*block_size_+j);

      for(int j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL(b.v[j], ref_bs.v[j]);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::generic);
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, float, std::uint64_t, 2, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::mkl);
//#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMaxElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMinElementTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMinElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMinElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_size_> ValueType;
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)));
        }
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ min = a.min_element();

      TEST_CHECK_EQUAL(min, DT_(0));
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, double, std::uint64_t, 3, PreferredBackend::generic);
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, float, std::uint64_t, 2, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, double, std::uint64_t, 3, PreferredBackend::mkl);
//#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedMinElementBlockedTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedMinElementBlockedTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedMinElementBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    typedef Tiny::Vector<DT_, block_size_> ValueType;
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    for (Index size(1) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)));
        }
      }

      //Adjacency::Permutation prm_rnd(a.size(), rng);
      //a.permute(prm_rnd);

      Tiny::Vector<DT_, block_size_> b(a.min_element_blocked());
      Tiny::Vector<DT_, block_size_> ref_bs;
      for (int j(0) ; j < block_size_ ; ++j)
        ref_bs.v[j] = DT_(j);

      for(int j(0); j < block_size_; ++j)
        TEST_CHECK_EQUAL(b.v[j], ref_bs.v[j]);
    }
  }
};
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::generic);
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, float, std::uint64_t, 2, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::mkl);
//#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, float, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, double, std::uint32_t, 2, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, float, std::uint64_t, 3, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedMinElementBlockedTest, double, std::uint64_t, 3, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_,
  int block_size_>
class DenseVectorBlockedPermuteTest
  : public UnitTest
{
public:
  explicit DenseVectorBlockedPermuteTest(PreferredBackend backend)
    : UnitTest("DenseVectorBlockedPermuteTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    typedef Tiny::Vector<DT_, block_size_> ValueType;

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(10) ; size < Index(1000) ; size *= 2)
    {
      DenseVectorBlocked<DT_, IT_, block_size_> a(size);
      {
        Memory::TypedView<ValueType> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          for (int j(0) ; j < block_size_ ; ++j)
            va[i][j]  = DT_(DT_((int(i) * block_size_ + j)) * (i%2 == 0 ? DT_(1) : DT_(-1)));
        }
      }

      DenseVectorBlocked<DT_, IT_, block_size_> ref(a.clone());

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      //allow identity permutation once per test due to rng...
      if(prm_rnd.is_identity())
      {
        TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), tol);
      }
      else
      {
        TEST_CHECK(a.max_rel_diff(ref) > tol);
      }

      auto prm_inv = prm_rnd.inverse();
      a.permute(prm_inv);

      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), tol);
    }
  }
};

SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, float, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, double, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, float, std::uint64_t, 3, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, double, std::uint64_t, 3, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, float, std::uint64_t, 2, PreferredBackend::mkl);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, double, std::uint64_t, 3, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, __float128, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, __float128, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, Half, std::uint32_t, 2, PreferredBackend::generic);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, Half, std::uint64_t, 3, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, float, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, double, std::uint32_t, 2, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, float, std::uint64_t, 3, PreferredBackend::cuda);
SPAWN_UNIT_TEST_3T_P(DenseVectorBlockedPermuteTest, double, std::uint64_t, 3, PreferredBackend::cuda);
#endif
