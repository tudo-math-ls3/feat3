// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/lafem/sparse_vector_factory.hpp>
#include <kernel/util/binary_stream.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse vector blocked class.
 *
 * \test test description missing
 *
 * \tparam DT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class SparseVectorBlockedTest
  : public UnitTest
{
public:
  explicit SparseVectorBlockedTest(PreferredBackend backend)
    : UnitTest("SparseVectorBlockedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    SparseVectorBlocked<DT_, IT_, 2> zero;
    TEST_CHECK(zero.empty());

    Tiny::Vector<DT_, 2> tv1, tv2;
    tv1[0] = DT_(17);
    tv1[1] = DT_(19);
    tv2[0] = DT_(41);
    tv2[1] = DT_(43);

    SparseVectorFactory<DT_, IT_, 2> a_fac(10);
    a_fac.add(3, tv1);
    a_fac.add(3, tv2);
    a_fac.add(6, tv1);
    a_fac.add(3, tv1);
    a_fac.add(1, tv1);
    a_fac.add(6, tv2);

    SparseVectorBlocked<DT_, IT_, 2> a(a_fac.make_svb());
    TEST_CHECK(!a.empty());
    TEST_CHECK_EQUAL(a.size(), Index(10));
    TEST_CHECK_EQUAL(a.num_nzes(), Index(3));
    TEST_CHECK_EQUAL(a.size_raw(), Index(20));
    TEST_CHECK_EQUAL(a.num_nzes_raw(), Index(6));

    {
      Memory::TypedView<IT_> vi(a.indices_view_r());
      TEST_CHECK_EQUAL(vi(0), IT_(1));
      TEST_CHECK_EQUAL(vi(1), IT_(3));
      TEST_CHECK_EQUAL(vi(2), IT_(6));
    }
    {
      Memory::TypedView<Tiny::Vector<DT_, 2>> vx(a.elements_view_r());
      TEST_CHECK_EQUAL(vx(0)[0], DT_(17));
      TEST_CHECK_EQUAL(vx(0)[1], DT_(19));
      TEST_CHECK_EQUAL(vx(1)[0], DT_(17));
      TEST_CHECK_EQUAL(vx(1)[1], DT_(19));
      TEST_CHECK_EQUAL(vx(2)[0], DT_(41));
      TEST_CHECK_EQUAL(vx(2)[1], DT_(43));
    }

    SparseVectorBlocked<DT_, IT_, 2> b(a.clone());
    TEST_CHECK_LESS_THAN(a.max_rel_diff(b), eps);

    a.elements_view_rw()[1] = tv2;
    TEST_CHECK(a.max_rel_diff(b) > DT_(0.4));

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    Adjacency::Permutation prm_rnd(a.size(), rng);
    SparseVectorBlocked<DT_, IT_, 2> ap(a.clone());
    ap.permute(prm_rnd);
    prm_rnd = prm_rnd.inverse();
    ap.permute(prm_rnd);
    TEST_CHECK_LESS_THAN(ap.max_rel_diff(a), eps);
    TEST_CHECK_EQUAL(ap.num_nzes(), Index(3));
  }
};
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseVectorBlockedSerializeTest
  : public UnitTest
{
public:
  explicit SparseVectorBlockedSerializeTest(PreferredBackend backend)
    : UnitTest("SparseVectorBlockedSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    Tiny::Vector<DT_, 2> tv1, tv2;
    tv1[0] = DT_(17);
    tv1[1] = DT_(19);
    tv2[0] = DT_(41);
    tv2[1] = DT_(43);

    SparseVectorFactory<DT_, IT_, 2> a_fac(10);
    a_fac.add(3, tv1);
    a_fac.add(3, tv2);
    a_fac.add(6, tv1);
    a_fac.add(3, tv1);
    a_fac.add(1, tv1);
    a_fac.add(6, tv2);

    SparseVectorBlocked<DT_, IT_, 2> a(a_fac.make_svb());

    BinaryStream bs;
    a.write_out(FileMode::fm_svb, bs);
    bs.seekg(0);
    SparseVectorBlocked<DT_, IT_, 2> bin(FileMode::fm_svb, bs);
    TEST_CHECK_LESS_THAN(bin.max_rel_diff(a), eps);

    auto op = a.serialize(LAFEM::SerialConfig(false, false));
    SparseVectorBlocked<DT_, IT_, 2> o(op);
    TEST_CHECK_LESS_THAN(a.max_rel_diff(o), eps);
#ifdef FEAT_HAVE_ZLIB
    auto zl = a.serialize(LAFEM::SerialConfig(true, false));
    SparseVectorBlocked<DT_, IT_, 2> zlib(zl);
    TEST_CHECK_LESS_THAN(zlib.max_rel_diff(a), eps);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = a.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseVectorBlocked<DT_, IT_, 2> zfp(zf);
    TEST_CHECK_LESS_THAN(bin.max_rel_diff(a), DT_(1e-7));
#endif
  }
};
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseVectorBlockedMaxRelDiffTest
  : public UnitTest
{
public:
  explicit SparseVectorBlockedMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseVectorBlockedMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    static constexpr int block_size = 2;
    const Index size = 100;
    SparseVectorFactory<DT_, IT_, block_size> a_fac(size);

    for(Index k = 0; k < 17; ++k)
      a_fac.add(((k+23)*131071) % size, DT_(k+1));

    SparseVectorBlocked<DT_, IT_, block_size> a(a_fac.make_svb());
    SparseVectorBlocked<DT_, IT_, block_size> b(a.clone());

    TEST_CHECK_LESS_THAN(a.max_rel_diff(b), eps*eps);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(a), eps*eps);

    b.elements_view_rw()[7][1] += DT_(0.314);
    TEST_CHECK_LESS_THAN(eps, b.max_rel_diff(a));
  }
};
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedMaxRelDiffTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseVectorBlockedSameLayoutTest
  : public UnitTest
{
public:
  explicit SparseVectorBlockedSameLayoutTest(PreferredBackend backend)
    : UnitTest("SparseVectorBlockedSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    static constexpr int block_size = 2;

    const Index size = 100;
    SparseVectorFactory<DT_, IT_, block_size> a_fac(size), z_fac(size+2);

    for(Index k = 0; k < 17; ++k)
    {
      a_fac.add(((k+23)*131071) % size, DT_(k+1));
      z_fac.add(((k+23)*131071) % size, DT_(k+1));
    }

    SparseVectorBlocked<DT_, IT_, block_size> a(a_fac.make_svb());

    // different sizes
    SparseVectorBlocked<DT_, IT_, block_size> z(z_fac.make_svb());
    TEST_CHECK(!a.same_layout(z));

    // weak copy
    auto b = a.clone(CloneMode::Weak);
    TEST_CHECK(a.same_layout(b));

    // shallow copy
    auto c = a.clone(CloneMode::Shallow);
    TEST_CHECK(a.same_layout(c));

    // different values at same position
    SparseVectorBlocked<DT_, IT_, block_size> d(a.clone());
    d.elements_view_rw()[7][1] += DT_(0.5);
    TEST_CHECK(a.same_layout(d));
  }
};

SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorBlockedSameLayoutTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
