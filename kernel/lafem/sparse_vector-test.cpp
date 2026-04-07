// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_vector.hpp>
#include <kernel/lafem/sparse_vector_factory.hpp>
#include <kernel/util/binary_stream.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse vector class.
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
class SparseVectorTest
  : public UnitTest
{
public:
  explicit SparseVectorTest(PreferredBackend backend)
    : UnitTest("SparseVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    SparseVector<DT_, IT_> zero;
    TEST_CHECK(zero.empty());

    SparseVectorFactory<DT_, IT_> a_fac(10);
    a_fac.add(3, DT_(7));
    a_fac.add(3, DT_(2));
    a_fac.add(6, DT_(1));
    a_fac.add(5, DT_(6));
    a_fac.add(6, DT_(8));

    SparseVector<DT_, IT_> a(a_fac.make_sv());
    TEST_CHECK(!a.empty());
    TEST_CHECK_EQUAL(a.size(), Index(10));
    TEST_CHECK_EQUAL(a.num_nzes(), Index(3));
    {
      Memory::TypedView<IT_> vi(a.indices_view_r());
      TEST_CHECK_EQUAL(vi(0), IT_(3));
      TEST_CHECK_EQUAL(vi(1), IT_(5));
      TEST_CHECK_EQUAL(vi(2), IT_(6));
    }
    {
      Memory::TypedView<DT_> va(a.elements_view_r());
      TEST_CHECK_EQUAL(va(0), DT_(2));
      TEST_CHECK_EQUAL(va(1), DT_(6));
      TEST_CHECK_EQUAL(va(2), DT_(8));
    }

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    Adjacency::Permutation prm_rnd(a.size(), rng);
    SparseVector<DT_, IT_> ap(a.clone());
    ap.permute(prm_rnd);
    prm_rnd = prm_rnd.inverse();
    ap.permute(prm_rnd);
    TEST_CHECK_LESS_THAN(ap.max_rel_diff(a), eps);
    TEST_CHECK_EQUAL(ap.num_nzes(), Index(3));

    SparseVector<DT_, IT_> b;
    b.convert(a);
    TEST_CHECK_LESS_THAN(a.max_rel_diff(b), eps);
    b.elements_view_rw()[2] = DT_(1);
    TEST_CHECK_LESS_THAN(eps, a.max_rel_diff(b));
    b.clone(a);
    b.elements_view_rw()[2] = DT_(3);
    TEST_CHECK_LESS_THAN(eps, a.max_rel_diff(b));
    TEST_CHECK_NOT_EQUAL(a.elements_arbiter(), b.elements_arbiter());
    TEST_CHECK_NOT_EQUAL(a.indices_arbiter(), b.indices_arbiter());
    b = a.clone();
    TEST_CHECK_NOT_EQUAL(a.elements_arbiter(), b.elements_arbiter());
    TEST_CHECK_NOT_EQUAL(a.indices_arbiter(), b.indices_arbiter());

    SparseVector<float, unsigned int> c;
    c.convert(a);
    SparseVector<float, unsigned int> d;
    d.clone(c);
    SparseVector<float, unsigned int> e;
    e.convert(a);
    TEST_CHECK_LESS_THAN(d.max_rel_diff(e), float(eps));
    c.elements_view_rw()[2] = DT_(1);
    TEST_CHECK_LESS_THAN(float(eps), c.max_rel_diff(e));

    a.format();
    TEST_CHECK_EQUAL(a.num_nzes(), Index(3));
    {
      Memory::TypedView<DT_> va(a.elements_view_r());
      TEST_CHECK_EQUAL(va(0), DT_(0));
      TEST_CHECK_EQUAL(va(1), DT_(0));
      TEST_CHECK_EQUAL(va(2), DT_(0));
    }
  }
};
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseVectorSerializeTest
  : public UnitTest
{
public:
  explicit SparseVectorSerializeTest(PreferredBackend backend)
    : UnitTest("SparseVectorSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    SparseVectorFactory<DT_, IT_> a_fac(10);
    a_fac.add(3, DT_(7));
    a_fac.add(3, DT_(2));
    a_fac.add(6, DT_(1));
    a_fac.add(5, DT_(6));
    a_fac.add(6, DT_(8));

    SparseVector<DT_, IT_> a(a_fac.make_sv());

    std::stringstream ts;
    a.write_out(FileMode::fm_mtx, ts);
    SparseVector<DT_, IT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_LESS_THAN(j.max_rel_diff(a), eps);

    BinaryStream bs;
    a.write_out(FileMode::fm_sv, bs);
    bs.seekg(0);
    SparseVector<DT_, IT_> bin(FileMode::fm_sv, bs);
    TEST_CHECK_LESS_THAN(bin.max_rel_diff(a), eps);

    auto op = a.serialize(LAFEM::SerialConfig(false, false));
    SparseVector<DT_, IT_> o(op);
    TEST_CHECK_LESS_THAN(o.max_rel_diff(a), eps);
#ifdef FEAT_HAVE_ZLIB
    auto zl = a.serialize(LAFEM::SerialConfig(true, false));
    SparseVector<DT_, IT_> zlib(zl);
    TEST_CHECK_LESS_THAN(zlib.max_rel_diff(a), eps);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = a.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseVector<DT_, IT_> zfp(zf);
    TEST_CHECK_LESS_THAN(zfp.max_rel_diff(a), DT_(1e-7));
#endif
  }
};
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseVectorMaxRelDiffTest
  : public UnitTest
{
public:
  explicit SparseVectorMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseVectorMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    const Index size = 100;
    SparseVectorFactory<DT_, IT_> a_fac(size);

    for(Index k = 0; k < 17; ++k)
      a_fac.add(((k+23)*131071) % size, DT_(k+1));

    SparseVector<DT_, IT_> a(a_fac.make_sv());
    SparseVector<DT_, IT_> b(a.clone());

    TEST_CHECK_LESS_THAN(a.max_rel_diff(b), eps*eps);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(a), eps*eps);

    b.elements_view_rw()[7] += DT_(0.314);
    TEST_CHECK_LESS_THAN(eps, b.max_rel_diff(a));
  }
};
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorMaxRelDiffTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseVectorSameLayoutTest
  : public UnitTest
{
public:
  explicit SparseVectorSameLayoutTest(PreferredBackend backend)
    : UnitTest("SparseVectorSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const Index size = 100;
    SparseVectorFactory<DT_, IT_> a_fac(size), z_fac(size+2);

    for(Index k = 0; k < 17; ++k)
    {
      a_fac.add(((k+23)*131071) % size, DT_(k+1));
      z_fac.add(((k+23)*131071) % size, DT_(k+1));
    }

    SparseVector<DT_, IT_> a(a_fac.make_sv());

    // different sizes
    SparseVector<DT_, IT_> z(z_fac.make_sv());
    TEST_CHECK(!a.same_layout(z));

    // weak copy
    auto b = a.clone(CloneMode::Weak);
    TEST_CHECK(a.same_layout(b));

    // shallow copy
    auto c = a.clone(CloneMode::Shallow);
    TEST_CHECK(a.same_layout(c));

    // different values at same position
    SparseVector<DT_, IT_> d(a.clone());
    d.elements_view_rw()[7] += DT_(0.5);
    TEST_CHECK(a.same_layout(d));
  }
};
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SparseVectorSameLayoutTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
