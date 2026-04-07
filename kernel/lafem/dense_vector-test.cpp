// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/type_traits.hpp>

#include <list>
#include <sstream>
#include <cstdio>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the dense vector class.
 *
 * \test test description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class DenseVectorTest
  : public UnitTest
{
public:
  explicit DenseVectorTest(PreferredBackend backend)
    : UnitTest("DenseVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    DenseVector<DT_, IT_> zero;
    TEST_CHECK(zero.empty());
    TEST_CHECK_EQUAL(zero.bytes(), 0);

    DenseVector<DT_, IT_> a(16, DT_(7)); //use multiple of 4 to circumnavigate memory padding
    TEST_CHECK(!a.empty());
    TEST_CHECK_EQUAL(a.size(), 16);
    TEST_CHECK_EQUAL(a.bytes(), 16 * sizeof(DT_));
    TEST_CHECK_EQUAL_WITHIN_EPS(a.elements_view_r()(5), DT_(7), tol);

    DenseVector<DT_, IT_> b(16, DT_(5));
    b.elements_view_rw()[7] = DT_(42);
    TEST_CHECK_EQUAL(b.elements_view_r()(7), DT_(42));
    TEST_CHECK_EQUAL(b.elements_view_r()(3), DT_(5));

    DenseVector<DT_, IT_> b_r(b, 5, 3);
    TEST_CHECK_EQUAL(b_r.elements_view_r()(0), b.elements_view_r()(0+3));
    TEST_CHECK_EQUAL(b_r.elements_view_r()(4), b.elements_view_r()(4+3));
    auto b_rc = b_r.clone();
    TEST_CHECK_EQUAL(b_rc.elements_view_r()(0), b.elements_view_r()(0+3));
    TEST_CHECK_EQUAL(b_rc.elements_view_r()(4), b.elements_view_r()(4+3));

    DenseVector<DT_, IT_> c(b.clone());
    TEST_CHECK_EQUAL(c.size(), b.size());
    {
      Memory::TypedView<DT_> cv = c.elements_view_r();
      Memory::TypedView<DT_> bv = b.elements_view_r();
      for (Index i(0) ; i < c.size() ; ++i)
        TEST_CHECK_EQUAL(cv(i), bv(i));
    }
    TEST_CHECK_LESS_THAN(c.max_rel_diff(b), tol);
    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c.elements_view_r()(7), b.elements_view_r()(7));
    TEST_CHECK_LESS_THAN(c.max_rel_diff(b), tol);
    DenseVector<float, unsigned int> d;
    d.convert(c);
    DenseVector<float, unsigned int> e;
    e.convert(b);
    TEST_CHECK_EQUAL(e.size(), d.size());
    TEST_CHECK_EQUAL(e.elements_view_r()(7), d.elements_view_r()(7));
    TEST_CHECK_LESS_THAN(e.max_rel_diff(d), float(tol));
    e.clone(a);
    {
      const Memory::TypedView<DT_> av = a.elements_view_r();
      const Memory::TypedView<float> ev = e.elements_view_r();
      for (Index i(0) ; i < a.size() ; ++i)
        TEST_CHECK_EQUAL(DT_(ev(i)), av(i));
    }

    b.clone(a);
    TEST_CHECK_NOT_EQUAL(b.elements_arbiter(), a.elements_arbiter());
    c.convert(a);
    TEST_CHECK_EQUAL(c.elements_arbiter(), a.elements_arbiter());
    TEST_CHECK_LESS_THAN(b.max_rel_diff(c), tol);

    DenseVector<DT_, IT_> g(b.size(), b.elements_arbiter().attach());
    TEST_CHECK_LESS_THAN(g.max_rel_diff(b), tol);
    TEST_CHECK_EQUAL(g.elements_arbiter(), b.elements_arbiter());

    DenseVector<DT_, IT_> ap(a.clone());
    Adjacency::Permutation prm_nil;
    ap.permute(prm_nil);
    Adjacency::Permutation prm_rnd(a.size(), rng);
    ap.permute(prm_rnd);
    prm_rnd = prm_rnd.inverse();
    ap.permute(prm_rnd);
    TEST_CHECK_LESS_THAN(ap.max_rel_diff(a), tol);

    // random constructor check
    DT_ rnd_range[2];
    IT_ rnd_size = 1234;
    rnd_range[0] = DT_(-10);
    rnd_range[1] = DT_(+10);
    DenseVector<DT_, IT_> rnd_vec(rng, rnd_size, rnd_range[0], rnd_range[1]);
    TEST_CHECK_EQUAL(rnd_vec.size(), rnd_size);
    DT_ rnd_max = rnd_vec.max_abs_element();
    TEST_CHECK_IN_RANGE(rnd_max, rnd_range[0], rnd_range[1]);
    rnd_vec.scale(rnd_vec, DT_(-1));
    DT_ rnd_min = -rnd_vec.max_abs_element();
    TEST_CHECK_IN_RANGE(rnd_min, rnd_range[0], rnd_range[1]);

    // new clone testing
    auto clone1 = a.clone(CloneMode::Deep);
    TEST_CHECK_LESS_THAN(clone1.max_rel_diff(a), tol);
    clone1.elements_view_rw()[7] = DT_(132);
    TEST_CHECK_LESS_THAN(tol, clone1.max_rel_diff(a));
    TEST_CHECK_NOT_EQUAL(clone1.elements_arbiter(), a.elements_arbiter());
    DenseVector<DT_, IT_> clone2 = clone1.clone(CloneMode::Layout);
    clone2.format(DT_(4713));
    TEST_CHECK_NOT_EQUAL(clone2.elements_view_r()(7), clone1.elements_view_r()(7));
    TEST_CHECK_NOT_EQUAL(clone2.elements_arbiter(), clone1.elements_arbiter());
    DenseVector<DT_, IT_> clone3 = clone1.clone(CloneMode::Weak);
    TEST_CHECK_LESS_THAN(clone3.max_rel_diff(clone1), tol);
    clone3.elements_view_rw()[7] = DT_(133);
    TEST_CHECK_LESS_THAN(tol, clone3.max_rel_diff(clone1));
    TEST_CHECK_NOT_EQUAL(clone3.elements_arbiter(), clone1.elements_arbiter());
    DenseVector<DT_, IT_> clone4 = clone1.clone(CloneMode::Shallow);
    TEST_CHECK_LESS_THAN(clone4.max_rel_diff(clone1), tol);
    clone4.elements_view_rw()[7] = DT_(134);
    TEST_CHECK_LESS_THAN(clone4.max_rel_diff(clone1), tol);
    TEST_CHECK_EQUAL(clone4.elements_arbiter(), clone1.elements_arbiter());
    auto clone5 = a.clone(CloneMode::Allocate);
    TEST_CHECK_NOT_EQUAL(clone5.elements_arbiter(), a.elements_arbiter());
    TEST_CHECK_EQUAL(clone5.size(), a.size());
  }
};
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMaxRelDiffTest
  : public UnitTest
{
public:
  explicit DenseVectorMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("DenseVectorMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(12); size < Index(1000); size *= 2)
    {
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);

      const Index off0 = (3*size) / 8;
      const Index off1 = (1*size) / 8;
      const Index off2 = (6*size) / 8;

      // a = i, b = i
      {
        Memory::TypedView<DT_> va = a.elements_view_w();
        Memory::TypedView<DT_> vb = b.elements_view_w();
        for (Index i(0); i < size; ++i)
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
      a.elements_view_rw()[off0] += delta_a0;
      b.elements_view_rw()[off0] -= delta_b0;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref0, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref0, tol);

      const DT_ delta1 = DT_(0.17);
      const DT_ ref1 = delta1 / (DT_(off0 - off1)*DT_(0.246) + delta1 + DT_(1));
      a.elements_view_rw()[off1] -= delta1;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref1, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref1, tol);

      const DT_ delta2 = DT_(0.73);
      const DT_ ref2 = delta2 / (DT_(off2 - off0)*DT_(0.246) + delta2 + DT_(1));
      b.elements_view_rw()[off2] += delta2;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref2, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref2, tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
// no MKL implementation available
//#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, float, std::uint64_t, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, double, std::uint64_t, PreferredBackend::mkl);
//#endif
// no CUDA implementation available
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::generic);
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::cuda);
//#endif
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, double, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxRelDiffTest, double, std::uint64_t, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorSerializeTest
  : public UnitTest
{
public:
  explicit DenseVectorSerializeTest(PreferredBackend backend)
    : UnitTest("DenseVectorSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
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

    {
      std::stringstream ioss;
      k.write_out(FileMode::fm_mtx, ioss);
      DenseVector<DT_, IT_> l(FileMode::fm_mtx, ioss);
      TEST_CHECK_LESS_THAN(l.max_rel_diff(k), tol);
    }

    {
      std::stringstream ioss;
      k.write_out(FileMode::fm_exp, ioss);
      DenseVector<DT_, IT_> m(FileMode::fm_exp, ioss);
      TEST_CHECK_LESS_THAN(m.max_rel_diff(k), tol);
    }

    {
      BinaryStream bs;
      k.write_out(FileMode::fm_dv, bs);
      bs.seekg(0);
      DenseVector<DT_, IT_> n(FileMode::fm_dv, bs);
      TEST_CHECK_LESS_THAN(n.max_rel_diff(k), tol);
    }

    {
      auto op = k.serialize(LAFEM::SerialConfig(false,false));
      DenseVector<DT_, IT_> o(op);
      TEST_CHECK_LESS_THAN(o.max_rel_diff(k), tol);
    }

    {
#ifdef FEAT_HAVE_ZLIB
      auto zb = k.serialize(LAFEM::SerialConfig(true,false));
      DenseVector<DT_, IT_> z(zb);
      TEST_CHECK_LESS_THAN(z.max_rel_diff(k), tol);
#endif
    }

    {
#ifdef FEAT_HAVE_ZFP
      auto zp = k.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-5)));
      DenseVector<DT_, IT_> y(zp);
      TEST_CHECK_LESS_THAN(y.max_rel_diff(k), DT_(2e-5));
#endif
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, double, std::uint64_t, PreferredBackend::generic);
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorAxpyTest
  : public UnitTest
{
public:
  explicit DenseVectorAxpyTest(PreferredBackend backend)
    : UnitTest("DenseVectorAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    const DT_ s(DT_(0.4711));
    Index max_size(1000);

#ifdef FEAT_HAVE_HALFMATH
    if (typeid(DT_) == typeid(Half))
      max_size = 129;
#endif

    for (Index size(1) ; size < max_size ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);
      DenseVector<DT_, IT_> ref_a(size);
      DenseVector<DT_, IT_> ref_b(size);

      {
        Memory::TypedView<DT_> va = a.elements_view_w();
        Memory::TypedView<DT_> vb = b.elements_view_w();
        Memory::TypedView<DT_> vra = ref_a.elements_view_w();
        Memory::TypedView<DT_> vrb = ref_b.elements_view_w();
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(i % 100) * DT_(1.234);
          vb[i] = DT_(2) - DT_(i % 42);
          vra[i] = s * va(i) + vb(i);
          vrb[i] = s * vb(i) + vb(i);
        }
      }

      // r != x
      a.scale(a, s);
      a.axpy(b);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref_a), tol);

      // r == x
      b.axpy(b, s);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref_b), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorAxpyTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorDotTest
  : public UnitTest
{
public:
  explicit DenseVectorDotTest(PreferredBackend backend)
    : UnitTest("DenseVectorDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(1) ; size < Index(1000) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);
      const DT_ den(DT_(1) / DT_(size));
      {
        Memory::TypedView<DT_> va = a.elements_view_w();
        Memory::TypedView<DT_> vb = b.elements_view_w();
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(i+1) * den;    // a[i] = (i+1) / n
          vb[i] = DT_(1) / DT_(i+1); // b[i] = 1 / (i+1)
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

SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorDotTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorTripleDotTest
  : public UnitTest
{
public:
  explicit DenseVectorTripleDotTest(PreferredBackend backend)
    : UnitTest("DenseVectorTripleDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(1) ; size < Index(1000) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);
      DenseVector<DT_, IT_> c(size);

      const DT_ den( DT_(1) / Math::sqrt(DT_(size)) );
      {
        Memory::TypedView<DT_> va = a.elements_view_w();
        Memory::TypedView<DT_> vb = b.elements_view_w();
        Memory::TypedView<DT_> vc = c.elements_view_w();
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(i+1) * den;    // a[i] = (i+1) / n
          vb[i] = DT_(1) / DT_(i+1); // b[i] = 1 / (i+1)
          vc[i] = den;
        }
      }

      // a^T diag(c) b = 1
      DT_ ref(DT_(1));

      DT_ res  = a.triple_dot(b,c);
      TEST_CHECK_RELATIVE(res, ref, tol);
      res  = a.triple_dot(c,b);
      TEST_CHECK_RELATIVE(res, ref, tol);

      res = b.triple_dot(a,c);
      TEST_CHECK_RELATIVE(res, ref, tol);
      res = b.triple_dot(c,a);
      TEST_CHECK_RELATIVE(res, ref, tol);

      res = c.triple_dot(a,b);
      TEST_CHECK_RELATIVE(res, ref, tol);
      res = c.triple_dot(b,a);
      TEST_CHECK_RELATIVE(res, ref, tol);

    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorTripleDotTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorComponentProductTest
  : public UnitTest
{
public:
  explicit DenseVectorComponentProductTest(PreferredBackend backend)
    : UnitTest("DenseVectorComponentProductTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(1) ; size < Index(1000) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);
      DenseVector<DT_, IT_> ref1(size);
      DenseVector<DT_, IT_> ref2(size);
      {
        Memory::TypedView<DT_> va = a.elements_view_w();
        Memory::TypedView<DT_> vb = b.elements_view_w();
        Memory::TypedView<DT_> vr1 = ref1.elements_view_w();
        Memory::TypedView<DT_> vr2 = ref2.elements_view_w();
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(i)/DT_(100) * DT_(1.234);
          vb[i] = DT_(size*2 - i);
          vr1[i] = va(i) * vb(i);
          vr2[i] = va(i) * va(i);
        }
      }

      DenseVector<DT_, IT_> c(size);
      c.component_product(a, b);
      TEST_CHECK_LESS_THAN(c.max_rel_diff(ref1), tol);

      b.component_product(a, b);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref1), tol);

      a.component_product(a, a);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref2), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentProductTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorComponentInvertTest
  : public UnitTest
{
public:
  explicit DenseVectorComponentInvertTest(PreferredBackend backend)
    : UnitTest("DenseVectorComponentInvertTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    const DT_ alpha(Math::pi<DT_>());

    for (Index size(1) ; size < Index(1000) ; size*=2)
    {
      // create a vector
      DenseVector<DT_, IT_>  vec(size);
      DenseVector<DT_, IT_>  ref2(size, alpha);
      DenseVector<DT_, IT_>  ref3(size);
      {
        Memory::TypedView<DT_> v(vec.elements_view_w());
        Memory::TypedView<DT_> r3(ref3.elements_view_w());
        for (Index i(0); i < size; ++i)
        {
          v[i] = DT_(7.63) * DT_(i % 3 + 1) - DT_(9.3);
          r3[i] = DT_(1) / v[i];
        }
      }

      DenseVector<DT_, IT_>  vec2(vec.clone());
      vec2.component_invert(vec2, alpha);
      vec2.component_product(vec2, vec);
      TEST_CHECK_LESS_THAN(vec2.max_rel_diff(ref2), tol);

      DenseVector<DT_, IT_>  vec3(size);
      vec3.component_invert(vec);
      TEST_CHECK_LESS_THAN(vec3.max_rel_diff(ref3), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, float, Index, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, double, Index, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, __float128, Index, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, float, Index, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, double, Index, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, Half, Index, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, Half, Index, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, float, Index, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorComponentInvertTest, double, Index, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorScaleTest
  : public UnitTest
{
public:
  explicit DenseVectorScaleTest(PreferredBackend backend)
    : UnitTest("DenseVectorScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));
    for (Index size(1) ; size < Index(1000) ; size*=2)
    {
      DT_ s(DT_(4.321));
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> ref(size);
      {
        Memory::TypedView<DT_> va = a.elements_view_w();
        Memory::TypedView<DT_> vr = ref.elements_view_w();
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(DT_(i) * DT_(1.234));
          vr[i] = va(i) * s;
        }
      }

      DenseVector<DT_, IT_> b(size);
      b.scale(a, s);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref), tol);

      a.scale(a, s);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorScaleTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorNorm2Test
  : public UnitTest
{
public:
  explicit DenseVectorNorm2Test(PreferredBackend backend)
    : UnitTest("DenseVectorNorm2Test", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(1) ; size < Index(1000) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      {
        Memory::TypedView<DT_> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
          va[i] = Math::pow(DT_(0.5), DT_(0.5) * DT_(i));
        }
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      const DT_ ref(Math::sqrt(DT_(2) - Math::pow(DT_(0.5), DT_(size-1))));

      DT_ c = a.norm2();
      TEST_CHECK_RELATIVE(c, ref, tol);

      c = a.norm2sqr();
      TEST_CHECK_RELATIVE(c, ref*ref, tol);
    }
  }
};
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorNorm2Test, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMaxAbsElementTest
  : public UnitTest
{
public:
  explicit DenseVectorMaxAbsElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorMaxAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(1) ; size < Index(1000) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      {
        Memory::TypedView<DT_> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(i) * (i%2 == 0 ? DT_(1) : DT_(-1));
        }
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ max = a.max_abs_element();

      TEST_CHECK_RELATIVE(max, DT_(size-1), tol);
    }
  }
};
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxAbsElementTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMinAbsElementTest
  : public UnitTest
{
public:
  explicit DenseVectorMinAbsElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorMinAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    for (Index size(1) ; size < Index(1000) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      {
        Memory::TypedView<DT_> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(i+1) * (i%2 == 0 ? DT_(1) : DT_(-1));
        }
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ min = a.min_abs_element();

      TEST_CHECK_RELATIVE(min, DT_(1), tol);
    }
  }
};
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinAbsElementTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMaxElementTest
  : public UnitTest
{
public:
  explicit DenseVectorMaxElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorMaxElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    for (Index size(5) ; size < Index(1000) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      {
        Memory::TypedView<DT_> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(i);
        }
        va[3] = DT_(-5);
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ max = a.max_element();

      TEST_CHECK_RELATIVE(max, DT_(size-1), tol);
    }
  }
};
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, Half, std::uint64_t, PreferredBackend::generic);
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, Half, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, Half, std::uint64_t, PreferredBackend::cuda);
//#endif
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, double, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMaxElementTest, double, std::uint64_t, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMinElementTest
  : public UnitTest
{
public:
  explicit DenseVectorMinElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorMinElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    for (Index size(5) ; size < Index(1000) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      {
        Memory::TypedView<DT_> va(a.elements_view_w());
        for (Index i(0) ; i < size ; ++i)
        {
          va[i] = DT_(i) - DT_(3);
        }
      }

      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ min = a.min_element();

      TEST_CHECK_RELATIVE(min, DT_(-3.0), tol);
    }
  }
};
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, Half, std::uint64_t, PreferredBackend::generic);
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, Half, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, Half, std::uint64_t, PreferredBackend::cuda);
//#endif
#endif
//#ifdef FEAT_HAVE_CUDA
//SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, float, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, double, std::uint32_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, float, std::uint64_t, PreferredBackend::cuda);
//SPAWN_UNIT_TEST_2T_P(DenseVectorMinElementTest, double, std::uint64_t, PreferredBackend::cuda);
//#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorSameLayoutTest
  : public UnitTest
{
public:
  explicit DenseVectorSameLayoutTest(PreferredBackend backend)
    : UnitTest("DenseVectorSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    for (Index size(2); size < Index(1000); size *= 2)
    {
      const Index diff_elem = size / 2;

      DenseVector<DT_, IT_> a(size);
      a.format();

      // a = i
      {
        Memory::TypedView<DT_> va(a.elements_view_w());
        for (Index i(0); i < size; ++i)
        {
          va[i] = DT_(i);
        }
      }

      // weak copy
      DenseVector<DT_, IT_> b = a.clone(CloneMode::Weak);
      TEST_CHECK(a.same_layout(b));

      // shallow copy
      DenseVector<DT_, IT_> c = a.clone(CloneMode::Shallow);
      TEST_CHECK(a.same_layout(c));

      // change one element
      c.elements_view_rw()[diff_elem] = DT_(0.5);
      TEST_CHECK(a.same_layout(c));

      // different sizes
      DenseVector<DT_, IT_> d(size + 2, DT_(10));
      DenseVector<DT_, IT_> e(size, DT_(10));
      TEST_CHECK(!d.same_layout(e));

      // one different element
      DenseVector<DT_, IT_> f(size);
      f.format();
      f.elements_view_rw()[diff_elem] = DT_(10);
      DenseVector<DT_, IT_> g(size);
      g.format();
      g.elements_view_rw()[diff_elem] = DT_(0.5);
      TEST_CHECK(f.same_layout(g));
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseVectorSameLayoutTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
