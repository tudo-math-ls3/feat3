// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
  DenseVectorTest(PreferredBackend backend)
    : UnitTest("DenseVectorTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorTest()
  {
  }


  virtual void run() const override
  {
    DenseVector<DT_, IT_> zero1;
    DenseVector<DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);
    zero2.convert(zero1);

    DenseVector<DT_, IT_> a(16, DT_(7)); //use multiple of 4 to circumanivate memory padding in MemoryPool
    TEST_CHECK_EQUAL(a.bytes(), 16 * sizeof(DT_) + 1 * sizeof(Index));

    TEST_CHECK_EQUAL(MemoryPool::allocated_memory(), a.bytes() - sizeof(Index));
    DenseVector<DT_, IT_> b(16, DT_(5));
    b(7, DT_(42));
    TEST_CHECK_EQUAL(b(7), DT_(42));
    TEST_CHECK_EQUAL(b(3), DT_(5));

    DenseVector<DT_, IT_> b_r(b, 5, 3);
    TEST_CHECK_EQUAL(b_r(0), b(0+3));
    TEST_CHECK_EQUAL(b_r(4), b(4+3));
    auto b_rc = b_r.clone();
    TEST_CHECK_EQUAL(b_rc(0), b(0+3));
    TEST_CHECK_EQUAL(b_rc(4), b(4+3));

    DenseVector<DT_, IT_> c(b.clone());
    TEST_CHECK_EQUAL(c.size(), b.size());
    for (Index i(0) ; i < c.size() ; ++i)
      TEST_CHECK_EQUAL(c(i), b(i));
    TEST_CHECK_EQUAL(c, b);
    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c(7), b(7));
    TEST_CHECK_EQUAL(c, b);
    DenseVector<float, unsigned int> d;
    d.convert(c);
    DenseVector<float, unsigned int> e;
    e.convert(b);
    TEST_CHECK_EQUAL(e.size(), d.size());
    TEST_CHECK_EQUAL(e(7), d(7));
    TEST_CHECK_EQUAL(e, d);
    e.clone(a);
    for (Index i(0) ; i < a.size() ; ++i)
      TEST_CHECK_EQUAL(DT_(e(i)), a(i));

    b.clone(a);
    TEST_CHECK_NOT_EQUAL((void*)b.elements(), (void*)a.elements());
    c.convert(a);
    TEST_CHECK_EQUAL((void*)c.elements(), (void*)a.elements());
    TEST_CHECK_EQUAL(b, c);

    DenseVector<DT_, IT_> g(b.size(), b.elements());
    TEST_CHECK_EQUAL(g, b);
    TEST_CHECK_EQUAL((void*)g.elements(), (void*)b.elements());

    DenseVector<DT_, IT_> ap(a.clone());
    Adjacency::Permutation prm_nil;
    ap.permute(prm_nil);
    Random::SeedType seed(Random::SeedType(time(nullptr)));
    std::cout << "seed: " << seed << std::endl;
    Random rng(seed);
    Adjacency::Permutation prm_rnd(a.size(), rng);
    ap.permute(prm_rnd);
    prm_rnd = prm_rnd.inverse();
    ap.permute(prm_rnd);
    TEST_CHECK_EQUAL(ap, a);

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
    TEST_CHECK_EQUAL(clone1, a);
    clone1(7, DT_(132));
    TEST_CHECK_NOT_EQUAL(clone1, a);
    TEST_CHECK_NOT_EQUAL((void*)clone1.elements(), (void*)a.elements());
    DenseVector<DT_, IT_> clone2 = clone1.clone(CloneMode::Layout);
    MemoryPool::set_memory(clone2.elements(), DT_(4713), clone2.size());
    TEST_CHECK_NOT_EQUAL(clone2(7), clone1(7));
    TEST_CHECK_NOT_EQUAL((void*)clone2.elements(), (void*)clone1.elements());
    DenseVector<DT_, IT_> clone3 = clone1.clone(CloneMode::Weak);
    TEST_CHECK_EQUAL(clone3, clone1);
    clone3(7, DT_(133));
    TEST_CHECK_NOT_EQUAL(clone3, clone1);
    TEST_CHECK_NOT_EQUAL((void*)clone3.elements(), (void*)clone1.elements());
    DenseVector<DT_, IT_> clone4 = clone1.clone(CloneMode::Shallow);
    TEST_CHECK_EQUAL(clone4, clone1);
    clone4(7, DT_(134));
    TEST_CHECK_EQUAL(clone4, clone1);
    TEST_CHECK_EQUAL((void*)clone4.elements(), (void*)clone1.elements());
    auto clone5 = a.clone(CloneMode::Allocate);
    TEST_CHECK_NOT_EQUAL((void*)clone5.elements(), (void*)a.elements());
    TEST_CHECK_EQUAL(clone5.size(), a.size());
  }
};
DenseVectorTest <float, std::uint32_t> dv_test_float_uint32(PreferredBackend::generic);
DenseVectorTest <double, std::uint32_t> dv_test_double_uint32(PreferredBackend::generic);
DenseVectorTest <float, std::uint64_t> dv_test_float_uint64(PreferredBackend::generic);
DenseVectorTest <double, std::uint64_t> dv_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorTest <__float128, std::uint32_t> dv_test_float128_uint32(PreferredBackend::generic);
DenseVectorTest <__float128, std::uint64_t> dv_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorTest <float, std::uint64_t> mkl_dv_test_float_uint64(PreferredBackend::mkl);
DenseVectorTest <double, std::uint64_t> mkl_dv_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorTest <Half, std::uint32_t> dv_test_half_uint32(PreferredBackend::generic);
DenseVectorTest <Half, std::uint64_t> dv_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorTest <Half, std::uint32_t> cuda_dv_test_half_uint32(PreferredBackend::cuda);
DenseVectorTest <Half, std::uint64_t> cuda_dv_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorTest <float, std::uint32_t> cuda_dv_test_float_uint32(PreferredBackend::cuda);
DenseVectorTest <double, std::uint32_t> cuda_dv_test_double_uint32(PreferredBackend::cuda);
DenseVectorTest <float, std::uint64_t> cuda_dv_test_float_uint64(PreferredBackend::cuda);
DenseVectorTest <double, std::uint64_t> cuda_dv_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorSerializeTest
  : public UnitTest
{
public:
  DenseVectorSerializeTest(PreferredBackend backend)
    : UnitTest("DenseVectorSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorSerializeTest()
  {
  }


  virtual void run() const override
  {
    Index io_vector_size = 1234;

    DenseVector<DT_, IT_> k(io_vector_size);
    for (Index i(0) ; i < k.size() ; ++i)
      k(i, DT_(i) / DT_(12));

    {
      std::stringstream mts;
      k.write_out(FileMode::fm_mtx, mts);
      DenseVector<DT_, IT_> l(FileMode::fm_mtx, mts);
      for (Index i(0) ; i < k.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(l(i), k(i), DT_(1e-4));
    }

    {
      std::stringstream ts;
      k.write_out(FileMode::fm_exp, ts);
      DenseVector<DT_, IT_> m(FileMode::fm_exp, ts);
      for (Index i(0) ; i < k.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(m(i), k(i), DT_(1e-4));
    }

    {
      BinaryStream bs;
      k.write_out(FileMode::fm_dv, bs);
      bs.seekg(0);
      DenseVector<DT_, IT_> n(FileMode::fm_dv, bs);
      for (Index i(0) ; i < k.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(n(i), k(i), DT_(1e-5));
    }

    {
      auto op = k.serialize(LAFEM::SerialConfig(false,false));
      DenseVector<DT_, IT_> o(op);
      for (Index i(0) ; i < k.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(o(i), k(i), DT_(1e-5));
#ifdef FEAT_HAVE_ZLIB
      auto zb = k.serialize(LAFEM::SerialConfig(true,false));
      DenseVector<DT_, IT_> zlib(zb);
      for (Index i(0) ; i < k.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(zlib(i), k(i), DT_(1e-5));
#endif
#ifdef FEAT_HAVE_ZFP
      auto zp = k.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-5)));
      DenseVector<DT_, IT_> zfp(zp);
      for (Index i(0) ; i < k.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(zfp(i), k(i), DT_(1e-5));
#endif
    }
  }
};

DenseVectorSerializeTest <float, std::uint32_t> dv_serialize_test_float_uint32(PreferredBackend::generic);
DenseVectorSerializeTest <double, std::uint32_t> dv_serialize_test_double_uint32(PreferredBackend::generic);
DenseVectorSerializeTest <float, std::uint64_t> dv_serialize_test_float_uint64(PreferredBackend::generic);
DenseVectorSerializeTest <double, std::uint64_t> dv_serialize_test_double_uint64(PreferredBackend::generic);
//#ifdef FEAT_HAVE_QUADMATH
//DenseVectorSerializeTest <__float128, std::uint32_t> dv_serialize_test_float128_uint32(PreferredBackend::generic);
//DenseVectorSerializeTest <__float128, std::uint64_t> dv_serialize_test_float128_uint64(PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_MKL
DenseVectorSerializeTest <float, std::uint64_t> mkl_dv_serialize_test_float_uint64(PreferredBackend::mkl);
DenseVectorSerializeTest <double, std::uint64_t> mkl_dv_serialize_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorSerializeTest <Half, std::uint32_t> dv_serialize_test_half_uint32(PreferredBackend::generic);
DenseVectorSerializeTest <Half, std::uint64_t> dv_serialize_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorSerializeTest <Half, std::uint32_t> cuda_dv_serialize_test_half_uint32(PreferredBackend::cuda);
DenseVectorSerializeTest <Half, std::uint64_t> cuda_dv_serialize_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorSerializeTest <float, std::uint32_t> cuda_dv_serialize_test_float_uint32(PreferredBackend::cuda);
DenseVectorSerializeTest <double, std::uint32_t> cuda_dv_serialize_test_double_uint32(PreferredBackend::cuda);
DenseVectorSerializeTest <float, std::uint64_t> cuda_dv_serialize_test_float_uint64(PreferredBackend::cuda);
DenseVectorSerializeTest <double, std::uint64_t> cuda_dv_serialize_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorAxpyTest
  : public UnitTest
{
public:
  DenseVectorAxpyTest(PreferredBackend backend)
    : UnitTest("DenseVectorAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorAxpyTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.4));

    DT_ s(DT_(47.11));
    Index max_size(1000);
#ifdef FEAT_HAVE_HALFMATH
    if (typeid(DT_) == typeid(Half))
      max_size = 129;
#endif

    for (Index size(1) ; size < max_size ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);
      DenseVector<DT_, IT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i % 100) * DT_(1.234));
        b(i, DT_(2) - DT_(i % 42));
        ref(i, s * a(i) + b(i));
      }

      DenseVector<DT_, IT_> c(size);
      c.axpy(a, b, s); //a != b != r
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(c(i), ref(i), DT_(eps));

      a.axpy(a, b, s); //r == a
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(a(i), ref(i), DT_(eps));

      for (Index i(0) ; i < size ; ++i)
      {
        ref(i, s * a(i) + b(i));
      }

      b.axpy(a, b, s); //r == b
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(b(i), ref(i), DT_(eps));
    }
  }
};
DenseVectorAxpyTest <float, std::uint32_t> dv_axpy_test_float_uint32(PreferredBackend::generic);
DenseVectorAxpyTest <double, std::uint32_t> dv_axpy_test_double_uint32(PreferredBackend::generic);
DenseVectorAxpyTest <float, std::uint64_t> dv_axpy_test_float_uint64(PreferredBackend::generic);
DenseVectorAxpyTest <double, std::uint64_t> dv_axpy_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorAxpyTest <__float128, std::uint32_t> dv_axpy_test_float128_uint32(PreferredBackend::generic);
DenseVectorAxpyTest <__float128, std::uint64_t> dv_axpy_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorAxpyTest <float, std::uint64_t> mkl_dv_axpy_test_float_uint64(PreferredBackend::mkl);
DenseVectorAxpyTest <double, std::uint64_t> mkl_dv_axpy_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorAxpyTest <Half, std::uint32_t> dv_axpy_test_half_uint32(PreferredBackend::generic);
DenseVectorAxpyTest <Half, std::uint64_t> dv_axpy_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorAxpyTest <Half, std::uint32_t> cuda_dv_axpy_test_half_uint32(PreferredBackend::cuda);
DenseVectorAxpyTest <Half, std::uint64_t> cuda_dv_axpy_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorAxpyTest <float, std::uint32_t> cuda_dv_axpy_test_float_uint32(PreferredBackend::cuda);
DenseVectorAxpyTest <double, std::uint32_t> cuda_dv_axpy_test_double_uint32(PreferredBackend::cuda);
DenseVectorAxpyTest <float, std::uint64_t> cuda_dv_axpy_test_float_uint64(PreferredBackend::cuda);
DenseVectorAxpyTest <double, std::uint64_t> cuda_dv_axpy_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorDotTest
  : public UnitTest
{
public:
  DenseVectorDotTest(PreferredBackend backend)
    : UnitTest("DenseVectorDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorDotTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);
      const DT_ den(DT_(1) / DT_(size));
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i+1) * den);    // a[i] = (i+1) / n
        b(i, DT_(1) / DT_(i+1)); // b[i] = 1 / (i+1)
      }

      // a*b = 1
      DT_ ref(DT_(1));
      DT_ c  = a.dot(b);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = b.dot(a);
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
      c = b.dot(b);
      ref = b.norm2();
      ref *= ref;
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);
    }
  }
};
DenseVectorDotTest <float, std::uint32_t> dv_dot_product_test_float_uint32(PreferredBackend::generic);
DenseVectorDotTest <double, std::uint32_t> dv_dot_product_test_double_uint32(PreferredBackend::generic);
DenseVectorDotTest <float, std::uint64_t> dv_dot_product_test_float_uint64(PreferredBackend::generic);
DenseVectorDotTest <double, std::uint64_t> dv_dot_product_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorDotTest <__float128, std::uint32_t> dv_dot_product_test_float128_uint32(PreferredBackend::generic);
DenseVectorDotTest <__float128, std::uint64_t> dv_dot_product_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorDotTest <float, std::uint64_t> mkl_dv_dot_product_test_float_uint64(PreferredBackend::mkl);
DenseVectorDotTest <double, std::uint64_t> mkl_dv_dot_product_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorDotTest <Half, std::uint32_t> dv_dot_product_test_half_uint32(PreferredBackend::generic);
DenseVectorDotTest <Half, std::uint64_t> dv_dot_product_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorDotTest <Half, std::uint32_t> cuda_dv_dot_product_test_half_uint32(PreferredBackend::cuda);
DenseVectorDotTest <Half, std::uint64_t> cuda_dv_dot_product_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorDotTest <float, std::uint32_t> cuda_dv_dot_product_test_float_uint32(PreferredBackend::cuda);
DenseVectorDotTest <double, std::uint32_t> cuda_dv_dot_product_test_double_uint32(PreferredBackend::cuda);
DenseVectorDotTest <float, std::uint64_t> cuda_dv_dot_product_test_float_uint64(PreferredBackend::cuda);
DenseVectorDotTest <double, std::uint64_t> cuda_dv_dot_product_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorTripleDotTest
  : public UnitTest
{
public:
  DenseVectorTripleDotTest(PreferredBackend backend)
    : UnitTest("DenseVectorTripleDotTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorTripleDotTest()
  {
  }


  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.7));
    if (Backend::get_preferred_backend() == PreferredBackend::cuda)
      eps = Math::pow(Math::eps<DT_>(), DT_(0.4));

    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);
      DenseVector<DT_, IT_> c(size);

      const DT_ den( DT_(1) / Math::sqrt(DT_(size)) );

      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i+1) * den);    // a[i] = (i+1) / n
        b(i, DT_(1) / DT_(i+1)); // b[i] = 1 / (i+1)
        c(i, den);
      }

      // a^T diag(c) b = 1
      DT_ ref(DT_(1));

      DT_ res  = a.triple_dot(b,c);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);
      res  = a.triple_dot(c,b);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);

      res = b.triple_dot(a,c);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);
      res = b.triple_dot(c,a);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);

      res = c.triple_dot(a,b);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);
      res = c.triple_dot(b,a);
      TEST_CHECK_EQUAL_WITHIN_EPS(res, ref, eps);

    }
  }
};
DenseVectorTripleDotTest <float, std::uint32_t> dv_triple_dot_product_test_float_uint32(PreferredBackend::generic);
DenseVectorTripleDotTest <double, std::uint32_t> dv_triple_dot_product_test_double_uint32(PreferredBackend::generic);
DenseVectorTripleDotTest <float, std::uint64_t> dv_triple_dot_product_test_float_uint64(PreferredBackend::generic);
DenseVectorTripleDotTest <double, std::uint64_t> dv_triple_dot_product_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorTripleDotTest <__float128, std::uint32_t> dv_triple_dot_product_test_float128_uint32(PreferredBackend::generic);
DenseVectorTripleDotTest <__float128, std::uint64_t> dv_triple_dot_product_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorTripleDotTest <float, std::uint64_t> mkl_dv_triple_dot_product_test_float_uint64(PreferredBackend::mkl);
DenseVectorTripleDotTest <double, std::uint64_t> mkl_dv_triple_dot_product_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorTripleDotTest <Half, std::uint32_t> dv_triple_dot_product_test_half_uint32(PreferredBackend::generic);
DenseVectorTripleDotTest <Half, std::uint64_t> dv_triple_dot_product_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorTripleDotTest <Half, std::uint32_t> cuda_dv_triple_dot_product_test_half_uint32(PreferredBackend::cuda);
DenseVectorTripleDotTest <Half, std::uint64_t> cuda_dv_triple_dot_product_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorTripleDotTest <float, std::uint32_t> cuda_dv_triple_dot_product_test_float_uint32(PreferredBackend::cuda);
DenseVectorTripleDotTest <double, std::uint32_t> cuda_dv_triple_dot_product_test_double_uint32(PreferredBackend::cuda);
DenseVectorTripleDotTest <float, std::uint64_t> cuda_dv_triple_dot_product_test_float_uint64(PreferredBackend::cuda);
DenseVectorTripleDotTest <double, std::uint64_t> cuda_dv_triple_dot_product_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorComponentProductTest
  : public UnitTest
{
public:
  DenseVectorComponentProductTest(PreferredBackend backend)
    : UnitTest("DenseVectorComponentProductTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorComponentProductTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.7));
      if (Backend::get_preferred_backend() == PreferredBackend::cuda)
        eps = Math::pow(Math::eps<DT_>(), DT_(0.2));

      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> b(size);
      DenseVector<DT_, IT_> ref(size);
      DenseVector<DT_, IT_> ref2(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(DT_(i)/DT_(100) * DT_(1.234)));
        b(i, DT_(size*2 - i));
        ref(i, a(i) * b(i));
        ref2(i, a(i) * a(i));
      }

      DenseVector<DT_, IT_> c(size);
      c.component_product(a, b);
      for (Index i(0); i < c.template size<Perspective::pod>(); ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(c.template elements<Perspective::pod>()[i], ref.template elements<Perspective::pod>()[i], eps);
      }
      //TEST_CHECK_EQUAL(c, ref);

      b.component_product(a, b);
      for (Index i(0); i < b.template size<Perspective::pod>(); ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(b.template elements<Perspective::pod>()[i], ref.template elements<Perspective::pod>()[i], eps);
      }
      //TEST_CHECK_EQUAL(b, ref);

      a.component_product(a, a);
      for (Index i(0); i < a.template size<Perspective::pod>(); ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(a.template elements<Perspective::pod>()[i], ref2.template elements<Perspective::pod>()[i], eps);
      }
      //TEST_CHECK_EQUAL(a, ref2);
    }
  }
};
DenseVectorComponentProductTest <float, std::uint32_t> dv_component_product_test_float_uint32(PreferredBackend::generic);
DenseVectorComponentProductTest <double, std::uint32_t> dv_component_product_test_double_uint32(PreferredBackend::generic);
DenseVectorComponentProductTest <float, std::uint64_t> dv_component_product_test_float_uint64(PreferredBackend::generic);
DenseVectorComponentProductTest <double, std::uint64_t> dv_component_product_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorComponentProductTest <__float128, std::uint32_t> dv_component_product_test_float128_uint32(PreferredBackend::generic);
DenseVectorComponentProductTest <__float128, std::uint64_t> dv_component_product_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorComponentProductTest <float, std::uint64_t> mkl_dv_component_product_test_float_uint64(PreferredBackend::mkl);
DenseVectorComponentProductTest <double, std::uint64_t> mkl_dv_component_product_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorComponentProductTest <Half, std::uint32_t> dv_component_product_test_half_uint32(PreferredBackend::generic);
DenseVectorComponentProductTest <Half, std::uint64_t> dv_component_product_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorComponentProductTest <Half, std::uint32_t> cuda_dv_component_product_test_half_uint32(PreferredBackend::cuda);
DenseVectorComponentProductTest <Half, std::uint64_t> cuda_dv_component_product_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorComponentProductTest <float, std::uint32_t> cuda_dv_component_product_test_float_uint32(PreferredBackend::cuda);
DenseVectorComponentProductTest <double, std::uint32_t> cuda_dv_component_product_test_double_uint32(PreferredBackend::cuda);
DenseVectorComponentProductTest <float, std::uint64_t> cuda_dv_component_product_test_float_uint64(PreferredBackend::cuda);
DenseVectorComponentProductTest <double, std::uint64_t> cuda_dv_component_product_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorScaleTest
  : public UnitTest
{
public:
  DenseVectorScaleTest(PreferredBackend backend)
    : UnitTest("DenseVectorScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorScaleTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DT_ s(DT_(4.321));
      DenseVector<DT_, IT_> a(size);
      DenseVector<DT_, IT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(DT_(i) * DT_(1.234)));
        ref(i, a(i) * s);
      }

      DenseVector<DT_, IT_> b(size);
      b.scale(a, s);
      TEST_CHECK_EQUAL(b, ref);

      a.scale(a, s);
      TEST_CHECK_EQUAL(a, ref);
    }
  }
};
DenseVectorScaleTest <float, std::uint32_t> dv_scale_test_float_uint32(PreferredBackend::generic);
DenseVectorScaleTest <double, std::uint32_t> dv_scale_test_double_uint32(PreferredBackend::generic);
DenseVectorScaleTest <float, std::uint64_t> dv_scale_test_float_uint64(PreferredBackend::generic);
DenseVectorScaleTest <double, std::uint64_t> dv_scale_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorScaleTest <__float128, std::uint32_t> dv_scale_test_float128_uint32(PreferredBackend::generic);
DenseVectorScaleTest <__float128, std::uint64_t> dv_scale_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorScaleTest <float, std::uint64_t> mkl_dv_scale_product_test_float_uint64(PreferredBackend::mkl);
DenseVectorScaleTest <double, std::uint64_t> mkl_dv_scale_product_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorScaleTest <Half, std::uint32_t> dv_scale_product_test_half_uint32(PreferredBackend::generic);
DenseVectorScaleTest <Half, std::uint64_t> dv_scale_product_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorScaleTest <Half, std::uint32_t> cuda_dv_scale_product_test_half_uint32(PreferredBackend::cuda);
DenseVectorScaleTest <Half, std::uint64_t> cuda_dv_scale_product_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorScaleTest <float, std::uint32_t> cuda_dv_scale_test_float_uint32(PreferredBackend::cuda);
DenseVectorScaleTest <double, std::uint32_t> cuda_dv_scale_test_double_uint32(PreferredBackend::cuda);
DenseVectorScaleTest <float, std::uint64_t> cuda_dv_scale_test_float_uint64(PreferredBackend::cuda);
DenseVectorScaleTest <double, std::uint64_t> cuda_dv_scale_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorNorm2Test
  : public UnitTest
{
public:
  DenseVectorNorm2Test(PreferredBackend backend)
    : UnitTest("DenseVectorNorm2Test", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorNorm2Test()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        // a[i] = 1/sqrt(2^i) = (1/2)^(i/2)
        a(i, Math::pow(DT_(0.5), DT_(0.5) * DT_(i)));
      }

      // ||a||_2 = sqrt(2 - 2^{1-n})
      const DT_ ref(Math::sqrt(DT_(2) - Math::pow(DT_(0.5), DT_(size-1))));

      DT_ c = a.norm2();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, eps);

      c = a.norm2sqr();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref*ref, eps);
    }
  }
};
DenseVectorNorm2Test <float, std::uint32_t> dv_norm2_test_float_uint32(PreferredBackend::generic);
DenseVectorNorm2Test <double, std::uint32_t> dv_norm2_test_double_uint32(PreferredBackend::generic);
DenseVectorNorm2Test <float, std::uint64_t> dv_norm2_test_float_uint64(PreferredBackend::generic);
DenseVectorNorm2Test <double, std::uint64_t> dv_norm2_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorNorm2Test <__float128, std::uint32_t> dv_norm2_test_float128_uint32(PreferredBackend::generic);
DenseVectorNorm2Test <__float128, std::uint64_t> dv_norm2_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorNorm2Test <float, std::uint64_t> mkl_dv_norm2_test_float_uint64(PreferredBackend::mkl);
DenseVectorNorm2Test <double, std::uint64_t> mkl_dv_norm2_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorNorm2Test <Half, std::uint32_t> dv_norm2_test_half_uint32(PreferredBackend::generic);
DenseVectorNorm2Test <Half, std::uint64_t> dv_norm2_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorNorm2Test <Half, std::uint32_t> cuda_dv_norm2_test_half_uint32(PreferredBackend::cuda);
DenseVectorNorm2Test <Half, std::uint64_t> cuda_dv_norm2_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorNorm2Test <float, std::uint32_t> cuda_dv_norm2_test_float_uint32(PreferredBackend::cuda);
DenseVectorNorm2Test <double, std::uint32_t> cuda_dv_norm2_test_double_uint32(PreferredBackend::cuda);
DenseVectorNorm2Test <float, std::uint64_t> cuda_dv_norm2_test_float_uint64(PreferredBackend::cuda);
DenseVectorNorm2Test <double, std::uint64_t> cuda_dv_norm2_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorComponentInvertTest
  : public UnitTest
{
public:
  DenseVectorComponentInvertTest(PreferredBackend backend)
    : UnitTest("DenseVectorComponentInvertTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorComponentInvertTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.7));
    if (Backend::get_preferred_backend() == PreferredBackend::cuda)
      eps = Math::pow(Math::eps<DT_>(), DT_(0.4));

    const DT_ alpha(Math::pi<DT_>());

    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      // create a vector
      DenseVector<DT_, IT_>  vec(size);
      for (Index i(0); i < size; ++i)
      {
        vec(i, DT_(7.63) * DT_(i % 3 + 1) - DT_(9.3));
      }

      DenseVector<DT_, IT_>  vec2(vec.clone());
      vec2.component_invert(vec2, alpha);
      vec2.component_product(vec2, vec);
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(vec2(i), alpha, eps);
      }

      DenseVector<DT_, IT_>  vec3(size);
      vec3.component_invert(vec);
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(vec3(i), DT_(1.0) / vec(i), eps);
      }
    }
  }
};

DenseVectorComponentInvertTest<float, Index> dv_component_invert_test_float(PreferredBackend::generic);
DenseVectorComponentInvertTest<double, Index> dv_component_invert_test_double(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorComponentInvertTest<__float128, Index> dv_component_invert_test_float128(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorComponentInvertTest<float, Index> mkl_dv_component_invert_test_float(PreferredBackend::mkl);
DenseVectorComponentInvertTest<double, Index> mkl_dv_component_invert_test_double(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorComponentInvertTest<Half, Index> dv_component_invert_test_half(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorComponentInvertTest<Half, Index> cuda_dv_component_invert_test_half(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorComponentInvertTest<float, Index> cuda_dv_component_invert_test_float(PreferredBackend::cuda);
DenseVectorComponentInvertTest<double, Index> cuda_dv_component_invert_test_double(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMaxAbsElementTest
  : public UnitTest
{
public:
  DenseVectorMaxAbsElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorMaxAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorMaxAbsElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i) * (i%2 == 0 ? DT_(1) : DT_(-1)));
      }

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ max = a.max_abs_element();

      TEST_CHECK_EQUAL(max, DT_(size-1));
    }
  }
};
DenseVectorMaxAbsElementTest <float, std::uint32_t> dv_max_abs_element_test_float_uint32(PreferredBackend::generic);
DenseVectorMaxAbsElementTest <double, std::uint32_t> dv_max_abs_element_test_double_uint32(PreferredBackend::generic);
DenseVectorMaxAbsElementTest <float, std::uint64_t> dv_max_abs_element_test_float_uint64(PreferredBackend::generic);
DenseVectorMaxAbsElementTest <double, std::uint64_t> dv_max_abs_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorMaxAbsElementTest <__float128, std::uint32_t> dv_max_abs_element_test_float128_uint32(PreferredBackend::generic);
DenseVectorMaxAbsElementTest <__float128, std::uint64_t> dv_max_abs_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorMaxAbsElementTest <float, std::uint64_t> mkl_dv_max_abs_test_float_uint64(PreferredBackend::mkl);
DenseVectorMaxAbsElementTest <double, std::uint64_t> mkl_dv_max_abs_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorMaxAbsElementTest <Half, std::uint32_t> dv_max_abs_test_half_uint32(PreferredBackend::generic);
DenseVectorMaxAbsElementTest <Half, std::uint64_t> dv_max_abs_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorMaxAbsElementTest <Half, std::uint32_t> cuda_dv_max_abs_test_half_uint32(PreferredBackend::cuda);
DenseVectorMaxAbsElementTest <Half, std::uint64_t> cuda_dv_max_abs_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorMaxAbsElementTest <float, std::uint32_t> cuda_dv_max_abs_element_test_float_uint32(PreferredBackend::cuda);
DenseVectorMaxAbsElementTest <double, std::uint32_t> cuda_dv_max_abs_element_test_double_uint32(PreferredBackend::cuda);
DenseVectorMaxAbsElementTest <float, std::uint64_t> cuda_dv_max_abs_element_test_float_uint64(PreferredBackend::cuda);
DenseVectorMaxAbsElementTest <double, std::uint64_t> cuda_dv_max_abs_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMinAbsElementTest
  : public UnitTest
{
public:
  DenseVectorMinAbsElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorMinAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorMinAbsElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i) * (i%2 == 0 ? DT_(1) : DT_(-1)));
      }

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ min = a.min_abs_element();

      TEST_CHECK_EQUAL(min, DT_(0));
    }
  }
};
DenseVectorMinAbsElementTest <float, std::uint32_t> dv_min_abs_element_test_float_uint32(PreferredBackend::generic);
DenseVectorMinAbsElementTest <double, std::uint32_t> dv_min_abs_element_test_double_uint32(PreferredBackend::generic);
DenseVectorMinAbsElementTest <float, std::uint64_t> dv_min_abs_element_test_float_uint64(PreferredBackend::generic);
DenseVectorMinAbsElementTest <double, std::uint64_t> dv_min_abs_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorMinAbsElementTest <__float128, std::uint32_t> dv_min_abs_element_test_float128_uint32(PreferredBackend::generic);
DenseVectorMinAbsElementTest <__float128, std::uint64_t> dv_min_abs_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorMinAbsElementTest <float, std::uint64_t> mkl_dv_min_abs_test_float_uint64(PreferredBackend::mkl);
DenseVectorMinAbsElementTest <double, std::uint64_t> mkl_dv_min_abs_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorMinAbsElementTest <Half, std::uint32_t> dv_min_abs_test_half_uint32(PreferredBackend::generic);
DenseVectorMinAbsElementTest <Half, std::uint64_t> dv_min_abs_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorMinAbsElementTest <Half, std::uint32_t> cuda_dv_min_abs_test_half_uint32(PreferredBackend::cuda);
DenseVectorMinAbsElementTest <Half, std::uint64_t> cuda_dv_min_abs_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorMinAbsElementTest <float, std::uint32_t> cuda_dv_min_abs_element_test_float_uint32(PreferredBackend::cuda);
DenseVectorMinAbsElementTest <double, std::uint32_t> cuda_dv_min_abs_element_test_double_uint32(PreferredBackend::cuda);
DenseVectorMinAbsElementTest <float, std::uint64_t> cuda_dv_min_abs_element_test_float_uint64(PreferredBackend::cuda);
DenseVectorMinAbsElementTest <double, std::uint64_t> cuda_dv_min_abs_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMaxElementTest
  : public UnitTest
{
public:
  DenseVectorMaxElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorMaxElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorMaxElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(5) ; size < Index(1e3) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(i));
      }
      a(0, DT_(-5));

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ max = a.max_element();

      TEST_CHECK_EQUAL(max, DT_(size-1));
    }
  }
};
DenseVectorMaxElementTest <float, std::uint32_t> dv_max_element_test_float_uint32(PreferredBackend::generic);
DenseVectorMaxElementTest <double, std::uint32_t> dv_max_element_test_double_uint32(PreferredBackend::generic);
DenseVectorMaxElementTest <float, std::uint64_t> dv_max_element_test_float_uint64(PreferredBackend::generic);
DenseVectorMaxElementTest <double, std::uint64_t> dv_max_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorMaxElementTest <__float128, std::uint32_t> dv_max_element_test_float128_uint32(PreferredBackend::generic);
DenseVectorMaxElementTest <__float128, std::uint64_t> dv_max_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorMaxElementTest <float, std::uint64_t> mkl_dv_max_test_float_uint64(PreferredBackend::mkl);
DenseVectorMaxElementTest <double, std::uint64_t> mkl_dv_max_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorMaxElementTest <Half, std::uint32_t> dv_max_test_half_uint32(PreferredBackend::generic);
DenseVectorMaxElementTest <Half, std::uint64_t> dv_max_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorMaxElementTest <Half, std::uint32_t> cuda_dv_max_test_half_uint32(PreferredBackend::cuda);
DenseVectorMaxElementTest <Half, std::uint64_t> cuda_dv_max_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorMaxElementTest <float, std::uint32_t> cuda_dv_max_element_test_float_uint32(PreferredBackend::cuda);
DenseVectorMaxElementTest <double, std::uint32_t> cuda_dv_max_element_test_double_uint32(PreferredBackend::cuda);
DenseVectorMaxElementTest <float, std::uint64_t> cuda_dv_max_element_test_float_uint64(PreferredBackend::cuda);
DenseVectorMaxElementTest <double, std::uint64_t> cuda_dv_max_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseVectorMinElementTest
  : public UnitTest
{
public:
  DenseVectorMinElementTest(PreferredBackend backend)
    : UnitTest("DenseVectorMinElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseVectorMinElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(5) ; size < Index(1e3) ; size*=2)
    {
      DenseVector<DT_, IT_> a(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a(i, DT_(DT_(i) - DT_(3)));
      }

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
      Adjacency::Permutation prm_rnd(a.size(), rng);
      a.permute(prm_rnd);

      DT_ min = a.min_element();

      TEST_CHECK_EQUAL(min, DT_(-3.));
    }
  }
};
DenseVectorMinElementTest <float, std::uint32_t> dv_min_element_test_float_uint32(PreferredBackend::generic);
DenseVectorMinElementTest <double, std::uint32_t> dv_min_element_test_double_uint32(PreferredBackend::generic);
DenseVectorMinElementTest <float, std::uint64_t> dv_min_element_test_float_uint64(PreferredBackend::generic);
DenseVectorMinElementTest <double, std::uint64_t> dv_min_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_QUADMATH
DenseVectorMinElementTest <__float128, std::uint32_t> dv_min_element_test_float128_uint32(PreferredBackend::generic);
DenseVectorMinElementTest <__float128, std::uint64_t> dv_min_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_MKL
DenseVectorMinElementTest <float, std::uint64_t> mkl_dv_min_test_float_uint64(PreferredBackend::mkl);
DenseVectorMinElementTest <double, std::uint64_t> mkl_dv_min_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseVectorMinElementTest <Half, std::uint32_t> dv_min_test_half_uint32(PreferredBackend::generic);
DenseVectorMinElementTest <Half, std::uint64_t> dv_min_test_half_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseVectorMinElementTest <Half, std::uint32_t> cuda_dv_min_test_half_uint32(PreferredBackend::cuda);
DenseVectorMinElementTest <Half, std::uint64_t> cuda_dv_min_test_half_uint64(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseVectorMinElementTest <float, std::uint32_t> cuda_dv_min_element_test_float_uint32(PreferredBackend::cuda);
DenseVectorMinElementTest <double, std::uint32_t> cuda_dv_min_element_test_double_uint32(PreferredBackend::cuda);
DenseVectorMinElementTest <float, std::uint64_t> cuda_dv_min_element_test_float_uint64(PreferredBackend::cuda);
DenseVectorMinElementTest <double, std::uint64_t> cuda_dv_min_element_test_double_uint64(PreferredBackend::cuda);
#endif
