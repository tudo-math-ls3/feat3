// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/random.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>

#include <sstream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse matrix csr class.
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
  class SparseMatrixCSRTest
  : public UnitTest
{
public:
  SparseMatrixCSRTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    SparseMatrixCSR<DT_, IT_> zero;
    TEST_CHECK(zero.empty());

    SparseMatrixCSR<DT_, IT_> zero3(10, 11, 12);
    TEST_CHECK_EQUAL(zero3.used_elements(), 12);
    TEST_CHECK_EQUAL(zero3.rows(), 10);
    TEST_CHECK_EQUAL(zero3.columns(), 11);
    TEST_CHECK_EQUAL(zero3.size(), 110);


    SparseMatrixCSR<DT_, IT_> empty1(11, 12, 111);
    SparseMatrixCSR<DT_, IT_> empty2;
    SparseMatrixCSR<DT_, IT_> empty3;
    empty2.convert(empty1);
    empty3.convert(empty2);
    TEST_CHECK_EQUAL(empty1.rows(), empty3.rows());
    TEST_CHECK_EQUAL(empty1.columns(), empty3.columns());
    TEST_CHECK_EQUAL(empty1.used_elements(), empty3.used_elements());

    SparseMatrixCSR<DT_, IT_> empty4(empty2.layout());
    SparseMatrixCSR<DT_, IT_> empty5(empty3.layout());
    empty4.convert(empty1);
    empty5.convert(empty4);
    TEST_CHECK_EQUAL(empty5.rows(), empty5.rows());
    TEST_CHECK_EQUAL(empty5.columns(), empty5.columns());
    TEST_CHECK_EQUAL(empty5.used_elements(), empty5.used_elements());
    empty5.convert(zero);
    TEST_CHECK_EQUAL(empty5.rows(), 0);
    TEST_CHECK_EQUAL(empty5.columns(), 0);
    TEST_CHECK_EQUAL(empty5.used_elements(), 0);

    SparseMatrixFactory<DT_, IT_> a01(IT_(10), IT_(10));
    SparseMatrixCSR<DT_, IT_> b01(a01.make_csr());
    TEST_CHECK_EQUAL(b01.size(), 100ul);
    TEST_CHECK_EQUAL(b01.rows(), 10ul);
    TEST_CHECK_EQUAL(b01.columns(), 10ul);
    TEST_CHECK_EQUAL(b01.used_elements(), 0ul);

    SparseMatrixFactory<DT_, IT_> a(IT_(10), IT_(10));
    a.add(IT_(1), IT_(2), DT_(7));
    a.add(IT_(5), IT_(5), DT_(2));
    a.add(IT_(5), IT_(7), DT_(3));
    a.add(IT_(5), IT_(2), DT_(4));
    SparseMatrixCSR<DT_, IT_> b(a.make_csr());
    TEST_CHECK(!b.empty());
    TEST_CHECK_EQUAL(b.used_elements(), a.used_elements());
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b(1, 2), DT_(7));
    TEST_CHECK_EQUAL(b(5, 5), DT_(2));
    TEST_CHECK_EQUAL(b(5, 7), DT_(3));
    TEST_CHECK_EQUAL(b(5, 2), DT_(4));
    TEST_CHECK_EQUAL(b(1, 1), DT_(0));

    Index bandw, bandw_idx;
    b.bandwidth_row(bandw, bandw_idx);
    TEST_CHECK_EQUAL(bandw, Index(6));
    TEST_CHECK_EQUAL(bandw_idx, Index(5));
    b.bandwidth_column(bandw, bandw_idx);
    TEST_CHECK_EQUAL(bandw, Index(5));
    TEST_CHECK_EQUAL(bandw_idx, Index(2));

    Index radius, radius_idx;
    b.radius_row(radius, radius_idx);
    TEST_CHECK_EQUAL(radius, Index(3));
    TEST_CHECK_EQUAL(radius_idx, Index(5));
    b.radius_column(radius, radius_idx);
    TEST_CHECK_EQUAL(radius, Index(3));
    TEST_CHECK_EQUAL(radius_idx, Index(2));

    SparseMatrixCSR<DT_, IT_> bl(b.layout());
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    bl = b.layout();
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    typename SparseLayout<IT_, SparseLayoutId::lt_csr>::template MatrixType<DT_> x(b.layout());
    TEST_CHECK_EQUAL((void*)x.row_ptr(), (void*)b.row_ptr());
    TEST_CHECK_NOT_EQUAL((void*)x.val(), (void*)b.val());
    /// \compilerhack icc 14.x and msvc do not understand the following single line, so we need a typedef detour here
#if defined(FEAT_COMPILER_MICROSOFT) || (defined(FEAT_COMPILER_INTEL) && __INTEL_COMPILER < 1500)
    typedef decltype(b.layout()) LayoutId;
    typename LayoutId::template MatrixType<DT_> y(b.layout());
#else
    typename decltype(b.layout())::template MatrixType<DT_> y(b.layout());
#endif
    TEST_CHECK_EQUAL((void*)y.row_ptr(), (void*)b.row_ptr());
    TEST_CHECK_NOT_EQUAL((void*)y.val(), (void*)b.val());


    SparseMatrixCSR<DT_, IT_> z;
    z.convert(b);
    TEST_CHECK_EQUAL(z.used_elements(), 4ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z(1, 2), DT_(7));
    TEST_CHECK_EQUAL(z(5, 5), DT_(2));

    SparseMatrixCSR<DT_, IT_> c;
    c.clone(b);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_EQUAL((void*)c.col_ind(), (void*)b.col_ind());
    c = b.clone(CloneMode::Deep);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_NOT_EQUAL((void*)c.col_ind(), (void*)b.col_ind());

    DenseVector<IT_, IT_> col_ind(c.used_elements(), c.col_ind());
    DenseVector<DT_, IT_> val(c.used_elements(), c.val());
    DenseVector<IT_, IT_> row_ptr(c.rows() + 1, c.row_ptr());
    SparseMatrixCSR<DT_, IT_> d(c.rows(), c.columns(), col_ind, val, row_ptr);
    TEST_CHECK_LESS_THAN(d.max_rel_diff(c), eps);

    SparseMatrixCSR<DT_, IT_> e;
    e.convert(c);
    TEST_CHECK_LESS_THAN(e.max_rel_diff(c), eps);
    e.copy(c);
    TEST_CHECK_LESS_THAN(e.max_rel_diff(c), eps);
    e.clone(c);
    b.clone(e);
    TEST_CHECK_LESS_THAN(b.max_rel_diff(c), eps);


    // new clone testing
    auto clone1 = b.clone(CloneMode::Deep);
    TEST_CHECK_LESS_THAN(clone1.max_rel_diff(b), eps);
    MemoryPool::set_memory(clone1.val() + 1, DT_(132));
    TEST_CHECK_LESS_THAN(eps, clone1.max_rel_diff(b));
    TEST_CHECK_NOT_EQUAL((void*)clone1.val(), (void*)b.val());
    TEST_CHECK_NOT_EQUAL((void*)clone1.row_ptr(), (void*)b.row_ptr());
    auto clone2 = clone1.clone(CloneMode::Layout);
    MemoryPool::set_memory(clone2.val(), DT_(4713), clone2.used_elements());
    TEST_CHECK_NOT_EQUAL(clone2(5, 5), clone1(5, 5));
    TEST_CHECK_NOT_EQUAL((void*)clone2.val(), (void*)clone1.val());
    TEST_CHECK_EQUAL((void*)clone2.row_ptr(), (void*)clone1.row_ptr());
    auto clone3 = clone1.clone(CloneMode::Weak);
    TEST_CHECK_LESS_THAN(clone3.max_rel_diff(clone1), eps);
    MemoryPool::set_memory(clone3.val() + 1, DT_(133));
    TEST_CHECK_LESS_THAN(eps, clone3.max_rel_diff(clone1));
    TEST_CHECK_NOT_EQUAL((void*)clone3.val(), (void*)clone1.val());
    TEST_CHECK_EQUAL((void*)clone3.row_ptr(), (void*)clone1.row_ptr());
    auto clone4 = clone1.clone(CloneMode::Shallow);
    TEST_CHECK_LESS_THAN(clone4.max_rel_diff(clone1), eps);
    MemoryPool::set_memory(clone4.val() + 1, DT_(134));
    TEST_CHECK_LESS_THAN(clone4.max_rel_diff(clone1), eps);
    TEST_CHECK_EQUAL((void*)clone4.val(), (void*)clone1.val());
    TEST_CHECK_EQUAL((void*)clone4.row_ptr(), (void*)clone1.row_ptr());

    SparseMatrixFactory<DT_, IT_> ffac(IT_(10), IT_(10));
    for (IT_ row(0); row < ffac.rows(); ++row)
    {
      for (IT_ col(0); col < ffac.columns(); ++col)
      {
        if (row == col)
          ffac.add(row, col, DT_(2));
        else if ((row == col + 1) || (row + 1 == col))
          ffac.add(row, col, DT_(-1));
      }
    }
    SparseMatrixCSR<DT_, IT_> f(ffac.make_csr());

    // shrink test
    SparseMatrixCSR<DT_, IT_> l(f.clone());
    l.shrink(DT_(1.9));
    TEST_CHECK_EQUAL(l.used_elements(), 10ul);
  }

};

SparseMatrixCSRTest <float, std::uint64_t> cpu_sparse_matrix_csr_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRTest <double, std::uint64_t> cpu_sparse_matrix_csr_test_double_uint64(PreferredBackend::generic);
SparseMatrixCSRTest <float, std::uint32_t> cpu_sparse_matrix_csr_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRTest <double, std::uint32_t> cpu_sparse_matrix_csr_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRTest <float, std::uint64_t> mkl_cpu_sparse_matrix_csr_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRTest <double, std::uint64_t> mkl_cpu_sparse_matrix_csr_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRTest <__float128, std::uint64_t> cpu_sparse_matrix_csr_test_float128_uint64(PreferredBackend::generic);
SparseMatrixCSRTest <__float128, std::uint32_t> cpu_sparse_matrix_csr_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRTest <Half, std::uint32_t> sparse_matrix_csr_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRTest <Half, std::uint64_t> sparse_matrix_csr_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRTest <float, std::uint64_t> cuda_sparse_matrix_csr_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRTest <double, std::uint64_t> cuda_sparse_matrix_csr_test_double_uint64(PreferredBackend::cuda);
SparseMatrixCSRTest <float, std::uint32_t> cuda_sparse_matrix_csr_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRTest <double, std::uint32_t> cuda_sparse_matrix_csr_test_double_uint32(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRSerializeTest
  : public UnitTest
{
public:
  SparseMatrixCSRSerializeTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRSerializeTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    SparseMatrixFactory<DT_, IT_> ffac(IT_(10), IT_(10));
    for (IT_ row(0); row < ffac.rows(); ++row)
    {
      for (IT_ col(0); col < ffac.columns(); ++col)
      {
        if (row == col)
          ffac.add(row, col, DT_(2));
        else if ((row == col + 1) || (row + 1 == col))
          ffac.add(row, col, DT_(-1));
      }
    }
    SparseMatrixCSR<DT_, IT_> f(ffac.make_csr());

    BinaryStream bs;
    f.write_out(FileMode::fm_csr, bs);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));
    bs.seekg(0);
    SparseMatrixCSR<DT_, IT_> g(FileMode::fm_csr, bs);
    TEST_CHECK_LESS_THAN(g.max_rel_diff(f), eps);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));

    std::stringstream ts;
    f.write_out(FileMode::fm_mtx, ts);
    SparseMatrixCSR<DT_, IT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_LESS_THAN(j.max_rel_diff(f), eps);

    std::stringstream ts2;
    f.write_out(FileMode::fm_mtx, ts2, true);
    SparseMatrixCSR<DT_, IT_> j2(FileMode::fm_mtx, ts2);
    TEST_CHECK_LESS_THAN(j2.max_rel_diff(f), eps);

    auto kp = f.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixCSR<DT_, IT_> k(kp);
    TEST_CHECK_LESS_THAN(k.max_rel_diff(f), eps);
#ifdef FEAT_HAVE_ZLIB
    auto zl = f.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixCSR<DT_, IT_> zlib(zl);
    TEST_CHECK_LESS_THAN(zlib.max_rel_diff(f), eps);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = f.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixCSR<DT_, IT_> zfp(zf);
    for (Index row(0); row < f.rows(); ++row)
    {
      for (Index col(0); col < f.columns(); ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(zfp(row, col), f(row, col), Math::pow(Math::eps<DT_>(), DT_(0.8)));
      }
    }
#endif
  }
};
SparseMatrixCSRSerializeTest <float, std::uint64_t> cpu_sparse_matrix_csr_serialize_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRSerializeTest <double, std::uint64_t> cpu_sparse_matrix_csr_serialize_test_double_uint64(PreferredBackend::generic);
SparseMatrixCSRSerializeTest <float, std::uint32_t> cpu_sparse_matrix_csr_serialize_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRSerializeTest <double, std::uint32_t> cpu_sparse_matrix_csr_serialize_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRSerializeTest <float, std::uint64_t> mkl_cpu_sparse_matrix_csr_serialize_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRSerializeTest <double, std::uint64_t> mkl_cpu_sparse_matrix_csr_serialize_test_double_uint64(PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SparseMatrixCSRSerializeTest <__float128, std::uint64_t> cpu_sparse_matrix_csr_serialize_test_float128_uint64(PreferredBackend::generic);
//SparseMatrixCSRSerializeTest <__float128, std::uint32_t> cpu_sparse_matrix_csr_serialize_test_float128_uint32(PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRSerializeTest <Half, std::uint32_t> sparse_matrix_csr_serialize_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRSerializeTest <Half, std::uint64_t> sparse_matrix_csr_serialize_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRSerializeTest <float, std::uint64_t> cuda_sparse_matrix_csr_serialize_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRSerializeTest <double, std::uint64_t> cuda_sparse_matrix_csr_serialize_test_double_uint64(PreferredBackend::cuda);
SparseMatrixCSRSerializeTest <float, std::uint32_t> cuda_sparse_matrix_csr_serialize_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRSerializeTest <double, std::uint32_t> cuda_sparse_matrix_csr_serialize_test_double_uint32(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRApplyTest
  : public UnitTest
{
public:
  SparseMatrixCSRApplyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRApplyTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.7));
    DT_ s(DT_(0.4711));
    for (IT_ size(1); size < IT_(1e3); size *= IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size);
      DenseVector<DT_, IT_> x(size);
      DenseVector<DT_, IT_> y(size);
      DenseVector<DT_, IT_> ref(size);
      DenseVector<DT_, IT_> ax(size);
      for (Index i(0); i < size; ++i)
      {
        x(i, DT_(i % 100) * DT_(1.234));
        y(i, DT_(2) - DT_(i % 42));
        if (i == 0 && size > 1)
          ax(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i + 1) % 100)));
        else if (i == size - 1 && size > 1)
          ax(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i - 1) % 100)));
        else if (size == 1 && i == 0)
          ax(i, DT_(0));
        else
          ax(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i - 1) % 100) - DT_((i + 1) % 100)));
      }

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
          }
          else if ((row == col + 1) || (row + 1 == col))
          {
            a_fac.add(row, col, DT_(-1));
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());

      DenseVector<DT_, IT_> r(size);

      // apply-test for alpha = 0.0
      a.apply(r, x, y, DT_(0.0));
      ref.copy(y);
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-test for alpha = -1.0
      a.apply(r, x, y, DT_(-1.0));
      a.apply(ref, x);
      ref.scale(ref, DT_(-1.0));
      ref.axpy(y);
      for(Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-test for alpha = -1.0 and &r==&y
      r.copy(y);
      a.apply(r, x, r, DT_(-1.0));
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-test for alpha = 4711.1
      //r.axpy(s, a, x, y);
      a.apply(r, x, y, s);
      //ref.product_matvec(a, x);
      a.apply(ref, x);
      ref.scale(ref, s);
      ref.axpy(y);
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-test for alpha = 4711.1 and &r==&y
      r.copy(y);
      a.apply(r, x, r, s);
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      a.apply(r, x);
      ref.copy(ax);
      for(Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(2)*eps);

      // transposed apply-test for alpha = -1
      a.apply_transposed(r, x, y,  DT_(-1.0));
      SparseMatrixCSR<DT_, IT_> at = a.transpose();
      at.apply(ref, x);
      ref.scale(ref,  DT_(-1.0));
      ref.axpy(y);

      for(Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);
    }
  }
};

SparseMatrixCSRApplyTest <float, std::uint64_t> sm_csr_apply_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRApplyTest <double, std::uint64_t> sm_csr_apply_test_double_uint64(PreferredBackend::generic);
SparseMatrixCSRApplyTest <float, std::uint32_t> sm_csr_apply_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRApplyTest <double, std::uint32_t> sm_csr_apply_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRApplyTest <float, std::uint64_t> mkl_sm_csr_apply_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRApplyTest <double, std::uint64_t> mkl_sm_csr_apply_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRApplyTest <__float128, std::uint64_t> sm_csr_apply_test_float128_uint64(PreferredBackend::generic);
SparseMatrixCSRApplyTest <__float128, std::uint32_t> sm_csr_apply_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRApplyTest <Half, std::uint32_t> sm_csr_apply_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRApplyTest <Half, std::uint64_t> sm_csr_apply_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRApplyTest <float, std::uint64_t> cuda_sm_csr_apply_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRApplyTest <double, std::uint64_t> cuda_sm_csr_apply_test_double_uint64(PreferredBackend::cuda);
SparseMatrixCSRApplyTest <float, std::uint32_t> cuda_sm_csr_apply_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRApplyTest <double, std::uint32_t> cuda_sm_csr_apply_test_double_uint32(PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRApplyTest <Half, std::uint32_t> cuda_sm_csr_apply_test_half_uint32(PreferredBackend::cuda);
SparseMatrixCSRApplyTest <Half, std::uint64_t> cuda_sm_csr_apply_test_half_uint64(PreferredBackend::cuda);
#endif
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRBApplyTest
  : public UnitTest
{
public:
  SparseMatrixCSRBApplyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRBApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRBApplyTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.7));
    for (IT_ size(1); size < IT_(1e3); size *= IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size);
      DenseVector<DT_, IT_> ref_x(size);
      DenseVector<DT_, IT_> aref_x(size);
      for (Index i(0); i < size; ++i)
      {
        ref_x(i, DT_(i % 100) * DT_(1.234));
        if (i == 0 && size > 1)
          aref_x(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i + 1) % 100)));
        else if (i == size - 1 && size > 1)
          aref_x(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i - 1) % 100)));
        else if (i == 0 && size == 1)
          aref_x(i, DT_(1.234) * DT_(2) * DT_(i % 100));
        else
          aref_x(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i - 1) % 100) - DT_((i + 1) % 100)));
      }

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
          }
          else if ((row == col + 1) || (row + 1 == col))
          {
            a_fac.add(row, col, DT_(-1));
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());


      DenseVectorBlocked<DT_, IT_, 3> x(size);
      for (Index i(0); i < size; ++i)
      {
        auto temp = x(i);
        temp[0] = ref_x(i);
        temp[1] = DT_(0.5) * ref_x(i);
        temp[2] = DT_(2.0) * ref_x(i);
        x(i, temp);
      }

      DenseVectorBlocked<DT_, IT_, 3> y(size);
      for (Index i(0); i < size; ++i)
      {
        auto temp = y(i);
        temp[0] = DT_(i);
        temp[1] = DT_(i * 2);
        temp[2] = DT_(i * 3);
        y(i, temp);
      }

      DenseVectorBlocked<DT_, IT_, 3> r(size);

      a.apply(r, x);
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[0], aref_x(i), eps);
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[1], aref_x(i) * DT_(0.5), eps);
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[2], aref_x(i) * DT_(2.0), eps*DT_(10));
      }

      a.apply(r, x, y, DT_(-1));
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_RELATIVE(r(i)[0], y(i)[0] - aref_x(i), eps);
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[1], y(i)[1] - aref_x(i) * DT_(0.5), eps);
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[2], y(i)[2] - aref_x(i) * DT_(2.0), eps);
      }

      DT_ alpha(0.75);
      a.apply(r, x, y, alpha);
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[0], y(i)[0] + alpha * aref_x(i), eps);
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[1], y(i)[1] + alpha * aref_x(i) * DT_(0.5), eps);
        TEST_CHECK_RELATIVE(r(i)[2], y(i)[2] + alpha * aref_x(i) * DT_(2.0), eps);
      }
    }
  }
};

SparseMatrixCSRBApplyTest <float, std::uint64_t> sm_csrib_apply_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRBApplyTest <double, std::uint64_t> sm_csrib_apply_test_double_uint64(PreferredBackend::generic);
SparseMatrixCSRBApplyTest <float, std::uint32_t> sm_csrib_apply_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRBApplyTest <double, std::uint32_t> sm_csrib_apply_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRBApplyTest <float, std::uint64_t> mkl_sm_csrib_apply_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRBApplyTest <double, std::uint64_t> mkl_sm_csrib_apply_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRBApplyTest <__float128, std::uint64_t> sm_csrib_apply_test_float128_uint64(PreferredBackend::generic);
SparseMatrixCSRBApplyTest <__float128, std::uint32_t> sm_csrib_apply_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRBApplyTest <Half, std::uint32_t> sm_csrib_apply_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRBApplyTest <Half, std::uint64_t> sm_csrib_apply_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRBApplyTest <float, std::uint64_t> cuda_sm_csrib_apply_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRBApplyTest <double, std::uint64_t> cuda_sm_csrib_apply_test_double_uint64(PreferredBackend::cuda);
SparseMatrixCSRBApplyTest <float, std::uint32_t> cuda_sm_csrib_apply_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRBApplyTest <double, std::uint32_t> cuda_sm_csrib_apply_test_double_uint32(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRScaleTest
  : public UnitTest
{
public:
  SparseMatrixCSRScaleTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRScaleTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRScaleTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    for (Index size(2); size < Index(3e2); size *= 2)
    {
      DT_ s(DT_(4.321));
      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_fac(IT_(size), IT_(size + 2));

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
            a_fac.add(row, col, DT_(2));
          else if ((row == col + 1) || (row + 1 == col))
            a_fac.add(row, col, DT_(-1));

          if (row == col)
            ref_fac.add(row, col, DT_(2) * s);
          else if ((row == col + 1) || (row + 1 == col))
            ref_fac.add(row, col, DT_(-1) * s);
        }
      }

      SparseMatrixCSR<DT_, IT_> ref(ref_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> b;
      b.clone(a);

      b.scale(a, s);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref), eps);

      a.scale(a, s);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), eps);
    }
  }
};

SparseMatrixCSRScaleTest <float, std::uint32_t> sm_csr_scale_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRScaleTest <double, std::uint32_t> sm_csr_scale_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRScaleTest <float, std::uint64_t> sm_csr_scale_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRScaleTest <double, std::uint64_t> sm_csr_scale_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRScaleTest <float, std::uint64_t> mkl_sm_csr_scale_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRScaleTest <double, std::uint64_t> mkl_sm_csr_scale_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRScaleTest <__float128, std::uint32_t> sm_csr_scale_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRScaleTest <__float128, std::uint64_t> sm_csr_scale_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRScaleTest <Half, std::uint32_t> sm_csr_scale_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRScaleTest <Half, std::uint64_t> sm_csr_scale_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRScaleTest <float, std::uint32_t> cuda_sm_csr_scale_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRScaleTest <double, std::uint32_t> cuda_sm_csr_scale_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRScaleTest <float, std::uint64_t> cuda_sm_csr_scale_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRScaleTest <double, std::uint64_t> cuda_sm_csr_scale_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRScaleRowColTest
  : public UnitTest
{
public:
  SparseMatrixCSRScaleRowColTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRScaleRowColTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRScaleRowColTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2); size < Index(3e2); size *= 3)
    {
      const DT_ pi(Math::pi<DT_>());
      const DT_ eps(Math::pow(Math::eps<DT_>(), DT_(0.8)));

      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
            a_fac.add(row, col, DT_(2));
          else if ((row == col + 1) || (row + 1 == col))
            a_fac.add(row, col, DT_(-1));
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> b(a.clone());

      // Scale rows
      DenseVector<DT_, IT_> s1(a.rows());
      for (Index i(0); i < s1.size(); ++i)
      {
        s1(i, pi * DT_(i % 3 + 1) - DT_(5.21) + DT_(i));
      }
      b.scale_rows(b, s1);
      for (Index row(0); row < a.rows(); ++row)
      {
        for (Index col(0); col < a.columns(); ++col)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(b(row, col), a(row, col) * s1(row), eps);
        }
      }

      // Scale cols
      DenseVector<DT_, IT_> s2(a.columns());
      for (Index i(0); i < s2.size(); ++i)
      {
        s2(i, pi * DT_(i % 3 + 1) - DT_(5.21) + DT_(i));
      }
      b.scale_cols(a, s2);
      for (Index row(0); row < a.rows(); ++row)
      {
        for (Index col(0); col < a.columns(); ++col)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(b(row, col), a(row, col) * s2(col), eps);
        }
      }
    }
  }
};

SparseMatrixCSRScaleRowColTest <float, std::uint32_t> sm_csr_scale_row_col_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest <double, std::uint32_t> sm_csr_scale_row_col_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest <float, std::uint64_t> sm_csr_scale_row_col_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest <double, std::uint64_t> sm_csr_scale_row_col_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRScaleRowColTest <float, std::uint64_t> mkl_sm_csr_scale_row_col_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRScaleRowColTest <double, std::uint64_t> mkl_sm_csr_scale_row_col_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRScaleRowColTest <__float128, std::uint32_t> sm_csr_scale_row_col_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest <__float128, std::uint64_t> sm_csr_scale_row_col_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRScaleRowColTest <Half, std::uint32_t> sm_csr_scale_row_col_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest <Half, std::uint64_t> sm_csr_scale_row_col_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRScaleRowColTest <float, std::uint32_t> cuda_sm_csr_scale_row_col_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRScaleRowColTest <double, std::uint32_t> cuda_sm_csr_scale_row_col_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRScaleRowColTest <float, std::uint64_t> cuda_sm_csr_scale_row_col_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRScaleRowColTest <double, std::uint64_t> cuda_sm_csr_scale_row_col_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRTranspositionTest
  : public UnitTest
{

public:
  SparseMatrixCSRTranspositionTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRTranspositionTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRTranspositionTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    for (Index size(2); size < Index(3e2); size *= 4)
    {
      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
            a_fac.add(row, col, DT_(2));
          else if (row == col + 1)
            a_fac.add(row, col, DT_(-1));
          else if (row + 1 == col)
            a_fac.add(row, col, DT_(-3));
        }
      }
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());


      SparseMatrixCSR<DT_, IT_> b;
      b.transpose(a);

      for (Index i(0); i < a.rows(); ++i)
      {
        for (Index j(0); j < a.columns(); ++j)
        {
          TEST_CHECK_EQUAL(b(j, i), a(i, j));
        }
      }

      b = b.transpose();

      TEST_CHECK_LESS_THAN(a.max_rel_diff(b), eps);
    }
  }
};

SparseMatrixCSRTranspositionTest <float, std::uint32_t> sm_csr_transposition_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest <double, std::uint32_t> sm_csr_transposition_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest <float, std::uint64_t> sm_csr_transposition_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest <double, std::uint64_t> sm_csr_transposition_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRTranspositionTest <float, std::uint64_t> mkl_sm_csr_transposition_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRTranspositionTest <double, std::uint64_t> mkl_sm_csr_transposition_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRTranspositionTest <__float128, std::uint32_t> sm_csr_transposition_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest <__float128, std::uint64_t> sm_csr_transposition_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRTranspositionTest <Half, std::uint32_t> sm_csr_transposition_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest <Half, std::uint64_t> sm_csr_transposition_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRTranspositionTest <float, std::uint32_t> cuda_sm_csr_transposition_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRTranspositionTest <double, std::uint32_t> cuda_sm_csr_transposition_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRTranspositionTest <float, std::uint64_t> cuda_sm_csr_transposition_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRTranspositionTest <double, std::uint64_t> cuda_sm_csr_transposition_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRPermuteTest
  : public UnitTest
{
public:
  SparseMatrixCSRPermuteTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRPermuteTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRPermuteTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    for (IT_ size(25); size < IT_(1e3); size *= IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size);
      DenseVector<DT_, IT_> x(size);
      for (Index i(0); i < size; ++i)
      {
        x(i, DT_(i % 100) * DT_(1.234));
      }

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
          }
          else if ((row == col + 7) || (row + 7 == col))
          {
            a_fac.add(row, col, DT_(-1));
          }
          else if ((row == col + 15) || (row + 15 == col))
          {
            a_fac.add(row, col, DT_(-2));
          }
          else if ((row == col + a_fac.columns() / 2) || (row + a_fac.rows() / 2 == col))
          {
            a_fac.add(row, col, DT_(1));
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());

      DenseVector<DT_, IT_> r(size);
      a.apply(r, x);
      DT_ ref_norm = r.norm2();

      auto a_backup = a.clone(CloneMode::Deep);

      Random rng;
      std::cout << "RNG Seed: " << rng.get_seed() << "\n";
      Adjacency::Permutation perm(a.rows(), rng);

      a.permute(perm, perm);
      x.permute(perm);

      a.apply(r, x);
      DT_ norm = r.norm2();
      DT_ deviation = norm / ref_norm;
      TEST_CHECK(deviation > DT_(0.99));
      TEST_CHECK(deviation < DT_(1.01));

      a = a_backup.clone(CloneMode::Deep);
      auto perm_inv = perm.inverse();
      a.permute(perm_inv, perm);
      a.permute(perm, perm_inv);
      TEST_CHECK_LESS_THAN(a.max_rel_diff(a_backup), eps);
    }
  }
};

SparseMatrixCSRPermuteTest <float, std::uint64_t> sm_csr_permute_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRPermuteTest <double, std::uint64_t> sm_csr_permute_test_double_uint64(PreferredBackend::generic);
SparseMatrixCSRPermuteTest <float, std::uint32_t> sm_csr_permute_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRPermuteTest <double, std::uint32_t> sm_csr_permute_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRPermuteTest <float, std::uint64_t> mkl_sm_csr_permute_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRPermuteTest <double, std::uint64_t> mkl_sm_csr_permute_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRPermuteTest <__float128, std::uint64_t> sm_csr_permute_test_float128_uint64(PreferredBackend::generic);
SparseMatrixCSRPermuteTest <__float128, std::uint32_t> sm_csr_permute_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRPermuteTest <Half, std::uint32_t> sm_csr_permute_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRPermuteTest <Half, std::uint64_t> sm_csr_permute_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRPermuteTest <float, std::uint64_t> cuda_sm_csr_permute_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRPermuteTest <double, std::uint64_t> cuda_sm_csr_permute_test_double_uint64(PreferredBackend::cuda);
SparseMatrixCSRPermuteTest <float, std::uint32_t> cuda_sm_csr_permute_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRPermuteTest <double, std::uint32_t> cuda_sm_csr_permute_test_double_uint32(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRDiagTest
  : public UnitTest
{
public:
  SparseMatrixCSRDiagTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRDiagTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRDiagTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    for (IT_ size(2); size < IT_(3e2); size *= IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size);
      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
            a_fac.add(row, col, DT_(DT_(col % 100) / DT_(2)));
          else if ((row == col + 1) || (row + 1 == col))
            a_fac.add(row, col, DT_(-1));
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());

      auto ref = a.create_vector_l();
      for (Index i(0); i < a_fac.rows(); ++i)
      {
        ref(i, DT_(DT_(i % 100) / DT_(2)));
      }

      auto diag = a.extract_diag();
      TEST_CHECK_LESS_THAN(diag.max_rel_diff(ref), eps);
    }

    SparseMatrixFactory<DT_, IT_> b_fac(16, 16);
    b_fac.add(0, 0, DT_(1.0));
    b_fac.add(2, 2, DT_(2.0));
    b_fac.add(3, 3, DT_(3.0));
    b_fac.add(7, 7, DT_(7.0));

    SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());

    auto ref = b.create_vector_l();
    ref.format(DT_(0.0));
    ref(0, DT_(1.0));
    ref(2, DT_(2.0));
    ref(3, DT_(3.0));
    ref(7, DT_(7.0));
    auto diag = b.extract_diag();
    TEST_CHECK_LESS_THAN(diag.max_rel_diff(ref), eps);
  }
};

SparseMatrixCSRDiagTest <float, std::uint32_t> sm_csr_diag_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRDiagTest <double, std::uint32_t> sm_csr_diag_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRDiagTest <float, std::uint64_t> sm_csr_diag_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRDiagTest <double, std::uint64_t> sm_csr_diag_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRDiagTest <float, std::uint64_t> mkl_sm_csr_diag_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRDiagTest <double, std::uint64_t> mkl_sm_csr_diag_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRDiagTest <__float128, std::uint32_t> sm_csr_diag_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRDiagTest <__float128, std::uint64_t> sm_csr_diag_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRDiagTest <Half, std::uint32_t> sm_csr_diag_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRDiagTest <Half, std::uint64_t> sm_csr_diag_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRDiagTest <float, std::uint32_t> cuda_sm_csr_diag_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRDiagTest <double, std::uint32_t> cuda_sm_csr_diag_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRDiagTest <float, std::uint64_t> cuda_sm_csr_diag_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRDiagTest <double, std::uint64_t> cuda_sm_csr_diag_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRAxpyTest
  : public UnitTest
{
public:
  SparseMatrixCSRAxpyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRAxpyTest()
  {
  }

  virtual void run() const override
  {
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    for (Index size(2); size < Index(3e2); size *= 2)
    {
      DT_ s(DT_(4.321));

      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> b_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_fac2(IT_(size), IT_(size + 2));

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
            b_fac.add(row, col, DT_((row + col) % 15));
            ref_fac.add(row, col, DT_(2) * s + DT_((row + col) % 15));
            ref_fac2.add(row, col, DT_((row + col) % 15) * s + DT_((row + col) % 15));
          }

          else if ((row == col + 1) || (row + 1 == col))
          {
            a_fac.add(row, col, DT_(-1));
            b_fac.add(row, col, DT_((row + col + 1) % 15));
            ref_fac.add(row, col, DT_(-1) * s + DT_((row + col + 1) % 15));
            ref_fac2.add(row, col, DT_((row + col + 1) % 15) * s + DT_((row + col + 1) % 15));
          }
        }
      }

      SparseMatrixCSR<DT_, IT_> ref(ref_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> ref2(ref_fac2.make_csr());
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());

      // r != x
      a.scale(a, s);
      a.axpy(b); /// \todo use axpby here
      TEST_CHECK_LESS_THAN(a.max_rel_diff(ref), eps);

      // r == x
      b.axpy(b, s);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(ref2), eps);
    }
  }
};

SparseMatrixCSRAxpyTest <float, std::uint32_t> sm_csr_axpy_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRAxpyTest <double, std::uint32_t> sm_csr_axpy_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRAxpyTest <float, std::uint64_t> sm_csr_axpy_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRAxpyTest <double, std::uint64_t> sm_csr_axpy_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRAxpyTest <float, std::uint64_t> mkl_sm_csr_axpy_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRAxpyTest <double, std::uint64_t> mkl_sm_csr_axpy_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRAxpyTest <__float128, std::uint32_t> sm_csr_axpy_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRAxpyTest <__float128, std::uint64_t> sm_csr_axpy_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRAxpyTest <Half, std::uint32_t> sm_csr_axpy_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRAxpyTest <Half, std::uint64_t> sm_csr_axpy_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRAxpyTest <float, std::uint32_t> cuda_sm_csr_axpy_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRAxpyTest <double, std::uint32_t> cuda_sm_csr_axpy_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRAxpyTest <float, std::uint64_t> cuda_sm_csr_axpy_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRAxpyTest <double, std::uint64_t> cuda_sm_csr_axpy_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRFrobeniusTest
  : public UnitTest
{
public:
  SparseMatrixCSRFrobeniusTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRFrobeniusTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRFrobeniusTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2); size < Index(3e2); size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
            a_fac.add(row, col, DT_(2));
          else if ((row == col + 1) || (row + 1 == col))
            a_fac.add(row, col, DT_(-1));
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());

      DenseVector<DT_, IT_> refv(a.used_elements(), a.val());
      DT_ ref = refv.norm2();
      DT_ c = a.norm_frobenius();
      TEST_CHECK_EQUAL(c, ref);
    }
  }
};

SparseMatrixCSRFrobeniusTest <float, std::uint32_t> sm_csr_frobenius_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest <double, std::uint32_t> sm_csr_frobenius_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest <float, std::uint64_t> sm_csr_frobenius_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest <double, std::uint64_t> sm_csr_frobenius_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRFrobeniusTest <float, std::uint64_t> mkl_sm_csr_frobenius_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRFrobeniusTest <double, std::uint64_t> mkl_sm_csr_frobenius_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRFrobeniusTest <__float128, std::uint32_t> sm_csr_frobenius_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest <__float128, std::uint64_t> sm_csr_frobenius_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRFrobeniusTest <Half, std::uint32_t> sm_csr_frobenius_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest <Half, std::uint64_t> sm_csr_frobenius_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRFrobeniusTest <float, std::uint32_t> cuda_sm_csr_frobenius_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRFrobeniusTest <double, std::uint32_t> cuda_sm_csr_frobenius_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRFrobeniusTest <float, std::uint64_t> cuda_sm_csr_frobenius_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRFrobeniusTest <double, std::uint64_t> cuda_sm_csr_frobenius_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRLumpTest
  : public UnitTest
{
public:
  SparseMatrixCSRLumpTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRLumpTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRLumpTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (IT_ size(2); size < IT_(3e2); size *= IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size);
      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
            a_fac.add(row, col, DT_(DT_(col % 100) / DT_(2)));
          else if ((row == col + 1) || (row + 1 == col))
            a_fac.add(row, col, DT_(-1));
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      auto lump = a.lump_rows();
      auto one = a.create_vector_r();
      auto res = a.create_vector_r();
      one.format(DT_(-1));
      a.apply(res, one, lump);

      TEST_CHECK_EQUAL_WITHIN_EPS(res.norm2(), DT_(0), tol);
    }
  }
};

SparseMatrixCSRLumpTest <float, std::uint32_t> sm_csr_lump_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRLumpTest <double, std::uint32_t> sm_csr_lump_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRLumpTest <float, std::uint64_t> sm_csr_lump_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRLumpTest <double, std::uint64_t> sm_csr_lump_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRLumpTest <float, std::uint64_t> mkl_sm_csr_lump_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRLumpTest <double, std::uint64_t> mkl_sm_csr_lump_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRLumpTest <__float128, std::uint32_t> sm_csr_lump_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRLumpTest <__float128, std::uint64_t> sm_csr_lump_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRLumpTest <Half, std::uint32_t> sm_csr_lump_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRLumpTest <Half, std::uint64_t> sm_csr_lump_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRLumpTest <float, std::uint32_t> cuda_sm_csr_lump_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRLumpTest <double, std::uint32_t> cuda_sm_csr_lump_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRLumpTest <float, std::uint64_t> cuda_sm_csr_lump_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRLumpTest <double, std::uint64_t> cuda_sm_csr_lump_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRCompressionTest
  : public UnitTest
{
public:
  SparseMatrixCSRCompressionTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRCompressionTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRCompressionTest()
  {
  }

  virtual void run() const override
  {
    Index mat_rows = 56;
    Index mat_cols = 56;
    DenseVector<IT_, IT_> col_ind(400);
    DenseVector<DT_, IT_> val(400);
    DenseVector<IT_, IT_> row_ptr(mat_rows + 1);
    for (Index i(0); i < 400; ++i)
    {
      col_ind(i, IT_(i / 40 + (3 * i) % 16));
      val(i, DT_(7) / DT_(3 * (i + 1)) - DT_(13) / DT_((i + 1) * (i + 1)));
    }
    row_ptr(0, 0);
    for (Index i(1); i < mat_rows; ++i)
    {
      row_ptr(i, IT_(i * Index(400 / mat_rows - 1)));
    }
    row_ptr(mat_rows, 400);
    SparseMatrixCSR<DT_, IT_> a(mat_rows, mat_cols, col_ind, val, row_ptr);

#ifdef FEAT_HAVE_ZLIB
#ifdef FEAT_HAVE_ZFP
    LAFEM::SerialConfig config(false, false);
    config.set_tolerance(FEAT::Real(1e-2));
    std::vector<char> uncompressed = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_zlib);
    std::vector<char> elements_compressed_zlib = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_zfp);
    std::vector<char> elements_compressed_zfp = a.serialize(config);
    config.set_indices_compression(LAFEM::CompressionModes::indices_zlib);
    std::vector<char> elements_indices_compressed_zfp = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_zlib);
    std::vector<char> elements_indices_compressed_zlib = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_off);
    std::vector<char> indices_compressed_zlib = a.serialize(config);
    config.set_elements_compression(LAFEM::CompressionModes::elements_zfp);
    config.set_indices_compression(LAFEM::CompressionModes::indices_off);
    config.set_tolerance(FEAT::Real(1e-4));
    std::vector<char> elements_compressed_zfp_e4 = a.serialize(config);
    config.set_tolerance(FEAT::Real(1e-6));
    std::vector<char> elements_compressed_zfp_e6 = a.serialize(config);

    XASSERTM(uncompressed.size() > elements_compressed_zlib.size(), "ele zlib is not smaller than uncomp");
    XASSERTM(elements_compressed_zlib.size() > elements_compressed_zfp.size(), "ele zfp is not smaleer than ele zlib");
    XASSERTM(elements_compressed_zfp.size() > elements_indices_compressed_zfp.size(), "ele + ind zfp is not smaller than ele zfp");
    XASSERTM(elements_compressed_zlib.size() > elements_indices_compressed_zlib.size(), "ele + ind zlib is not smaller than ele zlib");
    XASSERTM(uncompressed.size() > indices_compressed_zlib.size(), "ind zlib is not smaller than uncomp");
    XASSERTM(indices_compressed_zlib.size() > elements_indices_compressed_zlib.size(), "ele + ind zlib is not smaller than ind zlib");
    XASSERTM(elements_compressed_zfp_e4.size() > elements_compressed_zfp.size(), "e4 smaller than e2");
    XASSERTM(elements_compressed_zfp_e6.size() > elements_compressed_zfp_e4.size(), "e6 smaller than e4");
#endif
#endif
  }
};

SparseMatrixCSRCompressionTest <float, std::uint32_t> sm_csr_comp_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRCompressionTest <double, std::uint32_t> sm_csr_comp_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRCompressionTest <float, std::uint64_t> sm_csr_comp_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRCompressionTest <double, std::uint64_t> sm_csr_comp_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRCompressionTest <float, std::uint64_t> mkl_sm_csr_comp_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRCompressionTest <double, std::uint64_t> mkl_sm_csr_comp_test_double_uint64(PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SparseMatrixCSRCompressionTest <__float128, std::uint32_t> sm_csr_comp_test_float128_uint32(PreferredBackend::generic);
//SparseMatrixCSRCompressionTest <__float128, std::uint64_t> sm_csr_comp_test_float128_uint64(PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRCompressionTest <Half, std::uint32_t> sm_csr_comp_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRCompressionTest <Half, std::uint64_t> sm_csr_comp_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRCompressionTest <float, std::uint32_t> cuda_sm_csr_comp_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRCompressionTest <double, std::uint32_t> cuda_sm_csr_comp_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRCompressionTest <float, std::uint64_t> cuda_sm_csr_comp_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRCompressionTest <double, std::uint64_t> cuda_sm_csr_comp_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRMaxAbsElementTest
  : public UnitTest
{
public:
  SparseMatrixCSRMaxAbsElementTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRMaxAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRMaxAbsElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2); size < Index(3e2); size *= 2)
    {
      DT_ s(DT_(4.321));
      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> b_fac(IT_(size), IT_(size + 2));

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
            b_fac.add(row, col, DT_(s));
          }
          else if ((row == col + 1) || (row + 1 == col))
          {
            a_fac.add(row, col, DT_(-1));
            b_fac.add(row, col, DT_(-(s + s)));
          }
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      DT_ result;
      DT_ ref(2);
      result = a.max_abs_element();

      TEST_CHECK_EQUAL(result, ref);

      SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());
      ref = s + s;
      result = b.max_abs_element();

      TEST_CHECK_EQUAL(result, ref);
    }
  }
};

SparseMatrixCSRMaxAbsElementTest <float, std::uint32_t> sm_csr_max_abs_element_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest <double, std::uint32_t> sm_csr_max_abs_element_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest <float, std::uint64_t> sm_csr_max_abs_element_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest <double, std::uint64_t> sm_csr_max_abs_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMaxAbsElementTest <float, std::uint64_t> mkl_sm_csr_max_abs_element_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRMaxAbsElementTest <double, std::uint64_t> mkl_sm_csr_max_abs_element_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMaxAbsElementTest <__float128, std::uint32_t> sm_csr_max_abs_element_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest <__float128, std::uint64_t> sm_csr_max_abs_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMaxAbsElementTest <Half, std::uint32_t> sm_csr_max_abs_element_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest <Half, std::uint64_t> sm_csr_max_abs_element_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMaxAbsElementTest <float, std::uint32_t> cuda_sm_csr_max_abs_element_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRMaxAbsElementTest <double, std::uint32_t> cuda_sm_csr_max_abs_element_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRMaxAbsElementTest <float, std::uint64_t> cuda_sm_csr_max_abs_element_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRMaxAbsElementTest <double, std::uint64_t> cuda_sm_csr_max_abs_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRMinAbsElementTest
  : public UnitTest
{
public:
  SparseMatrixCSRMinAbsElementTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRMinAbsElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRMinAbsElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2); size < Index(3e2); size *= 2)
    {
      DT_ s(DT_(4.321));
      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> b_fac(IT_(size), IT_(size + 2));

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
            b_fac.add(row, col, DT_(s));
          }
          else if ((row == col + 1) || (row + 1 == col))
          {
            a_fac.add(row, col, DT_(-1));
            b_fac.add(row, col, DT_(-(s + s)));
          }
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      DT_ result;
      DT_ ref(1);
      result = a.min_abs_element();

      TEST_CHECK_EQUAL(result, ref);

      SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());
      ref = s;
      result = b.min_abs_element();

      TEST_CHECK_EQUAL(result, ref);
    }
  }
};

SparseMatrixCSRMinAbsElementTest <float, std::uint32_t> sm_csr_min_abs_element_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest <double, std::uint32_t> sm_csr_min_abs_element_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest <float, std::uint64_t> sm_csr_min_abs_element_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest <double, std::uint64_t> sm_csr_min_abs_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMinAbsElementTest <float, std::uint64_t> mkl_sm_csr_min_abs_element_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRMinAbsElementTest <double, std::uint64_t> mkl_sm_csr_min_abs_element_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMinAbsElementTest <__float128, std::uint32_t> sm_csr_min_abs_element_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest <__float128, std::uint64_t> sm_csr_min_abs_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMinAbsElementTest <Half, std::uint32_t> sm_csr_min_abs_element_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest <Half, std::uint64_t> sm_csr_min_abs_element_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMinAbsElementTest <float, std::uint32_t> cuda_sm_csr_min_abs_element_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRMinAbsElementTest <double, std::uint32_t> cuda_sm_csr_min_abs_element_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRMinAbsElementTest <float, std::uint64_t> cuda_sm_csr_min_abs_element_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRMinAbsElementTest <double, std::uint64_t> cuda_sm_csr_min_abs_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRMaxElementTest
  : public UnitTest
{
public:
  SparseMatrixCSRMaxElementTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRMaxElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRMaxElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2); size < Index(3e2); size *= 2)
    {
      DT_ s(DT_(4.321));
      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> b_fac(IT_(size), IT_(size + 2));

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
            b_fac.add(row, col, DT_(s));
          }
          else if ((row == col + 1) || (row + 1 == col))
          {
            a_fac.add(row, col, DT_(-1));
            b_fac.add(row, col, DT_(-(s + s)));
          }
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      DT_ result;
      DT_ ref(2);
      result = a.max_element();

      TEST_CHECK_EQUAL(result, ref);

      SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());
      ref = s;
      result = b.max_element();

      TEST_CHECK_EQUAL(result, ref);
    }
  }
};

SparseMatrixCSRMaxElementTest <float, std::uint32_t> sm_csr_max_element_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest <double, std::uint32_t> sm_csr_max_element_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest <float, std::uint64_t> sm_csr_max_element_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest <double, std::uint64_t> sm_csr_max_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMaxElementTest <float, std::uint64_t> mkl_sm_csr_max_element_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRMaxElementTest <double, std::uint64_t> mkl_sm_csr_max_element_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMaxElementTest <__float128, std::uint32_t> sm_csr_max_element_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest <__float128, std::uint64_t> sm_csr_max_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMaxElementTest <Half, std::uint32_t> sm_csr_max_element_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest <Half, std::uint64_t> sm_csr_max_element_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMaxElementTest <float, std::uint32_t> cuda_sm_csr_max_element_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRMaxElementTest <double, std::uint32_t> cuda_sm_csr_max_element_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRMaxElementTest <float, std::uint64_t> cuda_sm_csr_max_element_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRMaxElementTest <double, std::uint64_t> cuda_sm_csr_max_element_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class SparseMatrixCSRMinElementTest
  : public UnitTest
{
public:
  SparseMatrixCSRMinElementTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRMinElementTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRMinElementTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2); size < Index(3e2); size *= 2)
    {
      DT_ s(DT_(4.321));
      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> b_fac(IT_(size), IT_(size + 2));

      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
            b_fac.add(row, col, DT_(s));
          }
          else if ((row == col + 1) || (row + 1 == col))
          {
            a_fac.add(row, col, DT_(-1));
            b_fac.add(row, col, DT_(-(s + s)));
          }
        }
      }

      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      DT_ result;
      DT_ ref(-1);
      result = a.min_element();

      TEST_CHECK_EQUAL(result, ref);

      SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());
      ref = -(s + s);
      result = b.min_element();

      TEST_CHECK_EQUAL(result, ref);
    }
  }
};

SparseMatrixCSRMinElementTest <float, std::uint32_t> sm_csr_min_element_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRMinElementTest <double, std::uint32_t> sm_csr_min_element_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRMinElementTest <float, std::uint64_t> sm_csr_min_element_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRMinElementTest <double, std::uint64_t> sm_csr_min_element_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMinElementTest <float, std::uint64_t> mkl_sm_csr_min_element_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRMinElementTest <double, std::uint64_t> mkl_sm_csr_min_element_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMinElementTest <__float128, std::uint32_t> sm_csr_min_element_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRMinElementTest <__float128, std::uint64_t> sm_csr_min_element_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMinElementTest <Half, std::uint32_t> sm_csr_min_element_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRMinElementTest <Half, std::uint64_t> sm_csr_min_element_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMinElementTest <float, std::uint32_t> cuda_sm_csr_min_element_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRMinElementTest <double, std::uint32_t> cuda_sm_csr_min_element_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRMinElementTest <float, std::uint64_t> cuda_sm_csr_min_element_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRMinElementTest <double, std::uint64_t> cuda_sm_csr_min_element_test_double_uint64(PreferredBackend::cuda);
#endif


template<typename DT_, typename IT_>
class SparseMatrixCSRAddDoubleMatMultTest
  : public UnitTest
{
public:
  SparseMatrixCSRAddDoubleMatMultTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRAddDoubleMatMultTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSRAddDoubleMatMultTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.7));
    SparseMatrixFactory<DT_, IT_> b_fac(IT_(4), IT_(3));
    SparseMatrixFactory<DT_, IT_> d_fac(IT_(3), IT_(4));
    SparseMatrixFactory<DT_, IT_> s_fac(IT_(3), IT_(3));
    SparseMatrixFactory<DT_, IT_> s_fac2(IT_(3), IT_(3));

    b_fac.add(IT_(0), IT_(0), DT_(2));
    b_fac.add(IT_(0), IT_(1), DT_(3));
    b_fac.add(IT_(1), IT_(0), DT_(4));
    b_fac.add(IT_(1), IT_(2), DT_(5));
    b_fac.add(IT_(2), IT_(1), DT_(6));
    b_fac.add(IT_(2), IT_(2), DT_(7));
    b_fac.add(IT_(3), IT_(1), DT_(8));
    d_fac.add(IT_(0), IT_(0), DT_(5));
    d_fac.add(IT_(0), IT_(1), DT_(7));
    d_fac.add(IT_(1), IT_(0), DT_(3));
    d_fac.add(IT_(1), IT_(2), DT_(2));
    d_fac.add(IT_(2), IT_(1), DT_(8));
    d_fac.add(IT_(2), IT_(2), DT_(6));
    d_fac.add(IT_(2), IT_(3), DT_(4));
    s_fac.add(IT_(0), IT_(0), DT_(104));
    s_fac.add(IT_(0), IT_(1), DT_(30));
    s_fac.add(IT_(0), IT_(2), DT_(105));
    s_fac.add(IT_(1), IT_(0), DT_(12));
    s_fac.add(IT_(1), IT_(1), DT_(66));
    s_fac.add(IT_(1), IT_(2), DT_(56));
    s_fac.add(IT_(2), IT_(0), DT_(96));
    s_fac.add(IT_(2), IT_(1), DT_(304));
    s_fac.add(IT_(2), IT_(2), DT_(288));
    s_fac2.add(IT_(0), IT_(0), DT_(696));
    s_fac2.add(IT_(0), IT_(1), DT_(120));
    s_fac2.add(IT_(0), IT_(2), DT_(770));
    s_fac2.add(IT_(1), IT_(0), DT_(48));
    s_fac2.add(IT_(1), IT_(1), DT_(504));
    s_fac2.add(IT_(1), IT_(2), DT_(504));
    s_fac2.add(IT_(2), IT_(0), DT_(704));
    s_fac2.add(IT_(2), IT_(1), DT_(2896));
    s_fac2.add(IT_(2), IT_(2), DT_(2392));

    SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());
    SparseMatrixCSR<DT_, IT_> d(d_fac.make_csr());
    SparseMatrixCSR<DT_, IT_> s(s_fac.make_csr());

    DenseVector<DT_, IT_> a(4u);
    a.elements()[0] = DT_(2);
    a.elements()[1] = DT_(3);
    a.elements()[2] = DT_(4);
    a.elements()[3] = DT_(5);

    // perform double matrix product and check result
    s.add_double_mat_product(d, a, b, -DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(s.norm_frobenius(), DT_(0), tol);

    // convert matrices to blocked variants
    SparseMatrixBCSR<DT_, IT_, 2, 1> b2(b.layout());
    SparseMatrixBCSR<DT_, IT_, 1, 2> d2(d.layout());
    SparseMatrixCSR<DT_, IT_> s2(s_fac2.make_csr());
    for(Index i(0); i < b.used_elements();  ++i)
      b2.val()[i][0][0] = DT_(2) * (b2.val()[i][1][0] = b.val()[i]);
    for(Index i(0); i < d.used_elements();  ++i)
      d2.val()[i][0][0] = DT_(3) * (d2.val()[i][0][1] = d.val()[i]);

    DenseVectorBlocked<DT_, IT_, 2> a2(4u);
    a2.elements()[0][0] = DT_(1);
    a2.elements()[0][1] = DT_(2);
    a2.elements()[1][0] = DT_(3);
    a2.elements()[1][1] = DT_(4);
    a2.elements()[2][0] = DT_(5);
    a2.elements()[2][1] = DT_(6);
    a2.elements()[3][0] = DT_(7);
    a2.elements()[3][1] = DT_(8);

    s2.add_double_mat_product(d2, a2, b2, -DT_(1));
    TEST_CHECK_EQUAL_WITHIN_EPS(s2.norm_frobenius(), DT_(0), tol);
  }
};

SparseMatrixCSRAddDoubleMatMultTest <float, std::uint32_t> sm_csr_add_double_mat_mult_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRAddDoubleMatMultTest <double, std::uint32_t> sm_csr_add_double_mat_mult_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRAddDoubleMatMultTest <float, std::uint64_t> sm_csr_add_double_mat_mult_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRAddDoubleMatMultTest <double, std::uint64_t> sm_csr_add_double_mat_mult_test_double_uint64(PreferredBackend::generic);

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSRMaxRelDiffTest
  : public UnitTest
{
public:
  SparseMatrixCSRMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSRMaxRelDiffTest",
               Type::Traits<DT_>::name(),
               Type::Traits<IT_>::name(),
               backend)
  {
  }

  virtual ~SparseMatrixCSRMaxRelDiffTest() {}

  virtual void run() const override
  {
    const DT_ eps   = Math::pow(Math::eps<DT_>(), DT_(0.8));
    const DT_ delta = DT_(123.5);

    for (Index size(4); size < Index(128); size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> fac_a(size, size);
      for (IT_ i = 0; i < IT_(size); ++i)
        fac_a.add(i, i, DT_(i));
      SparseMatrixCSR<DT_, IT_> a(fac_a.make_csr());

      auto b = a.clone();

      // add delta
      a.val()[size / 2] += delta;

      const DT_ ref = delta / (DT_(size) + delta);

      // ||a-b||_inf
      const DT_ diff1 = a.max_rel_diff(b);
      TEST_CHECK_RELATIVE(diff1, ref, eps);

      // ||b-a||_inf
      const DT_ diff2 = b.max_rel_diff(a);
      TEST_CHECK_RELATIVE(diff2, ref, eps);
    }
  }
};
SparseMatrixCSRMaxRelDiffTest <float,  std::uint32_t> sm_csr_max_rel_diff_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxRelDiffTest <double, std::uint32_t> sm_csr_max_rel_diff_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxRelDiffTest <float,  std::uint64_t> sm_csr_max_rel_diff_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSRMaxRelDiffTest <double, std::uint64_t> sm_csr_max_rel_diff_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMaxRelDiffTest <float,  std::uint64_t> mkl_sm_csr_max_rel_diff_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSRMaxRelDiffTest <double, std::uint64_t> mkl_sm_csr_max_rel_diff_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMaxRelDiffTest <__float128, std::uint32_t> sm_csr_max_rel_diff_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxRelDiffTest <__float128, std::uint64_t> sm_csr_max_rel_diff_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMaxRelDiffTest <Half, std::uint32_t> sm_csr_max_rel_diff_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSRMaxRelDiffTest <Half, std::uint64_t> sm_csr_max_rel_diff_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMaxRelDiffTest <float,  std::uint32_t> cuda_sm_csr_max_rel_diff_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSRMaxRelDiffTest <double, std::uint32_t> cuda_sm_csr_max_rel_diff_test_double_uint32(PreferredBackend::cuda);
SparseMatrixCSRMaxRelDiffTest <float,  std::uint64_t> cuda_sm_csr_max_rel_diff_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSRMaxRelDiffTest <double, std::uint64_t> cuda_sm_csr_max_rel_diff_test_double_uint64(PreferredBackend::cuda);
#endif
