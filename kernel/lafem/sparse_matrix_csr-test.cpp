// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
    SparseMatrixCSR<DT_, IT_> zero1;
    SparseMatrixCSR<DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);
    zero2.convert(zero1);

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
    empty5.convert(zero1);
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
    TEST_CHECK_EQUAL(d, c);

    SparseMatrixCSR<DT_, IT_> e;
    e.convert(c);
    TEST_CHECK_EQUAL(e, c);
    e.copy(c);
    TEST_CHECK_EQUAL(e, c);
    e.clone(c);
    b.clone(e);
    TEST_CHECK_EQUAL(b, c);


    // new clone testing
    auto clone1 = b.clone(CloneMode::Deep);
    TEST_CHECK_EQUAL(clone1, b);
    MemoryPool::set_memory(clone1.val() + 1, DT_(132));
    TEST_CHECK_NOT_EQUAL(clone1, b);
    TEST_CHECK_NOT_EQUAL((void*)clone1.val(), (void*)b.val());
    TEST_CHECK_NOT_EQUAL((void*)clone1.row_ptr(), (void*)b.row_ptr());
    auto clone2 = clone1.clone(CloneMode::Layout);
    MemoryPool::set_memory(clone2.val(), DT_(4713), clone2.used_elements());
    TEST_CHECK_NOT_EQUAL(clone2(5, 5), clone1(5, 5));
    TEST_CHECK_NOT_EQUAL((void*)clone2.val(), (void*)clone1.val());
    TEST_CHECK_EQUAL((void*)clone2.row_ptr(), (void*)clone1.row_ptr());
    auto clone3 = clone1.clone(CloneMode::Weak);
    TEST_CHECK_EQUAL(clone3, clone1);
    MemoryPool::set_memory(clone3.val() + 1, DT_(133));
    TEST_CHECK_NOT_EQUAL(clone3, clone1);
    TEST_CHECK_NOT_EQUAL((void*)clone3.val(), (void*)clone1.val());
    TEST_CHECK_EQUAL((void*)clone3.row_ptr(), (void*)clone1.row_ptr());
    auto clone4 = clone1.clone(CloneMode::Shallow);
    TEST_CHECK_EQUAL(clone4, clone1);
    MemoryPool::set_memory(clone4.val() + 1, DT_(134));
    TEST_CHECK_EQUAL(clone4, clone1);
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

SparseMatrixCSRTest<float, unsigned long> cpu_sparse_matrix_csr_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRTest<double, unsigned long> cpu_sparse_matrix_csr_test_double_ulong(PreferredBackend::generic);
SparseMatrixCSRTest<float, unsigned int> cpu_sparse_matrix_csr_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRTest<double, unsigned int> cpu_sparse_matrix_csr_test_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRTest<float, unsigned long> mkl_cpu_sparse_matrix_csr_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRTest<double, unsigned long> mkl_cpu_sparse_matrix_csr_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRTest<__float128, unsigned long> cpu_sparse_matrix_csr_test_float128_ulong(PreferredBackend::generic);
SparseMatrixCSRTest<__float128, unsigned int> cpu_sparse_matrix_csr_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRTest<Half, unsigned int> sparse_matrix_csr_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRTest<Half, unsigned long> sparse_matrix_csr_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRTest<float, unsigned long> cuda_sparse_matrix_csr_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRTest<double, unsigned long> cuda_sparse_matrix_csr_test_double_ulong(PreferredBackend::cuda);
SparseMatrixCSRTest<float, unsigned int> cuda_sparse_matrix_csr_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRTest<double, unsigned int> cuda_sparse_matrix_csr_test_double_uint(PreferredBackend::cuda);
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
    TEST_CHECK_EQUAL(g, f);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));

    std::stringstream ts;
    f.write_out(FileMode::fm_mtx, ts);
    SparseMatrixCSR<DT_, IT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_EQUAL(j, f);

    std::stringstream ts2;
    f.write_out(FileMode::fm_mtx, ts2, true);
    SparseMatrixCSR<DT_, IT_> j2(FileMode::fm_mtx, ts2);
    TEST_CHECK_EQUAL(j2, f);

    auto kp = f.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixCSR<DT_, IT_> k(kp);
    TEST_CHECK_EQUAL(k, f);
#ifdef FEAT_HAVE_ZLIB
    auto zl = f.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixCSR<DT_, IT_> zlib(zl);
    TEST_CHECK_EQUAL(zlib, f);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = f.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixCSR<DT_, IT_> zfp(zf);
    for (Index row(0); row < f.rows(); ++row)
    {
      for (Index col(0); col < f.columns(); ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(zfp(row, col), f(row, col), DT_(1e-4));
      }
    }
#endif
  }
};
SparseMatrixCSRSerializeTest<float, unsigned long> cpu_sparse_matrix_csr_serialize_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRSerializeTest<double, unsigned long> cpu_sparse_matrix_csr_serialize_test_double_ulong(PreferredBackend::generic);
SparseMatrixCSRSerializeTest<float, unsigned int> cpu_sparse_matrix_csr_serialize_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRSerializeTest<double, unsigned int> cpu_sparse_matrix_csr_serialize_test_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRSerializeTest<float, unsigned long> mkl_cpu_sparse_matrix_csr_serialize_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRSerializeTest<double, unsigned long> mkl_cpu_sparse_matrix_csr_serialize_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRSerializeTest<__float128, unsigned long> cpu_sparse_matrix_csr_serialize_test_float128_ulong(PreferredBackend::generic);
SparseMatrixCSRSerializeTest<__float128, unsigned int> cpu_sparse_matrix_csr_serialize_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRSerializeTest<Half, unsigned int> sparse_matrix_csr_serialize_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRSerializeTest<Half, unsigned long> sparse_matrix_csr_serialize_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRSerializeTest<float, unsigned long> cuda_sparse_matrix_csr_serialize_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRSerializeTest<double, unsigned long> cuda_sparse_matrix_csr_serialize_test_double_ulong(PreferredBackend::cuda);
SparseMatrixCSRSerializeTest<float, unsigned int> cuda_sparse_matrix_csr_serialize_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRSerializeTest<double, unsigned int> cuda_sparse_matrix_csr_serialize_test_double_uint(PreferredBackend::cuda);
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
    DT_ s(DT_(4711.1));
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
        y(i, DT_(2 - DT_(i % 42)));
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
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-2));

      // apply-test for alpha = -1.0
      a.apply(r, x, y, DT_(-1.0));
      a.apply(ref, x);
      ref.scale(ref, DT_(-1.0));
      ref.axpy(ref, y);
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-2));

      // apply-test for alpha = -1.0 and &r==&y
      r.copy(y);
      a.apply(r, x, r, DT_(-1.0));
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-2));

      // apply-test for alpha = 4711.1
      //r.axpy(s, a, x, y);
      a.apply(r, x, y, s);
      //ref.product_matvec(a, x);
      a.apply(ref, x);
      ref.scale(ref, s);
      ref.axpy(ref, y);
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(5e-2));

      // apply-test for alpha = 4711.1 and &r==&y
      r.copy(y);
      a.apply(r, x, r, s);
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(5e-2));

      a.apply(r, x);
      ref.copy(ax);
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-2));

      // transposed apply-test for alpha = 4711.1
      if (!(Runtime::get_preferred_backend() == PreferredBackend::cuda && typeid(IT_) == typeid(unsigned long)))
      {
        a.apply(r, x, y, s, true);

        SparseMatrixCSR<DT_, IT_> at = a.transpose();
        at.apply(ref, x);
        ref.scale(ref, s);
        ref.axpy(ref, y);
        MemoryPool::synchronize();
        for (Index i(0); i < size; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(7e-2));
      }
    }
  }
};

SparseMatrixCSRApplyTest<float, unsigned long> sm_csr_apply_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRApplyTest<double, unsigned long> sm_csr_apply_test_double_ulong(PreferredBackend::generic);
SparseMatrixCSRApplyTest<float, unsigned int> sm_csr_apply_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRApplyTest<double, unsigned int> sm_csr_apply_test_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRApplyTest<float, unsigned long> mkl_sm_csr_apply_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRApplyTest<double, unsigned long> mkl_sm_csr_apply_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRApplyTest<__float128, unsigned long> sm_csr_apply_test_float128_ulong(PreferredBackend::generic);
SparseMatrixCSRApplyTest<__float128, unsigned int> sm_csr_apply_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRApplyTest<Half, unsigned int> sm_csr_apply_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRApplyTest<Half, unsigned long> sm_csr_apply_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRApplyTest<float, unsigned long> cuda_sm_csr_apply_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRApplyTest<double, unsigned long> cuda_sm_csr_apply_test_double_ulong(PreferredBackend::cuda);
SparseMatrixCSRApplyTest<float, unsigned int> cuda_sm_csr_apply_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRApplyTest<double, unsigned int> cuda_sm_csr_apply_test_double_uint(PreferredBackend::cuda);
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRApplyTest<Half, unsigned int> cuda_sm_csr_apply_test_half_uint(PreferredBackend::cuda);
SparseMatrixCSRApplyTest<Half, unsigned long> cuda_sm_csr_apply_test_half_ulong(PreferredBackend::cuda);
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
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[0], aref_x(i), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[1], aref_x(i) * DT_(0.5), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[2], aref_x(i) * DT_(2.0), DT_(1e-4));
      }

      a.apply(r, x, y, DT_(-1));
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[0], y(i)[0] - aref_x(i), DT_(1e-4));
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[1], y(i)[1] - aref_x(i) * DT_(0.5), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[2], y(i)[2] - aref_x(i) * DT_(2.0), DT_(1e-4));
      }

      DT_ alpha(0.75);
      a.apply(r, x, y, alpha);
      MemoryPool::synchronize();
      for (Index i(0); i < size; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[0], y(i)[0] + alpha * aref_x(i), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[1], y(i)[1] + alpha * aref_x(i) * DT_(0.5), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i)[2], y(i)[2] + alpha * aref_x(i) * DT_(2.0), DT_(1e-4));
      }
    }
  }
};

SparseMatrixCSRBApplyTest<float, unsigned long> sm_csrib_apply_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRBApplyTest<double, unsigned long> sm_csrib_apply_test_double_ulong(PreferredBackend::generic);
SparseMatrixCSRBApplyTest<float, unsigned int> sm_csrib_apply_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRBApplyTest<double, unsigned int> sm_csrib_apply_test_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRBApplyTest<float, unsigned long> mkl_sm_csrib_apply_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRBApplyTest<double, unsigned long> mkl_sm_csrib_apply_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRBApplyTest<__float128, unsigned long> sm_csrib_apply_test_float128_ulong(PreferredBackend::generic);
SparseMatrixCSRBApplyTest<__float128, unsigned int> sm_csrib_apply_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRBApplyTest<Half, unsigned int> sm_csrib_apply_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRBApplyTest<Half, unsigned long> sm_csrib_apply_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRBApplyTest<float, unsigned long> cuda_sm_csrib_apply_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRBApplyTest<double, unsigned long> cuda_sm_csrib_apply_test_double_ulong(PreferredBackend::cuda);
SparseMatrixCSRBApplyTest<float, unsigned int> cuda_sm_csrib_apply_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRBApplyTest<double, unsigned int> cuda_sm_csrib_apply_test_double_uint(PreferredBackend::cuda);
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
      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(b, ref);

      a.scale(a, s);
      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(a, ref);
    }
  }
};

SparseMatrixCSRScaleTest<float, unsigned int> sm_csr_scale_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRScaleTest<double, unsigned int> sm_csr_scale_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRScaleTest<float, unsigned long> sm_csr_scale_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRScaleTest<double, unsigned long> sm_csr_scale_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRScaleTest<float, unsigned long> mkl_sm_csr_scale_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRScaleTest<double, unsigned long> mkl_sm_csr_scale_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRScaleTest<__float128, unsigned int> sm_csr_scale_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRScaleTest<__float128, unsigned long> sm_csr_scale_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRScaleTest<Half, unsigned int> sm_csr_scale_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRScaleTest<Half, unsigned long> sm_csr_scale_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRScaleTest<float, unsigned int> cuda_sm_csr_scale_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRScaleTest<double, unsigned int> cuda_sm_csr_scale_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRScaleTest<float, unsigned long> cuda_sm_csr_scale_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRScaleTest<double, unsigned long> cuda_sm_csr_scale_test_double_ulong(PreferredBackend::cuda);
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
        s1(i, pi * (i % 3 + 1) - DT_(5.21) + DT_(i));
      }
      b.scale_rows(b, s1);
      MemoryPool::synchronize();
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
        s2(i, pi * (i % 3 + 1) - DT_(5.21) + DT_(i));
      }
      b.scale_cols(a, s2);
      MemoryPool::synchronize();
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

SparseMatrixCSRScaleRowColTest<float, unsigned int> sm_csr_scale_row_col_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest<double, unsigned int> sm_csr_scale_row_col_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest<float, unsigned long> sm_csr_scale_row_col_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest<double, unsigned long> sm_csr_scale_row_col_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRScaleRowColTest<float, unsigned long> mkl_sm_csr_scale_row_col_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRScaleRowColTest<double, unsigned long> mkl_sm_csr_scale_row_col_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRScaleRowColTest<__float128, unsigned int> sm_csr_scale_row_col_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest<__float128, unsigned long> sm_csr_scale_row_col_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRScaleRowColTest<Half, unsigned int> sm_csr_scale_row_col_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRScaleRowColTest<Half, unsigned long> sm_csr_scale_row_col_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRScaleRowColTest<float, unsigned int> cuda_sm_csr_scale_row_col_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRScaleRowColTest<double, unsigned int> cuda_sm_csr_scale_row_col_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRScaleRowColTest<float, unsigned long> cuda_sm_csr_scale_row_col_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRScaleRowColTest<double, unsigned long> cuda_sm_csr_scale_row_col_test_double_ulong(PreferredBackend::cuda);
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

      MemoryPool::synchronize();
      for (Index i(0); i < a.rows(); ++i)
      {
        for (Index j(0); j < a.columns(); ++j)
        {
          TEST_CHECK_EQUAL(b(j, i), a(i, j));
        }
      }

      b = b.transpose();

      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(a, b);
    }
  }
};

SparseMatrixCSRTranspositionTest<float, unsigned int> sm_csr_transposition_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest<double, unsigned int> sm_csr_transposition_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest<float, unsigned long> sm_csr_transposition_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest<double, unsigned long> sm_csr_transposition_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRTranspositionTest<float, unsigned long> mkl_sm_csr_transposition_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRTranspositionTest<double, unsigned long> mkl_sm_csr_transposition_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRTranspositionTest<__float128, unsigned int> sm_csr_transposition_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest<__float128, unsigned long> sm_csr_transposition_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRTranspositionTest<Half, unsigned int> sm_csr_transposition_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRTranspositionTest<Half, unsigned long> sm_csr_transposition_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRTranspositionTest<float, unsigned int> cuda_sm_csr_transposition_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRTranspositionTest<double, unsigned int> cuda_sm_csr_transposition_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRTranspositionTest<float, unsigned long> cuda_sm_csr_transposition_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRTranspositionTest<double, unsigned long> cuda_sm_csr_transposition_test_double_ulong(PreferredBackend::cuda);
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

      Random::SeedType seed(Random::SeedType(time(nullptr)));
      std::cout << "seed: " << seed << std::endl;
      Random rng(seed);
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
      TEST_CHECK_EQUAL(a, a_backup);
    }
  }
};

SparseMatrixCSRPermuteTest<float, unsigned long> sm_csr_permute_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRPermuteTest<double, unsigned long> sm_csr_permute_test_double_ulong(PreferredBackend::generic);
SparseMatrixCSRPermuteTest<float, unsigned int> sm_csr_permute_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRPermuteTest<double, unsigned int> sm_csr_permute_test_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRPermuteTest<float, unsigned long> mkl_sm_csr_permute_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRPermuteTest<double, unsigned long> mkl_sm_csr_permute_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRPermuteTest<__float128, unsigned long> sm_csr_permute_test_float128_ulong(PreferredBackend::generic);
SparseMatrixCSRPermuteTest<__float128, unsigned int> sm_csr_permute_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRPermuteTest<Half, unsigned int> sm_csr_permute_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRPermuteTest<Half, unsigned long> sm_csr_permute_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRPermuteTest<float, unsigned long> cuda_sm_csr_permute_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRPermuteTest<double, unsigned long> cuda_sm_csr_permute_test_double_ulong(PreferredBackend::cuda);
SparseMatrixCSRPermuteTest<float, unsigned int> cuda_sm_csr_permute_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRPermuteTest<double, unsigned int> cuda_sm_csr_permute_test_double_uint(PreferredBackend::cuda);
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
      TEST_CHECK_EQUAL(diag, ref);
    }
  }
};

SparseMatrixCSRDiagTest<float, unsigned int> sm_csr_diag_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRDiagTest<double, unsigned int> sm_csr_diag_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRDiagTest<float, unsigned long> sm_csr_diag_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRDiagTest<double, unsigned long> sm_csr_diag_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRDiagTest<float, unsigned long> mkl_sm_csr_diag_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRDiagTest<double, unsigned long> mkl_sm_csr_diag_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRDiagTest<__float128, unsigned int> sm_csr_diag_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRDiagTest<__float128, unsigned long> sm_csr_diag_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRDiagTest<Half, unsigned int> sm_csr_diag_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRDiagTest<Half, unsigned long> sm_csr_diag_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRDiagTest<float, unsigned int> cuda_sm_csr_diag_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRDiagTest<double, unsigned int> cuda_sm_csr_diag_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRDiagTest<float, unsigned long> cuda_sm_csr_diag_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRDiagTest<double, unsigned long> cuda_sm_csr_diag_test_double_ulong(PreferredBackend::cuda);
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
    for (Index size(2); size < Index(3e2); size *= 2)
    {
      DT_ s(DT_(4.321));

      SparseMatrixFactory<DT_, IT_> a_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> b_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_fac(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_fac_plus(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_fac_minus(IT_(size), IT_(size + 2));
      for (IT_ row(0); row < a_fac.rows(); ++row)
      {
        for (IT_ col(0); col < a_fac.columns(); ++col)
        {
          if (row == col)
          {
            a_fac.add(row, col, DT_(2));
            b_fac.add(row, col, DT_((row + col) % 15));
            ref_fac.add(row, col, DT_(2) * s + DT_((row + col) % 15));
            ref_fac_plus.add(row, col, DT_(2) + DT_((row + col) % 15));
            ref_fac_minus.add(row, col, DT_((row + col) % 15) - DT_(2));
          }

          else if ((row == col + 1) || (row + 1 == col))
          {
            a_fac.add(row, col, DT_(-1));
            b_fac.add(row, col, DT_((row + col) % 15));
            ref_fac.add(row, col, DT_(-1) * s + DT_((row + col) % 15));
            ref_fac_plus.add(row, col, DT_(-1) + DT_((row + col) % 15));
            ref_fac_minus.add(row, col, DT_((row + col) % 15) - DT_(-1));
          }
        }
      }

      SparseMatrixCSR<DT_, IT_> ref(ref_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> ref_plus(ref_fac_plus.make_csr());
      SparseMatrixCSR<DT_, IT_> ref_minus(ref_fac_minus.make_csr());
      SparseMatrixCSR<DT_, IT_> a(a_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> b(b_fac.make_csr());
      SparseMatrixCSR<DT_, IT_> c;

      c.clone(a);
      c.axpy(c, b, s);
      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(c, ref);

      c.clone(b);
      c.axpy(a, c, s);
      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(c, ref);

      c.axpy(a, b, s);
      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(c, ref);

      s = DT_(0);
      ref.clone(b);
      c.axpy(a, b, s);
      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(c, ref);

      s = DT_(1);
      c.axpy(a, b, s);
      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(c, ref_plus);

      s = DT_(-1);
      c.axpy(a, b, s);
      MemoryPool::synchronize();
      TEST_CHECK_EQUAL(c, ref_minus);
    }
  }
};

SparseMatrixCSRAxpyTest<float, unsigned int> sm_csr_axpy_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRAxpyTest<double, unsigned int> sm_csr_axpy_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRAxpyTest<float, unsigned long> sm_csr_axpy_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRAxpyTest<double, unsigned long> sm_csr_axpy_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRAxpyTest<float, unsigned long> mkl_sm_csr_axpy_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRAxpyTest<double, unsigned long> mkl_sm_csr_axpy_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRAxpyTest<__float128, unsigned int> sm_csr_axpy_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRAxpyTest<__float128, unsigned long> sm_csr_axpy_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRAxpyTest<Half, unsigned int> sm_csr_axpy_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRAxpyTest<Half, unsigned long> sm_csr_axpy_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRAxpyTest<float, unsigned int> cuda_sm_csr_axpy_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRAxpyTest<double, unsigned int> cuda_sm_csr_axpy_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRAxpyTest<float, unsigned long> cuda_sm_csr_axpy_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRAxpyTest<double, unsigned long> cuda_sm_csr_axpy_test_double_ulong(PreferredBackend::cuda);
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

SparseMatrixCSRFrobeniusTest<float, unsigned int> sm_csr_frobenius_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest<double, unsigned int> sm_csr_frobenius_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest<float, unsigned long> sm_csr_frobenius_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest<double, unsigned long> sm_csr_frobenius_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRFrobeniusTest<float, unsigned long> mkl_sm_csr_frobenius_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRFrobeniusTest<double, unsigned long> mkl_sm_csr_frobenius_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRFrobeniusTest<__float128, unsigned int> sm_csr_frobenius_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest<__float128, unsigned long> sm_csr_frobenius_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRFrobeniusTest<Half, unsigned int> sm_csr_frobenius_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRFrobeniusTest<Half, unsigned long> sm_csr_frobenius_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRFrobeniusTest<float, unsigned int> cuda_sm_csr_frobenius_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRFrobeniusTest<double, unsigned int> cuda_sm_csr_frobenius_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRFrobeniusTest<float, unsigned long> cuda_sm_csr_frobenius_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRFrobeniusTest<double, unsigned long> cuda_sm_csr_frobenius_test_double_ulong(PreferredBackend::cuda);
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

SparseMatrixCSRLumpTest<float, unsigned int> sm_csr_lump_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRLumpTest<double, unsigned int> sm_csr_lump_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRLumpTest<float, unsigned long> sm_csr_lump_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRLumpTest<double, unsigned long> sm_csr_lump_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRLumpTest<float, unsigned long> mkl_sm_csr_lump_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRLumpTest<double, unsigned long> mkl_sm_csr_lump_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRLumpTest<__float128, unsigned int> sm_csr_lump_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRLumpTest<__float128, unsigned long> sm_csr_lump_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRLumpTest<Half, unsigned int> sm_csr_lump_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRLumpTest<Half, unsigned long> sm_csr_lump_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRLumpTest<float, unsigned int> cuda_sm_csr_lump_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRLumpTest<double, unsigned int> cuda_sm_csr_lump_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRLumpTest<float, unsigned long> cuda_sm_csr_lump_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRLumpTest<double, unsigned long> cuda_sm_csr_lump_test_double_ulong(PreferredBackend::cuda);
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

SparseMatrixCSRCompressionTest<float, unsigned int> sm_csr_comp_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRCompressionTest<double, unsigned int> sm_csr_comp_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRCompressionTest<float, unsigned long> sm_csr_comp_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRCompressionTest<double, unsigned long> sm_csr_comp_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRCompressionTest<float, unsigned long> mkl_sm_csr_comp_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRCompressionTest<double, unsigned long> mkl_sm_csr_comp_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRCompressionTest<__float128, unsigned int> sm_csr_comp_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRCompressionTest<__float128, unsigned long> sm_csr_comp_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRCompressionTest<Half, unsigned int> sm_csr_comp_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRCompressionTest<Half, unsigned long> sm_csr_comp_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRCompressionTest<float, unsigned int> cuda_sm_csr_comp_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRCompressionTest<double, unsigned int> cuda_sm_csr_comp_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRCompressionTest<float, unsigned long> cuda_sm_csr_comp_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRCompressionTest<double, unsigned long> cuda_sm_csr_comp_test_double_ulong(PreferredBackend::cuda);
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

SparseMatrixCSRMaxAbsElementTest<float, unsigned int> sm_csr_max_abs_element_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest<double, unsigned int> sm_csr_max_abs_element_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest<float, unsigned long> sm_csr_max_abs_element_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest<double, unsigned long> sm_csr_max_abs_element_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMaxAbsElementTest<float, unsigned long> mkl_sm_csr_max_abs_element_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRMaxAbsElementTest<double, unsigned long> mkl_sm_csr_max_abs_element_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMaxAbsElementTest<__float128, unsigned int> sm_csr_max_abs_element_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest<__float128, unsigned long> sm_csr_max_abs_element_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMaxAbsElementTest<Half, unsigned int> sm_csr_max_abs_element_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRMaxAbsElementTest<Half, unsigned long> sm_csr_max_abs_element_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMaxAbsElementTest<float, unsigned int> cuda_sm_csr_max_abs_element_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRMaxAbsElementTest<double, unsigned int> cuda_sm_csr_max_abs_element_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRMaxAbsElementTest<float, unsigned long> cuda_sm_csr_max_abs_element_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRMaxAbsElementTest<double, unsigned long> cuda_sm_csr_max_abs_element_test_double_ulong(PreferredBackend::cuda);
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

SparseMatrixCSRMinAbsElementTest<float, unsigned int> sm_csr_min_abs_element_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest<double, unsigned int> sm_csr_min_abs_element_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest<float, unsigned long> sm_csr_min_abs_element_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest<double, unsigned long> sm_csr_min_abs_element_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMinAbsElementTest<float, unsigned long> mkl_sm_csr_min_abs_element_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRMinAbsElementTest<double, unsigned long> mkl_sm_csr_min_abs_element_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMinAbsElementTest<__float128, unsigned int> sm_csr_min_abs_element_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest<__float128, unsigned long> sm_csr_min_abs_element_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMinAbsElementTest<Half, unsigned int> sm_csr_min_abs_element_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRMinAbsElementTest<Half, unsigned long> sm_csr_min_abs_element_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMinAbsElementTest<float, unsigned int> cuda_sm_csr_min_abs_element_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRMinAbsElementTest<double, unsigned int> cuda_sm_csr_min_abs_element_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRMinAbsElementTest<float, unsigned long> cuda_sm_csr_min_abs_element_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRMinAbsElementTest<double, unsigned long> cuda_sm_csr_min_abs_element_test_double_ulong(PreferredBackend::cuda);
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

SparseMatrixCSRMaxElementTest<float, unsigned int> sm_csr_max_element_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest<double, unsigned int> sm_csr_max_element_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest<float, unsigned long> sm_csr_max_element_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest<double, unsigned long> sm_csr_max_element_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMaxElementTest<float, unsigned long> mkl_sm_csr_max_element_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRMaxElementTest<double, unsigned long> mkl_sm_csr_max_element_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMaxElementTest<__float128, unsigned int> sm_csr_max_element_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest<__float128, unsigned long> sm_csr_max_element_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMaxElementTest<Half, unsigned int> sm_csr_max_element_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRMaxElementTest<Half, unsigned long> sm_csr_max_element_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMaxElementTest<float, unsigned int> cuda_sm_csr_max_element_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRMaxElementTest<double, unsigned int> cuda_sm_csr_max_element_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRMaxElementTest<float, unsigned long> cuda_sm_csr_max_element_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRMaxElementTest<double, unsigned long> cuda_sm_csr_max_element_test_double_ulong(PreferredBackend::cuda);
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

SparseMatrixCSRMinElementTest<float, unsigned int> sm_csr_min_element_test_float_uint(PreferredBackend::generic);
SparseMatrixCSRMinElementTest<double, unsigned int> sm_csr_min_element_test_double_uint(PreferredBackend::generic);
SparseMatrixCSRMinElementTest<float, unsigned long> sm_csr_min_element_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSRMinElementTest<double, unsigned long> sm_csr_min_element_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSRMinElementTest<float, unsigned long> mkl_sm_csr_min_element_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSRMinElementTest<double, unsigned long> mkl_sm_csr_min_element_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRMinElementTest<__float128, unsigned int> sm_csr_min_element_test_float128_uint(PreferredBackend::generic);
SparseMatrixCSRMinElementTest<__float128, unsigned long> sm_csr_min_element_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSRMinElementTest<Half, unsigned int> sm_csr_min_element_test_half_uint(PreferredBackend::generic);
SparseMatrixCSRMinElementTest<Half, unsigned long> sm_csr_min_element_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRMinElementTest<float, unsigned int> cuda_sm_csr_min_element_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSRMinElementTest<double, unsigned int> cuda_sm_csr_min_element_test_double_uint(PreferredBackend::cuda);
SparseMatrixCSRMinElementTest<float, unsigned long> cuda_sm_csr_min_element_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSRMinElementTest<double, unsigned long> cuda_sm_csr_min_element_test_double_ulong(PreferredBackend::cuda);
#endif
