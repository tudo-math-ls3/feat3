// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
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
 * \tparam Mem_
 * description missing
 *
 * \tparam DT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRTest")
  {
  }

  virtual ~SparseMatrixCSRTest()
  {
  }

  virtual void run() const override
  {
    SparseMatrixCSR<Mem_, DT_, IT_> zero1;
    SparseMatrixCSR<Mem::Main, DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);
    zero2.convert(zero1);

    SparseMatrixCSR<Mem_, DT_, IT_> zero3(10, 11, 12);
    TEST_CHECK_EQUAL(zero3.used_elements(), 12);
    TEST_CHECK_EQUAL(zero3.rows(), 10);
    TEST_CHECK_EQUAL(zero3.columns(), 11);
    TEST_CHECK_EQUAL(zero3.size(), 110);


    SparseMatrixCSR<Mem_, DT_, IT_> empty1(11, 12, 111);
    SparseMatrixCSR<Mem::Main, DT_, IT_> empty2;
    SparseMatrixCSR<Mem_, DT_, IT_> empty3;
    empty2.convert(empty1);
    empty3.convert(empty2);
    TEST_CHECK_EQUAL(empty1.rows(), empty3.rows());
    TEST_CHECK_EQUAL(empty1.columns(), empty3.columns());
    TEST_CHECK_EQUAL(empty1.used_elements(), empty3.used_elements());

    SparseMatrixCSR<Mem::Main, DT_, IT_> empty4(empty2.layout());
    SparseMatrixCSR<Mem_, DT_, IT_> empty5(empty3.layout());
    empty4.convert(empty1);
    empty5.convert(empty4);
    TEST_CHECK_EQUAL(empty5.rows(), empty5.rows());
    TEST_CHECK_EQUAL(empty5.columns(), empty5.columns());
    TEST_CHECK_EQUAL(empty5.used_elements(), empty5.used_elements());
    empty5.convert(zero1);
    TEST_CHECK_EQUAL(empty5.rows(), 0);
    TEST_CHECK_EQUAL(empty5.columns(), 0);
    TEST_CHECK_EQUAL(empty5.used_elements(), 0);

    SparseMatrixCOO<Mem::Main, DT_, IT_> a00;
    SparseMatrixCSR<Mem_, DT_, IT_> b00(a00);
    TEST_CHECK_EQUAL(b00.size(), 0ul);
    TEST_CHECK_EQUAL(b00.rows(), 0ul);
    TEST_CHECK_EQUAL(b00.columns(), 0ul);
    TEST_CHECK_EQUAL(b00.used_elements(), 0ul);

    SparseMatrixFactory<DT_, IT_> a01(IT_(10), IT_(10));
    SparseMatrixCSR<Mem_, DT_, IT_> b01(a01.make_csr());
    TEST_CHECK_EQUAL(b01.size(), 100ul);
    TEST_CHECK_EQUAL(b01.rows(), 10ul);
    TEST_CHECK_EQUAL(b01.columns(), 10ul);
    TEST_CHECK_EQUAL(b01.used_elements(), 0ul);

    SparseMatrixFactory<DT_, IT_> a(IT_(10), IT_(10));
    a.add(IT_(1),IT_(2),DT_(7));
    a.add(IT_(5),IT_(5),DT_(2));
    a.add(IT_(5),IT_(7),DT_(3));
    a.add(IT_(5),IT_(2),DT_(4));
    SparseMatrixCSR<Mem_, DT_, IT_> b(a.make_csr());
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

    SparseMatrixCSR<Mem_, DT_, IT_> bl(b.layout());
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    bl = b.layout();
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    typename SparseLayout<Mem_, IT_, SparseLayoutId::lt_csr>::template MatrixType<DT_> x(b.layout());
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


    SparseMatrixCSR<Mem_, DT_, IT_> z;
    z.convert(b);
    TEST_CHECK_EQUAL(z.used_elements(), 4ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z(1, 2), DT_(7));
    TEST_CHECK_EQUAL(z(5, 5), DT_(2));

    SparseMatrixCSR<Mem_, DT_, IT_> c;
    c.clone(b);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_EQUAL((void*)c.col_ind(), (void*)b.col_ind());
    c = b.clone(CloneMode::Deep);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_NOT_EQUAL((void*)c.col_ind(), (void*)b.col_ind());

    DenseVector<Mem_, IT_, IT_> col_ind(c.used_elements(), c.col_ind());
    DenseVector<Mem_, DT_, IT_> val(c.used_elements(), c.val());
    DenseVector<Mem_, IT_, IT_> row_ptr(c.rows() + 1, c.row_ptr());
    SparseMatrixCSR<Mem_, DT_, IT_> d(c.rows(), c.columns(), col_ind, val, row_ptr);
    TEST_CHECK_EQUAL(d, c);

    SparseMatrixCSR<Mem::Main, DT_, IT_> e;
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
    MemoryPool<Mem_>::set_memory(clone1.val() + 1, DT_(132));
    TEST_CHECK_NOT_EQUAL(clone1, b);
    TEST_CHECK_NOT_EQUAL((void*)clone1.val(), (void*)b.val());
    TEST_CHECK_NOT_EQUAL((void*)clone1.row_ptr(), (void*)b.row_ptr());
    auto clone2 = clone1.clone(CloneMode::Layout);
    MemoryPool<Mem_>::set_memory(clone2.val(), DT_(4713), clone2.used_elements());
    TEST_CHECK_NOT_EQUAL(clone2(5, 5), clone1(5, 5));
    TEST_CHECK_NOT_EQUAL((void*)clone2.val(), (void*)clone1.val());
    TEST_CHECK_EQUAL((void*)clone2.row_ptr(), (void*)clone1.row_ptr());
    auto clone3 = clone1.clone(CloneMode::Weak);
    TEST_CHECK_EQUAL(clone3, clone1);
    MemoryPool<Mem_>::set_memory(clone3.val() + 1, DT_(133));
    TEST_CHECK_NOT_EQUAL(clone3, clone1);
    TEST_CHECK_NOT_EQUAL((void*)clone3.val(), (void*)clone1.val());
    TEST_CHECK_EQUAL((void*)clone3.row_ptr(), (void*)clone1.row_ptr());
    auto clone4 = clone1.clone(CloneMode::Shallow);
    TEST_CHECK_EQUAL(clone4, clone1);
    MemoryPool<Mem_>::set_memory(clone4.val() + 1, DT_(134));
    TEST_CHECK_EQUAL(clone4, clone1);
    TEST_CHECK_EQUAL((void*)clone4.val(), (void*)clone1.val());
    TEST_CHECK_EQUAL((void*)clone4.row_ptr(), (void*)clone1.row_ptr());

    SparseMatrixFactory<DT_, IT_> ffac(IT_(10), IT_(10));
    for (IT_ row(0) ; row < ffac.rows() ; ++row)
    {
      for (IT_ col(0) ; col < ffac.columns() ; ++col)
      {
        if(row == col)
          ffac.add(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          ffac.add(row, col, DT_(-1));
      }
    }
    SparseMatrixCSR<Mem_, DT_, IT_> f(ffac.make_csr());

    // shrink test
    SparseMatrixCSR<Mem_, DT_, IT_> l(f.clone());
    l.shrink(DT_(1.9));
    TEST_CHECK_EQUAL(l.used_elements(), 10ul);
  }

};

SparseMatrixCSRTest<Mem::Main, float, unsigned long> cpu_sparse_matrix_csr_test_float_ulong;
SparseMatrixCSRTest<Mem::Main, double, unsigned long> cpu_sparse_matrix_csr_test_double_ulong;
SparseMatrixCSRTest<Mem::Main, float, unsigned int> cpu_sparse_matrix_csr_test_float_uint;
SparseMatrixCSRTest<Mem::Main, double, unsigned int> cpu_sparse_matrix_csr_test_double_uint;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRTest<Mem::Main, __float128, unsigned long> cpu_sparse_matrix_csr_test_float128_ulong;
SparseMatrixCSRTest<Mem::Main, __float128, unsigned int> cpu_sparse_matrix_csr_test_float128_uint;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRTest<Mem::CUDA, float, unsigned long> cuda_sparse_matrix_csr_test_float_ulong;
SparseMatrixCSRTest<Mem::CUDA, double, unsigned long> cuda_sparse_matrix_csr_test_double_ulong;
SparseMatrixCSRTest<Mem::CUDA, float, unsigned int> cuda_sparse_matrix_csr_test_float_uint;
SparseMatrixCSRTest<Mem::CUDA, double, unsigned int> cuda_sparse_matrix_csr_test_double_uint;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRSerializeTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRSerializeTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRSerializeTest")
  {
  }

  virtual ~SparseMatrixCSRSerializeTest()
  {
  }

  virtual void run() const override
  {
    SparseMatrixFactory<DT_, IT_> ffac(IT_(10), IT_(10));
    for (IT_ row(0) ; row < ffac.rows() ; ++row)
    {
      for (IT_ col(0) ; col < ffac.columns() ; ++col)
      {
        if(row == col)
          ffac.add(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          ffac.add(row, col, DT_(-1));
      }
    }
    SparseMatrixCSR<Mem_, DT_, IT_> f(ffac.make_csr());

    BinaryStream bs;
    f.write_out(FileMode::fm_csr, bs);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));
    bs.seekg(0);
    SparseMatrixCSR<Mem_, DT_, IT_> g(FileMode::fm_csr, bs);
    TEST_CHECK_EQUAL(g, f);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));

    std::stringstream ts;
    f.write_out(FileMode::fm_mtx, ts);
    SparseMatrixCSR<Mem::Main, DT_, IT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_EQUAL(j, f);

    std::stringstream ts2;
    f.write_out(FileMode::fm_mtx, ts2, true);
    SparseMatrixCSR<Mem::Main, DT_, IT_> j2(FileMode::fm_mtx, ts2);
    TEST_CHECK_EQUAL(j2, f);

    auto kp = f.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixCSR<Mem_, DT_, IT_> k(kp);
    TEST_CHECK_EQUAL(k, f);
#ifdef FEAT_HAVE_ZLIB
    auto zl = f.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixCSR<Mem_, DT_, IT_> zlib(zl);
    TEST_CHECK_EQUAL(zlib, f);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = f.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixCSR<Mem_, DT_, IT_> zfp(zf);
    for(Index row(0) ; row < f.rows() ; ++row)
    {
      for(Index col(0) ; col < f.columns() ; ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(zfp(row, col), f(row, col), DT_(1e-4));
      }
    }
#endif
  }
};
SparseMatrixCSRSerializeTest<Mem::Main, float, unsigned long> cpu_sparse_matrix_csr_serialize_test_float_ulong;
SparseMatrixCSRSerializeTest<Mem::Main, double, unsigned long> cpu_sparse_matrix_csr_serialize_test_double_ulong;
SparseMatrixCSRSerializeTest<Mem::Main, float, unsigned int> cpu_sparse_matrix_csr_serialize_test_float_uint;
SparseMatrixCSRSerializeTest<Mem::Main, double, unsigned int> cpu_sparse_matrix_csr_serialize_test_double_uint;
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRSerializeTest<Mem::CUDA, float, unsigned long> cuda_sparse_matrix_csr_serialize_test_float_ulong;
SparseMatrixCSRSerializeTest<Mem::CUDA, double, unsigned long> cuda_sparse_matrix_csr_serialize_test_double_ulong;
SparseMatrixCSRSerializeTest<Mem::CUDA, float, unsigned int> cuda_sparse_matrix_csr_serialize_test_float_uint;
SparseMatrixCSRSerializeTest<Mem::CUDA, double, unsigned int> cuda_sparse_matrix_csr_serialize_test_double_uint;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRApplyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRApplyTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRApplyTest")
  {
  }

  virtual ~SparseMatrixCSRApplyTest()
  {
  }

  virtual void run() const override
  {
    DT_ s(DT_(4711.1));
    for (IT_ size(1) ; size < IT_(1e3) ; size*=IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_local(size, size);
      DenseVector<Mem::Main, DT_, IT_> x_local(size);
      DenseVector<Mem::Main, DT_, IT_> y_local(size);
      DenseVector<Mem::Main, DT_, IT_> ref_local(size);
      DenseVector<Mem_, DT_, IT_> ref(size);
      DenseVector<Mem::Main, DT_, IT_> result_local(size);
      DenseVector<Mem::Main, DT_, IT_> ax_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        x_local(i, DT_(i % 100) * DT_(1.234));
        y_local(i, DT_(2 - DT_(i % 42)));
        if (i == 0 && size > 1)
          ax_local(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i + 1) % 100)));
        else if (i == size - 1 && size > 1)
          ax_local(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i - 1) % 100)));
        else if (size == 1 && i == 0)
          ax_local(i, DT_(0));
        else
        ax_local(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i - 1) % 100) - DT_((i + 1) % 100)));
      }
      DenseVector<Mem_, DT_, IT_> x(size);
      x.copy(x_local);
      DenseVector<Mem_, DT_, IT_> y(size);
      y.copy(y_local);

      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
          {
            a_local.add(row, col, DT_(2));
          }
          else if((row == col+1) || (row+1 == col))
          {
            a_local.add(row, col, DT_(-1));
          }
        }
      }
      SparseMatrixCSR<Mem_,DT_, IT_> a(a_local.make_csr());

      DenseVector<Mem_, DT_, IT_> r(size);

      // apply-test for alpha = 0.0
      a.apply(r, x, y, DT_(0.0));
      result_local.copy(r);
      ref_local.copy(y);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(1e-2));

      // apply-test for alpha = -1.0
      a.apply(r, x, y, DT_(-1.0));
      result_local.copy(r);
      a.apply(ref, x);
      ref.scale(ref, DT_(-1.0));
      ref.axpy(ref, y);
      ref_local.copy(ref);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(1e-2));

      // apply-test for alpha = -1.0 and &r==&y
      r.copy(y);
      a.apply(r, x, r, DT_(-1.0));
      result_local.copy(r);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(1e-2));

      // apply-test for alpha = 4711.1
      //r.axpy(s, a, x, y);
      a.apply(r, x, y, s);
      result_local.copy(r);
      //ref.product_matvec(a, x);
      a.apply(ref, x);
      ref.scale(ref, s);
      ref.axpy(ref, y);
      ref_local.copy(ref);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(5e-2));

      // apply-test for alpha = 4711.1 and &r==&y
      r.copy(y);
      a.apply(r, x, r, s);
      result_local.copy(r);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(5e-2));

     a.apply(r, x);
      result_local.copy(r);
     // a_local.apply(ref_local, x_local);
      ref_local.copy(ax_local);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(1e-2));

      // transposed apply-test for alpha = 4711.1
      if (!(typeid(Mem::CUDA) == typeid(Mem_) && typeid(IT_) == typeid(unsigned long)))
      {
        a.apply(r, x, y, s, true);
        result_local.copy(r);

        SparseMatrixCSR<Mem_, DT_, IT_> at = a.transpose();
        at.apply(ref, x);
        ref.scale(ref, s);
        ref.axpy(ref, y);
        ref_local.copy(ref);
        for (Index i(0) ; i < size ; ++i)
          TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(7e-2));
      }
    }
  }
};

SparseMatrixCSRApplyTest<Mem::Main, float, unsigned long> sm_csr_apply_test_float_ulong;
SparseMatrixCSRApplyTest<Mem::Main, double, unsigned long> sm_csr_apply_test_double_ulong;
SparseMatrixCSRApplyTest<Mem::Main, float, unsigned int> sm_csr_apply_test_float_uint;
SparseMatrixCSRApplyTest<Mem::Main, double, unsigned int> sm_csr_apply_test_double_uint;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRApplyTest<Mem::Main, __float128, unsigned long> sm_csr_apply_test_float128_ulong;
SparseMatrixCSRApplyTest<Mem::Main, __float128, unsigned int> sm_csr_apply_test_float128_uint;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRApplyTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_apply_test_float_ulong;
SparseMatrixCSRApplyTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_apply_test_double_ulong;
SparseMatrixCSRApplyTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_apply_test_float_uint;
SparseMatrixCSRApplyTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_apply_test_double_uint;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRBApplyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRBApplyTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRBApplyTest")
  {
  }

  virtual ~SparseMatrixCSRBApplyTest()
  {
  }

  virtual void run() const override
  {
    for (IT_ size(1) ; size < IT_(1e3) ; size*=IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_local(size, size);
      DenseVector<Mem::Main, DT_, IT_> ref_x_local(size);
     // DenseVector<Mem::Main, DT_, IT_> ref_local(size);
      DenseVector<Mem::Main, DT_, IT_> aref_x_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        ref_x_local(i, DT_(i % 100) * DT_(1.234));
        if (i == 0 && size > 1)
          aref_x_local(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i + 1) % 100)));
        else if (i == size - 1 && size > 1)
          aref_x_local(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i - 1) % 100)));
        else if (i == 0 && size == 1)
          aref_x_local(i, DT_(1.234) * DT_(2) * DT_(i % 100));
        else
          aref_x_local(i, DT_(1.234) * (DT_(2) * DT_(i % 100) - DT_((i - 1) % 100) - DT_((i + 1) % 100)));
      }

      //a_local.apply(ref_local, ref_x_local);
      //ref_local = aref_x_local;

      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
          {
            a_local.add(row, col, DT_(2));
          }
          else if((row == col+1) || (row+1 == col))
          {
            a_local.add(row, col, DT_(-1));
          }
        }
      }
      SparseMatrixCSR<Mem_,DT_, IT_> a(a_local.make_csr());


      DenseVectorBlocked<Mem::Main, DT_, IT_, 3> x_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        auto temp = x_local(i);
        temp[0] = ref_x_local(i);
        temp[1] = DT_(0.5) * ref_x_local(i);
        temp[2] = DT_(2.0) * ref_x_local(i);
        x_local(i, temp);
      }
      DenseVectorBlocked<Mem_, DT_, IT_, 3> x;
      x.convert(x_local);

      DenseVectorBlocked<Mem::Main, DT_, IT_, 3> y_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        auto temp = y_local(i);
        temp[0] = DT_(i);
        temp[1] = DT_(i*2);
        temp[2] = DT_(i*3);
        y_local(i, temp);
      }
      DenseVectorBlocked<Mem_, DT_, IT_, 3> y;
      y.convert(y_local);

      DenseVectorBlocked<Mem_, DT_, IT_, 3> r(size);
      DenseVectorBlocked<Mem_, DT_, IT_, 3> r_local;

      a.apply(r, x);
      r_local.convert(r);
      for (Index i(0) ; i < size ; ++i)
      {
        //std::cout << i << std::endl;
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[0], aref_x_local(i), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[1], aref_x_local(i) * DT_(0.5), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[2], aref_x_local(i) * DT_(2.0), DT_(1e-4));
      }

      a.apply(r, x, y, DT_(-1));
      r_local.convert(r);
      for (Index i(0) ; i < size ; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[0], y_local(i)[0] - aref_x_local(i), DT_(1e-4));
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[1], y_local(i)[1] - aref_x_local(i) * DT_(0.5), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[2], y_local(i)[2] - aref_x_local(i) * DT_(2.0), DT_(1e-4));
      }

      DT_ alpha(0.75);
      a.apply(r, x, y, alpha);
      r_local.convert(r);
      for (Index i(0) ; i < size ; ++i)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[0], y_local(i)[0] + alpha * aref_x_local(i), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[1], y_local(i)[1] + alpha * aref_x_local(i) * DT_(0.5), DT_(1e-5));
        TEST_CHECK_EQUAL_WITHIN_EPS(r_local(i)[2], y_local(i)[2] + alpha * aref_x_local(i) * DT_(2.0), DT_(1e-4));
      }
    }
  }
};

SparseMatrixCSRBApplyTest<Mem::Main, float, unsigned long> sm_csrib_apply_test_float_ulong;
SparseMatrixCSRBApplyTest<Mem::Main, double, unsigned long> sm_csrib_apply_test_double_ulong;
SparseMatrixCSRBApplyTest<Mem::Main, float, unsigned int> sm_csrib_apply_test_float_uint;
SparseMatrixCSRBApplyTest<Mem::Main, double, unsigned int> sm_csrib_apply_test_double_uint;
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRBApplyTest<Mem::CUDA, float, unsigned long> cuda_sm_csrib_apply_test_float_ulong;
SparseMatrixCSRBApplyTest<Mem::CUDA, double, unsigned long> cuda_sm_csrib_apply_test_double_ulong;
SparseMatrixCSRBApplyTest<Mem::CUDA, float, unsigned int> cuda_sm_csrib_apply_test_float_uint;
SparseMatrixCSRBApplyTest<Mem::CUDA, double, unsigned int> cuda_sm_csrib_apply_test_double_uint;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRScaleTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRScaleTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRScaleTest")
  {
  }

  virtual ~SparseMatrixCSRScaleTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=2)
    {
      DT_ s(DT_(4.321));
      SparseMatrixFactory<DT_, IT_> a_local(IT_(size), IT_(size+2));
      SparseMatrixFactory<DT_, IT_> ref_local(IT_(size), IT_(size+2));

      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local.add(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local.add(row, col, DT_(-1));

          if(row == col)
            ref_local.add(row, col, DT_(2) * s);
          else if((row == col+1) || (row+1 == col))
            ref_local.add(row, col, DT_(-1) * s);
        }
      }

      SparseMatrixCSR<Mem_, DT_, IT_> ref(ref_local.make_csr());
      SparseMatrixCSR<Mem_, DT_, IT_> a(a_local.make_csr());
      SparseMatrixCSR<Mem_, DT_, IT_> b;
      b.clone(a);

      b.scale(a, s);
      TEST_CHECK_EQUAL(b, ref);

      a.scale(a, s);
      TEST_CHECK_EQUAL(a, ref);
    }
  }
};

SparseMatrixCSRScaleTest<Mem::Main, float, unsigned int> sm_csr_scale_test_float_uint;
SparseMatrixCSRScaleTest<Mem::Main, double, unsigned int> sm_csr_scale_test_double_uint;
SparseMatrixCSRScaleTest<Mem::Main, float, unsigned long> sm_csr_scale_test_float_ulong;
SparseMatrixCSRScaleTest<Mem::Main, double, unsigned long> sm_csr_scale_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRScaleTest<Mem::Main, __float128, unsigned int> sm_csr_scale_test_float128_uint;
SparseMatrixCSRScaleTest<Mem::Main, __float128, unsigned long> sm_csr_scale_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRScaleTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_scale_test_float_uint;
SparseMatrixCSRScaleTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_scale_test_double_uint;
SparseMatrixCSRScaleTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_scale_test_float_ulong;
SparseMatrixCSRScaleTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_scale_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRScaleRowColTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRScaleRowColTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRScaleRowColTest")
  {
  }

  virtual ~SparseMatrixCSRScaleRowColTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=3)
    {
      const DT_ pi(Math::pi<DT_>());
      const DT_ eps(Math::pow(Math::eps<DT_>(), DT_(0.8)));

      SparseMatrixFactory<DT_, IT_> a_local(IT_(size), IT_(size + 2));
      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local.add(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local.add(row, col, DT_(-1));
        }
      }

      SparseMatrixCSR<Mem_, DT_, IT_> a(a_local.make_csr());
      SparseMatrixCSR<Mem_, DT_, IT_> b(a.clone());

      // Scale rows
      DenseVector<Mem_, DT_, IT_> s1(a.rows());
      for (Index i(0); i < s1.size(); ++i)
      {
        s1(i, pi * (i % 3 + 1) - DT_(5.21) + DT_(i));
      }
      b.scale_rows(b, s1);
      for (Index row(0) ; row < a.rows() ; ++row)
      {
        for (Index col(0) ; col < a.columns() ; ++col)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(b(row, col), a(row, col) * s1(row), eps);
        }
      }

      // Scale rows
      DenseVector<Mem_, DT_, IT_> s2(a.columns());
      for (Index i(0); i < s2.size(); ++i)
      {
        s2(i, pi * (i % 3 + 1) - DT_(5.21) + DT_(i));
      }
      b.scale_cols(a, s2);
      for (Index row(0) ; row < a.rows() ; ++row)
      {
        for (Index col(0) ; col < a.columns() ; ++col)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(b(row, col), a(row, col) * s2(col), eps);
        }
      }
    }
  }
};

SparseMatrixCSRScaleRowColTest<Mem::Main, float, unsigned int> sm_csr_scale_row_col_test_float_uint;
SparseMatrixCSRScaleRowColTest<Mem::Main, double, unsigned int> sm_csr_scale_row_col_test_double_uint;
SparseMatrixCSRScaleRowColTest<Mem::Main, float, unsigned long> sm_csr_scale_row_col_test_float_ulong;
SparseMatrixCSRScaleRowColTest<Mem::Main, double, unsigned long> sm_csr_scale_row_col_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRScaleRowColTest<Mem::Main, __float128, unsigned int> sm_csr_scale_row_col_test_float128_uint;
SparseMatrixCSRScaleRowColTest<Mem::Main, __float128, unsigned long> sm_csr_scale_row_col_test_float128_ulong;
#endif
// #ifdef FEAT_HAVE_CUDA
// SparseMatrixCSRScaleRowColTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_scale_row_col_test_float_uint;
// SparseMatrixCSRScaleRowColTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_scale_row_col_test_double_uint;
// SparseMatrixCSRScaleRowColTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_scale_row_col_test_float_ulong;
// SparseMatrixCSRScaleRowColTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_scale_row_col_test_double_ulong;
// #endif


/**
 * \brief Test class for the transposition of a SparseMatrixCSR
 *
 * \test test description missing
 *
 * \tparam Mem_
 * description missing
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \author Christoph Lohmann
 */
template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRTranspositionTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{

public:
  typedef SparseMatrixCSR<Mem_, DT_, IT_> MatrixType;

   SparseMatrixCSRTranspositionTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRTranspositionTest")
  {
  }

  virtual ~SparseMatrixCSRTranspositionTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=4)
    {
      SparseMatrixFactory<DT_, IT_> a_local(IT_(size), IT_(size + 2));

      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local.add(row, col, DT_(2));
          else if(row == col+1)
            a_local.add(row, col, DT_(-1));
          else if(row+1 == col)
            a_local.add(row, col, DT_(-3));
        }
      }
      MatrixType a(a_local.make_csr());


      MatrixType b;
      b.transpose(a);

      for (Index i(0) ; i < a.rows() ; ++i)
      {
        for (Index j(0) ; j < a.columns() ; ++j)
        {
          TEST_CHECK_EQUAL(b(j, i), a(i, j));
        }
      }

      b = b.transpose();

      TEST_CHECK_EQUAL(a, b);
    }
  }
};

SparseMatrixCSRTranspositionTest<Mem::Main, float, unsigned int> sm_csr_transposition_test_float_uint;
SparseMatrixCSRTranspositionTest<Mem::Main, double, unsigned int> sm_csr_transposition_test_double_uint;
SparseMatrixCSRTranspositionTest<Mem::Main, float, unsigned long> sm_csr_transposition_test_float_ulong;
SparseMatrixCSRTranspositionTest<Mem::Main, double, unsigned long> sm_csr_transposition_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRTranspositionTest<Mem::Main, __float128, unsigned int> sm_csr_transposition_test_float128_uint;
SparseMatrixCSRTranspositionTest<Mem::Main, __float128, unsigned long> sm_csr_transposition_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRTranspositionTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_transposition_test_float_uint;
SparseMatrixCSRTranspositionTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_transposition_test_double_uint;
SparseMatrixCSRTranspositionTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_transposition_test_float_ulong;
SparseMatrixCSRTranspositionTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_transposition_test_double_ulong;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRPermuteTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRPermuteTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRPermuteTest")
  {
  }

  virtual ~SparseMatrixCSRPermuteTest()
  {
  }

  virtual void run() const override
  {
    for (IT_ size(25) ; size < IT_(1e3) ; size*=IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_local(size, size);
      DenseVector<Mem::Main, DT_, IT_> x_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        x_local(i, DT_(i % 100) * DT_(1.234));
      }
      DenseVector<Mem_, DT_, IT_> x(size);
      x.copy(x_local);

      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
          {
            a_local.add(row, col, DT_(2));
          }
          else if((row == col+7) || (row+7 == col))
          {
            a_local.add(row, col, DT_(-1));
          }
          else if((row == col+15) || (row+15 == col))
          {
            a_local.add(row, col, DT_(-2));
          }
          else if((row == col+a_local.columns()/2) || (row+a_local.rows()/2 == col))
          {
            a_local.add(row, col, DT_(1));
          }
        }
      }
      SparseMatrixCSR<Mem_, DT_, IT_> a(a_local.make_csr());

      DenseVector<Mem_, DT_, IT_> r(size);
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

SparseMatrixCSRPermuteTest<Mem::Main, float, unsigned long> sm_csr_permute_test_float_ulong;
SparseMatrixCSRPermuteTest<Mem::Main, double, unsigned long> sm_csr_permute_test_double_ulong;
SparseMatrixCSRPermuteTest<Mem::Main, float, unsigned int> sm_csr_permute_test_float_uint;
SparseMatrixCSRPermuteTest<Mem::Main, double, unsigned int> sm_csr_permute_test_double_uint;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRPermuteTest<Mem::Main, __float128, unsigned long> sm_csr_permute_test_float128_ulong;
SparseMatrixCSRPermuteTest<Mem::Main, __float128, unsigned int> sm_csr_permute_test_float128_uint;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRPermuteTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_permute_test_float_ulong;
SparseMatrixCSRPermuteTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_permute_test_double_ulong;
SparseMatrixCSRPermuteTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_permute_test_float_uint;
SparseMatrixCSRPermuteTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_permute_test_double_uint;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRDiagTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRDiagTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRDiagTest")
  {
  }

  virtual ~SparseMatrixCSRDiagTest()
  {
  }

  virtual void run() const override
  {
    for (IT_ size(2) ; size < IT_(3e2) ; size*=IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_local(size, size);
      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local.add(row, col, DT_(DT_(col % 100) / DT_(2)));
          else if((row == col+1) || (row+1 == col))
            a_local.add(row, col, DT_(-1));
        }
      }

      SparseMatrixCSR<Mem_, DT_, IT_> a(a_local.make_csr());

      auto ref = a.create_vector_l();
      //auto ref_local = a_local.create_vector_l();
      DenseVector<Mem::Main, DT_, IT_> ref_local(a_local.rows());
      for (Index i(0) ; i < a_local.rows() ; ++i)
      {
        ref_local(i, DT_(DT_(i % 100) / DT_(2)));
      }
      ref.convert(ref_local);

      auto diag = a.extract_diag();
      TEST_CHECK_EQUAL(diag, ref);
    }
  }
};

SparseMatrixCSRDiagTest<Mem::Main, float, unsigned int> sm_csr_diag_test_float_uint;
SparseMatrixCSRDiagTest<Mem::Main, double, unsigned int> sm_csr_diag_test_double_uint;
SparseMatrixCSRDiagTest<Mem::Main, float, unsigned long> sm_csr_diag_test_float_ulong;
SparseMatrixCSRDiagTest<Mem::Main, double, unsigned long> sm_csr_diag_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRDiagTest<Mem::Main, __float128, unsigned int> sm_csr_diag_test_float128_uint;
SparseMatrixCSRDiagTest<Mem::Main, __float128, unsigned long> sm_csr_diag_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRDiagTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_diag_test_float_uint;
SparseMatrixCSRDiagTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_diag_test_double_uint;
SparseMatrixCSRDiagTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_diag_test_float_ulong;
SparseMatrixCSRDiagTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_diag_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRAxpyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRAxpyTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRAxpyTest")
  {
  }

  virtual ~SparseMatrixCSRAxpyTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=2)
    {
      DT_ s(DT_(4.321));

      SparseMatrixFactory<DT_, IT_> a_local(IT_(size), IT_(size+2));
      SparseMatrixFactory<DT_, IT_> b_local(IT_(size), IT_(size+2));
      SparseMatrixFactory<DT_, IT_> ref_local(IT_(size), IT_(size+2));
      SparseMatrixFactory<DT_, IT_> ref_local_plus(IT_(size), IT_(size + 2));
      SparseMatrixFactory<DT_, IT_> ref_local_minus(IT_(size), IT_(size + 2));
      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if (row == col)
          {
            a_local.add(row, col, DT_(2));
            b_local.add(row, col, DT_((row + col) % 15));
            ref_local.add(row, col, DT_(2) * s + DT_((row + col) % 15));
            ref_local_plus.add(row, col, DT_(2)  + DT_((row + col) % 15));
            ref_local_minus.add(row, col,   DT_((row + col) % 15)- DT_(2));
          }

          else if ((row == col + 1) || (row + 1 == col))
          {
            a_local.add(row, col, DT_(-1));
            b_local.add(row, col, DT_((row + col) % 15));
            ref_local.add(row, col, DT_(-1) * s + DT_((row + col) % 15));
            ref_local_plus.add(row, col, DT_(-1) + DT_((row + col) % 15));
            ref_local_minus.add(row, col, DT_((row + col) % 15) - DT_(-1));
          }
        }
      }

      SparseMatrixCSR<Mem_, DT_, IT_> ref(ref_local.make_csr());
      SparseMatrixCSR<Mem_, DT_, IT_> ref_plus(ref_local_plus.make_csr());
      SparseMatrixCSR<Mem_, DT_, IT_> ref_minus(ref_local_minus.make_csr());
      SparseMatrixCSR<Mem_, DT_, IT_> a(a_local.make_csr());
      SparseMatrixCSR<Mem_, DT_, IT_> b(b_local.make_csr());
      SparseMatrixCSR<Mem_, DT_, IT_> c;

      c.clone(a);
      c.axpy(c, b, s);
      TEST_CHECK_EQUAL(c, ref);

      c.clone(b);
      c.axpy(a, c, s);
      TEST_CHECK_EQUAL(c, ref);

      c.axpy(a, b, s);
      TEST_CHECK_EQUAL(c, ref);

      s = DT_(0);
      ref.clone(b);
      c.axpy(a, b, s);
      TEST_CHECK_EQUAL(c, ref);

      s = DT_(1);
      c.axpy(a, b, s);
      TEST_CHECK_EQUAL(c, ref_plus);

      s = DT_(-1);
      c.axpy(a, b, s);
      TEST_CHECK_EQUAL(c, ref_minus);
    }
  }
};

SparseMatrixCSRAxpyTest<Mem::Main, float, unsigned int> sm_csr_axpy_test_float_uint;
SparseMatrixCSRAxpyTest<Mem::Main, double, unsigned int> sm_csr_axpy_test_double_uint;
SparseMatrixCSRAxpyTest<Mem::Main, float, unsigned long> sm_csr_axpy_test_float_ulong;
SparseMatrixCSRAxpyTest<Mem::Main, double, unsigned long> sm_csr_axpy_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRAxpyTest<Mem::Main, __float128, unsigned int> sm_csr_axpy_test_float128_uint;
SparseMatrixCSRAxpyTest<Mem::Main, __float128, unsigned long> sm_csr_axpy_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRAxpyTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_axpy_test_float_uint;
SparseMatrixCSRAxpyTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_axpy_test_double_uint;
SparseMatrixCSRAxpyTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_axpy_test_float_ulong;
SparseMatrixCSRAxpyTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_axpy_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRFrobeniusTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRFrobeniusTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRFrobeniusTest")
  {
  }

  virtual ~SparseMatrixCSRFrobeniusTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=2)
    {
      SparseMatrixFactory<DT_, IT_> a_local(IT_(size), IT_(size + 2));
      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local.add(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local.add(row, col, DT_(-1));
        }
      }

      SparseMatrixCSR<Mem_, DT_, IT_> a(a_local.make_csr());

      DenseVector<Mem_, DT_, IT_> refv(a.used_elements(), a.val());
      DT_ ref = refv.norm2();
      DT_ c = a.norm_frobenius();
      TEST_CHECK_EQUAL(c, ref);
    }
  }
};

SparseMatrixCSRFrobeniusTest<Mem::Main, float, unsigned int> sm_csr_frobenius_test_float_uint;
SparseMatrixCSRFrobeniusTest<Mem::Main, double, unsigned int> sm_csr_frobenius_test_double_uint;
SparseMatrixCSRFrobeniusTest<Mem::Main, float, unsigned long> sm_csr_frobenius_test_float_ulong;
SparseMatrixCSRFrobeniusTest<Mem::Main, double, unsigned long> sm_csr_frobenius_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRFrobeniusTest<Mem::Main, __float128, unsigned int> sm_csr_frobenius_test_float128_uint;
SparseMatrixCSRFrobeniusTest<Mem::Main, __float128, unsigned long> sm_csr_frobenius_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRFrobeniusTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_frobenius_test_float_uint;
SparseMatrixCSRFrobeniusTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_frobenius_test_double_uint;
SparseMatrixCSRFrobeniusTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_frobenius_test_float_ulong;
SparseMatrixCSRFrobeniusTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_frobenius_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRLumpTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRLumpTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRLumpTest")
  {
  }

  virtual ~SparseMatrixCSRLumpTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (IT_ size(2) ; size < IT_(3e2) ; size*=IT_(2))
    {
      SparseMatrixFactory<DT_, IT_> a_local(size, size);
      for (IT_ row(0) ; row < a_local.rows() ; ++row)
      {
        for (IT_ col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local.add(row, col, DT_(DT_(col % 100) / DT_(2)));
          else if((row == col+1) || (row+1 == col))
            a_local.add(row, col, DT_(-1));
        }
      }

      SparseMatrixCSR<Mem_, DT_, IT_> a(a_local.make_csr());
      auto lump = a.lump_rows();
      auto one = a.create_vector_r();
      auto res = a.create_vector_r();
      one.format(DT_(-1));
      a.apply(res, one, lump);

      TEST_CHECK_EQUAL_WITHIN_EPS(res.norm2(), DT_(0), tol);
    }
  }
};

SparseMatrixCSRLumpTest<Mem::Main, float, unsigned int> sm_csr_lump_test_float_uint;
SparseMatrixCSRLumpTest<Mem::Main, double, unsigned int> sm_csr_lump_test_double_uint;
SparseMatrixCSRLumpTest<Mem::Main, float, unsigned long> sm_csr_lump_test_float_ulong;
SparseMatrixCSRLumpTest<Mem::Main, double, unsigned long> sm_csr_lump_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSRLumpTest<Mem::Main, __float128, unsigned int> sm_csr_lump_test_float128_uint;
SparseMatrixCSRLumpTest<Mem::Main, __float128, unsigned long> sm_csr_lump_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSRLumpTest<Mem::CUDA, float, unsigned int> cuda_sm_csr_lump_test_float_uint;
SparseMatrixCSRLumpTest<Mem::CUDA, double, unsigned int> cuda_sm_csr_lump_test_double_uint;
SparseMatrixCSRLumpTest<Mem::CUDA, float, unsigned long> cuda_sm_csr_lump_test_float_ulong;
SparseMatrixCSRLumpTest<Mem::CUDA, double, unsigned long> cuda_sm_csr_lump_test_double_ulong;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRCompressionTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSRCompressionTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSRCompressionTest")
  {
  }

  virtual ~SparseMatrixCSRCompressionTest()
  {
  }

  virtual void run() const override
  {
    Index mat_rows = 56;
    Index mat_cols = 56;
    DenseVector<Mem_, IT_, IT_> col_ind(400);
    DenseVector<Mem_, DT_, IT_> val(400);
    DenseVector<Mem_, IT_, IT_> row_ptr(mat_rows+1);
    for(Index i(0) ; i < 400 ; ++i)
    {
      col_ind(i, IT_(i/40 + (3*i)%16));
      val(i, DT_(7)/DT_(3*(i+1)) - DT_(13)/DT_((i+1)*(i+1)));
    }
    row_ptr(0, 0);
    for(Index i(1) ; i < mat_rows ; ++i)
    {
      row_ptr(i, IT_(i * Index(400/mat_rows -1)));
    }
    row_ptr(mat_rows, 400);
    SparseMatrixCSR<Mem_, DT_, IT_> a(mat_rows, mat_cols, col_ind, val, row_ptr);
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
SparseMatrixCSRCompressionTest<Mem::Main, float, unsigned int> sm_csr_comp_test_float_uint;
SparseMatrixCSRCompressionTest<Mem::Main, double, unsigned int> sm_csr_comp_test_double_uint;
SparseMatrixCSRCompressionTest<Mem::Main, float, unsigned long> sm_csr_comp_test_float_ulong;
SparseMatrixCSRCompressionTest<Mem::Main, double, unsigned long> sm_csr_comp_test_double_ulong;
