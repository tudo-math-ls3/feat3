// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_cscr.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/random.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>
//#include <kernel/adjacency/cuthill_mckee.hpp>

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
class SparseMatrixCSCRTest
  : public UnitTest
{
public:
   SparseMatrixCSCRTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSCRTest()
  {
  }

  virtual void run() const override
  {
    SparseMatrixCSCR<DT_, IT_> zero1;
    SparseMatrixCSCR<DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);
    zero2.convert(zero1);

    SparseMatrixCSCR<DT_, IT_> zero3(10, 11, 12, 3);
    TEST_CHECK_EQUAL(zero3.used_elements(), 12);
    TEST_CHECK_EQUAL(zero3.rows(), 10);
    TEST_CHECK_EQUAL(zero3.columns(), 11);
    TEST_CHECK_EQUAL(zero3.size(), 110);
    TEST_CHECK_EQUAL(zero3.used_rows(), 3);


    SparseMatrixCSCR<DT_, IT_> empty1(11, 12, 111, 3);
    SparseMatrixCSCR<DT_, IT_> empty2;
    SparseMatrixCSCR<DT_, IT_> empty3;
    empty2.convert(empty1);
    empty3.convert(empty2);
    TEST_CHECK_EQUAL(empty1.rows(), empty3.rows());
    TEST_CHECK_EQUAL(empty1.columns(), empty3.columns());
    TEST_CHECK_EQUAL(empty1.used_elements(), empty3.used_elements());
    TEST_CHECK_EQUAL(empty1.used_rows(), empty3.used_rows());

    SparseMatrixCSCR<DT_, IT_> empty4(empty2.layout());
    SparseMatrixCSCR<DT_, IT_> empty5(empty3.layout());
    empty4.convert(empty1);
    empty5.convert(empty4);
    TEST_CHECK_EQUAL(empty5.rows(), empty5.rows());
    TEST_CHECK_EQUAL(empty5.columns(), empty5.columns());
    TEST_CHECK_EQUAL(empty5.used_elements(), empty5.used_elements());
    TEST_CHECK_EQUAL(empty5.used_rows(), empty5.used_rows());
    empty5.convert(zero1);
    TEST_CHECK_EQUAL(empty5.rows(), 0);
    TEST_CHECK_EQUAL(empty5.columns(), 0);
    TEST_CHECK_EQUAL(empty5.used_elements(), 0);
    TEST_CHECK_EQUAL(empty5.used_rows(), 0);

    DenseVector<DT_, IT_> val(7);
    DenseVector<IT_, IT_> col_ind(7);
    for (Index i(0) ; i < val.size() ; ++i)
    {
      val(i, DT_(i+1));
      col_ind(i, IT_(i + 3));
    }
    DenseVector<IT_, IT_> row_ptr(4);
    row_ptr(0, IT_(0));
    row_ptr(1, IT_(2));
    row_ptr(2, IT_(5));
    row_ptr(3, IT_(6));
    DenseVector<IT_, IT_> row_numbers(3);
    row_numbers(0, IT_(0));
    row_numbers(1, IT_(4));
    row_numbers(2, IT_(8));

    SparseMatrixCSCR<DT_, IT_> a(10, 10, col_ind, val, row_ptr, row_numbers);
    TEST_CHECK_EQUAL(a.used_elements(), 7);
    TEST_CHECK_EQUAL(a.size(), 100);
    TEST_CHECK_EQUAL(a.rows(), 10);
    TEST_CHECK_EQUAL(a.columns(), 10);
    TEST_CHECK_EQUAL(a.used_rows(), 3);
    TEST_CHECK_EQUAL(a(0,0), DT_(0));
    TEST_CHECK_EQUAL(a(0,3), DT_(1));
    TEST_CHECK_EQUAL(a(4,4), DT_(0));
    TEST_CHECK_EQUAL(a(4,6), DT_(4));
    TEST_CHECK_EQUAL(a(9,9), DT_(0));

    SparseMatrixCSCR<DT_, IT_> b;
    b.convert(a);
    TEST_CHECK_EQUAL(b, a);
    SparseMatrixCSCR<DT_, IT_> c;
    c.convert(b);
    TEST_CHECK_EQUAL(c, b);
    SparseMatrixCSCR<DT_, IT_> d(a.clone());
    TEST_CHECK_EQUAL(d, c);


    SparseMatrixFactory<DT_, IT_> fac(IT_(10), IT_(10));
    fac.add(0, 3, DT_(1));
    fac.add(0, 4, DT_(2));
    fac.add(1, 5, DT_(7));
    fac.add(3, 3, DT_(8));
    fac.add(4, 5, DT_(3));
    fac.add(4, 6, DT_(4));
    fac.add(4, 7, DT_(5));
    fac.add(6, 6, DT_(9));
    fac.add(8, 8, DT_(6));
    SparseMatrixCSR<DT_, IT_> csr(10, 10);
    csr.convert(fac.make_csr());
    VectorMirror<DT_, IT_> mirror_main(10, 3);
    mirror_main.indices()[0] = 0;
    mirror_main.indices()[1] = 4;
    mirror_main.indices()[2] = 8;
    VectorMirror<DT_, IT_> mirror;
    mirror.convert(mirror_main);
    SparseMatrixCSCR<DT_, IT_> f(csr, mirror);
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(f(row, col), a(row, col), Math::eps<DT_>());
      }
    }

  }
};

SparseMatrixCSCRTest<float, unsigned long> cpu_sparse_matrix_cscr_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSCRTest<double, unsigned long> cpu_sparse_matrix_cscr_test_double_ulong(PreferredBackend::generic);
SparseMatrixCSCRTest<float, unsigned int> cpu_sparse_matrix_cscr_test_float_uint(PreferredBackend::generic);
SparseMatrixCSCRTest<double, unsigned int> cpu_sparse_matrix_cscr_test_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSCRTest<float, unsigned long> mkl_sparse_matrix_cscr_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSCRTest<double, unsigned long> mkl_sparse_matrix_cscr_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRTest<__float128, unsigned long> sparse_matrix_cscr_test_float128_ulong(PreferredBackend::generic);
SparseMatrixCSCRTest<__float128, unsigned int> sparse_matrix_cscr_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSCRTest<Half, unsigned int> sparse_matrix_cscr_test_half_uint(PreferredBackend::generic);
SparseMatrixCSCRTest<Half, unsigned long> sparse_matrix_cscr_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRTest<float, unsigned long> cuda_sparse_matrix_cscr_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSCRTest<double, unsigned long> cuda_sparse_matrix_cscr_test_double_ulong(PreferredBackend::cuda);
SparseMatrixCSCRTest<float, unsigned int> cuda_sparse_matrix_cscr_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSCRTest<double, unsigned int> cuda_sparse_matrix_cscr_test_double_uint(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRSerializeTest
  : public UnitTest
{
public:
   SparseMatrixCSCRSerializeTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSCRSerializeTest()
  {
  }

  virtual void run() const override
  {
    DenseVector<DT_, IT_> val(7);
    DenseVector<IT_, IT_> col_ind(7);
    for (Index i(0) ; i < val.size() ; ++i)
    {
      val(i, DT_(i+1));
      col_ind(i, IT_(i + 3));
    }
    DenseVector<IT_, IT_> row_ptr(4);
    row_ptr(0, IT_(0));
    row_ptr(1, IT_(2));
    row_ptr(2, IT_(5));
    row_ptr(3, IT_(6));
    DenseVector<IT_, IT_> row_numbers(3);
    row_numbers(0, IT_(0));
    row_numbers(1, IT_(4));
    row_numbers(2, IT_(8));

    SparseMatrixCSCR<DT_, IT_> a(10, 10, col_ind, val, row_ptr, row_numbers);

    BinaryStream bs;
    a.write_out(FileMode::fm_cscr, bs);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));
    bs.seekg(0);
    SparseMatrixCSCR<DT_, IT_> e(FileMode::fm_cscr, bs);
    TEST_CHECK_EQUAL(e.rows(), 10);
    TEST_CHECK_EQUAL(e.columns(), 10);
    TEST_CHECK_EQUAL(e.size(), 100);
    TEST_CHECK_EQUAL(e.used_rows(), 3);
    for (Index row(0) ; row < e.rows() ; ++row)
    {
      for (Index col(0) ; col < e.columns() ; ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(e(row, col), a(row, col), Math::eps<DT_>());
      }
    }
    auto l = a.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixCSCR<DT_, IT_> tst(l);
    TEST_CHECK_EQUAL(tst.rows(), a.rows());
    TEST_CHECK_EQUAL(tst.columns(), a.columns());
    TEST_CHECK_EQUAL(tst.size(), a.size());
    TEST_CHECK_EQUAL(tst.used_rows(), a.used_rows());
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(tst(row, col), a(row, col), Math::eps<DT_>());
      }
    }
#ifdef FEAT_HAVE_ZLIB
    auto zl = a.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixCSCR<DT_, IT_> zlib(zl);
    TEST_CHECK_EQUAL(zlib.rows(), a.rows());
    TEST_CHECK_EQUAL(zlib.columns(), a.columns());
    TEST_CHECK_EQUAL(zlib.size(), a.size());
    TEST_CHECK_EQUAL(zlib.used_rows(), a.used_rows());
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(zlib(row, col), a(row, col), Math::eps<DT_>());
      }
    }
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = a.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixCSCR<DT_, IT_> zfp(zf);
    TEST_CHECK_EQUAL(zfp.rows(), a.rows());
    TEST_CHECK_EQUAL(zfp.columns(), a.columns());
    TEST_CHECK_EQUAL(zfp.size(), a.size());
    TEST_CHECK_EQUAL(zfp.used_rows(), a.used_rows());
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(zfp(row, col), a(row, col), DT_(1e-4));
      }
    }
#endif

    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));

  }
};
SparseMatrixCSCRSerializeTest<float, unsigned long> cpu_sparse_matrix_cscr_serialize_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest<double, unsigned long> cpu_sparse_matrix_cscr_serialize_test_double_ulong(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest<float, unsigned int> cpu_sparse_matrix_cscr_serialize_test_float_uint(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest<double, unsigned int> cpu_sparse_matrix_cscr_serialize_test_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSCRSerializeTest<float, unsigned long> mkl_sparse_matrix_cscr_serialize_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSCRSerializeTest<double, unsigned long> mkl_sparse_matrix_cscr_serialize_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRSerializeTest<__float128, unsigned long> sparse_matrix_cscr_serialize_test_float128_ulong(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest<__float128, unsigned int> sparse_matrix_cscr_serialize_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSCRSerializeTest<Half, unsigned int> sparse_matrix_cscr_serialize_test_half_uint(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest<Half, unsigned long> sparse_matrix_cscr_serialize_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRSerializeTest<float, unsigned long> cuda_sparse_matrix_cscr_serialize_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSCRSerializeTest<double, unsigned long> cuda_sparse_matrix_cscr_serialize_test_double_ulong(PreferredBackend::cuda);
SparseMatrixCSCRSerializeTest<float, unsigned int> cuda_sparse_matrix_cscr_serialize_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSCRSerializeTest<double, unsigned int> cuda_sparse_matrix_cscr_serialize_test_double_uint(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRApplyTest
  : public UnitTest
{
public:
   SparseMatrixCSCRApplyTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSCRApplyTest()
  {
  }

  virtual void run() const override
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size);
      DenseVector<DT_, IT_> x(size);
      DenseVector<DT_, IT_> y(size);
      DenseVector<DT_, IT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        x(i, DT_(i % 100 * DT_(1.234)));
        y(i, DT_(2 - DT_(i % 42)));
      }

      for (Index row(0) ; row < a_fac.rows() ; ++row)
      {
        for (Index col(0) ; col < a_fac.columns() ; ++col)
        {
          if(row == col)
          {
            a_fac.add(row, col, DT_(2));
          }
          else if((row == col+1) || (row+1 == col))
          {
            a_fac.add(row, col, DT_(-1));
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> acsr(a_fac.make_csr());
      DenseVector<DT_, IT_> r(size);

      DenseVector<DT_, IT_> val(acsr.used_elements(), acsr.val());
      DenseVector<IT_, IT_> col_ind(acsr.used_elements(), acsr.col_ind());
      DenseVector<IT_, IT_> row_ptr(acsr.rows() + 1, acsr.row_ptr());
      DenseVector<IT_, IT_> row_numbers(acsr.rows());
      for (Index i(0) ; i < row_numbers.size() ; ++i)
      {
        row_numbers(i, IT_(i));
      }
      SparseMatrixCSCR<DT_, IT_> a(acsr.rows(), acsr.columns(), col_ind, val, row_ptr, row_numbers);

      // apply-test for alpha = 0.0
      a.apply(r, x, y, DT_(0.0));
      ref.copy(y);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-2));

      // apply-test for alpha = -1.0
      a.apply(r, x, y, DT_(-1.0));
      a.apply(ref, x);
      ref.scale(ref, DT_(-1.0));
      ref.axpy(ref, y);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-2));

      // apply-test for alpha = -1.0 and &r==&y
      r.copy(y);
      a.apply(r, x, r, DT_(-1.0));
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-2));

      // apply-test for alpha = 4711.1
      //r.axpy(s, a, x, y);
      a.apply(r, x, y, s);
      //ref.product_matvec(a, x);
      a.apply(ref, x);
      ref.scale(ref, s);
      ref.axpy(ref, y);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(5e-2));

      // apply-test for alpha = 4711.1 and &r==&y
      r.copy(y);
      a.apply(r, x, r, s);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(5e-2));

      a.apply(r, x);
      a.apply(ref, x);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), DT_(1e-2));
    }
  }
};

SparseMatrixCSCRApplyTest<float, unsigned long> sm_cscr_apply_test_float_ulong(PreferredBackend::generic);
SparseMatrixCSCRApplyTest<double, unsigned long> sm_cscr_apply_test_double_ulong(PreferredBackend::generic);
SparseMatrixCSCRApplyTest<float, unsigned int> sm_cscr_apply_test_float_uint(PreferredBackend::generic);
SparseMatrixCSCRApplyTest<double, unsigned int> sm_cscr_apply_test_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSCRApplyTest<float, unsigned long> mkl_sm_cscr_apply_test_float_ulong(PreferredBackend::mkl);
SparseMatrixCSCRApplyTest<double, unsigned long> mkl_sm_cscr_apply_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRApplyTest<__float128, unsigned long> sm_cscr_apply_test_float128_ulong(PreferredBackend::generic);
SparseMatrixCSCRApplyTest<__float128, unsigned int> sm_cscr_apply_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSCRApplyTest<Half, unsigned int> sm_cscr_apply_test_half_uint(PreferredBackend::generic);
SparseMatrixCSCRApplyTest<Half, unsigned long> sm_cscr_apply_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRApplyTest<float, unsigned long> cuda_sm_cscr_apply_test_float_ulong(PreferredBackend::cuda);
SparseMatrixCSCRApplyTest<double, unsigned long> cuda_sm_cscr_apply_test_double_ulong(PreferredBackend::cuda);
SparseMatrixCSCRApplyTest<float, unsigned int> cuda_sm_cscr_apply_test_float_uint(PreferredBackend::cuda);
SparseMatrixCSCRApplyTest<double, unsigned int> cuda_sm_cscr_apply_test_double_uint(PreferredBackend::cuda);
#endif
