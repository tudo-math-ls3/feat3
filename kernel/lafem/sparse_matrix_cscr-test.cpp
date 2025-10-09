// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_cscr.hpp>
#include <kernel/util/binary_stream.hpp>
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
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    SparseMatrixCSCR<DT_, IT_> zero;
    TEST_CHECK(zero.empty());

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
    empty5.convert(zero);
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
    TEST_CHECK(!a.empty());
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
    TEST_CHECK_LESS_THAN(b.max_rel_diff(a), eps);
    SparseMatrixCSCR<DT_, IT_> c;
    c.convert(b);
    TEST_CHECK_LESS_THAN(c.max_rel_diff(b), eps);
    SparseMatrixCSCR<DT_, IT_> d(a.clone());
    TEST_CHECK_LESS_THAN(d.max_rel_diff(c), eps);


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

SparseMatrixCSCRTest <float, std::uint64_t> cpu_sparse_matrix_cscr_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSCRTest <double, std::uint64_t> cpu_sparse_matrix_cscr_test_double_uint64(PreferredBackend::generic);
SparseMatrixCSCRTest <float, std::uint32_t> cpu_sparse_matrix_cscr_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSCRTest <double, std::uint32_t> cpu_sparse_matrix_cscr_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSCRTest <float, std::uint64_t> mkl_sparse_matrix_cscr_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSCRTest <double, std::uint64_t> mkl_sparse_matrix_cscr_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRTest <__float128, std::uint64_t> sparse_matrix_cscr_test_float128_uint64(PreferredBackend::generic);
SparseMatrixCSCRTest <__float128, std::uint32_t> sparse_matrix_cscr_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSCRTest <Half, std::uint32_t> sparse_matrix_cscr_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSCRTest <Half, std::uint64_t> sparse_matrix_cscr_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRTest <float, std::uint64_t> cuda_sparse_matrix_cscr_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSCRTest <double, std::uint64_t> cuda_sparse_matrix_cscr_test_double_uint64(PreferredBackend::cuda);
SparseMatrixCSCRTest <float, std::uint32_t> cuda_sparse_matrix_cscr_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSCRTest <double, std::uint32_t> cuda_sparse_matrix_cscr_test_double_uint32(PreferredBackend::cuda);
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
        TEST_CHECK_EQUAL_WITHIN_EPS(zfp(row, col), a(row, col), Math::eps<DT_>());
      }
    }
#endif

    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));

  }
};
SparseMatrixCSCRSerializeTest <float, std::uint64_t> cpu_sparse_matrix_cscr_serialize_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest <double, std::uint64_t> cpu_sparse_matrix_cscr_serialize_test_double_uint64(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest <float, std::uint32_t> cpu_sparse_matrix_cscr_serialize_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest <double, std::uint32_t> cpu_sparse_matrix_cscr_serialize_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSCRSerializeTest <float, std::uint64_t> mkl_sparse_matrix_cscr_serialize_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSCRSerializeTest <double, std::uint64_t> mkl_sparse_matrix_cscr_serialize_test_double_uint64(PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SparseMatrixCSCRSerializeTest <__float128, std::uint64_t> sparse_matrix_cscr_serialize_test_float128_uint64(PreferredBackend::generic);
//SparseMatrixCSCRSerializeTest <__float128, std::uint32_t> sparse_matrix_cscr_serialize_test_float128_uint32(PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSCRSerializeTest <Half, std::uint32_t> sparse_matrix_cscr_serialize_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSCRSerializeTest <Half, std::uint64_t> sparse_matrix_cscr_serialize_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRSerializeTest <float, std::uint64_t> cuda_sparse_matrix_cscr_serialize_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSCRSerializeTest <double, std::uint64_t> cuda_sparse_matrix_cscr_serialize_test_double_uint64(PreferredBackend::cuda);
SparseMatrixCSCRSerializeTest <float, std::uint32_t> cuda_sparse_matrix_cscr_serialize_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSCRSerializeTest <double, std::uint32_t> cuda_sparse_matrix_cscr_serialize_test_double_uint32(PreferredBackend::cuda);
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
    DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    DT_ s(DT_(4.7111));
    for (Index size(1) ; size < Index(1e3) ; size*=2)
    {
      SparseMatrixFactory<DT_, IT_> a_fac(size, size);
      DenseVector<DT_, IT_> x(size);
      DenseVector<DT_, IT_> y(size);
      DenseVector<DT_, IT_> ref(size);
      for (Index i(0) ; i < size ; ++i)
      {
        x(i, DT_(i % 100) * DT_(1.234));
        y(i, DT_(2) - DT_(i % 42));
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
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-transposed-test for alpha = 0.0
      a.apply_transposed(r, x, y, DT_(0.0));
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-test for alpha = -1.0
      a.apply(r, x, y, DT_(-1.0));
      a.apply(ref, x);
      ref.scale(ref, DT_(-1.0));
      ref.axpy(y);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-transposed-test for alpha = -1.0
      a.apply_transposed(r, x, y, DT_(-1.0));
      a.apply_transposed(ref, x);
      ref.scale(ref, DT_(-1.0));
      ref.axpy(y);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_RELATIVE(r(i), ref(i), eps);

      // apply-test for alpha = -1.0 and &r==&y
      r.copy(y);
      a.apply(r, x, r, DT_(-1.0));
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-transposed-test for alpha = -1.0 and &r==&y
      r.copy(y);
      a.apply_transposed(r, x, r, DT_(-1.0));
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_RELATIVE(r(i), ref(i), eps);

      // apply-test for alpha = 4.7111
      //r.axpy(s, a, x, y);
      a.apply(r, x, y, s);
      //ref.product_matvec(a, x);
      a.apply(ref, x);
      ref.scale(ref, s);
      ref.axpy(y);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      // apply-transposed-test for alpha = 4.7111
      a.apply_transposed(r, x, y, s);
      a.apply_transposed(ref, x);
      ref.scale(ref, s);
      ref.axpy(y);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_RELATIVE(r(i), ref(i), eps);

      // apply-test for alpha = 4.7111 and &r==&y
      r.copy(y);
      a.apply(r, x, r, s);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);

      a.apply(r, x);
      a.apply(ref, x);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(r(i), ref(i), eps);
    }
  }
};

SparseMatrixCSCRApplyTest <float, std::uint64_t> sm_cscr_apply_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSCRApplyTest <double, std::uint64_t> sm_cscr_apply_test_double_uint64(PreferredBackend::generic);
SparseMatrixCSCRApplyTest <float, std::uint32_t> sm_cscr_apply_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSCRApplyTest <double, std::uint32_t> sm_cscr_apply_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSCRApplyTest <float, std::uint64_t> mkl_sm_cscr_apply_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSCRApplyTest <double, std::uint64_t> mkl_sm_cscr_apply_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRApplyTest <__float128, std::uint64_t> sm_cscr_apply_test_float128_uint64(PreferredBackend::generic);
SparseMatrixCSCRApplyTest <__float128, std::uint32_t> sm_cscr_apply_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSCRApplyTest <Half, std::uint32_t> sm_cscr_apply_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSCRApplyTest <Half, std::uint64_t> sm_cscr_apply_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRApplyTest <float, std::uint64_t> cuda_sm_cscr_apply_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSCRApplyTest <double, std::uint64_t> cuda_sm_cscr_apply_test_double_uint64(PreferredBackend::cuda);
SparseMatrixCSCRApplyTest <float, std::uint32_t> cuda_sm_cscr_apply_test_float_uint32(PreferredBackend::cuda);
SparseMatrixCSCRApplyTest <double, std::uint32_t> cuda_sm_cscr_apply_test_double_uint32(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRMaxRelDiffTest
  : public UnitTest
{
public:
  SparseMatrixCSCRMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSCRMaxRelDiffTest()
  {
  }

  virtual void run() const override
  {
    const DT_ eps = Math::pow(Math::eps<DT_>(), DT_(0.8));
    const DT_ delta = DT_(123.5);

    const Index size = 10;
    const Index diff_row = 4;
    const Index diff_col = 6;
    const DT_ initial_value = DT_(10.0);

    // ref matrix b
    SparseMatrixFactory<DT_, IT_> fac_b(size, size);

    // b(4, 6) = 10
    fac_b.add(diff_row, diff_col, initial_value);

    // convert to SparseMatrixCSCR
    SparseMatrixCSR<DT_, IT_> b_csr(size, size);
    b_csr.convert(fac_b.make_csr());
    DenseVector<DT_, IT_> val_b(b_csr.used_elements(), b_csr.val());
    DenseVector<IT_, IT_> col_ind_b(b_csr.used_elements(), b_csr.col_ind());
    DenseVector<IT_, IT_> row_ptr_b(b_csr.rows() + 1, b_csr.row_ptr());
    DenseVector<IT_, IT_> row_numbers_b(b_csr.rows());
    SparseMatrixCSCR<DT_, IT_> b(b_csr.rows(), b_csr.columns(), col_ind_b, val_b, row_ptr_b, row_numbers_b);

    // copy b into a
    SparseMatrixCSCR<DT_, IT_> a = b.clone();

    // delta matrix with delta(4, 6) = 123.5
    SparseMatrixFactory<DT_, IT_> fac_delta(size, size);
    fac_delta.add(diff_row, diff_col, delta);
    SparseMatrixCSR<DT_, IT_> delta_mat_csr(size, size);
    delta_mat_csr.convert(fac_delta.make_csr());
    DenseVector<DT_, IT_> val_delta(delta_mat_csr.used_elements(), delta_mat_csr.val());
    DenseVector<IT_, IT_> col_ind_delta(delta_mat_csr.used_elements(), delta_mat_csr.col_ind());
    DenseVector<IT_, IT_> row_ptr_delta(delta_mat_csr.rows() + 1, delta_mat_csr.row_ptr());
    DenseVector<IT_, IT_> row_numbers_delta(delta_mat_csr.rows());
    SparseMatrixCSCR<DT_, IT_> delta_mat(delta_mat_csr.rows(), delta_mat_csr.columns(), col_ind_delta, val_delta, row_ptr_delta, row_numbers_delta);


    // a = a + 1.0 * delta_mat
    a.axpy(delta_mat, DT_(1.0));

    // reference value
    const DT_ ref = delta / (DT_(2) * initial_value + delta);

    // test ||a-b||_infty
    const DT_ diff_1 = a.max_rel_diff(b);
    TEST_CHECK_RELATIVE(diff_1, ref, eps);

    // test ||b-a||_infty
    const DT_ diff_2 = b.max_rel_diff(a);
    TEST_CHECK_RELATIVE(diff_2, ref, eps);
  }
};
SparseMatrixCSCRMaxRelDiffTest <float, std::uint32_t> sm_cscr_max_rel_diff_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSCRMaxRelDiffTest <double, std::uint32_t> sm_cscr_max_rel_diff_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSCRMaxRelDiffTest <float, std::uint64_t> sm_cscr_max_rel_diff_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSCRMaxRelDiffTest <double, std::uint64_t> sm_cscr_max_rel_diff_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSCRMaxRelDiffTest <float, std::uint64_t> mkl_sm_cscr_max_rel_diff_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSCRMaxRelDiffTest <double, std::uint64_t> mkl_sm_cscr_max_rel_diff_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRMaxRelDiffTest <__float128, std::uint32_t> sm_cscr_max_rel_diff_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSCRMaxRelDiffTest <__float128, std::uint64_t> sm_cscr_max_rel_diff_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSCRMaxRelDiffTest <Half, std::uint32_t> sm_cscr_max_rel_diff_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSCRMaxRelDiffTest <Half, std::uint64_t> sm_cscr_max_rel_diff_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRMaxRelDiffTest <float, std::uint64_t> cuda_sm_cscr_max_rel_diff_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSCRMaxRelDiffTest <double, std::uint64_t> cuda_sm_cscr_max_rel_diff_test_double_uint64(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class SparseMatrixCSCRSameLayoutTest
  : public UnitTest
{
public:
  SparseMatrixCSCRSameLayoutTest(PreferredBackend backend)
    : UnitTest("SparseMatrixCSCRSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixCSCRSameLayoutTest()
  {
  }

  virtual void run() const override
  {
    const Index size = 10;
    const Index diff_row = 4;
    const Index diff_col = 6;
    const DT_ initial_value = DT_(10.0);

    // ref matrix a
    SparseMatrixFactory<DT_, IT_> fac_a(size, size);
    fac_a.add(diff_row, diff_col, initial_value);

    // convert to SparseMatrixCSCR
    SparseMatrixCSR<DT_, IT_> a_csr(size, size);
    a_csr.convert(fac_a.make_csr());
    DenseVector<DT_, IT_> val_a(a_csr.used_elements(), a_csr.val());
    DenseVector<IT_, IT_> col_ind_a(a_csr.used_elements(), a_csr.col_ind());
    DenseVector<IT_, IT_> row_ptr_a(a_csr.rows() + 1, a_csr.row_ptr());
    DenseVector<IT_, IT_> row_numbers_a(a_csr.rows());
    SparseMatrixCSCR<DT_, IT_> a(a_csr.rows(), a_csr.columns(), col_ind_a, val_a, row_ptr_a, row_numbers_a);

    // weak copy
    auto b = a.clone(CloneMode::Weak);
    TEST_CHECK(a.same_layout(b));

    // shallow copy
    auto c = a.clone(CloneMode::Shallow);
    TEST_CHECK(a.same_layout(c));

    // different values at same position
    SparseMatrixFactory<DT_, IT_> fac_d(size, size);
    fac_d.add(diff_row, diff_col, DT_(0.5));
    SparseMatrixCSR<DT_, IT_> d_csr(size, size);
    d_csr.convert(fac_d.make_csr());
    DenseVector<DT_, IT_> val_d(d_csr.used_elements(), d_csr.val());
    DenseVector<IT_, IT_> col_ind_d(d_csr.used_elements(), d_csr.col_ind());
    DenseVector<IT_, IT_> row_ptr_d(d_csr.rows() + 1, d_csr.row_ptr());
    DenseVector<IT_, IT_> row_numbers_d(d_csr.rows());
    SparseMatrixCSCR<DT_, IT_> d(d_csr.rows(), d_csr.columns(), col_ind_d, val_d, row_ptr_d, row_numbers_d);
    TEST_CHECK(a.same_layout(d));

    // value at different position
    SparseMatrixFactory<DT_, IT_> fac_e(size, size);
    fac_e.add(diff_row + 1, diff_col, initial_value);
    SparseMatrixCSR<DT_, IT_> e_csr(size, size);
    e_csr.convert(fac_e.make_csr());
    DenseVector<DT_, IT_> val_e(e_csr.used_elements(), e_csr.val());
    DenseVector<IT_, IT_> col_ind_e(e_csr.used_elements(), e_csr.col_ind());
    DenseVector<IT_, IT_> row_ptr_e(e_csr.rows() + 1, e_csr.row_ptr());
    DenseVector<IT_, IT_> row_numbers_e(e_csr.rows());
    SparseMatrixCSCR<DT_, IT_> e(e_csr.rows(), e_csr.columns(),col_ind_e, val_e, row_ptr_e, row_numbers_e);
    TEST_CHECK(!a.same_layout(e));

    // different sizes
    SparseMatrixFactory<DT_, IT_> fac_f(size + 2, size);
    fac_f.add(diff_row, diff_col, initial_value);
    SparseMatrixCSR<DT_, IT_> f_csr(size + 2, size);
    f_csr.convert(fac_f.make_csr());
    DenseVector<DT_, IT_> val_f(f_csr.used_elements(), f_csr.val());
    DenseVector<IT_, IT_> col_ind_f(f_csr.used_elements(), f_csr.col_ind());
    DenseVector<IT_, IT_> row_ptr_f(f_csr.rows() + 1, f_csr.row_ptr());
    DenseVector<IT_, IT_> row_numbers_f(f_csr.rows());
    SparseMatrixCSCR<DT_, IT_> f(f_csr.rows(), f_csr.columns(), col_ind_f, val_f, row_ptr_f, row_numbers_f);
    TEST_CHECK(!a.same_layout(f));
  }
};
SparseMatrixCSCRSameLayoutTest <float, std::uint32_t> sm_cscr_same_layout_test_float_uint32(PreferredBackend::generic);
SparseMatrixCSCRSameLayoutTest <double, std::uint32_t> sm_cscr_same_layout_test_double_uint32(PreferredBackend::generic);
SparseMatrixCSCRSameLayoutTest <float, std::uint64_t> sm_cscr_same_layout_test_float_uint64(PreferredBackend::generic);
SparseMatrixCSCRSameLayoutTest <double, std::uint64_t> sm_cscr_same_layout_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixCSCRSameLayoutTest <float, std::uint64_t> mkl_sm_cscr_same_layout_test_float_uint64(PreferredBackend::mkl);
SparseMatrixCSCRSameLayoutTest <double, std::uint64_t> mkl_sm_cscr_same_layout_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRSameLayoutTest <__float128, std::uint32_t> sm_cscr_same_layout_test_float128_uint32(PreferredBackend::generic);
SparseMatrixCSCRSameLayoutTest <__float128, std::uint64_t> sm_cscr_same_layout_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixCSCRSameLayoutTest <Half, std::uint32_t> sm_cscr_same_layout_test_half_uint32(PreferredBackend::generic);
SparseMatrixCSCRSameLayoutTest <Half, std::uint64_t> sm_cscr_same_layout_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRSameLayoutTest <float, std::uint64_t> cuda_sm_cscr_same_layout_test_float_uint64(PreferredBackend::cuda);
SparseMatrixCSCRSameLayoutTest <double, std::uint64_t> cuda_sm_cscr_same_layout_test_double_uint64(PreferredBackend::cuda);
#endif
