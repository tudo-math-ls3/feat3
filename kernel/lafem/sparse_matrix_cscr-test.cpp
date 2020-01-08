// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_cscr.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/random.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>

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
class SparseMatrixCSCRTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSCRTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSCRTest")
  {
  }

  virtual ~SparseMatrixCSCRTest()
  {
  }

  virtual void run() const override
  {
    SparseMatrixCSCR<Mem_, DT_, IT_> zero1;
    SparseMatrixCSCR<Mem::Main, DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);
    zero2.convert(zero1);

    SparseMatrixCSCR<Mem_, DT_, IT_> zero3(10, 11, 12, 3);
    TEST_CHECK_EQUAL(zero3.used_elements(), 12);
    TEST_CHECK_EQUAL(zero3.rows(), 10);
    TEST_CHECK_EQUAL(zero3.columns(), 11);
    TEST_CHECK_EQUAL(zero3.size(), 110);
    TEST_CHECK_EQUAL(zero3.used_rows(), 3);


    SparseMatrixCSCR<Mem_, DT_, IT_> empty1(11, 12, 111, 3);
    SparseMatrixCSCR<Mem::Main, DT_, IT_> empty2;
    SparseMatrixCSCR<Mem_, DT_, IT_> empty3;
    empty2.convert(empty1);
    empty3.convert(empty2);
    TEST_CHECK_EQUAL(empty1.rows(), empty3.rows());
    TEST_CHECK_EQUAL(empty1.columns(), empty3.columns());
    TEST_CHECK_EQUAL(empty1.used_elements(), empty3.used_elements());
    TEST_CHECK_EQUAL(empty1.used_rows(), empty3.used_rows());

    SparseMatrixCSCR<Mem::Main, DT_, IT_> empty4(empty2.layout());
    SparseMatrixCSCR<Mem_, DT_, IT_> empty5(empty3.layout());
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

    DenseVector<Mem_, DT_, IT_> val(7);
    DenseVector<Mem_, IT_, IT_> col_ind(7);
    for (Index i(0) ; i < val.size() ; ++i)
    {
      val(i, DT_(i+1));
      col_ind(i, IT_(i + 3));
    }
    DenseVector<Mem_, IT_, IT_> row_ptr(4);
    row_ptr(0, IT_(0));
    row_ptr(1, IT_(2));
    row_ptr(2, IT_(5));
    row_ptr(3, IT_(6));
    DenseVector<Mem_, IT_, IT_> row_numbers(3);
    row_numbers(0, IT_(0));
    row_numbers(1, IT_(4));
    row_numbers(2, IT_(8));

    SparseMatrixCSCR<Mem_, DT_, IT_> a(10, 10, col_ind, val, row_ptr, row_numbers);
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

    SparseMatrixCSCR<Mem::Main, DT_, IT_> b;
    b.convert(a);
    TEST_CHECK_EQUAL(b, a);
    SparseMatrixCSCR<Mem_, DT_, IT_> c;
    c.convert(b);
    TEST_CHECK_EQUAL(c, b);
    SparseMatrixCSCR<Mem_, DT_, IT_> d(a.clone());
    TEST_CHECK_EQUAL(d, c);

    BinaryStream bs;
    a.write_out(FileMode::fm_cscr, bs);
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));
    bs.seekg(0);
    SparseMatrixCSCR<Mem_, DT_, IT_> e(FileMode::fm_cscr, bs);
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
    //TEST_CHECK_EQUAL(bs.tellg(), std::streampos(696));

    SparseMatrixCOO<Mem_, DT_, IT_> coo(10, 10);
    coo(0, 3, DT_(1));
    coo(0, 4, DT_(2));
    coo(1, 5, DT_(7));
    coo(3, 3, DT_(8));
    coo(4, 5, DT_(3));
    coo(4, 6, DT_(4));
    coo(4, 7, DT_(5));
    coo(6, 6, DT_(9));
    coo(8, 8, DT_(6));
    SparseMatrixCSR<Mem_, DT_, IT_> csr(10, 10);
    csr.convert(coo);
    VectorMirror<Mem::Main, DT_, IT_> mirror_main(10, 3);
    mirror_main.indices()[0] = 0;
    mirror_main.indices()[1] = 4;
    mirror_main.indices()[2] = 8;
    VectorMirror<Mem_, DT_, IT_> mirror;
    mirror.convert(mirror_main);
    SparseMatrixCSCR<Mem_, DT_, IT_> f(csr, mirror);
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(f(row, col), a(row, col), Math::eps<DT_>());
      }
    }

  }
};

SparseMatrixCSCRTest<Mem::Main, float, unsigned long> cpu_sparse_matrix_cscr_test_float_ulong;
SparseMatrixCSCRTest<Mem::Main, double, unsigned long> cpu_sparse_matrix_cscr_test_double_ulong;
SparseMatrixCSCRTest<Mem::Main, float, unsigned int> cpu_sparse_matrix_cscr_test_float_uint;
SparseMatrixCSCRTest<Mem::Main, double, unsigned int> cpu_sparse_matrix_cscr_test_double_uint;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRTest<Mem::Main, __float128, unsigned long> cpu_sparse_matrix_cscr_test_float128_ulong;
SparseMatrixCSCRTest<Mem::Main, __float128, unsigned int> cpu_sparse_matrix_cscr_test_float128_uint;
#endif
/*#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRTest<Mem::CUDA, float, unsigned long> cuda_sparse_matrix_cscr_test_float_ulong;
SparseMatrixCSCRTest<Mem::CUDA, double, unsigned long> cuda_sparse_matrix_cscr_test_double_ulong;
SparseMatrixCSCRTest<Mem::CUDA, float, unsigned int> cuda_sparse_matrix_cscr_test_float_uint;
SparseMatrixCSCRTest<Mem::CUDA, double, unsigned int> cuda_sparse_matrix_cscr_test_double_uint;
#endif*/

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixCSCRApplyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixCSCRApplyTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixCSCRApplyTest")
  {
  }

  virtual ~SparseMatrixCSCRApplyTest()
  {
  }

  virtual void run() const override
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_, IT_> a_local(size, size);
      DenseVector<Mem::Main, DT_, IT_> x_local(size);
      DenseVector<Mem::Main, DT_, IT_> y_local(size);
      DenseVector<Mem::Main, DT_, IT_> ref_local(size);
      DenseVector<Mem_, DT_, IT_> ref(size);
      DenseVector<Mem::Main, DT_, IT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        x_local(i, DT_(i % 100 * DT_(1.234)));
        y_local(i, DT_(2 - DT_(i % 42)));
      }
      DenseVector<Mem_, DT_, IT_> x(size);
      x.copy(x_local);
      DenseVector<Mem_, DT_, IT_> y(size);
      y.copy(y_local);

      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
          {
            a_local(row, col, DT_(2));
          }
          else if((row == col+1) || (row+1 == col))
          {
            a_local(row, col, DT_(-1));
          }
        }
      }
      SparseMatrixCSR<Mem_, DT_, IT_> acsr(a_local);
      DenseVector<Mem_, DT_, IT_> r(size);

      DenseVector<Mem_, DT_, IT_> val(acsr.used_elements(), acsr.val());
      DenseVector<Mem_, IT_, IT_> col_ind(acsr.used_elements(), acsr.col_ind());
      DenseVector<Mem_, IT_, IT_> row_ptr(acsr.rows() + 1, acsr.row_ptr());
      DenseVector<Mem::Main, IT_, IT_> row_numbers_local(acsr.rows());
      for (Index i(0) ; i < row_numbers_local.size() ; ++i)
      {
        row_numbers_local(i, IT_(i));
      }
      DenseVector<Mem_, IT_, IT_> row_numbers;
      row_numbers.convert(row_numbers_local);
      SparseMatrixCSCR<Mem_, DT_, IT_> a(acsr.rows(), acsr.columns(), col_ind, val, row_ptr, row_numbers);

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
      a_local.apply(ref_local, x_local);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(1e-2));
    }
  }
};

SparseMatrixCSCRApplyTest<Mem::Main, float, unsigned long> sm_cscr_apply_test_float_ulong;
SparseMatrixCSCRApplyTest<Mem::Main, double, unsigned long> sm_cscr_apply_test_double_ulong;
SparseMatrixCSCRApplyTest<Mem::Main, float, unsigned int> sm_cscr_apply_test_float_uint;
SparseMatrixCSCRApplyTest<Mem::Main, double, unsigned int> sm_cscr_apply_test_double_uint;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixCSCRApplyTest<Mem::Main, __float128, unsigned long> sm_cscr_apply_test_float128_ulong;
SparseMatrixCSCRApplyTest<Mem::Main, __float128, unsigned int> sm_cscr_apply_test_float128_uint;
#endif
/*#ifdef FEAT_HAVE_CUDA
SparseMatrixCSCRApplyTest<Mem::CUDA, float, unsigned long> cuda_sm_cscr_apply_test_float_ulong;
SparseMatrixCSCRApplyTest<Mem::CUDA, double, unsigned long> cuda_sm_cscr_apply_test_double_ulong;
SparseMatrixCSCRApplyTest<Mem::CUDA, float, unsigned int> cuda_sm_cscr_apply_test_float_uint;
SparseMatrixCSCRApplyTest<Mem::CUDA, double, unsigned int> cuda_sm_cscr_apply_test_double_uint;
#endif*/
