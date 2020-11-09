// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/random.hpp>

#include <cstdio>
#include <sstream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the sparse matrix ell class.
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
class SparseMatrixELLTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLTest")
  {
  }

  virtual ~SparseMatrixELLTest()
  {
  }

  virtual void run() const override
  {
    SparseMatrixELL<Mem_, DT_, IT_> zero1;
    SparseMatrixELL<Mem::Main, DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);
    zero2.convert(zero1);

    SparseMatrixCOO<Mem::Main, DT_, IT_> a(10, 12);
    a(1,2,7);
    a.format();
    a(1,2,7);
    a(5,5,2);
    SparseMatrixELL<Mem_, DT_, IT_> b(a);
    TEST_CHECK_EQUAL(b.used_elements(), a.used_elements());
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(b(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(b(5, 5), a(5, 5));

    SparseMatrixELL<Mem_, DT_, IT_> bl(b.layout());
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    bl = b.layout();
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    SparseMatrixELL<Mem_, DT_, IT_> z;
    z.convert(b);
    TEST_CHECK_EQUAL(z.used_elements(), 2ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z.C(), b.C());
    TEST_CHECK_EQUAL(z.num_of_chunks(), b.num_of_chunks());
    TEST_CHECK_EQUAL(z.val_size(), b.val_size());
    TEST_CHECK_EQUAL(z(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(z(5, 5), a(5, 5));
    TEST_CHECK_EQUAL(z(1, 3), a(1, 3));

    SparseMatrixELL<Mem::Main, DT_, IT_> e;
    e.convert(b);
    TEST_CHECK_EQUAL(e, b);
    e.copy(b);
    TEST_CHECK_EQUAL(e, b);

    SparseMatrixELL<Mem_, DT_, IT_> c;
    c.clone(b);
    TEST_CHECK_EQUAL(c, b);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_EQUAL((void*)c.col_ind(), (void*)b.col_ind());
    c = b.clone(CloneMode::Deep);
    TEST_CHECK_EQUAL(c, b);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_NOT_EQUAL((void*)c.col_ind(), (void*)b.col_ind());

    decltype(b) y(b.layout());
    TEST_CHECK_EQUAL((void*)y.cs(), (void*)b.cs());
    TEST_CHECK_EQUAL((void*)y.cl(), (void*)b.cl());
    TEST_CHECK_EQUAL((void*)y.rl(), (void*)b.rl());
  }
};

SparseMatrixELLTest<Mem::Main, float, unsigned long> cpu_sparse_matrix_ell_test_float_ulong;
SparseMatrixELLTest<Mem::Main, double, unsigned long> cpu_sparse_matrix_ell_test_double_ulong;
SparseMatrixELLTest<Mem::Main, float, unsigned int> cpu_sparse_matrix_ell_test_float_uint;
SparseMatrixELLTest<Mem::Main, double, unsigned int> cpu_sparse_matrix_ell_test_double_uint;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLTest<Mem::Main, __float128, unsigned long> cpu_sparse_matrix_ell_test_float128_ulong;
SparseMatrixELLTest<Mem::Main, __float128, unsigned int> cpu_sparse_matrix_ell_test_float128_uint;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLTest<Mem::CUDA, float, unsigned long> cuda_sparse_matrix_ell_test_float_ulong;
SparseMatrixELLTest<Mem::CUDA, double, unsigned long> cuda_sparse_matrix_ell_test_double_ulong;
SparseMatrixELLTest<Mem::CUDA, float, unsigned int> cuda_sparse_matrix_ell_test_float_uint;
SparseMatrixELLTest<Mem::CUDA, double, unsigned int> cuda_sparse_matrix_ell_test_double_uint;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixELLSerializeTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLSerializeTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLSerializeTest")
  {
  }

  virtual ~SparseMatrixELLSerializeTest()
  {
  }

  virtual void run() const override
  {
    SparseMatrixCOO<Mem::Main, DT_, IT_> fcoo(10, 10);
    for (Index row(0) ; row < fcoo.rows() ; ++row)
    {
      for (Index col(0) ; col < fcoo.columns() ; ++col)
      {
        if(row == col)
          fcoo(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          fcoo(row, col, DT_(-1));
      }
    }
    SparseMatrixELL<Mem_, DT_, IT_> f(fcoo);

    BinaryStream bs;
    f.write_out(FileMode::fm_ell, bs);
    bs.seekg(0);
    SparseMatrixELL<Mem_, DT_, IT_> g(FileMode::fm_ell, bs);
    TEST_CHECK_EQUAL(g, f);

    std::stringstream ts;
    f.write_out(FileMode::fm_mtx, ts);
    SparseMatrixELL<Mem::Main, DT_, IT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_EQUAL(j, f);

    auto kp = f.serialize(LAFEM::SerialConfig(false, false));
    SparseMatrixELL<Mem_, DT_, IT_> k(kp);
    TEST_CHECK_EQUAL(k, f);
#ifdef FEAT_HAVE_ZLIB
    auto zl = f.serialize(LAFEM::SerialConfig(true, false));
    SparseMatrixELL<Mem_, DT_, IT_> zlib(zl);
    TEST_CHECK_EQUAL(zlib, f);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zf = f.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-7)));
    SparseMatrixELL<Mem_, DT_, IT_> zfp(zf);
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
SparseMatrixELLSerializeTest<Mem::Main, float, unsigned long> cpu_sparse_matrix_ell_serialize_test_float_ulong;
SparseMatrixELLSerializeTest<Mem::Main, double, unsigned long> cpu_sparse_matrix_ell_serialize_test_double_ulong;
SparseMatrixELLSerializeTest<Mem::Main, float, unsigned int> cpu_sparse_matrix_ell_serialize_test_float_uint;
SparseMatrixELLSerializeTest<Mem::Main, double, unsigned int> cpu_sparse_matrix_ell_serialize_test_double_uint;
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLSerializeTest<Mem::CUDA, float, unsigned long> cuda_sparse_matrix_ell_serialize_test_float_ulong;
SparseMatrixELLSerializeTest<Mem::CUDA, double, unsigned long> cuda_sparse_matrix_ell_serialize_test_double_ulong;
SparseMatrixELLSerializeTest<Mem::CUDA, float, unsigned int> cuda_sparse_matrix_ell_serialize_test_float_uint;
SparseMatrixELLSerializeTest<Mem::CUDA, double, unsigned int> cuda_sparse_matrix_ell_serialize_test_double_uint;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixELLApplyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLApplyTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLApplyTest")
  {
  }

  virtual ~SparseMatrixELLApplyTest()
  {
  }

  virtual void run() const override
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < Index(1e3) ; size*=2)
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
      SparseMatrixELL<Mem_,DT_, IT_> a(a_local);

      DenseVector<Mem_, DT_, IT_> r(size);

      // apply-test for alpha = 0.0
      a.apply(r, x, y, DT_(0.0));
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(y(i), r(i), DT_(1e-2));

      // apply-test for alpha = -1.0
      a.apply(r, x, y, DT_(-1.0));
      a.apply(ref, x);
      ref.scale(ref, DT_(-1.0));
      ref.axpy(ref, y);

      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(ref(i), r(i), DT_(1e-2));

      // apply-test for alpha = 4711.1
      //r.gaxpyg(s, a, x, y);
      a.apply(r, x, y, s);
      result_local.copy(r);

      //ref.gproduct_matvecg(a, x);
      a.apply(ref, x);
      ref.scale(ref, s);
      ref.axpy(ref, y);
      ref_local.copy(ref);

      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(1e-2));

      a.apply(r, x);
      result_local.copy(r);
      a_local.apply(ref_local, x_local);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(1e-2));
    }
  }
};

SparseMatrixELLApplyTest<Mem::Main, float, unsigned long> sm_ell_apply_test_float_ulong;
SparseMatrixELLApplyTest<Mem::Main, double, unsigned long> sm_ell_apply_test_double_ulong;
SparseMatrixELLApplyTest<Mem::Main, float, unsigned int> sm_ell_apply_test_float_uint;
SparseMatrixELLApplyTest<Mem::Main, double, unsigned int> sm_ell_apply_test_double_uint;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLApplyTest<Mem::Main, __float128, unsigned long> sm_ell_apply_test_float128_ulong;
SparseMatrixELLApplyTest<Mem::Main, __float128, unsigned int> sm_ell_apply_test_float128_uint;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLApplyTest<Mem::CUDA, float, unsigned long> cuda_sm_ell_apply_test_float_ulong;
SparseMatrixELLApplyTest<Mem::CUDA, double, unsigned long> cuda_sm_ell_apply_test_double_ulong;
SparseMatrixELLApplyTest<Mem::CUDA, float, unsigned int> cuda_sm_ell_apply_test_float_uint;
SparseMatrixELLApplyTest<Mem::CUDA, double, unsigned int> cuda_sm_ell_apply_test_double_uint;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixELLScaleTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLScaleTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLScaleTest")
  {
  }

  virtual ~SparseMatrixELLScaleTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=2)
    {
      DT_ s(DT_(4.321));

      SparseMatrixCOO<Mem::Main, DT_, IT_> a_local(size, size + 2);
      SparseMatrixCOO<Mem::Main, DT_, IT_> ref_local(size, size + 2);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));

          if(row == col)
            ref_local(row, col, DT_(2) * s);
          else if((row == col+1) || (row+1 == col))
            ref_local(row, col, DT_(-1) * s);
        }
      }

      SparseMatrixELL<Mem_, DT_, IT_> a(a_local);
      SparseMatrixELL<Mem_, DT_, IT_> b(a.clone());

      b.scale(a, s);
      SparseMatrixCOO<Mem::Main, DT_, IT_> b_local(b);
      TEST_CHECK_EQUAL(b_local, ref_local);

      a.scale(a, s);
      SparseMatrixCOO<Mem_, DT_, IT_> a_coo(a);
      a_local.convert(a_coo);
      TEST_CHECK_EQUAL(a_local, ref_local);
    }
  }
};

SparseMatrixELLScaleTest<Mem::Main, float, unsigned int> sm_ell_scale_test_float_uint;
SparseMatrixELLScaleTest<Mem::Main, double, unsigned int> sm_ell_scale_test_double_uint;
SparseMatrixELLScaleTest<Mem::Main, float, unsigned long> sm_ell_scale_test_float_ulong;
SparseMatrixELLScaleTest<Mem::Main, double, unsigned long> sm_ell_scale_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLScaleTest<Mem::Main, __float128, unsigned int> sm_ell_scale_test_float128_uint;
SparseMatrixELLScaleTest<Mem::Main, __float128, unsigned long> sm_ell_scale_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLScaleTest<Mem::CUDA, float, unsigned int> cuda_sm_ell_scale_test_float_uint;
SparseMatrixELLScaleTest<Mem::CUDA, double, unsigned int> cuda_sm_ell_scale_test_double_uint;
SparseMatrixELLScaleTest<Mem::CUDA, float, unsigned long> cuda_sm_ell_scale_test_float_ulong;
SparseMatrixELLScaleTest<Mem::CUDA, double, unsigned long> cuda_sm_ell_scale_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixELLScaleRowColTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLScaleRowColTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLScaleRowColTest")
  {
  }

  virtual ~SparseMatrixELLScaleRowColTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(2e2) ; size*=3)
    {
      const DT_ pi(Math::pi<DT_>());
      const DT_ eps(Math::pow(Math::eps<DT_>(), DT_(0.8)));

      SparseMatrixCOO<Mem::Main, DT_, IT_> a_local(size, size + 2);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
        }
      }

      SparseMatrixELL<Mem_, DT_, IT_> a(a_local);
      SparseMatrixELL<Mem_, DT_, IT_> b(a.clone());

      // Scale rows
      DenseVector<Mem_, DT_, IT_> s1(a.rows());
      for (Index i(0); i < s1.size(); ++i)
      {
        s1(i, pi * (i % 3 + 1) - DT_(5.21) + DT_(i));
      }
      b.scale_rows(b, s1);

      SparseMatrixELL<Mem::Main, DT_, IT_> b_local(b.clone());
      for (Index row(0) ; row < a.rows() ; ++row)
      {
        for (Index col(0) ; col < a.columns() ; ++col)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(b_local(row, col), a_local(row, col) * s1(row), eps);
        }
      }

      // Scale rows
      DenseVector<Mem_, DT_, IT_> s2(a.columns());
      for (Index i(0); i < s2.size(); ++i)
      {
        s2(i, pi * (i % 3 + 1) - DT_(5.21) + DT_(i));
      }
      b.scale_cols(a, s2);

      b_local.convert(b);
      for (Index row(0) ; row < a.rows() ; ++row)
      {
        for (Index col(0) ; col < a.columns() ; ++col)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(b_local(row, col), a_local(row, col) * s2(col), eps);
        }
      }
    }
  }
};

SparseMatrixELLScaleRowColTest<Mem::Main, float, unsigned int> sm_ell_scale_row_col_test_float_uint;
SparseMatrixELLScaleRowColTest<Mem::Main, double, unsigned int> sm_ell_scale_row_col_test_double_uint;
SparseMatrixELLScaleRowColTest<Mem::Main, float, unsigned long> sm_ell_scale_row_col_test_float_ulong;
SparseMatrixELLScaleRowColTest<Mem::Main, double, unsigned long> sm_ell_scale_row_col_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLScaleRowColTest<Mem::Main, __float128, unsigned int> sm_ell_scale_row_col_test_float128_uint;
SparseMatrixELLScaleRowColTest<Mem::Main, __float128, unsigned long> sm_ell_scale_row_col_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLScaleRowColTest<Mem::CUDA, float, unsigned int> cuda_sm_ell_scale_row_col_test_float_uint;
SparseMatrixELLScaleRowColTest<Mem::CUDA, double, unsigned int> cuda_sm_ell_scale_row_col_test_double_uint;
SparseMatrixELLScaleRowColTest<Mem::CUDA, float, unsigned long> cuda_sm_ell_scale_row_col_test_float_ulong;
SparseMatrixELLScaleRowColTest<Mem::CUDA, double, unsigned long> cuda_sm_ell_scale_row_col_test_double_ulong;
#endif


/**
 * \brief Test class for the transposition of a SparseMatrixELL
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
class SparseMatrixELLTranspositionTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{

public:
  typedef SparseMatrixELL<Mem_, DT_, IT_> MatrixType;

   SparseMatrixELLTranspositionTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLTranspositionTest")
  {
  }

  virtual ~SparseMatrixELLTranspositionTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=4)
    {
      SparseMatrixCOO<Mem::Main, DT_, IT_> a_local(size, size + 2);

      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if(row == col+1)
            a_local(row, col, DT_(-1));
          else if(row+1 == col)
            a_local(row, col, DT_(-3));
        }
      }
      MatrixType a;
      a.convert(a_local);

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

SparseMatrixELLTranspositionTest<Mem::Main, float, unsigned int> sm_ell_transposition_test_float_uint;
SparseMatrixELLTranspositionTest<Mem::Main, double, unsigned int> sm_ell_transposition_test_double_uint;
SparseMatrixELLTranspositionTest<Mem::Main, float, unsigned long> sm_ell_transposition_test_float_ulong;
SparseMatrixELLTranspositionTest<Mem::Main, double, unsigned long> sm_ell_transposition_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLTranspositionTest<Mem::Main, __float128, unsigned int> sm_ell_transposition_test_float128_uint;
SparseMatrixELLTranspositionTest<Mem::Main, __float128, unsigned long> sm_ell_transposition_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLTranspositionTest<Mem::CUDA, float, unsigned int> cuda_sm_ell_transposition_test_float_uint;
SparseMatrixELLTranspositionTest<Mem::CUDA, double, unsigned int> cuda_sm_ell_transposition_test_double_uint;
SparseMatrixELLTranspositionTest<Mem::CUDA, float, unsigned long> cuda_sm_ell_transposition_test_float_ulong;
SparseMatrixELLTranspositionTest<Mem::CUDA, double, unsigned long> cuda_sm_ell_transposition_test_double_ulong;
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixELLDiagTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLDiagTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLDiagTest")
  {
  }

  virtual ~SparseMatrixELLDiagTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_, IT_> a_local(size, size);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(DT_(col % 100) / DT_(2)));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
        }
      }

      SparseMatrixELL<Mem_, DT_, IT_> a(a_local);

      auto ref = a.create_vector_l();
      auto ref_local = a_local.create_vector_l();
      for (Index i(0) ; i < a_local.rows() ; ++i)
      {
        ref_local(i, a_local(i, i));
      }
      ref.convert(ref_local);

      auto diag = a.extract_diag();
      TEST_CHECK_EQUAL(diag, ref);
    }
  }
};

SparseMatrixELLDiagTest<Mem::Main, float, unsigned int> sm_ell_diag_test_float_uint;
SparseMatrixELLDiagTest<Mem::Main, double, unsigned int> sm_ell_diag_test_double_uint;
SparseMatrixELLDiagTest<Mem::Main, float, unsigned long> sm_ell_diag_test_float_ulong;
SparseMatrixELLDiagTest<Mem::Main, double, unsigned long> sm_ell_diag_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLDiagTest<Mem::Main, __float128, unsigned int> sm_ell_diag_test_float128_uint;
SparseMatrixELLDiagTest<Mem::Main, __float128, unsigned long> sm_ell_diag_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLDiagTest<Mem::CUDA, float, unsigned int> cuda_sm_ell_diag_test_float_uint;
SparseMatrixELLDiagTest<Mem::CUDA, double, unsigned int> cuda_sm_ell_diag_test_double_uint;
SparseMatrixELLDiagTest<Mem::CUDA, float, unsigned long> cuda_sm_ell_diag_test_float_ulong;
SparseMatrixELLDiagTest<Mem::CUDA, double, unsigned long> cuda_sm_ell_diag_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixELLAxpyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLAxpyTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLAxpyTest")
  {
  }

  virtual ~SparseMatrixELLAxpyTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=2)
    {
      DT_ s(DT_(4.321));

      SparseMatrixCOO<Mem::Main, DT_, IT_> a_local(size, size + 2);
      SparseMatrixCOO<Mem::Main, DT_, IT_> b_local(size, size + 2);
      SparseMatrixCOO<Mem::Main, DT_, IT_> ref_local(size, size + 2);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));

          if((row == col+1) || (row+1 == col) || row==col)
          {
            b_local(row, col, DT_((row+col) % 15));
            ref_local(row, col, a_local(row, col) * s + b_local(row, col));
          }
        }
      }

      SparseMatrixELL<Mem_, DT_, IT_> ref(ref_local);
      SparseMatrixELL<Mem_, DT_, IT_> a(a_local);
      SparseMatrixELL<Mem_, DT_, IT_> b(b_local);
      SparseMatrixELL<Mem_, DT_, IT_> c;

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
      for (Index i(0) ; i < c.used_elements() ; ++i)
      {
        ref_local.val()[i] = a_local.val()[i] + b_local.val()[i];
      }
      ref.convert(ref_local);
      c.axpy(a, b, s);
      TEST_CHECK_EQUAL(c, ref);

      s = DT_(-1);
      for (Index i(0) ; i < c.used_elements() ; ++i)
      {
        ref_local.val()[i] = b_local.val()[i] - a_local.val()[i];
      }
      ref.convert(ref_local);
      c.axpy(a, b, s);
      TEST_CHECK_EQUAL(c, ref);
    }
  }
};

SparseMatrixELLAxpyTest<Mem::Main, float, unsigned int> sm_ell_axpy_test_float_uint;
SparseMatrixELLAxpyTest<Mem::Main, double, unsigned int> sm_ell_axpy_test_double_uint;
SparseMatrixELLAxpyTest<Mem::Main, float, unsigned long> sm_ell_axpy_test_float_ulong;
SparseMatrixELLAxpyTest<Mem::Main, double, unsigned long> sm_ell_axpy_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLAxpyTest<Mem::Main, __float128, unsigned int> sm_ell_axpy_test_float128_uint;
SparseMatrixELLAxpyTest<Mem::Main, __float128, unsigned long> sm_ell_axpy_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLAxpyTest<Mem::CUDA, float, unsigned int> cuda_sm_ell_axpy_test_float_uint;
SparseMatrixELLAxpyTest<Mem::CUDA, double, unsigned int> cuda_sm_ell_axpy_test_double_uint;
SparseMatrixELLAxpyTest<Mem::CUDA, float, unsigned long> cuda_sm_ell_axpy_test_float_ulong;
SparseMatrixELLAxpyTest<Mem::CUDA, double, unsigned long> cuda_sm_ell_axpy_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixELLFrobeniusTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLFrobeniusTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLFrobeniusTest")
  {
  }

  virtual ~SparseMatrixELLFrobeniusTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(2) ; size < Index(3e2) ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_, IT_> a_local(size, size + 2);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
        }
      }
      DenseVector<Mem::Main, DT_, IT_> refv(a_local.used_elements(), a_local.val());
      DT_ ref = refv.norm2();

      SparseMatrixELL<Mem_, DT_, IT_> a(a_local);

      DT_ c = a.norm_frobenius();
      TEST_CHECK_EQUAL_WITHIN_EPS(c, ref, DT_(1e-5));
    }
  }
};

SparseMatrixELLFrobeniusTest<Mem::Main, float, unsigned int> sm_ell_frobenius_test_float_uint;
SparseMatrixELLFrobeniusTest<Mem::Main, double, unsigned int> sm_ell_frobenius_test_double_uint;
SparseMatrixELLFrobeniusTest<Mem::Main, float, unsigned long> sm_ell_frobenius_test_float_ulong;
SparseMatrixELLFrobeniusTest<Mem::Main, double, unsigned long> sm_ell_frobenius_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLFrobeniusTest<Mem::Main, __float128, unsigned int> sm_ell_frobenius_test_float128_uint;
SparseMatrixELLFrobeniusTest<Mem::Main, __float128, unsigned long> sm_ell_frobenius_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLFrobeniusTest<Mem::CUDA, float, unsigned int> cuda_sm_ell_frobenius_test_float_uint;
SparseMatrixELLFrobeniusTest<Mem::CUDA, double, unsigned int> cuda_sm_ell_frobenius_test_double_uint;
SparseMatrixELLFrobeniusTest<Mem::CUDA, float, unsigned long> cuda_sm_ell_frobenius_test_float_ulong;
SparseMatrixELLFrobeniusTest<Mem::CUDA, double, unsigned long> cuda_sm_ell_frobenius_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixELLLumpTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   SparseMatrixELLLumpTest()
    : FullTaggedTest<Mem_, DT_, IT_>("SparseMatrixELLLumpTest")
  {
  }

  virtual ~SparseMatrixELLLumpTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.8));

    for (Index size(2) ; size < Index(3e2) ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_, IT_> a_local(size, size);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(DT_(col % 100) / DT_(2)));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
        }
      }

      SparseMatrixELL<Mem_, DT_, IT_> a(a_local);
      auto lump = a.lump_rows();
      auto one = a.create_vector_r();
      auto res = a.create_vector_r();
      one.format(DT_(-1));
      a.apply(res, one, lump);

      TEST_CHECK_EQUAL_WITHIN_EPS(res.norm2(), DT_(0), tol);
    }
  }
};

SparseMatrixELLLumpTest<Mem::Main, float, unsigned int> sm_ell_lump_test_float_uint;
SparseMatrixELLLumpTest<Mem::Main, double, unsigned int> sm_ell_lump_test_double_uint;
SparseMatrixELLLumpTest<Mem::Main, float, unsigned long> sm_ell_lump_test_float_ulong;
SparseMatrixELLLumpTest<Mem::Main, double, unsigned long> sm_ell_lump_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixELLLumpTest<Mem::Main, __float128, unsigned int> sm_ell_lump_test_float128_uint;
SparseMatrixELLLumpTest<Mem::Main, __float128, unsigned long> sm_ell_lump_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixELLLumpTest<Mem::CUDA, float, unsigned int> cuda_sm_ell_lump_test_float_uint;
SparseMatrixELLLumpTest<Mem::CUDA, double, unsigned int> cuda_sm_ell_lump_test_double_uint;
SparseMatrixELLLumpTest<Mem::CUDA, float, unsigned long> cuda_sm_ell_lump_test_float_ulong;
SparseMatrixELLLumpTest<Mem::CUDA, double, unsigned long> cuda_sm_ell_lump_test_double_ulong;
#endif
