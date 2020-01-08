// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_matrix.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the dense matrix class.
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
class DenseMatrixTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{

public:

  DenseMatrixTest()
    : FullTaggedTest<Mem_, DT_, IT_>("dense_matrix_test")
  {
  }

  virtual ~DenseMatrixTest()
  {
  }

  virtual void run() const override
  {
    DenseMatrix<Mem_, DT_, IT_> zero1;
    DenseMatrix<Mem::Main, DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseMatrix<Mem_, DT_, IT_> a(10, 11);
    TEST_CHECK_EQUAL(a.rows(), 10);
    TEST_CHECK_EQUAL(a.columns(), 11);
    TEST_CHECK_EQUAL(a.size() , 110);
    DenseMatrix<Mem_, DT_, IT_> b(10, 10, 5.);
    b(7, 6, DT_(42));
    DenseMatrix<Mem_, DT_, IT_> c;
    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c.rows(), b.rows());
    TEST_CHECK_EQUAL(c(7,6), b(7,6));
    TEST_CHECK_EQUAL(c, b);

    DenseMatrix<Mem_, DT_, IT_> e(11, 12, 5.);
    TEST_CHECK_EQUAL(e.rows(), 11ul);
    TEST_CHECK_EQUAL(e.columns(), 12ul);

    DenseMatrix<Mem_, DT_, IT_> h;
    h.clone(c);
    TEST_CHECK_EQUAL(h, c);
    h(1,2,3);
    TEST_CHECK_NOT_EQUAL(h, c);
    TEST_CHECK_NOT_EQUAL((void *)h.elements(), (void *)c.elements());

    auto kp = b.serialise();
    DenseMatrix<Mem_, DT_, IT_> k(kp);
    TEST_CHECK_EQUAL(k, b);

    // Lehmer matrix inverse test
    DenseMatrix<Mem_, DT_, IT_> l(11, 11);
    for(Index i(0) ; i < l.rows() ; ++i)
    {
      for(Index j(0) ; j < l.rows() ; ++j)
      {
        l(i,j, DT_(Math::min(i, j) + 1) / DT_(Math::max(i, j) + 1));
      }
    }
    auto m = l.inverse();
    m.invert();
    for(Index i(0) ; i < l.rows() ; ++i)
    {
      for(Index j(0) ; j < l.rows() ; ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(m(i, j), l(i, j), Math::template eps<DT_>() * DT_(100));
      }
    }

   // DenseMatrix write_out test
    DenseMatrix<Mem_, DT_, IT_> u(11, 11);
    for(Index i(0) ; i < u.rows() ; ++i)
    {
      for(Index j(0) ; j < u.columns() ; ++j)
      {
        u(i, j, DT_(i + j + 1));
      }
    }
    //Binary Test
    {
      BinaryStream bs;
      u.write_out(FileMode::fm_dm, bs);
      bs.seekg(0);
      DenseMatrix<Mem_, DT_, IT_> test(FileMode::fm_dm, bs);
      for(Index i(0) ; i < u.rows() ; ++i)
      {
        for(Index j(0) ; j < u.columns() ; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(u(i, j), test(i, j), Math::template eps<DT_>() * DT_(100));
        }
      }
    }
    //Mtx Test
    {
      std::stringstream ts;
      u.write_out(FileMode::fm_mtx, ts);
      DenseMatrix<Mem_, DT_, IT_> test2(FileMode::fm_mtx, ts);
      for(Index i(0) ; i < u.rows() ; ++i)
      {
        for(Index j(0) ; j < u.columns() ; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(u(i, j), test2(i, j), Math::template eps<DT_>() * DT_(100));
        }
      }
    }
    /*
    //FileTest-> for now... problem if write rights arent given...
    {
      String filename = "test_dense_matrix_file_bin.txt";
      std::ofstream (filename.c_str());
      u.write_out(FileMode::fm_dm, filename);
      DenseMatrix<Mem_, DT_, IT_> test3(FileMode::fm_dm, filename);
      for(Index i(0) ; i < u.rows() ; ++i)
      {
        for(Index j(0) ; j < u.columns() ; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(u(i, j), test3(i, j), Math::template eps<DT_>() * DT_(100));
        }
      }
      std::remove(filename.c_str());
    }
    {
      String filename = "test_dense_matrix_file_mtx.txt";
      std::ofstream (filename.c_str());
      u.write_out(FileMode::fm_mtx, filename);
      DenseMatrix<Mem_, DT_, IT_> test4(FileMode::fm_mtx, filename);
      for(Index i(0) ; i < u.rows() ; ++i)
      {
        for(Index j(0) ; j < u.columns() ; ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(u(i, j), test4(i, j), Math::template eps<DT_>() * DT_(100));
        }
      }
      std::remove(filename.c_str());
    }
    */
  }
};
DenseMatrixTest<Mem::Main, float, unsigned int> cpu_dense_matrix_test_float_uint;
DenseMatrixTest<Mem::Main, double, unsigned int> cpu_dense_matrix_test_double_uint;
DenseMatrixTest<Mem::Main, float, unsigned long> cpu_dense_matrix_test_float_ulong;
DenseMatrixTest<Mem::Main, double, unsigned long> cpu_dense_matrix_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixTest<Mem::Main, __float128, unsigned int> cpu_dense_matrix_test_float128_uint;
DenseMatrixTest<Mem::Main, __float128, unsigned long> cpu_dense_matrix_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixTest<Mem::CUDA, float, unsigned int> cuda_dense_matrix_test_float_uint;
DenseMatrixTest<Mem::CUDA, double, unsigned int> cuda_dense_matrix_test_double_uint;
DenseMatrixTest<Mem::CUDA, float, unsigned long> cuda_dense_matrix_test_float_ulong;
DenseMatrixTest<Mem::CUDA, double, unsigned long> cuda_dense_matrix_test_double_ulong;
#endif


template<
  typename Mem_,
  typename DT_,
  typename IT_>
class DenseMatrixApplyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  double _eps;

  explicit DenseMatrixApplyTest(double eps)
    : FullTaggedTest<Mem_, DT_, IT_>("DenseMatrixApplyTest"),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixApplyTest()
  {
  }

  virtual void run() const override
  {
    DT_ s(DT_(0.123));
    for (Index size(1) ; size < 100 ; size*=2)
    {
      DenseMatrix<Mem::Main, DT_, IT_> a_local(size, size + 1, DT_(0));
      DenseVector<Mem::Main, DT_, IT_> x_local(size + 1);

      DenseVector<Mem::Main, DT_, IT_> y_local(size);
      DenseVector<Mem::Main, DT_, IT_> ref_local(size);
      DenseVector<Mem_, DT_, IT_> ref(size);
      DenseVector<Mem::Main, DT_, IT_> result_local(size);

      for (Index i(0) ; i < x_local.size() ; ++i)
      {
        x_local(i, DT_(DT_(1) / (DT_(1) + DT_(i % 100) * DT_(1.234))));
      }
      for (Index i(0) ; i < y_local.size() ; ++i)
      {
        y_local(i, DT_(DT_(2) - DT_(i % 42)));
      }
      DenseVector<Mem_, DT_, IT_> x(x_local.size());
      x.copy(x_local);
      DenseVector<Mem_, DT_, IT_> y(y_local.size());
      y.copy(y_local);

      for (Index i(0) ; i < a_local.size() ; ++i)
      {
        a_local.elements()[i] = DT_(DT_(i % 100) * DT_(1.234));
      }
      DenseMatrix<Mem_, DT_, IT_> a(a_local.rows(), a_local.columns());
      a.copy(a_local);

      DenseVector<Mem_, DT_, IT_> r(result_local.size());

      // apply-test for alpha = 0.0
      a.apply(r, x, y, DT_(0.0));
      result_local.copy(r);
      for (Index i(0) ; i < r.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(y_local(i), result_local(i), DT_(_eps));

      a.apply(r, x);
      result_local.copy(r);
      for (Index i(0) ; i < a_local.rows() ; ++i)
      {
        DT_ sum(0);
        for (Index j(0) ; j < a_local.columns() ; ++j)
        {
          sum += a_local(i, j) * x_local(j);
        }
        ref_local(i, sum);
      }
      for (Index i(0) ; i < r.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(_eps));

      // apply-test for alpha = -1.0
      a.apply(r, x, y, DT_(-1.0));
      result_local.copy(r);
      for (Index i(0) ; i < a_local.rows() ; ++i)
      {
        DT_ sum(0);
        for (Index j(0) ; j < a_local.columns() ; ++j)
        {
          sum += a_local(i, j) * x_local(j);
        }
        ref_local(i, y_local(i) - sum);
      }

      for (Index i(0) ; i < r.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(_eps));

      // apply-test for s = 0.123
      a.apply(r, x, y, s);
      result_local.copy(r);

      a.apply(ref, x);
      ref.axpy(ref, y, s);
      ref_local.copy(ref);

      for (Index i(0) ; i < r.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), DT_(_eps));
    }
  }
};

DenseMatrixApplyTest<Mem::Main, float, unsigned int> dense_matrix_apply_test_float_uint(1e-3);
DenseMatrixApplyTest<Mem::Main, double, unsigned int> dense_matrix_apply_test_double_uint(1e-5);
DenseMatrixApplyTest<Mem::Main, float, unsigned long> dense_matrix_apply_test_float_ulong(1e-3);
DenseMatrixApplyTest<Mem::Main, double, unsigned long> dense_matrix_apply_test_double_ulong(1e-5);
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixApplyTest<Mem::Main, __float128, unsigned int> dense_matrix_apply_test_float128_uint(1e-6);
DenseMatrixApplyTest<Mem::Main, __float128, unsigned long> dense_matrix_apply_test_float128_ulong(1e-6);
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixApplyTest<Mem::CUDA, float, unsigned int> cuda_dense_matrix_apply_test_float_uint(1e-2);
DenseMatrixApplyTest<Mem::CUDA, double, unsigned int> cuda_dense_matrix_apply_test_double_uint(1e-4);
DenseMatrixApplyTest<Mem::CUDA, float, unsigned long> cuda_dense_matrix_apply_test_float_ulong(1e-2);
DenseMatrixApplyTest<Mem::CUDA, double, unsigned long> cuda_dense_matrix_apply_test_double_ulong(1e-4);
#endif

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class DenseMatrixMultiplyTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
   double _eps;

  explicit DenseMatrixMultiplyTest(double eps)
    : FullTaggedTest<Mem_, DT_, IT_>("DenseMatrixMultiplyTest"),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixMultiplyTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1) ; size < 100 ; size*=2)
    {
      DenseMatrix<Mem::Main, DT_, IT_> x_local(size, size+2, DT_(0));
      DenseMatrix<Mem::Main, DT_, IT_> y_local(size+2, size+1, DT_(0));
      DenseMatrix<Mem::Main, DT_, IT_> ref_local(size, size+1, DT_(4711));
      DenseMatrix<Mem::Main, DT_, IT_> result_local(size, size+1, DT_(0));
      DenseMatrix<Mem_, DT_, IT_> result(size, size+1, DT_(1234));

      for (Index i(0) ; i < x_local.size() ; ++i)
      {
        x_local.elements()[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
      }
      for (Index i(0) ; i < y_local.size() ; ++i)
      {
        y_local.elements()[i] = DT_(1) / (DT_(1) + DT_(i % 42));
      }
      DenseMatrix<Mem_, DT_, IT_> x(x_local.rows(), x_local.columns());
      DenseMatrix<Mem_, DT_, IT_> y(y_local.rows(), y_local.columns());
      x.copy(x_local);
      y.copy(y_local);

      for (Index i(0) ; i < result_local.rows() ; ++i)
      {
        for (Index k(0) ; k < result_local.columns() ; ++k)
        {
          DT_ sum(0);
          for (Index j(0) ; j < x_local.columns() ; ++j)
          {
            sum += x_local(i, j) * y_local(j, k);
          }
          ref_local(i, k, sum);
        }
      }

      result.multiply(x, y);
      result_local.copy(result);

      for (Index i(0) ; i < result.rows() ; ++i)
      {
        for (Index j(0) ; j < result.columns() ; ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i, j), ref_local(i, j), DT_(_eps));
      }
    }
  }
};

DenseMatrixMultiplyTest<Mem::Main, float, unsigned int> dense_matrix_multiply_test_float_uint(1e-3);
DenseMatrixMultiplyTest<Mem::Main, double, unsigned int> dense_matrix_multiply_test_double_uint(1e-6);
DenseMatrixMultiplyTest<Mem::Main, float, unsigned long> dense_matrix_multiply_test_float_ulong(1e-3);
DenseMatrixMultiplyTest<Mem::Main, double, unsigned long> dense_matrix_multiply_test_double_ulong(1e-6);
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixMultiplyTest<Mem::Main, __float128, unsigned int> dense_matrix_multiply_test_float128_uint(1e-6);
DenseMatrixMultiplyTest<Mem::Main, __float128, unsigned long> dense_matrix_multiply_test_float128_ulong(1e-6);
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixMultiplyTest<Mem::CUDA, float, unsigned int> cuda_dense_matrix_multiply_test_float_uint(1e-1);
DenseMatrixMultiplyTest<Mem::CUDA, double, unsigned int> cuda_dense_matrix_multiply_test_double_uint(1e-6);
DenseMatrixMultiplyTest<Mem::CUDA, float, unsigned long> cuda_dense_matrix_multiply_test_float_ulong(1e-1);
DenseMatrixMultiplyTest<Mem::CUDA, double, unsigned long> cuda_dense_matrix_multiply_test_double_ulong(1e-6);
#endif
