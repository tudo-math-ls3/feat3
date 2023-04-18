// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the dense matrix class.
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
  class DenseMatrixTest
  : public UnitTest
{

public:

  DenseMatrixTest(PreferredBackend backend)
    : UnitTest("DenseMatrixTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseMatrixTest()
  {
  }

  virtual void run() const override
  {
    DenseMatrix<DT_, IT_> zero1;
    DenseMatrix<DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    DenseMatrix<DT_, IT_> a(10, 11);
    TEST_CHECK_EQUAL(a.rows(), 10);
    TEST_CHECK_EQUAL(a.columns(), 11);
    TEST_CHECK_EQUAL(a.size(), 110);
    DenseMatrix<DT_, IT_> b(10, 10, 5.);
    b(7, 6, DT_(42));
    DenseMatrix<DT_, IT_> c;
    c.convert(b);
    TEST_CHECK_EQUAL(c.size(), b.size());
    TEST_CHECK_EQUAL(c.rows(), b.rows());
    TEST_CHECK_EQUAL(c(7, 6), b(7, 6));
    TEST_CHECK_EQUAL(c, b);

    DenseMatrix<DT_, IT_> e(11, 12, 5.);
    TEST_CHECK_EQUAL(e.rows(), 11ul);
    TEST_CHECK_EQUAL(e.columns(), 12ul);

    DenseMatrix<DT_, IT_> h;
    h.clone(c);
    TEST_CHECK_EQUAL(h, c);
    h(1, 2, 3);
    TEST_CHECK_NOT_EQUAL(h, c);
    TEST_CHECK_NOT_EQUAL((void*)h.elements(), (void*)c.elements());

    // Lehmer matrix inverse test
    DenseMatrix<DT_, IT_> l(11, 11);
    for (Index i(0); i < l.rows(); ++i)
    {
      for (Index j(0); j < l.rows(); ++j)
      {
        l(i, j, DT_(Math::min(i, j) + 1) / DT_(Math::max(i, j) + 1));
      }
    }
    auto m = l.inverse();
    m.invert();
    for (Index i(0); i < l.rows(); ++i)
    {
      for (Index j(0); j < l.rows(); ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(m(i, j), l(i, j), Math::template eps<DT_>() * DT_(100));
      }
    }

  }
};
DenseMatrixTest<float, unsigned int> cpu_dense_matrix_test_float_uint(PreferredBackend::generic);
DenseMatrixTest<double, unsigned int> cpu_dense_matrix_test_double_uint(PreferredBackend::generic);
DenseMatrixTest<float, unsigned long> cpu_dense_matrix_test_float_ulong(PreferredBackend::generic);
DenseMatrixTest<double, unsigned long> cpu_dense_matrix_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseMatrixTest<float, unsigned long> mkl_cpu_dense_matrix_test_float_ulong(PreferredBackend::mkl);
DenseMatrixTest<double, unsigned long> mkl_cpu_dense_matrix_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixTest<__float128, unsigned int> cpu_dense_matrix_test_float128_uint(PreferredBackend::generic);
DenseMatrixTest<__float128, unsigned long> cpu_dense_matrix_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseMatrixTest<Half, unsigned int> cpu_dense_matrix_test_half_uint(PreferredBackend::generic);
DenseMatrixTest<Half, unsigned long> cpu_dense_matrix_test_half_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
DenseMatrixTest<Half, unsigned int> cuda_dense_matrix_test_half_uint(PreferredBackend::cuda);
DenseMatrixTest<Half, unsigned long> cuda_dense_matrix_test_half_ulong(PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixTest<float, unsigned int> cuda_dense_matrix_test_float_uint(PreferredBackend::cuda);
DenseMatrixTest<double, unsigned int> cuda_dense_matrix_test_double_uint(PreferredBackend::cuda);
DenseMatrixTest<float, unsigned long> cuda_dense_matrix_test_float_ulong(PreferredBackend::cuda);
DenseMatrixTest<double, unsigned long> cuda_dense_matrix_test_double_ulong(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixSerializeTest
  : public UnitTest
{

public:

  DenseMatrixSerializeTest(PreferredBackend backend)
    : UnitTest("DenseMatrixSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~DenseMatrixSerializeTest()
  {
  }

  virtual void run() const override
  {
    DenseMatrix<DT_, IT_> test_mat(40, 40);
    for (Index i(0); i < test_mat.rows(); ++i)
    {
      for (Index j(0); j < test_mat.rows(); ++j)
      {
        test_mat(i, j, DT_(Math::min(i, j) + 1) / DT_(Math::max(i, j) + 1));
      }
    }
    auto kp = test_mat.serialize(LAFEM::SerialConfig(false, false));
    DenseMatrix<DT_, IT_> k(kp);
    TEST_CHECK_EQUAL(k, test_mat);
#ifdef FEAT_HAVE_ZLIB
    auto zl = test_mat.serialize(LAFEM::SerialConfig(true, false));
    DenseMatrix<DT_, IT_> k1(zl);
    TEST_CHECK_EQUAL(k1, test_mat);
#endif
#ifdef FEAT_HAVE_ZFP
    auto zfp = test_mat.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-6)));
    DenseMatrix<DT_, IT_> k2(zfp);
    for (Index i(0); i < k2.rows(); ++i)
    {
      for (Index j(0); j < k2.rows(); ++j)
      {
        TEST_CHECK_EQUAL_WITHIN_EPS(k2(i, j), test_mat(i, j), (DT_)1e-6);
      }
    }
#endif

    // DenseMatrix write_out test
    DenseMatrix<DT_, IT_> u(11, 11);
    for (Index i(0); i < u.rows(); ++i)
    {
      for (Index j(0); j < u.columns(); ++j)
      {
        u(i, j, DT_(i + j + 1));
      }
    }
    //Binary Test
    {
      BinaryStream bs;
      u.write_out(FileMode::fm_dm, bs);
      bs.seekg(0);
      DenseMatrix<DT_, IT_> test(FileMode::fm_dm, bs);
      for (Index i(0); i < u.rows(); ++i)
      {
        for (Index j(0); j < u.columns(); ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(u(i, j), test(i, j), Math::template eps<DT_>() * DT_(100));
        }
      }
    }
    //Mtx Test
    {
      std::stringstream ts;
      u.write_out(FileMode::fm_mtx, ts);
      DenseMatrix<DT_, IT_> test2(FileMode::fm_mtx, ts);
      for (Index i(0); i < u.rows(); ++i)
      {
        for (Index j(0); j < u.columns(); ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(u(i, j), test2(i, j), Math::template eps<DT_>() * DT_(100));
        }
      }
    }
    //*
    //FileTest-> for now... problem if write rights arent given...
    {
      String filename = "test_dense_matrix_file_bin.txt";
      std::ofstream(filename.c_str());
      u.write_out(FileMode::fm_dm, filename);
      DenseMatrix<DT_, IT_> test3(FileMode::fm_dm, filename);
      for (Index i(0); i < u.rows(); ++i)
      {
        for (Index j(0); j < u.columns(); ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(u(i, j), test3(i, j), Math::template eps<DT_>() * DT_(100));
        }
      }
      std::remove(filename.c_str());
    }
    {
      String filename = "test_dense_matrix_file_mtx.txt";
      std::ofstream(filename.c_str());
      u.write_out(FileMode::fm_mtx, filename);
      DenseMatrix<DT_, IT_> test4(FileMode::fm_mtx, filename);
      for (Index i(0); i < u.rows(); ++i)
      {
        for (Index j(0); j < u.columns(); ++j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(u(i, j), test4(i, j), Math::template eps<DT_>() * DT_(100));
        }
      }
      std::remove(filename.c_str());
    }
    //*/
  }
};
DenseMatrixSerializeTest<float, unsigned int> cpu_dense_matrix_serialize_test_float_uint(PreferredBackend::generic);
DenseMatrixSerializeTest<double, unsigned int> cpu_dense_matrix_serialize_test_double_uint(PreferredBackend::generic);
DenseMatrixSerializeTest<float, unsigned long> cpu_dense_matrix_serialize_test_float_ulong(PreferredBackend::generic);
DenseMatrixSerializeTest<double, unsigned long> cpu_dense_matrix_serialize_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
DenseMatrixSerializeTest<float, unsigned long> mkl_cpu_dense_matrix_serialize_test_float_ulong(PreferredBackend::mkl);
DenseMatrixSerializeTest<double, unsigned long> mkl_cpu_dense_matrix_serialize_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixSerializeTest<__float128, unsigned int> cpu_dense_matrix_serialize_test_float128_uint(PreferredBackend::generic);
DenseMatrixSerializeTest<__float128, unsigned long> cpu_dense_matrix_serialize_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
//DenseMatrixSerializeTest<__half, unsigned int> cpu_dense_matrix_serialize_test_half_uint(PreferredBackend::generic);
//DenseMatrixSerializeTest<__half, unsigned long> cpu_dense_matrix_serialize_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixSerializeTest<float, unsigned int> cuda_dense_matrix_serialize_test_float_uint(PreferredBackend::cuda);
DenseMatrixSerializeTest<double, unsigned int> cuda_dense_matrix_serialize_test_double_uint(PreferredBackend::cuda);
DenseMatrixSerializeTest<float, unsigned long> cuda_dense_matrix_serialize_test_float_ulong(PreferredBackend::cuda);
DenseMatrixSerializeTest<double, unsigned long> cuda_dense_matrix_serialize_test_double_ulong(PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixApplyTest
  : public UnitTest
{
public:
  double _eps;

  explicit DenseMatrixApplyTest(PreferredBackend backend, double eps)
    : UnitTest("DenseMatrixApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixApplyTest()
  {
  }

  virtual void run() const override
  {
    DT_ s(DT_(0.123));
    for (Index size(16); size < 100; size *= 2)
    {
      DenseMatrix<DT_, IT_> a(size, size + 1, DT_(0));
      DenseVector<DT_, IT_> x(size + 1);

      DenseVector<DT_, IT_> y(size);
      DenseVector<DT_, IT_> ref(size);
      DenseVector<DT_, IT_> result(size, DT_(4711));

      for (Index i(0); i < x.size(); ++i)
      {
        x(i, DT_(DT_(1) / (DT_(1) + DT_(i % 100) * DT_(1.234))));
      }
      for (Index i(0); i < y.size(); ++i)
      {
        y(i, DT_(DT_(2) - DT_(i % 42)));
      }

      for (Index i(0); i < a.size(); ++i)
      {
        a.elements()[i] = DT_(DT_(i % 100) * DT_(1.234));
      }

      // apply-test for alpha = 0.0
      a.apply(result, x, y, DT_(0.0));
      for (Index i(0); i < result.size(); ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(y(i), result(i), DT_(_eps));

      //apply-test for reduced apply call
      a.apply(result, x);
      for (Index i(0); i < a.rows(); ++i)
      {
        DT_ sum(0);
        for (Index j(0); j < a.columns(); ++j)
        {
          sum += a(i, j) * x(j);
        }
        ref(i, sum);
      }
      for (Index i(0); i < result.size(); ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result(i), ref(i), DT_(_eps));

      //apply-test for alpha = -1.0
      a.apply(result, x, y, DT_(-1.0));
      for (Index i(0); i < a.rows(); ++i)
      {
        DT_ sum(0);
        for (Index j(0); j < a.columns(); ++j)
        {
          sum += a(i, j) * x(j);
        }
        ref(i, y(i) - sum);
      }

      for (Index i(0); i < result.size(); ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result(i), ref(i), DT_(_eps));

      // apply-test for s = 0.123
      a.apply(result, x, y, s);

      Runtime::set_preferred_backend(PreferredBackend::generic);
      a.apply(ref, x);
      ref.axpy(ref, y, s);

      for (Index i(0); i < result.size(); ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result(i), ref(i), DT_(_eps));
    }
  }
};

DenseMatrixApplyTest<float, unsigned int> dense_matrix_apply_test_float_uint(PreferredBackend::generic, 1e-3);
DenseMatrixApplyTest<double, unsigned int> dense_matrix_apply_test_double_uint(PreferredBackend::generic, 1e-5);
DenseMatrixApplyTest<float, unsigned long> dense_matrix_apply_test_float_ulong(PreferredBackend::generic, 1e-3);
DenseMatrixApplyTest<double, unsigned long> dense_matrix_apply_test_double_ulong(PreferredBackend::generic, 1e-5);
#ifdef FEAT_HAVE_MKL
DenseMatrixApplyTest<float, unsigned long> mkl_dense_matrix_apply_test_float_ulong(PreferredBackend::mkl, 1e-3);
DenseMatrixApplyTest<double, unsigned long> mkl_dense_matrix_apply_test_double_ulong(PreferredBackend::mkl, 1e-5);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixApplyTest<__float128, unsigned int> dense_matrix_apply_test_float128_uint(PreferredBackend::generic, 1e-6);
DenseMatrixApplyTest<__float128, unsigned long> dense_matrix_apply_test_float128_ulong(PreferredBackend::generic, 1e-6);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseMatrixApplyTest<Half, unsigned int> dense_matrix_apply_test_half_uint(PreferredBackend::generic, 5e-2);
DenseMatrixApplyTest<Half, unsigned long> dense_matrix_apply_test_half_ulong(PreferredBackend::generic, 5e-2);
#ifdef FEAT_HAVE_CUDA
DenseMatrixApplyTest<Half, unsigned int> cuda_dense_matrix_apply_test_half_uint(PreferredBackend::cuda, 5e-1);
DenseMatrixApplyTest<Half, unsigned long> cuda_dense_matrix_apply_test_half_ulong(PreferredBackend::cuda, 5e-1);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixApplyTest<float, unsigned int> cuda_dense_matrix_apply_test_float_uint(PreferredBackend::cuda, 1e-2);
DenseMatrixApplyTest<double, unsigned int> cuda_dense_matrix_apply_test_double_uint(PreferredBackend::cuda, 1e-4);
DenseMatrixApplyTest<float, unsigned long> cuda_dense_matrix_apply_test_float_ulong(PreferredBackend::cuda, 1e-2);
DenseMatrixApplyTest<double, unsigned long> cuda_dense_matrix_apply_test_double_ulong(PreferredBackend::cuda, 1e-4);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixAxpyTest
  : public UnitTest
{
public:
  double _eps;

  explicit DenseMatrixAxpyTest(PreferredBackend backend, double eps)
    : UnitTest("DenseMatrixAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixAxpyTest()
  {
  }

  virtual void run() const override
  {
    double eps(_eps);

    for (Index size(1); size < 65; size *= 2)
    {
      DenseMatrix<DT_, IT_> x(size, size + 2, DT_(0.));
      DenseMatrix<DT_, IT_> y(size, size + 2, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 2, DT_(4711.));
      DenseMatrix<DT_, IT_> result(size, size + 2, DT_(1234.));
      DT_ alpha(DT_(1.5));

      for (Index i(0); i < x.size(); ++i)
      {
        x.elements()[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
        y.elements()[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
      }

      for (Index i(0); i < ref.rows(); ++i)
      {
        for (Index k(0); k < ref.columns(); ++k)
        {
          ref(i,k, alpha * x(i,k) + y(i,k));
        }
      }

      result.axpy(x, y, alpha);

      for (Index i(0); i < result.rows(); ++i)
      {
        for (Index j(0); j < result.columns(); ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result(i, j), ref(i, j), DT_(eps));
      }
#ifdef FEAT_HAVE_HALFMATH
      if (typeid(DT_) == typeid(Half))
        eps *= 4.;
#endif
    }
  }
};

DenseMatrixAxpyTest<float, unsigned int> dense_matrix_axpy_test_float_uint(PreferredBackend::generic, 1e-3);
DenseMatrixAxpyTest<double, unsigned int> dense_matrix_axpy_test_double_uint(PreferredBackend::generic, 1e-6);
DenseMatrixAxpyTest<float, unsigned long> dense_matrix_axpy_test_float_ulong(PreferredBackend::generic, 1e-3);
DenseMatrixAxpyTest<double, unsigned long> dense_matrix_axpy_test_double_ulong(PreferredBackend::generic, 1e-6);
#ifdef FEAT_HAVE_MKL
DenseMatrixAxpyTest<float, unsigned long> mkl_dense_matrix_axpy_test_float_ulong(PreferredBackend::mkl, 1e-3);
DenseMatrixAxpyTest<double, unsigned long> mkl_dense_matrix_axpy_test_double_ulong(PreferredBackend::mkl, 1e-6);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixAxpyTest<__float128, unsigned int> dense_matrix_axpy_test_float128_uint(PreferredBackend::generic, 1e-6);
DenseMatrixAxpyTest<__float128, unsigned long> dense_matrix_axpy_test_float128_ulong(PreferredBackend::generic, 1e-6);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseMatrixAxpyTest<Half, unsigned int> dense_matrix_axpy_test_half_uint(PreferredBackend::generic, 1e-2);
DenseMatrixAxpyTest<Half, unsigned long> dense_matrix_axpy_test_half_ulong(PreferredBackend::generic, 1e-2);
#ifdef FEAT_HAVE_CUDA
DenseMatrixAxpyTest<Half, unsigned int> cuda_dense_matrix_axpy_test_half_uint(PreferredBackend::cuda, 1e-2);
DenseMatrixAxpyTest<Half, unsigned long> cuda_dense_matrix_axpy_test_half_ulong(PreferredBackend::cuda, 1e-2);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixAxpyTest<float, unsigned int> cuda_dense_matrix_axpy_test_float_uint(PreferredBackend::cuda, 1e-3);
DenseMatrixAxpyTest<double, unsigned int> cuda_dense_matrix_axpy_test_double_uint(PreferredBackend::cuda, 1e-6);
DenseMatrixAxpyTest<float, unsigned long> cuda_dense_matrix_axpy_test_float_ulong(PreferredBackend::cuda, 1e-3);
DenseMatrixAxpyTest<double, unsigned long> cuda_dense_matrix_axpy_test_double_ulong(PreferredBackend::cuda, 1e-6);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixMultiplyTest
  : public UnitTest
{
public:
  double _eps;

  explicit DenseMatrixMultiplyTest(PreferredBackend backend, double eps)
    : UnitTest("DenseMatrixMultiplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixMultiplyTest()
  {
  }

  virtual void run() const override
  {
    double eps(_eps);

    for (Index size(1); size < 32; size *= 2)
    {
      DenseMatrix<DT_, IT_> x(size, size + 2, DT_(0.));
      DenseMatrix<DT_, IT_> y(size + 2, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 1, DT_(4711.));
      DenseMatrix<DT_, IT_> result(size, size + 1, DT_(1234.));

      for (Index i(0); i < x.size(); ++i)
      {
        x.elements()[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
      }
      for (Index i(0); i < y.size(); ++i)
      {
        y.elements()[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
      }

      for (Index i(0); i < ref.rows(); ++i)
      {
        for (Index k(0); k < ref.columns(); ++k)
        {
          DT_ sum(0.);
          for (Index j(0); j < x.columns(); ++j)
          {
            sum = sum + x(i, j) * y(j, k);
          }
          ref(i, k, sum);
        }
      }

      result.multiply(x, y);

      for (Index i(0); i < result.rows(); ++i)
      {
        for (Index j(0); j < result.columns(); ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result(i, j), ref(i, j), DT_(eps));
      }
#ifdef FEAT_HAVE_HALFMATH
      if (typeid(DT_) == typeid(Half))
        eps *= 4.;
#endif
    }
  }
};

DenseMatrixMultiplyTest<float, unsigned int> dense_matrix_multiply_test_float_uint(PreferredBackend::generic, 1e-3);
DenseMatrixMultiplyTest<double, unsigned int> dense_matrix_multiply_test_double_uint(PreferredBackend::generic, 1e-6);
DenseMatrixMultiplyTest<float, unsigned long> dense_matrix_multiply_test_float_ulong(PreferredBackend::generic, 1e-3);
DenseMatrixMultiplyTest<double, unsigned long> dense_matrix_multiply_test_double_ulong(PreferredBackend::generic, 1e-6);
#ifdef FEAT_HAVE_MKL
DenseMatrixMultiplyTest<float, unsigned long> mkl_dense_matrix_multiply_test_float_ulong(PreferredBackend::mkl, 1e-3);
DenseMatrixMultiplyTest<double, unsigned long> mkl_dense_matrix_multiply_test_double_ulong(PreferredBackend::mkl, 1e-6);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixMultiplyTest<__float128, unsigned int> dense_matrix_multiply_test_float128_uint(PreferredBackend::generic, 1e-6);
DenseMatrixMultiplyTest<__float128, unsigned long> dense_matrix_multiply_test_float128_ulong(PreferredBackend::generic, 1e-6);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseMatrixMultiplyTest<Half, unsigned int> dense_matrix_multiply_test_half_uint(PreferredBackend::generic, 1e-2);
DenseMatrixMultiplyTest<Half, unsigned long> dense_matrix_multiply_test_half_ulong(PreferredBackend::generic, 1e-2);
#ifdef FEAT_HAVE_CUDA
DenseMatrixMultiplyTest<Half, unsigned int> cuda_dense_matrix_multiply_test_half_uint(PreferredBackend::cuda, 1e-2);
DenseMatrixMultiplyTest<Half, unsigned long> cuda_dense_matrix_multiply_test_half_ulong(PreferredBackend::cuda, 1e-2);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixMultiplyTest<float, unsigned int> cuda_dense_matrix_multiply_test_float_uint(PreferredBackend::cuda, 1e-3);
DenseMatrixMultiplyTest<double, unsigned int> cuda_dense_matrix_multiply_test_double_uint(PreferredBackend::cuda, 1e-6);
DenseMatrixMultiplyTest<float, unsigned long> cuda_dense_matrix_multiply_test_float_ulong(PreferredBackend::cuda, 1e-3);
DenseMatrixMultiplyTest<double, unsigned long> cuda_dense_matrix_multiply_test_double_ulong(PreferredBackend::cuda, 1e-6);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixCSRMultiplyTest
  : public UnitTest
{
public:
  double _eps;

  explicit DenseMatrixCSRMultiplyTest(PreferredBackend backend, double eps)
    : UnitTest("DenseMatrixCSRMultiplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixCSRMultiplyTest()
  {
  }

  virtual void run() const override
  {
    double eps(_eps);

    for (Index size(1); size < 32; size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> xfac(size, size + 2);
      DenseMatrix<DT_, IT_> x_dense(size, size + 2, DT_(0.));
      for (IT_ row(0); row < xfac.rows(); ++row)
      {
        for (IT_ col(0); col < xfac.columns(); ++col)
        {
          if (row == col)
          {
            xfac.add(row, col, DT_(2));
            x_dense.elements()[row * xfac.columns() + col] = DT_(2);
          }
          else if ((row == col + 1) || (row + 1 == col))
          {
            xfac.add(row, col, DT_(-1));
            x_dense.elements()[row * xfac.columns() + col] = DT_(-1);
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> x(xfac.make_csr());
      DenseMatrix<DT_, IT_> y(size + 2, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 1, DT_(4711.));
      DenseMatrix<DT_, IT_> result(size, size + 1, DT_(1234.));
      for (Index i(0); i < y.size(); ++i)
      {
        y.elements()[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
      }

      ref.multiply(x_dense, y);
      result.multiply(x, y);
      for (Index i(0); i < result.rows(); ++i)
      {
        for (Index j(0); j < result.columns(); ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result(i, j), ref(i, j), DT_(eps));
      }
#ifdef FEAT_HAVE_HALFMATH
      if (typeid(DT_) == typeid(Half))
        eps *= 4.;
#endif
    }
  }
};

DenseMatrixCSRMultiplyTest<float, unsigned int> dense_matrix_csr_multiply_test_float_uint(PreferredBackend::generic, 1e-3);
DenseMatrixCSRMultiplyTest<double, unsigned int> dense_matrix_csr_multiply_test_double_uint(PreferredBackend::generic, 1e-6);
DenseMatrixCSRMultiplyTest<float, unsigned long> dense_matrix_csr_multiply_test_float_ulong(PreferredBackend::generic, 1e-3);
DenseMatrixCSRMultiplyTest<double, unsigned long> dense_matrix_csr_multiply_test_double_ulong(PreferredBackend::generic, 1e-6);
#ifdef FEAT_HAVE_MKL
DenseMatrixCSRMultiplyTest<float, unsigned long> mkl_dense_matrix_csr_multiply_test_float_ulong(PreferredBackend::mkl, 1e-3);
DenseMatrixCSRMultiplyTest<double, unsigned long> mkl_dense_matrix_csr_multiply_test_double_ulong(PreferredBackend::mkl, 1e-6);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixCSRMultiplyTest<__float128, unsigned int> dense_csr_matrix_multiply_test_float128_uint(PreferredBackend::generic, 1e-6);
DenseMatrixCSRMultiplyTest<__float128, unsigned long> dense_csr_matrix_multiply_test_float128_ulong(PreferredBackend::generic, 1e-6);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseMatrixCSRMultiplyTest<Half, unsigned int> dense_matrix_csr_multiply_test_half_uint(PreferredBackend::generic, 1e-2);
DenseMatrixCSRMultiplyTest<Half, unsigned long> dense_matrix_csr_multiply_test_half_ulong(PreferredBackend::generic, 1e-2);
#ifdef FEAT_HAVE_CUDA
DenseMatrixCSRMultiplyTest<Half, unsigned int> cuda_dense_matrix_csr_multiply_test_half_uint(PreferredBackend::cuda, 1e-2);
DenseMatrixCSRMultiplyTest<Half, unsigned long> cuda_dense_matrix_csr_multiply_test_half_ulong(PreferredBackend::cuda, 1e-2);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixCSRMultiplyTest<float, unsigned int> cuda_dense_matrix_csr_multiply_test_float_uint(PreferredBackend::cuda, 1e-3);
DenseMatrixCSRMultiplyTest<double, unsigned int> cuda_dense_matrix_csr_multiply_test_double_uint(PreferredBackend::cuda, 1e-6);
DenseMatrixCSRMultiplyTest<float, unsigned long> cuda_dense_matrix_csr_multiply_test_float_ulong(PreferredBackend::cuda, 1e-3);
DenseMatrixCSRMultiplyTest<double, unsigned long> cuda_dense_matrix_csr_multiply_test_double_ulong(PreferredBackend::cuda, 1e-6);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixMultiply2Test
  : public UnitTest
{
public:
  double _eps;

  explicit DenseMatrixMultiply2Test(PreferredBackend backend, double eps)
    : UnitTest("DenseMatrixMultiply2Test", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixMultiply2Test()
  {
  }

  virtual void run() const override
  {
    double eps(_eps);

    for (Index size(1); size < 65; size *= 2)
    {
      DenseMatrix<DT_, IT_> x(size, size + 2, DT_(0.));
      DenseMatrix<DT_, IT_> y(size + 2, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> z(size, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 1, DT_(4711.));
      DenseMatrix<DT_, IT_> result(size, size + 1, DT_(1234.));
      DT_ alpha = DT_(3);
      DT_ beta = DT_(1.5);

      for (Index i(0); i < x.size(); ++i)
      {
        x.elements()[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
      }
      for (Index i(0); i < y.size(); ++i)
      {
        y.elements()[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
      }
      for (Index i(0); i < z.size(); ++i)
      {
        z.elements()[i] = DT_(1.) / (DT_(1.) + DT_(i % 37));
      }

      for (Index i(0); i < ref.rows(); ++i)
      {
        for (Index k(0); k < ref.columns(); ++k)
        {
          DT_ sum(0.);
          for (Index j(0); j < x.columns(); ++j)
          {
            sum = sum + alpha * x(i, j) * y(j, k);
          }
          sum += beta * z(i,k);
          ref(i, k, sum);
        }
      }

      result.multiply(x, y, z, alpha, beta);

      for (Index i(0); i < result.rows(); ++i)
      {
        for (Index j(0); j < result.columns(); ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(result(i, j), ref(i, j), DT_(eps));
      }
#ifdef FEAT_HAVE_HALFMATH
      if (typeid(DT_) == typeid(Half))
        eps *= 4.;
#endif
    }
  }
};

DenseMatrixMultiply2Test<float, unsigned int> dense_matrix_multiply_2_test_float_uint(PreferredBackend::generic, 1e-3);
DenseMatrixMultiply2Test<double, unsigned int> dense_matrix_multiply_2_test_double_uint(PreferredBackend::generic, 1e-6);
DenseMatrixMultiply2Test<float, unsigned long> dense_matrix_multiply_2_test_float_ulong(PreferredBackend::generic, 1e-3);
DenseMatrixMultiply2Test<double, unsigned long> dense_matrix_multiply_2_test_double_ulong(PreferredBackend::generic, 1e-6);
#ifdef FEAT_HAVE_MKL
//DenseMatrixMultiply2Test<float, unsigned long> mkl_dense_matrix_multiply_2_test_float_ulong(PreferredBackend::mkl, 1e-3);
//DenseMatrixMultiply2Test<double, unsigned long> mkl_dense_matrix_multiply_2_test_double_ulong(PreferredBackend::mkl, 1e-6);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixMultiply2Test<__float128, unsigned int> dense_matrix_multiply_2_test_float128_uint(PreferredBackend::generic, 1e-6);
DenseMatrixMultiply2Test<__float128, unsigned long> dense_matrix_multiply_2_test_float128_ulong(PreferredBackend::generic, 1e-6);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseMatrixMultiply2Test<Half, unsigned int> dense_matrix_multiply_2_test_half_uint(PreferredBackend::generic, 1e-2);
DenseMatrixMultiply2Test<Half, unsigned long> dense_matrix_multiply_2_test_half_ulong(PreferredBackend::generic, 1e-2);
#ifdef FEAT_HAVE_CUDA
DenseMatrixMultiply2Test<Half, unsigned int> cuda_dense_matrix_multiply_2_test_half_uint(PreferredBackend::cuda, 1e-2);
DenseMatrixMultiply2Test<Half, unsigned long> cuda_dense_matrix_multiply_2_test_half_ulong(PreferredBackend::cuda, 1e-2);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixMultiply2Test<float, unsigned int> cuda_dense_matrix_multiply_2_test_float_uint(PreferredBackend::cuda, 1e-3);
DenseMatrixMultiply2Test<double, unsigned int> cuda_dense_matrix_multiply_2_test_double_uint(PreferredBackend::cuda, 1e-6);
DenseMatrixMultiply2Test<float, unsigned long> cuda_dense_matrix_multiply_2_test_float_ulong(PreferredBackend::cuda, 1e-3);
DenseMatrixMultiply2Test<double, unsigned long> cuda_dense_matrix_multiply_2_test_double_ulong(PreferredBackend::cuda, 1e-6);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixCSRMultiply2Test
  : public UnitTest
{
public:
  double _eps;

  explicit DenseMatrixCSRMultiply2Test(PreferredBackend backend, double eps)
    : UnitTest("DenseMatrixCSRMultiply2Test", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixCSRMultiply2Test()
  {
  }

  virtual void run() const override
  {
    double eps(_eps);

    for (Index size(1); size < 65; size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> xfac(size, size + 2);
      DenseMatrix<DT_, IT_> x_dense(size, size + 2, DT_(0.));
      for (IT_ row(0); row < xfac.rows(); ++row)
      {
        for (IT_ col(0); col < xfac.columns(); ++col)
        {
          if (row == col)
          {
            xfac.add(row, col, DT_(2));
            x_dense.elements()[row * xfac.columns() + col] = DT_(2);
          }
          else if ((row == col + 1) || (row + 1 == col))
          {
            xfac.add(row, col, DT_(-1));
            x_dense.elements()[row * xfac.columns() + col] = DT_(-1);
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> x(xfac.make_csr());
      DenseMatrix<DT_, IT_> y(size + 2, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> z(size, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 1, DT_(4711.));
      DT_ alpha = DT_(3);
      DT_ beta = DT_(1.5);

      for (Index i(0); i < y.size(); ++i)
      {
        y.elements()[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
      }
      for (Index i(0); i < z.size(); ++i)
      {
        z.elements()[i] = DT_(1.) / (DT_(1.) + DT_(i % 37));
      }

      for (Index i(0); i < ref.rows(); ++i)
      {
        for (Index k(0); k < ref.columns(); ++k)
        {
          DT_ sum(0.);
          for (Index j(0); j < x_dense.columns(); ++j)
          {
            sum = sum + alpha * x_dense(i, j) * y(j, k);
          }
          sum += beta * z(i,k);
          ref(i, k, sum);
        }
      }

      z.multiply(x, y, alpha, beta);

      for (Index i(0); i < z.rows(); ++i)
      {
        for (Index j(0); j < z.columns(); ++j)
          TEST_CHECK_EQUAL_WITHIN_EPS(z(i, j), ref(i, j), DT_(eps));
      }
#ifdef FEAT_HAVE_HALFMATH
      if (typeid(DT_) == typeid(Half))
        eps *= 4.;
#endif
    }
  }
};

DenseMatrixCSRMultiply2Test<float, unsigned int> dense_matrix_csr_multiply_2_test_float_uint(PreferredBackend::generic, 1e-3);
DenseMatrixCSRMultiply2Test<double, unsigned int> dense_matrix_csr_multiply_2_test_double_uint(PreferredBackend::generic, 1e-6);
DenseMatrixCSRMultiply2Test<float, unsigned long> dense_matrix_csr_multiply_2_test_float_ulong(PreferredBackend::generic, 1e-3);
DenseMatrixCSRMultiply2Test<double, unsigned long> dense_matrix_csr_multiply_2_test_double_ulong(PreferredBackend::generic, 1e-6);
#ifdef FEAT_HAVE_MKL
DenseMatrixCSRMultiply2Test<float, unsigned long> mkl_dense_matrix_csr_multiply_2_test_float_ulong(PreferredBackend::mkl, 1e-3);
DenseMatrixCSRMultiply2Test<double, unsigned long> mkl_dense_matrix_csr_multiply_2_test_double_ulong(PreferredBackend::mkl, 1e-6);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixCSRMultiply2Test<__float128, unsigned int> dense_matrix_csr_multiply_2_test_float128_uint(PreferredBackend::generic, 1e-6);
DenseMatrixCSRMultiply2Test<__float128, unsigned long> dense_matrix_csr_multiply_2_test_float128_ulong(PreferredBackend::generic, 1e-6);
#endif
#ifdef FEAT_HAVE_HALFMATH
DenseMatrixCSRMultiply2Test<Half, unsigned int> dense_matrix_csr_multiply_2_test_half_uint(PreferredBackend::generic, 1e-2);
DenseMatrixCSRMultiply2Test<Half, unsigned long> dense_matrix_csr_multiply_2_test_half_ulong(PreferredBackend::generic, 1e-2);
#ifdef FEAT_HAVE_CUDA
DenseMatrixCSRMultiply2Test<Half, unsigned int> cuda_dense_matrix_csr_multiply_2_test_half_uint(PreferredBackend::cuda, 1e-2);
DenseMatrixCSRMultiply2Test<Half, unsigned long> cuda_dense_matrix_csr_multiply_2_test_half_ulong(PreferredBackend::cuda, 1e-2);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixCSRMultiply2Test<float, unsigned int> cuda_dense_matrix_csr_multiply_2_test_float_uint(PreferredBackend::cuda, 1e-3);
DenseMatrixCSRMultiply2Test<double, unsigned int> cuda_dense_matrix_csr_multiply_2_test_double_uint(PreferredBackend::cuda, 1e-6);
DenseMatrixCSRMultiply2Test<float, unsigned long> cuda_dense_matrix_csr_multiply_2_test_float_ulong(PreferredBackend::cuda, 1e-3);
DenseMatrixCSRMultiply2Test<double, unsigned long> cuda_dense_matrix_csr_multiply_2_test_double_ulong(PreferredBackend::cuda, 1e-6);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixTranposeTest
  : public UnitTest
{
public:
  double _eps;

  explicit DenseMatrixTranposeTest(PreferredBackend backend, double eps)
    : UnitTest("DenseMatrixTranposeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend),
    _eps(eps)
  {
  }

  virtual ~DenseMatrixTranposeTest()
  {
  }

  virtual void run() const override
  {
    for (Index size(1); size < 100; size *= 2)
    {
      DenseMatrix<DT_, IT_> x(size, size + 2, DT_(0));
      DenseMatrix<DT_, IT_> result;

      for (Index i(0); i < x.size(); ++i)
      {
        x.elements()[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
      }

      result.transpose(x);

      for (Index i(0); i < result.rows(); ++i)
      {
        for (Index j(0); j < result.columns(); ++j)
        {
          TEST_CHECK_EQUAL(result(i, j), x(j, i));
        }
      }

      if (Runtime::get_preferred_backend() != PreferredBackend::cuda)
      {
        result.transpose_inplace();

        for (Index i(0); i < result.rows(); ++i)
        {
          for (Index j(0); j < result.columns(); ++j)
          {
            TEST_CHECK_EQUAL(result(i, j), x(i, j));
          }
        }
      }
    }
  }
};

DenseMatrixTranposeTest<float, unsigned int> dense_matrix_transpose_test_float_uint(PreferredBackend::generic, 1e-3);
DenseMatrixTranposeTest<double, unsigned int> dense_matrix_transpose_test_double_uint(PreferredBackend::generic, 1e-6);
DenseMatrixTranposeTest<float, unsigned long> dense_matrix_transpose_test_float_ulong(PreferredBackend::generic, 1e-3);
DenseMatrixTranposeTest<double, unsigned long> dense_matrix_transpose_test_double_ulong(PreferredBackend::generic, 1e-6);
#ifdef FEAT_HAVE_MKL
DenseMatrixTranposeTest<float, unsigned long> mkl_dense_matrix_transpose_test_float_ulong(PreferredBackend::mkl, 1e-3);
DenseMatrixTranposeTest<double, unsigned long> mkl_dense_matrix_transpose_test_double_ulong(PreferredBackend::mkl, 1e-6);
#endif
#ifdef FEAT_HAVE_QUADMATH
DenseMatrixTranposeTest<__float128, unsigned int> dense_matrix_transpose_test_float128_uint(PreferredBackend::generic, 1e-6);
DenseMatrixTranposeTest<__float128, unsigned long> dense_matrix_transpose_test_float128_ulong(PreferredBackend::generic, 1e-6);
#endif
#ifdef FEAT_HAVE_HALFMATH
//DenseMatrixTranposeTest<__half, unsigned int> dense_matrix_transpose_test_half_uint(PreferredBackend::generic, 1e-6);
//DenseMatrixTranposeTest<__half, unsigned long> dense_matrix_transpose_test_half_ulong(PreferredBackend::generic, 1e-6);
#endif
#ifdef FEAT_HAVE_CUDA
DenseMatrixTranposeTest<float, unsigned int> cuda_dense_matrix_transpose_test_float_uint(PreferredBackend::cuda, 1e-1);
DenseMatrixTranposeTest<double, unsigned int> cuda_dense_matrix_transpose_test_double_uint(PreferredBackend::cuda, 1e-6);
DenseMatrixTranposeTest<float, unsigned long> cuda_dense_matrix_transpose_test_float_ulong(PreferredBackend::cuda, 1e-1);
DenseMatrixTranposeTest<double, unsigned long> cuda_dense_matrix_transpose_test_double_ulong(PreferredBackend::cuda, 1e-6);
#endif
