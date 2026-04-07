// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>
#include <kernel/util/binary_stream.hpp>

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
  explicit DenseMatrixTest(PreferredBackend backend)
    : UnitTest("DenseMatrixTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();
    DenseMatrix<DT_, IT_> zero;
    TEST_CHECK(zero.empty());

    DenseMatrix<DT_, IT_> a(10, 11);
    TEST_CHECK(!a.empty());
    TEST_CHECK_EQUAL(a.num_rows(), 10);
    TEST_CHECK_EQUAL(a.num_cols(), 11);
    TEST_CHECK_EQUAL(a.num_nzes(), 110);
    TEST_CHECK_EQUAL(a.size(), 110);
    DenseMatrix<DT_, IT_> b(10, 10, DT_(5));
    b.elements_view_rw()[7*11 + 6] = DT_(42);
    DenseMatrix<DT_, IT_> c;
    c.convert(b);
    TEST_CHECK_EQUAL(c.num_rows(), b.num_rows());
    TEST_CHECK_EQUAL(b.elements_view_r()(7*11 + 5), c.elements_view_r()(7*11 + 5));
    TEST_CHECK_EQUAL(b.elements_view_r()(7*11 + 6), c.elements_view_r()(7*11 + 6));
    TEST_CHECK_LESS_THAN(c.max_rel_diff(b), eps);

    DenseMatrix<DT_, IT_> e(11, 12, DT_(5));
    TEST_CHECK_EQUAL(e.num_rows(), 11ul);
    TEST_CHECK_EQUAL(e.num_cols(), 12ul);

    DenseMatrix<DT_, IT_> h;
    h.clone(c);
    TEST_CHECK_LESS_THAN(h.max_rel_diff(c), eps);
    h.elements_view_rw()[1*11+2] = DT_(3);
    TEST_CHECK_LESS_THAN(eps, h.max_rel_diff(c));
    TEST_CHECK_NOT_EQUAL(h.elements_arbiter(), c.elements_arbiter());
  }
};
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixInverseTest
  : public UnitTest
{
public:
  explicit DenseMatrixInverseTest(PreferredBackend backend)
    : UnitTest("DenseMatrixInverseTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    // Lehmer matrix inverse test
    for(Index n = 2; n < 5; ++n)
    {
      DenseMatrix<DT_, IT_> leh(n, n), ref(n,n);
      {
        Memory::TypedView<DT_> vl(leh.elements_view_w());
        for (Index i(0); i < n; ++i)
        {
          for (Index j(0); j < n; ++j)
          {
            vl[i*n+j] = DT_(Math::min(i, j) + 1) / DT_(Math::max(i, j) + 1);
          }
        }
      }

      ref.format();
      {
        Memory::TypedView<DT_> vr(ref.elements_view_rw());
        DT_ b = -DT_(2) / DT_(3);
        vr[0*n+0] = DT_(4) / DT_(3);
        vr[0*n+1] = b;
        for(Index i(1) ; i < n - 1 ; ++i)
        {
          vr[i*n+i-1] = b;
          vr[i*n+i  ] = DT_(4*Math::cub(i+1)) / DT_(4*Math::sqr(i+1) - 1);
          vr[i*n+i+1] = b = -DT_((i+2)*(i+1)) / DT_(2*i + 3);
        }
        vr[(n-1)*n + n-2] = b;
        vr[(n-1)*n + n-1] = DT_(n*n) / DT_(2*n - 1);
      }

      DenseMatrix<DT_, IT_> inv = leh.inverse();
      {
        Memory::TypedView<DT_> vx(inv.elements_view_r());
        Memory::TypedView<DT_> vr(ref.elements_view_r());
        for(Index i = 0; i < n*n; ++i)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(vx(i), vr(i), eps);
        }
      }
    }
  }
};
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixInverseTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixSerializeTest
  : public UnitTest
{
public:
  explicit DenseMatrixSerializeTest(PreferredBackend backend)
    : UnitTest("DenseMatrixSerializeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    const Index n = 16;
    DenseMatrix<DT_, IT_> test_mat(n, n);
    {
      Memory::TypedView<DT_> vl(test_mat.elements_view_w());
      for (Index i(0); i < n; ++i)
      {
        for (Index j(0); j < n; ++j)
        {
          vl[i*n+j] = DT_(Math::min(i, j) + 1) / DT_(Math::max(i, j) + 1);
        }
      }
    }

    auto kp = test_mat.serialize(LAFEM::SerialConfig(false, false));
    DenseMatrix<DT_, IT_> k(kp);
    TEST_CHECK_LESS_THAN(k.max_rel_diff(test_mat), eps);

#ifdef FEAT_HAVE_ZLIB
    auto zl = test_mat.serialize(LAFEM::SerialConfig(true, false));
    DenseMatrix<DT_, IT_> k1(zl);
    TEST_CHECK_LESS_THAN(k1.max_rel_diff(test_mat), eps);
#endif

#ifdef FEAT_HAVE_ZFP
    auto zfp = test_mat.serialize(LAFEM::SerialConfig(false, true, FEAT::Real(1e-6)));
    DenseMatrix<DT_, IT_> k2(zfp);
    TEST_CHECK_LESS_THAN(k2.max_rel_diff(test_mat), DT_(1e-6));
#endif

    //Binary Test
    {
      BinaryStream bs;
      test_mat.write_out(FileMode::fm_dm, bs);
      bs.seekg(0);
      DenseMatrix<DT_, IT_> test(FileMode::fm_dm, bs);
      TEST_CHECK_LESS_THAN(test.max_rel_diff(test_mat), eps);
    }

    //Mtx Test
    {
      std::stringstream ts;
      test_mat.write_out(FileMode::fm_mtx, ts);
      DenseMatrix<DT_, IT_> test(FileMode::fm_mtx, ts);
      TEST_CHECK_LESS_THAN(test.max_rel_diff(test_mat), eps);
    }

    //*
    //FileTest-> for now... problem if write rights arent given...
    {
      String filename = "test_dense_matrix_file_bin.txt";
      std::ofstream(filename.c_str());
      test_mat.write_out(FileMode::fm_dm, filename);
      DenseMatrix<DT_, IT_> test(FileMode::fm_dm, filename);
      TEST_CHECK_LESS_THAN(test.max_rel_diff(test_mat), eps);
      std::remove(filename.c_str());
    }
    {
      String filename = "test_dense_matrix_file_mtx.txt";
      std::ofstream(filename.c_str());
      test_mat.write_out(FileMode::fm_mtx, filename);
      DenseMatrix<DT_, IT_> test(FileMode::fm_mtx, filename);
      TEST_CHECK_LESS_THAN(test.max_rel_diff(test_mat), eps);
      std::remove(filename.c_str());
    }
    //*/
  }
};
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, __float128, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, __float128, std::uint64_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, __half, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, __half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSerializeTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixApplyTest
  : public UnitTest
{
public:
  explicit DenseMatrixApplyTest(PreferredBackend backend)
    : UnitTest("DenseMatrixApplyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();
    DT_ s(DT_(0.123));
    for (Index size(16); size < 100; size *= 2)
    {
      DenseMatrix<DT_, IT_> a(size, size + 1, DT_(0));
      DenseVector<DT_, IT_> x(size + 1);

      DenseVector<DT_, IT_> y(size);
      DenseVector<DT_, IT_> ref(size);
      DenseVector<DT_, IT_> result(size, DT_(4711));

      {
        Memory::TypedView<DT_> vx(x.elements_view_w());
        for (Index i(0); i < x.size(); ++i)
        {
          vx[i] = DT_(DT_(1) / (DT_(1) + DT_(i % 100) * DT_(1.234)));
        }
      }
      {
        Memory::TypedView<DT_> vy(y.elements_view_w());
        for (Index i(0); i < y.size(); ++i)
        {
          vy[i] = DT_(DT_(2) - DT_(i % 42));
        }
      }
      {
        Memory::TypedView<DT_> va(a.elements_view_w());
        for (Index i(0); i < a.size(); ++i)
        {
          va[i] = DT_(DT_(i % 100) * DT_(1.234));
        }
      }

      // apply-test for alpha = 0.0
      a.apply(result, x, y, DT_(0.0));
      TEST_CHECK_LESS_THAN(y.max_rel_diff(result), eps);

      //apply-test for reduced apply call
      a.apply(result, x);
      {
        Memory::TypedView<DT_> va(a.elements_view_r());
        Memory::TypedView<DT_> vx(x.elements_view_r());
        Memory::TypedView<DT_> vr(ref.elements_view_w());
        for (Index i(0); i < a.num_rows(); ++i)
        {
          DT_ sum(0);
          for (Index j(0); j < a.num_cols(); ++j)
          {
            sum += va(i * a.num_cols() + j) * vx(j);
          }
          vr[i] = sum;
        }
      }
      TEST_CHECK_LESS_THAN(ref.max_rel_diff(result), eps);

      //apply-test for alpha = -1.0
      a.apply(result, x, y, DT_(-1.0));
      {
        Memory::TypedView<DT_> va(a.elements_view_r());
        Memory::TypedView<DT_> vx(x.elements_view_r());
        Memory::TypedView<DT_> vy(y.elements_view_r());
        Memory::TypedView<DT_> vr(ref.elements_view_w());
        for (Index i(0); i < a.num_rows(); ++i)
        {
          DT_ sum(0);
          for (Index j(0); j < a.num_cols(); ++j)
          {
            sum += va(i * a.num_cols() + j) * vx(j);
          }
          vr[i] = vy(i) - sum;
        }
      }
      TEST_CHECK_LESS_THAN(ref.max_rel_diff(result), eps);

      // apply-test for s = 0.123
      a.apply(result, x, y, s);

      Backend::set_preferred_backend(PreferredBackend::generic);
      a.apply(ref, x);
      ref.scale(ref, s);
      ref.axpy(y);
      TEST_CHECK_LESS_THAN(ref.max_rel_diff(result), eps);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixApplyTransposedTest
  : public UnitTest
{
public:
  explicit DenseMatrixApplyTransposedTest(PreferredBackend backend)
    : UnitTest("DenseMatrixApplyTransposedTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();
    DT_ s(DT_(0.123));
    for (Index size(16); size < 100; size *= 2)
    {
      DenseMatrix<DT_, IT_> a(size + 1, size, DT_(0));
      DenseVector<DT_, IT_> x(size + 1);

      DenseVector<DT_, IT_> y(size);
      DenseVector<DT_, IT_> ref(size);
      DenseVector<DT_, IT_> result(size, DT_(4711));

      {
        Memory::TypedView<DT_> vx(x.elements_view_w());
        for (Index i(0); i < x.size(); ++i)
        {
          vx[i] = DT_(DT_(1) / (DT_(1) + DT_(i % 100) * DT_(1.234)));
        }
      }
      {
        Memory::TypedView<DT_> vy(y.elements_view_w());
        for (Index i(0); i < y.size(); ++i)
        {
          vy[i] = DT_(DT_(2) - DT_(i % 42));
        }
      }
      {
        Memory::TypedView<DT_> va(a.elements_view_w());
        for (Index i(0); i < a.size(); ++i)
        {
          va[i] = DT_(DT_(i % 100) * DT_(1.234));
        }
      }

      // apply-test for alpha = 0.0
      a.apply_transposed(result, x, y, DT_(0.0));
      TEST_CHECK_LESS_THAN(y.max_rel_diff(result), eps);

      //apply-test for reduced apply call
      a.apply_transposed(result, x);
      {
        Memory::TypedView<DT_> va(a.elements_view_r());
        Memory::TypedView<DT_> vx(x.elements_view_r());
        Memory::TypedView<DT_> vr(ref.elements_view_w());
        for (Index i(0); i < a.num_cols(); ++i)
        {
          DT_ sum(0);
          for (Index j(0); j < a.num_rows(); ++j)
          {
            sum += va(j * a.num_cols() + i) * vx(j);
          }
          vr[i] = sum;
        }
      }
      TEST_CHECK_LESS_THAN(ref.max_rel_diff(result), eps);

      //apply-test for alpha = -1.2
      a.apply_transposed(result, x, y, DT_(-1.2));
      {
        Memory::TypedView<DT_> va(a.elements_view_r());
        Memory::TypedView<DT_> vx(x.elements_view_r());
        Memory::TypedView<DT_> vy(y.elements_view_r());
        Memory::TypedView<DT_> vr(ref.elements_view_w());
        for (Index i(0); i < a.num_cols(); ++i)
        {
          DT_ sum(0);
          for (Index j(0); j < a.num_rows(); ++j)
          {
            sum += va(j * a.num_cols() + i) * vx(j);
          }
          vr[i] = vy(i) - DT_(1.2)*sum;
        }
      }

      TEST_CHECK_LESS_THAN(ref.max_rel_diff(result), eps);

      // apply-test for s = 0.123
      a.apply_transposed(result, x, y, s);

      Backend::set_preferred_backend(PreferredBackend::generic);
      a.apply_transposed(ref, x);
      ref.scale(ref, s);
      ref.axpy(y);

      TEST_CHECK_LESS_THAN(ref.max_rel_diff(result), eps);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixApplyTransposedTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixAxpyTest
  : public UnitTest
{
public:
  explicit DenseMatrixAxpyTest(PreferredBackend backend)
    : UnitTest("DenseMatrixAxpyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    for (Index size(1); size < 33; size *= 2)
    {
      DenseMatrix<DT_, IT_> x(size, size + 2, DT_(0.));
      DenseMatrix<DT_, IT_> y(size, size + 2, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 2, DT_(4711.));
      DenseMatrix<DT_, IT_> ref2(size, size + 2, DT_(4711.));
      DenseMatrix<DT_, IT_> result(size, size + 2, DT_(1234.));
      DT_ alpha(DT_(1.5));

      {
        Memory::TypedView<DT_> vx = x.elements_view_w();
        Memory::TypedView<DT_> vy = y.elements_view_w();
        for (Index i(0); i < x.size(); ++i)
        {
          vx[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
          vy[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
        }
      }

      {
        Index n = ref.num_cols();
        Memory::TypedView<DT_> vx = x.elements_view_r();
        Memory::TypedView<DT_> vy = y.elements_view_r();
        Memory::TypedView<DT_> vr1 = ref.elements_view_w();
        Memory::TypedView<DT_> vr2 = ref2.elements_view_w();
        for (Index i(0); i < ref.num_rows(); ++i)
        {
          for (Index k(0); k < ref.num_cols(); ++k)
          {
            vr1[i*n+k] = alpha * vx(i*n+k) + vy(i*n+k);
            vr2[i*n+k] = alpha * vy(i*n+k) + vy(i*n+k);
          }
        }
      }
      // r == x
      result.copy(y);
      result.axpy(result, alpha);
      TEST_CHECK_LESS_THAN(result.max_rel_diff(ref2), eps);

      // r != x
      result.copy(y);
      result.axpy(x, alpha);

      TEST_CHECK_LESS_THAN(result.max_rel_diff(ref), eps);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAxpyTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixSetProductTest
  : public UnitTest
{
public:
  explicit DenseMatrixSetProductTest(PreferredBackend backend) :
    UnitTest("DenseMatrixSetProductTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    for (Index size(1); size < 33; size *= 2)
    {
      DenseMatrix<DT_, IT_> x(size, size + 2, DT_(0.));
      DenseMatrix<DT_, IT_> y(size + 2, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 1, DT_(4711.));
      DenseMatrix<DT_, IT_> result(size, size + 1, DT_(1234.));

      {
        Memory::TypedView<DT_> vx = x.elements_view_w();
        for (Index i(0); i < x.size(); ++i)
        {
          vx[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
        }
      }
      {
        Memory::TypedView<DT_> vy = y.elements_view_w();
        for (Index i(0); i < y.size(); ++i)
        {
          vy[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
        }
      }

      {
        Memory::TypedView<DT_> vx = x.elements_view_r();
        Memory::TypedView<DT_> vy = y.elements_view_r();
        Memory::TypedView<DT_> vr = ref.elements_view_w();
        const Index ncx = x.num_cols();
        const Index ncy = y.num_cols();
        const Index ncr = ref.num_cols();
        for (Index i(0); i < ref.num_rows(); ++i)
        {
          for (Index k(0); k < ref.num_cols(); ++k)
          {
            DT_ sum(0.);
            for (Index j(0); j < x.num_cols(); ++j)
            {
              sum = sum + vx(i*ncx + j) * vy(j*ncy + k);
            }
            vr[i*ncr + k] = sum;
          }
        }
      }

      result.set_product(x, y);

      TEST_CHECK_LESS_THAN(result.max_rel_diff(ref), eps);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSetProductTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixCSRSetProductTest
  : public UnitTest
{
public:
  explicit DenseMatrixCSRSetProductTest(PreferredBackend backend) :
    UnitTest("DenseMatrixCSRSetProductTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    for (Index size(1); size < 33; size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> xfac(size, size + 2);
      DenseMatrix<DT_, IT_> x_dense(size, size + 2, DT_(0.));
      {
        Memory::TypedView<DT_> vx = x_dense.elements_view_rw();
        for (IT_ row(0); row < xfac.num_rows(); ++row)
        {
          for (IT_ col(0); col < xfac.num_cols(); ++col)
          {
            if (row == col)
            {
              xfac.add(row, col, DT_(2));
              vx[row * xfac.num_cols() + col] = DT_(2);
            }
            else if ((row == col + 1) || (row + 1 == col))
            {
              xfac.add(row, col, DT_(-1));
              vx[row * xfac.num_cols() + col] = DT_(-1);
            }
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> x(xfac.make_csr());
      DenseMatrix<DT_, IT_> y(size + 2, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 1, DT_(4711.));
      DenseMatrix<DT_, IT_> result(size, size + 1, DT_(1234.));
      {
        Memory::TypedView<DT_> vy = y.elements_view_w();
        Index mn = y.num_rows() * y.num_cols();
        for (Index i(0); i < mn; ++i)
        {
          vy[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
        }
      }

      ref.set_product(x_dense, y);
      result.set_product(x, y);
      TEST_CHECK_LESS_THAN(result.max_rel_diff(ref), eps);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRSetProductTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixAddProductTest
  : public UnitTest
{
public:
  explicit DenseMatrixAddProductTest(PreferredBackend backend)
    : UnitTest("DenseMatrixAddProductTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    for (Index size(1); size < 33; size *= 2)
    {
      DenseMatrix<DT_, IT_> x(size, size + 2, DT_(0.));
      DenseMatrix<DT_, IT_> y(size + 2, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> z(size, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 1, DT_(4711.));
      DT_ alpha = DT_(3);

      {
        Memory::TypedView<DT_> vx = x.elements_view_w();
        for (Index i(0); i < x.size(); ++i)
        {
          vx[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
        }
      }
      {
        Memory::TypedView<DT_> vy = y.elements_view_w();
        for (Index i(0); i < y.size(); ++i)
        {
          vy[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
        }
      }
      {
        Memory::TypedView<DT_> vz = z.elements_view_w();
        for (Index i(0); i < z.size(); ++i)
        {
          vz[i] = DT_(1.) / (DT_(1.) + DT_(i % 37));
        }
      }

      {
        Memory::TypedView<DT_> vx = x.elements_view_r();
        Memory::TypedView<DT_> vy = y.elements_view_r();
        Memory::TypedView<DT_> vz = z.elements_view_r();
        Memory::TypedView<DT_> vr = ref.elements_view_w();
        const Index nr = ref.num_rows();
        const Index nc = ref.num_cols();
        const Index ncx = x.num_cols();
        const Index ncy = y.num_cols();
        for (Index i(0); i < nr; ++i)
        {
          for (Index k(0); k < nc; ++k)
          {
            DT_ sum(0.);
            for (Index j(0); j < ncx; ++j)
            {
              sum = sum + vx(i * ncx + j) * vy(j * ncy + k);
            }
            vr[i * nc + k] = vz(i * nc + k) + alpha*sum;
          }
        }
      }

      z.add_product(x, y, alpha);

      TEST_CHECK_LESS_THAN(z.max_rel_diff(ref), eps);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
//SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, float, std::uint64_t, PreferredBackend::mkl);
//SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixAddProductTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
  class DenseMatrixCSRAddProductTest
  : public UnitTest
{
public:
  explicit DenseMatrixCSRAddProductTest(PreferredBackend backend)
    : UnitTest("DenseMatrixCSRAddProductTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ eps = TestSystem::tol<DT_>();

    for (Index size(1); size < 33; size *= 2)
    {
      SparseMatrixFactory<DT_, IT_> xfac(size, size + 2);
      DenseMatrix<DT_, IT_> x_dense(size, size + 2, DT_(0.));
      {
        Memory::TypedView<DT_> vx = x_dense.elements_view_rw();
        for (IT_ row(0); row < xfac.num_rows(); ++row)
        {
          for (IT_ col(0); col < xfac.num_cols(); ++col)
          {
            if (row == col)
            {
              xfac.add(row, col, DT_(2));
              vx[row * xfac.num_cols() + col] = DT_(2);
            }
            else if ((row == col + 1) || (row + 1 == col))
            {
              xfac.add(row, col, DT_(-1));
              vx[row * xfac.num_cols() + col] = DT_(-1);
            }
          }
        }
      }
      SparseMatrixCSR<DT_, IT_> x(xfac.make_csr());
      DenseMatrix<DT_, IT_> y(size + 2, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> z(size, size + 1, DT_(0.));
      DenseMatrix<DT_, IT_> ref(size, size + 1, DT_(4711.));
      DT_ alpha = DT_(3);

      {
        Memory::TypedView<DT_> vy = y.elements_view_w();
        Index mn = y.num_rows() * y.num_cols();
        for (Index i(0); i < mn; ++i)
        {
          vy[i] = DT_(1.) / (DT_(1.) + DT_(i % 42));
        }
      }
      {
        Memory::TypedView<DT_> vz = z.elements_view_w();
        Index mn = z.num_rows() * z.num_cols();
        for (Index i(0); i < mn; ++i)
        {
          vz[i] = DT_(1.) / (DT_(1.) + DT_(i % 37));
        }
      }

      {
        Memory::TypedView<DT_> vx = x_dense.elements_view_r();
        Memory::TypedView<DT_> vy = y.elements_view_r();
        Memory::TypedView<DT_> vz = z.elements_view_r();
        Memory::TypedView<DT_> vr = ref.elements_view_w();
        for (Index i(0); i < ref.num_rows(); ++i)
        {
          for (Index k(0); k < ref.num_cols(); ++k)
          {
            DT_ sum(0.);
            for (Index j(0); j < x_dense.num_cols(); ++j)
            {
              sum = sum + vx(i * x.num_cols()+ j) * vy(j * y.num_cols() + k);
            }
            vr[i * ref.num_cols() + k] = vz(i * z.num_cols() + k) + alpha*sum;
          }
        }
      }

      z.add_product(x, y, alpha);

      {
        Memory::TypedView<DT_> vz = z.elements_view_r();
        Memory::TypedView<DT_> vr = ref.elements_view_r();
        for (Index i(0); i < z.num_rows(); ++i)
        {
          for (Index j(0); j < z.num_cols(); ++j)
            TEST_CHECK_EQUAL_WITHIN_EPS(vz(i * z.num_cols() + j), vr(i * ref.num_cols() + j), eps);
        }
      }
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixCSRAddProductTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixTranposeTest
  : public UnitTest
{
public:
  explicit DenseMatrixTranposeTest(PreferredBackend backend)
    : UnitTest("DenseMatrixTranposeTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    for (Index size(1); size < 33; size *= 2)
    {
      const Index m = size;
      const Index n = size + 2;

      DenseMatrix<DT_, IT_> x(m, n, DT_(0));

      {
        Memory::TypedView<DT_> v = x.elements_view_w();
        for (Index i(0); i < m*n; ++i)
        {
          v[i] = DT_(DT_(1 + i % 100) * DT_(1.234));
        }
      }

      DenseMatrix<DT_, IT_> result = x.transpose();

      {
        Memory::TypedView<DT_> v = x.elements_view_r();
        Memory::TypedView<DT_> w = result.elements_view_r();
        for (Index i(0); i < m; ++i)
        {
          for (Index j(0); j < n; ++j)
          {
            TEST_CHECK_EQUAL(v(i*n+j), w(j*m + i));
          }
        }
      }
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, __half, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, __half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixTranposeTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixMaxRelDiffTest
  : public UnitTest
{
public:
  explicit DenseMatrixMaxRelDiffTest(PreferredBackend backend)
    : UnitTest("DenseMatrixMaxRelDiffTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    for (Index size(12); size < Index(70); size *= 2)
    {
      DenseMatrix<DT_, IT_> a(size, size+2);
      DenseMatrix<DT_, IT_> b(size, size+2);

      Index nzes = a.num_nzes();

      const Index off0 = (3*nzes) / 8;
      const Index off1 = (1*nzes) / 8;
      const Index off2 = (6*nzes) / 8;

      // a = i, b = i
      {
        Memory::TypedView<DT_> va = a.elements_view_w();
        Memory::TypedView<DT_> vb = b.elements_view_w();
        for (Index i(0); i < nzes; ++i)
        {
          va[i] = vb[i] = DT_(int(i) - int(off0)) * DT_(0.123);
        }
      }

      // identical vectors, result should be zero
      TEST_CHECK_LESS_THAN(a.max_rel_diff(b), tol*tol);
      TEST_CHECK_LESS_THAN(b.max_rel_diff(a), tol*tol);

      // two values close to zero
      const DT_ delta_a0(Math::sqrt(Math::eps<DT_>()));
      const DT_ delta_b0(Math::sqr(Math::eps<DT_>()));
      const DT_ ref0 = (delta_a0 + delta_b0) / (DT_(1) + delta_a0 + delta_b0);
      a.elements_view_rw()[off0] += delta_a0;
      b.elements_view_rw()[off0] -= delta_b0;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref0, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref0, tol);

      const DT_ delta1 = DT_(0.17);
      const DT_ ref1 = delta1 / (DT_(off0 - off1)*DT_(0.246) + delta1 + DT_(1));
      a.elements_view_rw()[off1] -= delta1;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref1, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref1, tol);

      const DT_ delta2 = DT_(0.73);
      const DT_ ref2 = delta2 / (DT_(off2 - off0)*DT_(0.246) + delta2 + DT_(1));
      b.elements_view_rw()[off2] += delta2;
      TEST_CHECK_RELATIVE(a.max_rel_diff(b), ref2, tol);
      TEST_CHECK_RELATIVE(b.max_rel_diff(a), ref2, tol);
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixMaxRelDiffTest, double, std::uint64_t, PreferredBackend::cuda);
#endif

template<
  typename DT_,
  typename IT_>
class DenseMatrixSameLayoutTest
  : public UnitTest
{
public:
  explicit DenseMatrixSameLayoutTest(PreferredBackend backend)
    : UnitTest("DenseMatrixSameLayoutTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    for (Index size(2); size < 100; size *= 2)
    {
      const Index diff_row = size / 2;
      const Index diff_col = size - 1;
      const DT_ initial_value = DT_(10);

      // only initial value in a
      DenseMatrix<DT_, IT_> a(size, size + 2, initial_value);

      // weak copy
      DenseMatrix<DT_, IT_> b = a.clone(CloneMode::Weak);
      TEST_CHECK(a.same_layout(b));

      // shallow copy
      DenseMatrix<DT_, IT_> c = a.clone(CloneMode::Shallow);
      TEST_CHECK(a.same_layout(c));

      // change one element
      c.elements_view_rw()[diff_row*size+ diff_col] = DT_(0.5);
      TEST_CHECK(a.same_layout(c));

      // different sizes
      DenseMatrix<DT_, IT_> d(size, size + 2, initial_value);
      DenseMatrix<DT_, IT_> e(size, size, initial_value);
      TEST_CHECK(!d.same_layout(e));

      // one different element
      DenseMatrix<DT_, IT_> f(size, size + 2, DT_(0));
      f.elements_view_rw()[diff_row*size+ diff_col] = initial_value;
      DenseMatrix<DT_, IT_> g(size, size + 2, DT_(0));
      g.elements_view_rw()[diff_row*size+ diff_col] = DT_(0.5);
      TEST_CHECK(f.same_layout(g));
    }
  }
};

SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, Half, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, Half, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, Half, std::uint64_t, PreferredBackend::cuda);
#endif
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(DenseMatrixSameLayoutTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
