#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/axpy.hpp>
#include <kernel/lafem/algorithm.hpp>
#include <kernel/lafem/product_matvec.hpp>
#include <kernel/lafem/sum.hpp>
#include <kernel/lafem/scale.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVAxpyTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVAxpyTest()
    : TaggedTest<Arch_, DT_>("dv_axpy_test")
  {
  }

  virtual void run() const
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i % 100 * DT_(1.234)));
        b_local(i, DT_(2 - i % 42));
        ref(i, s * a_local(i) + b_local(i));
      }
      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);
      copy(b, b_local);

      DenseVector<Arch_, DT_> c(size);
      Axpy<Algo_>::value(c, s, a, b);
      copy(result_local, c);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      Axpy<Algo_>::value(a, s, a, b);
      copy(result_local, a);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      copy(a, a_local);
      Axpy<Algo_>::value(b, s, a, b);
      copy(result_local, b);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);
    }
  }
};
DVAxpyTest<Mem::Main, Algo::Generic, float> dv_axpy_test_float;
DVAxpyTest<Mem::Main, Algo::Generic, double> dv_axpy_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVAxpyTest<Mem::CUDA, Algo::CUDA, float> gpu_dv_axpy_test_float;
DVAxpyTest<Mem::CUDA, Algo::CUDA, double> gpu_dv_axpy_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVAxpyVTest
  : public TaggedTest<Arch_, DT_>
{

public:

  DVAxpyVTest()
    : TaggedTest<Arch_, DT_>("dv_axpy_v_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      DenseVector<Mem::Main, DT_> c_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i % 100 * DT_(1.234)));
        b_local(i, DT_(2 - i % 42));
        c_local(i, DT_(1 - i % 23));
        ref(i, c_local(i) * a_local(i) + b_local(i));
      }
      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);
      copy(b, b_local);
      DenseVector<Arch_, DT_> c(size);
      copy(c, c_local);

      DenseVector<Arch_, DT_> d(size);
      Axpy<Algo_>::value(d, c, a, b);
      copy(result_local, d);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      Axpy<Algo_>::value(a, c, a, b);
      copy(result_local, a);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      copy(a, a_local);
      Axpy<Algo_>::value(b, c, a, b);
      copy(result_local, b);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      copy(b, b_local);
      Axpy<Algo_>::value(c, c, a, b);
      copy(result_local, c);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);
    }
  }
};
DVAxpyVTest<Mem::Main, Algo::Generic, float> dv_axpy_v_test_float;
DVAxpyVTest<Mem::Main, Algo::Generic, double> dv_axpy_v_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVAxpyVTest<Mem::CUDA, Algo::CUDA, float> gpu_dv_axpy_v_test_float;
DVAxpyVTest<Mem::CUDA, Algo::CUDA, double> gpu_dv_axpy_v_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVAxpyMVTest
  : public TaggedTest<Arch_, DT_>
{

public:
  DVAxpyMVTest()
    : TaggedTest<Arch_, DT_>("dv_axpy_mv_test")
  {
  }

  virtual void run() const
  {
    DT_ s(DT_(4711.1));
    for (Index size(1) ; size < 1e3 ; size*=2)
    {
      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size);
      DenseVector<Mem::Main, DT_> x_local(size);
      DenseVector<Mem::Main, DT_> y_local(size);
      DenseVector<Mem::Main, DT_> ref_local(size);
      DenseVector<Arch_, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        x_local(i, DT_(i % 100 * DT_(1.234)));
        y_local(i, DT_(2 - i % 42));
      }
      DenseVector<Arch_, DT_> x(size);
      copy(x, x_local);
      DenseVector<Arch_, DT_> y(size);
      copy(y, y_local);

      for (unsigned long row(0) ; row < a_local.rows() ; ++row)
      {
        for (unsigned long col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
            a_local(row, col, DT_(2));
          else if((row == col+1) || (row+1 == col))
            a_local(row, col, DT_(-1));
        }
      }
      SparseMatrixCSR<Arch_, DT_> a(a_local);

      DenseVector<Arch_, DT_> r(size);
      Axpy<Algo_>::value(r, s, a, x, y);
      copy(result_local, r);

      ProductMatVec<Algo_>::value(ref, a, x);
      Scale<Algo_>::value(ref, ref, s);
      Sum<Algo_>::value(ref, ref, y);
      copy(ref_local, ref);

      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), 1e-2);
    }
  }
};
DVAxpyMVTest<Mem::Main, Algo::Generic, float> dv_axpy_mv_test_float;
DVAxpyMVTest<Mem::Main, Algo::Generic, double> dv_axpy_mv_test_double;
#ifdef FEAST_BACKENDS_CUDA
DVAxpyMVTest<Mem::CUDA, Algo::CUDA, float> gpu_dv_axpy_mv_test_float;
DVAxpyMVTest<Mem::CUDA, Algo::CUDA, double> gpu_dv_axpy_mv_test_double;
#endif
