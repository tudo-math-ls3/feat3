#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/algorithm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVAxpyTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVAxpyTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_axpy_test")
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
      c.template axpy<Algo_>(s, a, b);
      copy(result_local, c);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      a.template axpy<Algo_>(s, a, b);
      copy(result_local, a);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      copy(a, a_local);
      b.template axpy<Algo_>(s, a, b);
      copy(result_local, b);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);
    }
  }
};
DVAxpyTest<Mem::Main, Algo::Generic, float> dv_axpy_test_float;
DVAxpyTest<Mem::Main, Algo::Generic, double> dv_axpy_test_double;
#ifdef FEAST_GMP
DVAxpyTest<Mem::Main, Algo::Generic, mpf_class> dv_axpy_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_MKL
DVAxpyTest<Mem::Main, Algo::MKL, float> mkl_dv_axpy_test_float;
DVAxpyTest<Mem::Main, Algo::MKL, double> mkl_dv_axpy_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVAxpyTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_axpy_test_float;
DVAxpyTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_axpy_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVAxpyVTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVAxpyVTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_axpy_v_test")
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
      d.template axpy<Algo_>(c, a, b);
      copy(result_local, d);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      a.template axpy<Algo_>(c, a, b);
      copy(result_local, a);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      copy(a, a_local);
      b.template axpy<Algo_>(c, a, b);
      copy(result_local, b);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);

      copy(b, b_local);
      c.template axpy<Algo_>(c, a, b);
      copy(result_local, c);
      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref(i), 1e-2);
    }
  }
};
DVAxpyVTest<Mem::Main, Algo::Generic, float> dv_axpy_v_test_float;
DVAxpyVTest<Mem::Main, Algo::Generic, double> dv_axpy_v_test_double;
#ifdef FEAST_GMP
DVAxpyVTest<Mem::Main, Algo::Generic, mpf_class> dv_axpy_v_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVAxpyVTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_axpy_v_test_float;
DVAxpyVTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_axpy_v_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_,
  typename SM_>
class DVAxpyMVTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:
  DVAxpyMVTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_axpy_mv_test: " + SM_::type_name())
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
      SM_ a(a_local);

      DenseVector<Arch_, DT_> r(size);
      r.template axpy<Algo_>(s, a, x, y);
      copy(result_local, r);

      ref.template product_matvec<Algo_>(a, x);
      ref.template scale<Algo_>(ref, s);
      ref.template sum<Algo_>(ref, y);
      copy(ref_local, ref);

      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), 1e-2);
    }
  }
};
DVAxpyMVTest<Mem::Main, Algo::Generic, float, SparseMatrixCSR<Mem::Main, float> > dv_axpy_mv_csr_test_float;
DVAxpyMVTest<Mem::Main, Algo::Generic, double, SparseMatrixCSR<Mem::Main, double> > dv_axpy_mv_csr_test_double;
#ifdef FEAST_GMP
DVAxpyMVTest<Mem::Main, Algo::Generic, mpf_class, SparseMatrixCSR<Mem::Main, mpf_class> > dv_axpy_mv_csr_test_mpf_class;
#endif
DVAxpyMVTest<Mem::Main, Algo::Generic, float, SparseMatrixELL<Mem::Main, float> > dv_axpy_mv_ell_test_float;
DVAxpyMVTest<Mem::Main, Algo::Generic, double, SparseMatrixELL<Mem::Main, double> > dv_axpy_mv_ell_test_double;
#ifdef FEAST_GMP
DVAxpyMVTest<Mem::Main, Algo::Generic, mpf_class, SparseMatrixELL<Mem::Main, mpf_class> > dv_axpy_mv_ell_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVAxpyMVTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixCSR<Mem::CUDA, float> > cuda_dv_axpy_mv_csr_test_float;
DVAxpyMVTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixCSR<Mem::CUDA, double> > cuda_dv_axpy_mv_csr_test_double;
DVAxpyMVTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixELL<Mem::CUDA, float> > cuda_dv_axpy_mv_ell_test_float;
DVAxpyMVTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixELL<Mem::CUDA, double> > cuda_dv_axpy_mv_ell_test_double;
#endif
