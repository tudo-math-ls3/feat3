#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/util/binary_stream.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
* \brief Test class for the sparse matrix coo class.
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
  typename DT_>
class SparseMatrixCOOTest
  : public TaggedTest<Mem_, DT_>
{
public:
  SparseMatrixCOOTest()
    : TaggedTest<Mem_, DT_>("SparseMatrixCOOTest")
  {
  }

  virtual void run() const
  {
    SparseMatrixCOO<Mem_, DT_> zero1;
    SparseMatrixCOO<Mem::Main, DT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    SparseMatrixCOO<Mem_, DT_> x;
    SparseMatrixCOO<Mem_, DT_> a(10, 10);
    a(5,5,5);
    a(1,2,7);
    a(5,5,2);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);

    a.format();
    a(1,2,7);
    a(5,5,8);
    a(5,5,2);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);

    a.format();
    a(1,2,8);
    a(5,5,2);
    a(1,2,7);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);

    SparseMatrixCOO<Mem_, DT_> b;
    b.convert(a);
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(a(1,2), b(1,2));
    TEST_CHECK_EQUAL(a(0,2), b(0,2));
    TEST_CHECK_EQUAL(a, b);

    SparseMatrixCOO<Mem_, DT_> c(b.clone());
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_NOT_EQUAL((void*)c.row_indices(), (void*)b.row_indices());
    c(1,2,3);
    TEST_CHECK_NOT_EQUAL(c, b);
    c.clone(b);
    TEST_CHECK_EQUAL(c, b);
    c(1,2,3);
    TEST_CHECK_NOT_EQUAL(c, b);

    c.copy(b);
    TEST_CHECK_EQUAL(c, b);
    TEST_CHECK_NOT_EQUAL((void*)c.row_indices(), (void*)b.row_indices());
    c(1,2,3);
    TEST_CHECK_NOT_EQUAL(c, b);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());

    decltype(c) d(c.layout());
    TEST_CHECK_NOT_EQUAL((void*)d.row_indices(), (void*)c.row_indices());

    SparseMatrixCOO<Mem_, DT_> f(10, 10);
    for (Index row(0) ; row < f.rows() ; ++row)
    {
      for (Index col(0) ; col < f.columns() ; ++col)
      {
        if(row == col)
          f(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          f(row, col, DT_(-1));
      }
    }

    BinaryStream bs;
    f.write_out(FileMode::fm_coo, bs);
    bs.seekg(0);
    SparseMatrixCOO<Mem_, DT_> g(FileMode::fm_coo, bs);
    TEST_CHECK_EQUAL(g, f);

    std::stringstream ts;
    f.write_out(FileMode::fm_m, ts);
    SparseMatrixCOO<Mem_, DT_> i(FileMode::fm_m, ts);
    TEST_CHECK_EQUAL(i, f);

    std::stringstream ms;
    f.write_out(FileMode::fm_mtx, ms);
    SparseMatrixCOO<Mem_, DT_> j(FileMode::fm_mtx, ms);
    TEST_CHECK_EQUAL(j, f);
  }
};
SparseMatrixCOOTest<Mem::Main, float> sparse_matrix_coo_test_float;
SparseMatrixCOOTest<Mem::Main, double> sparse_matrix_coo_test_double;
#ifdef FEAST_HAVE_QUADMATH
SparseMatrixCOOTest<Mem::Main, float> sparse_matrix_coo_test_float128;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCOOTest<Mem::CUDA, float> cuda_sparse_matrix_coo_test_float;
SparseMatrixCOOTest<Mem::CUDA, double> cuda_sparse_matrix_coo_test_double;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class SparseMatrixCOOApplyTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  SparseMatrixCOOApplyTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("SparseMatrixCOOApplyTest")
  {
  }

  virtual void run() const
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

      Index ue(0);
      for (Index row(0) ; row < a_local.rows() ; ++row)
      {
        for (Index col(0) ; col < a_local.columns() ; ++col)
        {
          if(row == col)
          {
            a_local(row, col, DT_(2));
            ++ue;
          }
          else if((row == col+1) || (row+1 == col))
          {
            a_local(row, col, DT_(-1));
            ++ue;
          }
        }
      }
      SparseMatrixCOO<Mem_,DT_, IT_> a;
      a.convert(a_local);

      DenseVector<Mem_, DT_, IT_> r(size);
      //r.template axpy<Algo_>(s, a, x, y);
      a.template apply<Algo_>(r, x, y, s);
      result_local.copy(r);

      //ref.template product_matvec<Algo_>(a, x);
      a.template apply<Algo_>(ref, x);
      ref.template scale<Algo_>(ref, s);
      ref.template axpy<Algo_>(ref, y);
      ref_local.copy(ref);

      for (Index i(0) ; i < size ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(result_local(i), ref_local(i), 1e-2);
    }
  }
};

SparseMatrixCOOApplyTest<Mem::Main, Algo::Generic, float, unsigned int> sm_coo_apply_test_float_uint;
SparseMatrixCOOApplyTest<Mem::Main, Algo::Generic, double, unsigned int> sm_coo_apply_test_double_uint;
SparseMatrixCOOApplyTest<Mem::Main, Algo::Generic, float, unsigned long> sm_coo_apply_test_float_ulong;
SparseMatrixCOOApplyTest<Mem::Main, Algo::Generic, double, unsigned long> sm_coo_apply_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
SparseMatrixCOOApplyTest<Mem::Main, Algo::Generic, __float128, unsigned int> sm_coo_apply_test_float128_uint;
SparseMatrixCOOApplyTest<Mem::Main, Algo::Generic, __float128, unsigned long> sm_coo_apply_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_MKL
SparseMatrixCOOApplyTest<Mem::Main, Algo::MKL, float, unsigned long> mkl_sm_coo_apply_test_float_ulong;
SparseMatrixCOOApplyTest<Mem::Main, Algo::MKL, double, unsigned long> mkl_sm_coo_apply_test_double_ulong;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class SparseMatrixCOOScaleTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  SparseMatrixCOOScaleTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("SparseMatrixCOOScaleTest")
  {
  }

  virtual void run() const
  {
    for (Index size(2) ; size < 3e2 ; size*=2)
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

      SparseMatrixCOO<Mem_, DT_, IT_> a;
      a.convert(a_local);
      SparseMatrixCOO<Mem_, DT_, IT_> b;
      b.clone(a);

      b.template scale<Algo_>(a, s);
      TEST_CHECK_EQUAL(b, ref_local);

      a.template scale<Algo_>(a, s);
      TEST_CHECK_EQUAL(a, ref_local);
    }
  }
};
SparseMatrixCOOScaleTest<Mem::Main, Algo::Generic, float, unsigned int> sm_coo_scale_test_float_uint;
SparseMatrixCOOScaleTest<Mem::Main, Algo::Generic, double, unsigned int> sm_coo_scale_test_double_uint;
SparseMatrixCOOScaleTest<Mem::Main, Algo::Generic, float, unsigned long> sm_coo_scale_test_float_ulong;
SparseMatrixCOOScaleTest<Mem::Main, Algo::Generic, double, unsigned long> sm_coo_scale_test_double_ulong;
#ifdef FEAST_HAVE_QUADMATH
SparseMatrixCOOScaleTest<Mem::Main, Algo::Generic, __float128, unsigned int> sm_coo_scale_test_float128_uint;
SparseMatrixCOOScaleTest<Mem::Main, Algo::Generic, __float128, unsigned long> sm_coo_scale_test_float128_ulong;
#endif
#ifdef FEAST_BACKENDS_MKL
SparseMatrixCOOScaleTest<Mem::Main, Algo::MKL, float, unsigned int> mkl_sm_coo_scale_test_float_uint;
SparseMatrixCOOScaleTest<Mem::Main, Algo::MKL, double, unsigned int> mkl_sm_coo_scale_test_double_uint;
SparseMatrixCOOScaleTest<Mem::Main, Algo::MKL, float, unsigned long> mkl_sm_coo_scale_test_float_ulong;
SparseMatrixCOOScaleTest<Mem::Main, Algo::MKL, double, unsigned long> mkl_sm_coo_scale_test_double_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCOOScaleTest<Mem::CUDA, Algo::CUDA, float, unsigned int> cuda_sm_coo_scale_test_float_uint;
SparseMatrixCOOScaleTest<Mem::CUDA, Algo::CUDA, double, unsigned int> cuda_sm_coo_scale_test_double_uint;
SparseMatrixCOOScaleTest<Mem::CUDA, Algo::CUDA, float, unsigned long> cuda_sm_coo_scale_test_float_ulong;
SparseMatrixCOOScaleTest<Mem::CUDA, Algo::CUDA, double, unsigned long> cuda_sm_coo_scale_test_double_ulong;
#endif
