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
    SparseMatrixCOO<Mem_, DT_> x;
    SparseMatrixCOO<Mem_, DT_> a(10, 10);
    a(5,5,5);
    a(1,2,7);
    a(5,5,2);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);

    a.clear();
    a(1,2,7);
    a(5,5,8);
    a(5,5,2);
    TEST_CHECK_EQUAL(a.used_elements(), 2ul);
    TEST_CHECK_EQUAL(a(1, 2), 7.);
    TEST_CHECK_EQUAL(a(5, 5), 2.);

    a.clear();
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

    SparseMatrixCOO<Mem_, DT_> c;
    c.clone(b);
    TEST_CHECK_EQUAL(c, b);
    c(1,2,3);
    TEST_CHECK_NOT_EQUAL(c, b);

    c.copy(b);
    TEST_CHECK_EQUAL(c, b);
    c(1,2,3);
    TEST_CHECK_NOT_EQUAL(c, b);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());

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
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCOOTest<Mem::CUDA, float> cuda_sparse_matrix_coo_test_float;
SparseMatrixCOOTest<Mem::CUDA, double> cuda_sparse_matrix_coo_test_double;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class SparseMatrixCOOApplyTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  SparseMatrixCOOApplyTest()
    : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixCOOApplyTest")
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
      DenseVector<Mem_, DT_> ref(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        x_local(i, DT_(i % 100 * DT_(1.234)));
        y_local(i, DT_(2 - DT_(i % 42)));
      }
      DenseVector<Mem_, DT_> x(size);
      x.copy(x_local);
      DenseVector<Mem_, DT_> y(size);
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
      SparseMatrixCOO<Mem_,DT_> a;
      a.convert(a_local);

      DenseVector<Mem_, DT_> r(size);
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

SparseMatrixCOOApplyTest<Mem::Main, Algo::Generic, float> sm_coo_apply_test_float;
SparseMatrixCOOApplyTest<Mem::Main, Algo::Generic, double> sm_coo_apply_test_double;
#ifdef FEAST_BACKENDS_MKL
SparseMatrixCOOApplyTest<Mem::Main, Algo::MKL, float> mkl_sm_coo_apply_test_float;
SparseMatrixCOOApplyTest<Mem::Main, Algo::MKL, double> mkl_sm_coo_apply_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class SparseMatrixCOOScaleTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  SparseMatrixCOOScaleTest()
    : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixCOOScaleTest")
  {
  }

  virtual void run() const
  {
    for (Index size(2) ; size < 3e2 ; size*=2)
    {
      DT_ s(DT_(4.321));

      SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
      SparseMatrixCOO<Mem::Main, DT_> ref_local(size, size + 2);
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

      SparseMatrixCOO<Mem_, DT_> a;
      a.convert(a_local);
      SparseMatrixCOO<Mem_, DT_> b;
      b.clone(a);

      b.template scale<Algo_>(a, s);
      TEST_CHECK_EQUAL(b, ref_local);

      a.template scale<Algo_>(a, s);
      TEST_CHECK_EQUAL(a, ref_local);
    }
  }
};
SparseMatrixCOOScaleTest<Mem::Main, Algo::Generic, float> sm_coo_scale_test_float;
SparseMatrixCOOScaleTest<Mem::Main, Algo::Generic, double> sm_coo_scale_test_double;
#ifdef FEAST_BACKENDS_MKL
SparseMatrixCOOScaleTest<Mem::Main, Algo::MKL, float> mkl_sm_coo_scale_test_float;
SparseMatrixCOOScaleTest<Mem::Main, Algo::MKL, double> mkl_sm_coo_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCOOScaleTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_coo_scale_test_float;
SparseMatrixCOOScaleTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_coo_scale_test_double;
#endif
