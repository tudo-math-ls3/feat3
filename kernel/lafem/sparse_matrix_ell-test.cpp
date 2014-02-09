#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/util/binary_stream.hpp>

#include <cstdio>
#include <sstream>

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
class SparseMatrixELLTest
  : public TaggedTest<Mem_, DT_>
{
public:
  SparseMatrixELLTest()
    : TaggedTest<Mem_, DT_>("SparseMatrixELLTest")
  {
  }

  virtual void run() const
  {
    SparseMatrixCOO<Mem::Main, DT_> a(10, 12);
    a(1,2,7);
    a.clear();
    a(1,2,7);
    a(5,5,2);
    SparseMatrixELL<Mem_, DT_> b(a);
    TEST_CHECK_EQUAL(b.used_elements(), a.used_elements());
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(b(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(b(5, 5), a(5, 5));

    SparseMatrixELL<Mem_, DT_> bl(b.layout());
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    bl = b.layout();
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    //SparseMatrixCSR<Mem_, DT_> b2(b);
    SparseMatrixCSR<Mem_, DT_> b2(a);
    SparseMatrixELL<Mem_, DT_> b3(b);
    TEST_CHECK_EQUAL(b3, b);

    SparseMatrixELL<Mem_, DT_> z(b);
    TEST_CHECK_EQUAL(z.used_elements(), 2ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z.stride(), b.stride());
    TEST_CHECK_EQUAL(z.num_cols_per_row(), b.num_cols_per_row());
    TEST_CHECK_EQUAL(z(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(z(5, 5), a(5, 5));
    TEST_CHECK_EQUAL(z(1, 3), a(1, 3));

    SparseMatrixELL<Mem_, DT_> c;
    c = b;
    TEST_CHECK_EQUAL(c.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(c(0,2), b(0,2));
    TEST_CHECK_EQUAL(c(1,2), b(1,2));
    TEST_CHECK_EQUAL(c, b);

    SparseMatrixELL<Mem::Main, DT_> e(c);
    TEST_CHECK_EQUAL(e, c);
    e = c;
    TEST_CHECK_EQUAL(e, c);
    e = c.clone();
    TEST_CHECK_EQUAL(e, c);

    SparseMatrixCOO<Mem::Main, DT_> fcoo(10, 10);
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
    SparseMatrixELL<Mem_, DT_> f(fcoo);

    BinaryStream bs;
    f.write_out(fm_ell, bs);
    bs.seekg(0);
    SparseMatrixELL<Mem_, DT_> g(bs);
    TEST_CHECK_EQUAL(g, f);

    std::stringstream ts;
    f.write_out(fm_m, ts);
    SparseMatrixCOO<Mem::Main, DT_> i(fm_m, ts);
    TEST_CHECK_EQUAL(i, fcoo);

    f.write_out(fm_mtx, ts);
    SparseMatrixCOO<Mem::Main, DT_> j(fm_mtx, ts);
    TEST_CHECK_EQUAL(j, fcoo);
  }
};
SparseMatrixELLTest<Mem::Main, float> cpu_sparse_matrix_ell_test_float;
SparseMatrixELLTest<Mem::Main, double> cpu_sparse_matrix_ell_test_double;
#ifdef FEAST_GMP
SparseMatrixELLTest<Mem::Main, mpf_class> cpu_sparse_matrix_ell_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixELLTest<Mem::CUDA, float> cuda_sparse_matrix_ell_test_float;
SparseMatrixELLTest<Mem::CUDA, double> cuda_sparse_matrix_ell_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class SparseMatrixELLApplyTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  SparseMatrixELLApplyTest()
    : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixELLApplyTest")
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
      std::cout<<"real ue: "<<ue<<std::endl;
      SparseMatrixELL<Mem_,DT_> a(a_local);

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

SparseMatrixELLApplyTest<Mem::Main, Algo::Generic, float> sm_ell_apply_test_float;
SparseMatrixELLApplyTest<Mem::Main, Algo::Generic, double> sm_ell_apply_test_double;
#ifdef FEAST_GMP
SparseMatrixELLApplyTest<Mem::Main, Algo::Generic, mpf_class> sm_ell_apply_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixELLApplyTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_ell_apply_test_float;
SparseMatrixELLApplyTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_ell_apply_test_double;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class SparseMatrixELLScaleTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  SparseMatrixELLScaleTest()
    : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixELLScaleTest")
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

      SparseMatrixELL<Mem_, DT_> a(a_local);
      SparseMatrixELL<Mem_, DT_> b(a.clone());

      b.template scale<Algo_>(a, s);
      SparseMatrixCOO<Mem::Main, DT_> b_local(b);
      TEST_CHECK_EQUAL(b_local, ref_local);

      a.template scale<Algo_>(a, s);
      SparseMatrixCOO<Mem_, DT_> a_coo(a);
      a_local = a_coo;
      TEST_CHECK_EQUAL(a_local, ref_local);
    }
  }
};
SparseMatrixELLScaleTest<Mem::Main, Algo::Generic, float> sm_ell_scale_test_float;
SparseMatrixELLScaleTest<Mem::Main, Algo::Generic, double> sm_ell_scale_test_double;
#ifdef FEAST_GMP
SparseMatrixELLScaleTest<Mem::Main, Algo::Generic, mpf_class> sm_ell_scale_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_MKL
SparseMatrixELLScaleTest<Mem::Main, Algo::MKL, float> mkl_sm_ell_scale_test_float;
SparseMatrixELLScaleTest<Mem::Main, Algo::MKL, double> mkl_sm_ell_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixELLScaleTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_ell_scale_test_float;
SparseMatrixELLScaleTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_ell_scale_test_double;
#endif
