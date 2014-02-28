#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/util/binary_stream.hpp>

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
class SparseMatrixCSRTest
  : public TaggedTest<Mem_, DT_>
{
public:
  SparseMatrixCSRTest()
    : TaggedTest<Mem_, DT_>("SparseMatrixCSRTest")
  {
  }

  virtual void run() const
  {
    SparseMatrixCOO<Mem::Main, DT_> a(10, 10);
    a(1,2,7);
    a.clear();
    a(1,2,7);
    a(5,5,2);
    SparseMatrixCSR<Mem_, DT_> b(a);
    TEST_CHECK_EQUAL(b.used_elements(), a.used_elements());
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(b(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(b(5, 5), a(5, 5));

    SparseMatrixCSR<Mem_, DT_> bl(b.layout());
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    bl = b.layout();
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    SparseMatrixCSR<Mem_, DT_> z(b);
    TEST_CHECK_EQUAL(z.used_elements(), 2ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(z(5, 5), a(5, 5));

    SparseMatrixCSR<Mem_, DT_> c;
    c = b;
    TEST_CHECK_EQUAL(c.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(c(0,2), b(0,2));
    TEST_CHECK_EQUAL(c(1,2), b(1,2));
    TEST_CHECK_EQUAL(c, b);

    DenseVector<Mem_, Index> col_ind(c.used_elements(), c.col_ind());
    DenseVector<Mem_, DT_> val(c.used_elements(), c.val());
    DenseVector<Mem_, Index> row_ptr(c.rows() + 1, c.row_ptr());
    SparseMatrixCSR<Mem_, DT_> d(c.rows(), c.columns(), col_ind, val, row_ptr);
    TEST_CHECK_EQUAL(d, c);

    SparseMatrixCSR<Mem::Main, DT_> e(c);
    TEST_CHECK_EQUAL(e, c);
    e = c;
    TEST_CHECK_EQUAL(e, c);
    e = c.clone();
    TEST_CHECK_EQUAL(e, c);
    e.copy(c);
    TEST_CHECK_EQUAL(e, c);

    TEST_CHECK_NOT_EQUAL((void*)e.val(), (void*)c.val());

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
    SparseMatrixCSR<Mem_, DT_> f(fcoo);

    BinaryStream bs;
    f.write_out(FileMode::fm_csr, bs);
    bs.seekg(0);
    SparseMatrixCSR<Mem_, DT_> g(FileMode::fm_csr, bs);
    TEST_CHECK_EQUAL(g, f);

    std::stringstream ts;
    f.write_out(FileMode::fm_m, ts);
    SparseMatrixCSR<Mem::Main, DT_> i(FileMode::fm_m, ts);
    TEST_CHECK_EQUAL(i, f);

    f.write_out(FileMode::fm_mtx, ts);
    SparseMatrixCSR<Mem::Main, DT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_EQUAL(j, f);
  }
};
SparseMatrixCSRTest<Mem::Main, float> cpu_sparse_matrix_csr_test_float;
SparseMatrixCSRTest<Mem::Main, double> cpu_sparse_matrix_csr_test_double;
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCSRTest<Mem::CUDA, float> cuda_sparse_matrix_csr_test_float;
SparseMatrixCSRTest<Mem::CUDA, double> cuda_sparse_matrix_csr_test_double;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class SparseMatrixCSRApplyTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  SparseMatrixCSRApplyTest()
    : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixCSRApplyTest")
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
      SparseMatrixCSR<Mem_,DT_> a(a_local);

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

SparseMatrixCSRApplyTest<Mem::Main, Algo::Generic, float> sm_csr_apply_test_float;
SparseMatrixCSRApplyTest<Mem::Main, Algo::Generic, double> sm_csr_apply_test_double;
#ifdef HONEI_BACKENDS_MKL
SparseMatrixCSRApplyTest<Mem::Main, Algo::MKL, float> mkl_sm_csr_apply_test_float;
SparseMatrixCSRApplyTest<Mem::Main, Algo::MKL, double> mkl_sm_csr_apply_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCSRApplyTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_csr_apply_test_float;
SparseMatrixCSRApplyTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_csr_apply_test_double;
#endif

template<
  typename Mem_,
  typename Algo_,
  typename DT_>
class SparseMatrixCSRScaleTest
  : public TaggedTest<Mem_, DT_, Algo_>
{
public:
  SparseMatrixCSRScaleTest()
    : TaggedTest<Mem_, DT_, Algo_>("SparseMatrixCSRScaleTest")
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

      SparseMatrixCSR<Mem_, DT_> a(a_local);
      SparseMatrixCSR<Mem_, DT_> b(a.clone());

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
SparseMatrixCSRScaleTest<Mem::Main, Algo::Generic, float> sm_csr_scale_test_float;
SparseMatrixCSRScaleTest<Mem::Main, Algo::Generic, double> sm_csr_scale_test_double;
#ifdef FEAST_BACKENDS_MKL
SparseMatrixCSRScaleTest<Mem::Main, Algo::MKL, float> mkl_sm_csr_scale_test_float;
SparseMatrixCSRScaleTest<Mem::Main, Algo::MKL, double> mkl_sm_csr_scale_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCSRScaleTest<Mem::CUDA, Algo::CUDA, float> cuda_sm_csr_scale_test_float;
SparseMatrixCSRScaleTest<Mem::CUDA, Algo::CUDA, double> cuda_sm_csr_scale_test_double;
#endif
