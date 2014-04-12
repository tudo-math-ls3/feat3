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
* \brief Test class for the sparse matrix csr class.
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
  typename Algo_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  SparseMatrixCSRTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("SparseMatrixCSRTest")
  {
  }

  virtual void run() const
  {
    SparseMatrixCSR<Mem_, DT_, IT_> zero1;
    SparseMatrixCSR<Mem::Main, DT_, IT_> zero2;
    TEST_CHECK_EQUAL(zero1, zero2);

    SparseMatrixCOO<Mem::Main, DT_, IT_> a(10, 10);
    a(1,2,7);
    a.format();
    a(1,2,7);
    a(5,5,2);
    SparseMatrixCSR<Mem_, DT_, IT_> b(a);
    TEST_CHECK_EQUAL(b.used_elements(), a.used_elements());
    TEST_CHECK_EQUAL(b.size(), a.size());
    TEST_CHECK_EQUAL(b.rows(), a.rows());
    TEST_CHECK_EQUAL(b.columns(), a.columns());
    TEST_CHECK_EQUAL(b(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(b(5, 5), a(5, 5));

    SparseMatrixCSR<Mem_, DT_, IT_> bl(b.layout());
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    bl = b.layout();
    TEST_CHECK_EQUAL(bl.used_elements(), b.used_elements());
    TEST_CHECK_EQUAL(bl.size(), b.size());
    TEST_CHECK_EQUAL(bl.rows(), b.rows());
    TEST_CHECK_EQUAL(bl.columns(), b.columns());

    typename SparseLayout<Mem_, IT_, SparseLayoutId::lt_csr>::template MatrixType<DT_> x(b.layout());
    // icc 14.0.2 does not understand the following line, so we need a typedef hier
    //typename decltype(b.layout())::template MatrixType<DT_> y(b.layout());
    typedef decltype(b.layout()) LayoutId;
    typename LayoutId::template MatrixType<DT_> y(b.layout());
    TEST_CHECK_EQUAL((void*)x.row_ptr(), (void*)b.row_ptr());


    SparseMatrixCSR<Mem_, DT_, IT_> z;
    z.convert(b);
    TEST_CHECK_EQUAL(z.used_elements(), 2ul);
    TEST_CHECK_EQUAL(z.size(), a.size());
    TEST_CHECK_EQUAL(z.rows(), a.rows());
    TEST_CHECK_EQUAL(z.columns(), a.columns());
    TEST_CHECK_EQUAL(z(1, 2), a(1, 2));
    TEST_CHECK_EQUAL(z(5, 5), a(5, 5));

    SparseMatrixCSR<Mem_, DT_, IT_> c;
    c.clone(b);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_EQUAL((void*)c.col_ind(), (void*)b.col_ind());
    c = b.clone(true);
    TEST_CHECK_NOT_EQUAL((void*)c.val(), (void*)b.val());
    TEST_CHECK_NOT_EQUAL((void*)c.col_ind(), (void*)b.col_ind());

    DenseVector<Mem_, IT_, IT_> col_ind(c.used_elements(), c.col_ind());
    DenseVector<Mem_, DT_, IT_> val(c.used_elements(), c.val());
    DenseVector<Mem_, IT_, IT_> row_ptr(c.rows() + 1, c.row_ptr());
    SparseMatrixCSR<Mem_, DT_, IT_> d(c.rows(), c.columns(), col_ind, val, row_ptr);
    TEST_CHECK_EQUAL(d, c);

    SparseMatrixCSR<Mem::Main, DT_, IT_> e;
    e.convert(c);
    TEST_CHECK_EQUAL(e, c);
    e.copy(c);
    TEST_CHECK_EQUAL(e, c);

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
    SparseMatrixCSR<Mem_, DT_, IT_> f(fcoo);

    BinaryStream bs;
    f.write_out(FileMode::fm_csr, bs);
    bs.seekg(0);
    SparseMatrixCSR<Mem_, DT_, IT_> g(FileMode::fm_csr, bs);
    TEST_CHECK_EQUAL(g, f);

    std::stringstream ts;
    f.write_out(FileMode::fm_m, ts);
    SparseMatrixCSR<Mem::Main, DT_, IT_> i(FileMode::fm_m, ts);
    TEST_CHECK_EQUAL(i, f);

    f.write_out(FileMode::fm_mtx, ts);
    SparseMatrixCSR<Mem::Main, DT_, IT_> j(FileMode::fm_mtx, ts);
    TEST_CHECK_EQUAL(j, f);
  }
};
SparseMatrixCSRTest<Mem::Main, NotSet, float, unsigned long> cpu_sparse_matrix_csr_test_float_ulong;
SparseMatrixCSRTest<Mem::Main, NotSet, double, unsigned long> cpu_sparse_matrix_csr_test_double_ulong;
SparseMatrixCSRTest<Mem::Main, NotSet, float, unsigned int> cpu_sparse_matrix_csr_test_float_uint;
SparseMatrixCSRTest<Mem::Main, NotSet, double, unsigned int> cpu_sparse_matrix_csr_test_double_uint;
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCSRTest<Mem::CUDA, NotSet, float, unsigned long> cuda_sparse_matrix_csr_test_float_ulong;
SparseMatrixCSRTest<Mem::CUDA, NotSet, double, unsigned long> cuda_sparse_matrix_csr_test_double_ulong;
SparseMatrixCSRTest<Mem::CUDA, NotSet, float, unsigned int> cuda_sparse_matrix_csr_test_float_uint;
SparseMatrixCSRTest<Mem::CUDA, NotSet, double, unsigned int> cuda_sparse_matrix_csr_test_double_uint;
#endif


template<
  typename Mem_,
  typename Algo_,
  typename DT_,
  typename IT_>
class SparseMatrixCSRApplyTest
  : public FullTaggedTest<Mem_, Algo_, DT_, IT_>
{
public:
  SparseMatrixCSRApplyTest()
    : FullTaggedTest<Mem_, Algo_, DT_, IT_>("SparseMatrixCSRApplyTest")
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
      SparseMatrixCSR<Mem_,DT_, IT_> a(a_local);

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

SparseMatrixCSRApplyTest<Mem::Main, Algo::Generic, float, unsigned long> sm_csr_apply_test_float_ulong;
SparseMatrixCSRApplyTest<Mem::Main, Algo::Generic, double, unsigned long> sm_csr_apply_test_double_ulong;
SparseMatrixCSRApplyTest<Mem::Main, Algo::Generic, float, unsigned int> sm_csr_apply_test_float_uint;
SparseMatrixCSRApplyTest<Mem::Main, Algo::Generic, double, unsigned int> sm_csr_apply_test_double_uint;
#ifdef HONEI_BACKENDS_MKL
SparseMatrixCSRApplyTest<Mem::Main, Algo::MKL, float, unsigned long> mkl_sm_csr_apply_test_float_ulong;
SparseMatrixCSRApplyTest<Mem::Main, Algo::MKL, double, unsigned long> mkl_sm_csr_apply_test_double_ulong;
#endif
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCSRApplyTest<Mem::CUDA, Algo::CUDA, float, unsigned long> cuda_sm_csr_apply_test_float_ulong;
SparseMatrixCSRApplyTest<Mem::CUDA, Algo::CUDA, double, unsigned long> cuda_sm_csr_apply_test_double_ulong;
SparseMatrixCSRApplyTest<Mem::CUDA, Algo::CUDA, float, unsigned int> cuda_sm_csr_apply_test_float_uint;
SparseMatrixCSRApplyTest<Mem::CUDA, Algo::CUDA, double, unsigned int> cuda_sm_csr_apply_test_double_uint;
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

      SparseMatrixCSR<Mem_, DT_> ref(ref_local);
      SparseMatrixCSR<Mem_, DT_> a(a_local);
      SparseMatrixCSR<Mem_, DT_> b;
      b.clone(a);

      b.template scale<Algo_>(a, s);
      TEST_CHECK_EQUAL(b, ref);

      a.template scale<Algo_>(a, s);
      TEST_CHECK_EQUAL(a, ref);
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
