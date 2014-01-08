#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/algorithm.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

template<
  typename Arch_,
  typename Algo_,
  typename DT_>
class DVSumTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  DVSumTest()
    : TaggedTest<Arch_, DT_, Algo_>("dv_sum_test")
  {
  }

  virtual void run() const
  {
    for (Index size(1) ; size < 1e5 ; size*=2)
    {
      DenseVector<Mem::Main, DT_> a_local(size);
      DenseVector<Mem::Main, DT_> b_local(size);
      DenseVector<Mem::Main, DT_> ref(size);
      DenseVector<Mem::Main, DT_> ref2(size);
      DenseVector<Mem::Main, DT_> result_local(size);
      for (Index i(0) ; i < size ; ++i)
      {
        a_local(i, DT_(i * DT_(1.234)));
        b_local(i, DT_(size*2 - i));
        ref(i, a_local(i) + b_local(i));
        ref2(i, a_local(i) + a_local(i));
      }

      DenseVector<Arch_, DT_> a(size);
      copy(a, a_local);
      DenseVector<Arch_, DT_> b(size);
      copy(b, b_local);
      DenseVector<Arch_, DT_> c(size);

      c.template sum<Algo_>(a, b);
      copy(result_local, c);
      TEST_CHECK_EQUAL(result_local, ref);

      a.template sum<Algo_>(a, b);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(a, a_local);
      b.template sum<Algo_>(a, b);
      copy(result_local, b);
      TEST_CHECK_EQUAL(result_local, ref);

      copy(b, b_local);
      a.template sum<Algo_>(a, a);
      copy(result_local, a);
      TEST_CHECK_EQUAL(result_local, ref2);
    }
  }
};
DVSumTest<Mem::Main, Algo::Generic, float> dv_sum_test_float;
DVSumTest<Mem::Main, Algo::Generic, double> dv_sum_test_double;
#ifdef FEAST_GMP
DVSumTest<Mem::Main, Algo::Generic, mpf_class> dv_sum_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_MKL
DVSumTest<Mem::Main, Algo::MKL, float> mkl_dv_sum_test_float;
DVSumTest<Mem::Main, Algo::MKL, double> mkl_dv_sum_test_double;
#endif
#ifdef FEAST_BACKENDS_CUDA
DVSumTest<Mem::CUDA, Algo::CUDA, float> cuda_dv_sum_test_float;
DVSumTest<Mem::CUDA, Algo::CUDA, double> cuda_dv_sum_test_double;
#endif

template<
  typename Arch_,
  typename Algo_,
  typename DT_,
  typename SM_>
class SMSumTest
  : public TaggedTest<Arch_, DT_, Algo_>
{

public:

  SMSumTest()
    : TaggedTest<Arch_, DT_, Algo_>("smcsr_sum_test")
  {
  }

  virtual void run() const
  {
    Index size(123);
    SparseMatrixCOO<Mem::Main, DT_> a_local(size, size + 2);
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
    SM_ b(a.clone());
    b. template scale<Algo_>(b, DT_(2));

    SM_ r(a.clone());

    r.template sum<Algo_>(a, b);

    for (Index i(0) ; i < a_local.rows() ; ++i)
    {
      for (Index j(0) ; j < a_local.columns() ; ++j)
      {
        TEST_CHECK_EQUAL(r(i, j), a(i, j) + b(i, j));
      }
    }
  }
};
SMSumTest<Mem::Main, Algo::Generic, float, SparseMatrixCSR<Mem::Main, float> > smcsr_sum_test_float;
SMSumTest<Mem::Main, Algo::Generic, double, SparseMatrixCSR<Mem::Main, double> > smcsr_sum_test_double;
#ifdef FEAST_BACKENDS_MKL
SMSumTest<Mem::Main, Algo::MKL, float, SparseMatrixCSR<Mem::Main, float> > mkl_smcsr_sum_test_float;
SMSumTest<Mem::Main, Algo::MKL, double, SparseMatrixCSR<Mem::Main, double> > mkl_smcsr_sum_test_double;
#endif
#ifdef FEAST_GMP
SMSumTest<Mem::Main, Algo::Generic, mpf_class, SparseMatrixCSR<Mem::Main, mpf_class> > smcsr_sum_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
SMSumTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixCSR<Mem::CUDA, float> > cuda_smcsr_sum_test_float;
SMSumTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixCSR<Mem::CUDA, double> > cuda_smcsr_sum_test_double;
#endif

SMSumTest<Mem::Main, Algo::Generic, float, SparseMatrixELL<Mem::Main, float> > smell_sum_test_float;
SMSumTest<Mem::Main, Algo::Generic, double, SparseMatrixELL<Mem::Main, double> > smell_sum_test_double;
#ifdef FEAST_BACKENDS_MKL
SMSumTest<Mem::Main, Algo::MKL, float, SparseMatrixELL<Mem::Main, float> > mkl_smell_sum_test_float;
SMSumTest<Mem::Main, Algo::MKL, double, SparseMatrixELL<Mem::Main, double> > mkl_smell_sum_test_double;
#endif
#ifdef FEAST_GMP
SMSumTest<Mem::Main, Algo::Generic, mpf_class, SparseMatrixELL<Mem::Main, mpf_class> > smell_sum_test_mpf_class;
#endif
#ifdef FEAST_BACKENDS_CUDA
SMSumTest<Mem::CUDA, Algo::CUDA, float, SparseMatrixELL<Mem::CUDA, float> > cuda_smell_sum_test_float;
SMSumTest<Mem::CUDA, Algo::CUDA, double, SparseMatrixELL<Mem::CUDA, double> > cuda_smell_sum_test_double;
#endif
